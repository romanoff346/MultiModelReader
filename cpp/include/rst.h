// rst.h — структуры и API быстрого ридера файлов результатов ANSYS MAPDL (.rst/.rnn)
//
// Формат: FORTRAN unformatted sequential, little-endian, double precision.
// Все смещения в заголовках указаны в 4-байтных словах (word = 4 байта).
// Запись на диске: [int32 reclen_words][body][int32 reclen_words].
//
// Поддерживаются как НЕсжатые RST (MAPDL OUTRES,...,NOCOMPRESS / COMPRESS,OFF),
// так и СЖАТЫЕ записи формата ANSCOMP (OUTRES,...,,ANSYS; по умолчанию в MAPDL
// >= 2023 R1): сжатая FORTRAN-запись помечается магическим словом 'ANS-' в
// начале тела и распаковывается встроенным DEFLATE-декодером (zlite.cpp),
// внешняя библиотека zlib не требуется.
#pragma once

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <map>
#include <mutex>
#include <string>
#include <vector>

// ---- кроссплатформенный позиционный ввод (SSD-friendly, без общего seek) ----
#if defined(_WIN32)
  #define WIN32_LEAN_AND_MEAN
  #include <windows.h>
#else
  #include <fcntl.h>
  #include <unistd.h>
#endif

namespace rst {

constexpr int kWord = 4;

// Магическое слово сжатой FORTRAN-записи ANSYS (ANSCOMP): первые 8 байт
// тела — 'ANS-' + номер версии формата ('1'..'3'). reclen при этом = длина
// СЖАТОГО кадра.
inline bool is_comp_magic(const void* p) {
    const uint8_t* b = (const uint8_t*)p;
    return std::memcmp(b, "ANS-", 4) == 0 && b[4] >= '1' && b[4] <= '9';
}

// Нативный файловый дескриптор для pread-подобного чтения
using NativeFd =
#if defined(_WIN32)
    HANDLE;
#else
    int;
#endif

inline NativeFd fd_open(const std::string& p) {
#if defined(_WIN32)
    return ::CreateFileA(p.c_str(), GENERIC_READ,
                         FILE_SHARE_READ | FILE_SHARE_WRITE,
                         nullptr, OPEN_EXISTING,
                         FILE_ATTRIBUTE_NORMAL | FILE_FLAG_SEQUENTIAL_SCAN,
                         nullptr);
#else
    return ::open(p.c_str(), O_RDONLY);
#endif
}
inline void fd_close(NativeFd fd) {
#if defined(_WIN32)
    if (fd != INVALID_HANDLE_VALUE) ::CloseHandle(fd);
#else
    if (fd >= 0) ::close(fd);
#endif
}
inline bool fd_valid(NativeFd fd) {
#if defined(_WIN32)
    return fd != INVALID_HANDLE_VALUE;
#else
    return fd >= 0;
#endif
}
// Чтение count байт по абсолютному смещению без изменения позиции файла.
// Возвращает true, если прочитано ровно count байт.
inline bool fd_pread(NativeFd fd, void* dst, size_t count, int64_t offset) {
#if defined(_WIN32)
    OVERLAPPED ov{};
    LARGE_INTEGER li; li.QuadPart = offset;
    ov.Offset = li.u.LowPart;
    ov.OffsetHigh = li.u.HighPart;
    DWORD got = 0;
    if (!::ReadFile(fd, dst, (DWORD)count, &got, &ov)) return false;
    return got == (DWORD)count;
#else
    size_t done = 0;
    while (done < count) {
        ssize_t n = ::pread(fd, (char*)dst + done, count - done, offset + (int64_t)done);
        if (n <= 0) return false;
        done += (size_t)n;
    }
    return true;
#endif
}

inline int64_t rd_i64(const uint8_t* p) {
    int64_t v; std::memcpy(&v, p, 8); return v;
}
inline int32_t rd_i32(const uint8_t* p) {
    int32_t v; std::memcpy(&v, p, 4); return v;
}
inline double rd_f64(const uint8_t* p) {
    double v; std::memcpy(&v, p, 8); return v;
}

// Одна запись результата (set): шаг/подшаг, время, указатели
struct ResultSet {
    int32_t step = 0;
    int32_t substep = 0;
    double  time = 0.0;
    int64_t rpointer_words = 0;   // абсолютное смещение заголовка решения (words)
};

struct Component {
    std::string name;             // name selection (имя блока *CM)
    std::vector<int32_t> nodes;   // заполняется для узловых компонентов
    std::vector<int32_t> elems;   // заполняется для элементных компонентов
    bool is_node = false;         // тип: узловой блок?
    bool is_elem = false;         // тип: элементный блок?
};

// Таблица элементов: id ->Connectivity (1-based ids, ноды 1-based)
struct ElemTable {
    int32_t max_eid = 0;
    std::vector<int32_t> e2n;     // [eid*8 .. eid*8+7] — ноды элемента (0 = нет)
    std::vector<int32_t> ncomp;   // кол-во узлов на элемент
};

class File {
public:
    bool open(const std::string& path, std::string& err);
    void close();

    // ---- метаданные (заполняются в open()) ----
    int      nsets = 0;
    int64_t  n_nodes = 0;         // узлов в наборе решения (Glbnnod для split-файлов)
    int64_t  n_elem_ptr = 0;      // размер таблицы указателей элементов
    std::vector<ResultSet> sets;
    std::vector<Component> components;

    // Смещения таблиц геометрии (words)
    int64_t ptr_geo_words = 0;
    int64_t ptr_elem_words = 0;   // таблица указателей элементов
    int64_t ptr_comp_words = 0;
    int32_t maxcomp = 0;

    // Таблица связности элементов (заполняется лениво, load_elements())
    ElemTable elem_tab;
    bool elements_loaded = false;
    bool load_elements(std::string& err);           // читает ELMAP/ELCONN
    // Развернуть элементный компонент в уникальные отсортированные узлы
    bool expand_elem_component(const Component& c, std::vector<int32_t>& out_nodes,
                               std::string& err);

private:
    void parse_components(const uint8_t* p, int64_t total);
    int32_t ncms_hdr = 0, nnam_hdr = 0;             // заголовок блока PTRCMS
    int64_t file_bytes_ = 0;                        // размер файла для границ обхода
    int64_t nelem_hint_ = 0;                        // оценка числа элементов

    // ---- поддержка сжатых записей (ANSCOMP) ----
    bool compressed_ = false;                       // файл содержит сжатые записи
    // Кэш распакованных записей: word_off -> тело без padding. Нужен потому,
    // что смещения указателей в ANSCOMP-файлах — это адреса СЖАТЫХ кадров,
    // а таблицы внутри набора склеиваются из нескольких записей подряд
    // ("off += words + 2" работает только при известной распакованной длине).
    std::map<int64_t, std::vector<uint8_t>> unc_cache_;
    mutable std::mutex unc_mtx_;

public:
    FILE* fp = nullptr;
    std::string path;
    bool compressed() const { return compressed_; }

    // Прочитать одну FORTRAN-запись по словесному смещению.
    // Возвращает тело записи в buf (без маркеров длины); words_out — длина в словах.
    bool read_record(int64_t word_off, std::vector<uint8_t>& buf,
                     int64_t& words_out, std::string& err);

    // То же, но через потокобезопасный позиционный ввод (для reader.cpp),
    // с прозрачной декомпрессией ANSCOMP-записей и кэшем распакованных тел.
    // out — тело без padding; возвращает длину тела в словах (-1 при ошибке).
    int64_t read_pread(NativeFd fd, int64_t word_off, std::vector<uint8_t>& out,
                       std::string& err);

    // Найти запись результата по номеру substep (1-based) или по времени.
    int find_set_by_substep(int substep) const;
    int find_set_by_time(double time, double tol = 1e-9) const;
};

} // namespace rst
