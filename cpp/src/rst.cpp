// rst.cpp — реализация заголовочного парсера и параллельного чтения .rst
#include "rst.h"

#include <algorithm>
// (strcasecmp не используется; избегаем POSIX-специфичных заголовков ради MSVC)

namespace rst {

// Прочитать одну FORTRAN-запись по словесному смещению (body без маркеров длины)
bool File::read_record(int64_t word_off, std::vector<uint8_t>& buf,
                       int64_t& words_out, std::string& err) {
    const int64_t byte_off = word_off * kWord;
    int32_t len_pre = 0;
    if (std::fseek(fp, (long)byte_off, SEEK_SET) != 0) { err = "fseek fail"; return false; }
    if (std::fread(&len_pre, 4, 1, fp) != 1) { err = "read reclen fail"; return false; }
    if (len_pre <= 0) { err = "bad record length " + std::to_string(len_pre); return false; }
    const int64_t body_bytes = (int64_t)len_pre * kWord;
    buf.resize((size_t)body_bytes);
    if (std::fread(buf.data(), 1, (size_t)body_bytes, fp) != (size_t)body_bytes) {
        err = "read body fail"; return false;
    }
    int32_t len_post = 0;
    if (std::fread(&len_post, 4, 1, fp) != 1) { err = "read reclen2 fail"; return false; }
    if (len_post != len_pre) { err = "record markers mismatch"; return false; }
    words_out = len_pre;
    return true;
}

void File::close() {
    if (fp) { std::fclose(fp); fp = nullptr; }
}

bool File::open(const std::string& p, std::string& err) {
    path = p;
    fp = std::fopen(p.c_str(), "rb");
    if (!fp) { err = "cannot open " + p; return false; }

    // ---- Запись 1: указатели на таблицы заголовка ----
    std::vector<uint8_t> b; int64_t w = 0;
    if (!read_record(0, b, w, err)) return false;
    auto ptr = [&](int idx) -> int64_t {           // idx — 1-й элемент Фортрана
        return (int64_t)llround(rd_f64(b.data() + (idx - 1) * 8));
    };
    // Layout: [1]=PTRPTR, [2]=PTRDAT, [3]=PRSOL, [4]=PTRGEO, [5]=PTRELM, [6]=PTRCMS
    const int64_t ptr_ptr      = ptr(1); (void)ptr_ptr;   // таблица указателей (манифест)
    const int64_t ptr_data     = ptr(2);   // таблица данных
    const int64_t ptr_solution = ptr(3);   // таблица решений (RSPEC)
    const int64_t ptr_geom     = ptr(4);   // геометрия/узлы
    const int64_t ptr_elem     = ptr(5);   // элементы
    const int64_t ptr_comp     = ptr(6);   // name selections (компоненты)
    ptr_geo_words  = ptr_geom;
    ptr_comp_words = ptr_comp;

    // ---- Таблица данных: Glbnnod, Maxcomp, ... ----
    if (!read_record(ptr_data, b, w, err)) return false;
    n_nodes  = (int64_t)llround(rd_f64(b.data() + 0 * 8));
    maxcomp  = (int32_t)llround(rd_f64(b.data() + 2 * 8));

    // ---- Таблица решений: Nsets, RPNTR, TPRNT3D ----
    if (!read_record(ptr_solution, b, w, err)) return false;
    nsets = (int32_t)llround(rd_f64(b.data() + 0 * 8));
    const int64_t rpntr_abs  = (int64_t)llround(rd_f64(b.data() + 1 * 8));
    const int64_t tprnt3d_w  = (int64_t)llround(rd_f64(b.data() + 2 * 8));  // слова от начала файла
    if (nsets <= 0 || nsets > 10000000) { err = "bad NSETS"; return false; }

    // ---- RPNTR: абсолютные словесные смещения заголовков решений ----
    if (!read_record(rpntr_abs, b, w, err)) return false;
    sets.resize(nsets);
    for (int i = 0; i < nsets; ++i)
        sets[i].rpointer_words = (int64_t)llround(rd_f64(b.data() + i * 8));

    // ---- TPRNT3D: шаг/подшаг/время каждого набора ----
    // Таблица может быть разбита на несколько FORTRAN-записей.
    int64_t off = tprnt3d_w;
    for (int i = 0; i < nsets; ) {
        if (!read_record(off, b, w, err)) return false;
        int64_t nrows = w / 3;                        // 3 слова на строку
        for (int64_t r = 0; r < nrows && i < nsets; ++r, ++i) {
            sets[i].step    = rd_i32(b.data() + r * 12 + 0);
            sets[i].substep = rd_i32(b.data() + r * 12 + 4);
            sets[i].time    = (double)rd_i32(b.data() + r * 12 + 8);
        }
        off += w + 2;                                 // пропуск маркеров длины
    }

    // ---- Компоненты (name selections): имена узловых групп ----
    // Layout записи PTRCMS: [NCMS,NNAM,NVAL][names NCMS*NNAM байт][counts NCMS*i32]
    //                       [node lists int32...][elem lists ...]
    if (ptr_comp != 0) {
        if (read_record(ptr_comp, b, w, err)) {
            const uint8_t* p = b.data();
            int32_t ncms = rd_i32(p + 0);
            int32_t nnam = rd_i32(p + 4);
            // Автовыравнивание блока имён: счётчики идут сразу за именами,
            // поэтому сдвиг блока counts подбираем по условию cnt <= Glbnnod.
            int64_t names_bytes = (int64_t)ncms * (nnam + 1);
            int64_t counts_off = -1;
            for (int64_t pad : {0LL, 1LL, 2LL, 3LL}) {
                int64_t co = 12 + names_bytes + pad;
                bool ok = true;
                for (int32_t c = 0; c < ncms && ok; ++c) {
                    int32_t cnt = rd_i32(p + co + c * 4);
                    if (cnt < 0 || (int64_t)cnt > n_nodes) ok = false;
                }
                if (ok) { counts_off = co; break; }
            }
            if (counts_off < 0) counts_off = 12 + names_bytes;
            components.clear();
            int64_t nodes_off = counts_off + (int64_t)ncms * 4;
            int64_t cur = nodes_off;
            for (int32_t c = 0; c < ncms; ++c) {
                std::string nm((const char*)p + 12 + (int64_t)c * (nnam + 1), nnam);
                auto end = std::find_if(nm.rbegin(), nm.rend(),
                                        [](char ch){ return ch != ' ' && ch != '\0'; });
                nm.resize((size_t)(end - nm.rbegin()));
                if (nm.empty()) continue;
                Component comp; comp.name = nm;
                int32_t cnt = rd_i32(p + counts_off + c * 4);
                comp.nodes.resize(cnt);
                for (int32_t k = 0; k < cnt; ++k)
                    comp.nodes[k] = rd_i32(p + cur + k * 4);
                cur += (int64_t)cnt * 4;
                components.push_back(std::move(comp));
            }
            err.clear();                              // компоненты — опциональны
        }
    }

    // ---- Смещение таблицы указателей элементов (для фильтрации по элементам) ----
    if (ptr_elem != 0 && read_record(ptr_elem, b, w, err)) {
        n_elem_ptr = rd_i64(b.data());
    }
    err.clear();
    return true;
}

int File::find_set_by_substep(int substep) const {
    for (size_t i = 0; i < sets.size(); ++i)
        if (sets[i].substep == substep) return (int)i;
    return -1;
}

int File::find_set_by_time(double t, double tol) const {
    int best = -1; double bd = tol;
    for (size_t i = 0; i < sets.size(); ++i) {
        double d = std::fabs(sets[i].time - t);
        if (d <= bd) { bd = d; best = (int)i; }
    }
    return best;
}

} // namespace rst
