// reader.cpp — параллельное чтение наборов результатов из .rst
//
// Модель параллелизма: каждый поток берёт из очереди свой набор (set) и
// читает его целиком через независимый файловый дескриптор + POSIX pread()
// (файл на SSD, случайные позиционные чтения масштабируются почти линейно).
#include "reader.h"

#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <mutex>
#include <queue>

namespace rst {

// ---- вспомогательный ридер FORTRAN-записей поверх pread (без общего seek) ----
namespace {

struct PreadCtx {
    int fd = -1;
    std::vector<uint8_t> buf;     // сырой кадр: [4][body][4]
    int64_t words = 0;            // длина тела в словах
    const uint8_t* body() const { return buf.data() + 4; }

    bool read(int64_t word_off, std::string& err) {
        int32_t len = 0;
        ssize_t n = ::pread(fd, &len, 4, word_off * kWord);
        if (n != 4) { err = "pread reclen"; return false; }
        if (len <= 0 || len > (1 << 28)) { err = "bad reclen " + std::to_string(len); return false; }
        size_t total = 8 + (size_t)len * kWord;
        buf.resize(total);
        memcpy(buf.data(), &len, 4);
        n = ::pread(fd, buf.data() + 4, total - 8, word_off * kWord + 4);
        if ((size_t)n != total - 8) { err = "pread body"; return false; }
        int32_t len2; memcpy(&len2, buf.data() + 4 + (size_t)len * kWord, 4);
        if (len2 != len) { err = "reclen mismatch"; return false; }
        words = len;
        return true;
    }
};

inline double f64(const uint8_t* p, int64_t idx) { return rd_f64(p + idx * 8); }
inline int32_t i32(const uint8_t* p, int64_t idx) { return rd_i32(p + idx * 4); }
// 64-битное смещение (в словах), хранящееся как FORTRAN REAL*8
inline int64_t f64_as_word(const uint8_t* p, int64_t idx) {
    return (int64_t)llround(rd_f64(p + idx * 8));
}

} // namespace

SetData read_set(File& meta, int set_index, const std::vector<int32_t>* node_filter) {
    SetData out; out.set_index = set_index;
    PreadCtx c;
    c.fd = ::open(meta.path.c_str(), O_RDONLY);
    if (c.fd < 0) { out.err = "open failed"; return out; }
    std::string err;

    const ResultSet& rs = meta.sets[set_index];

    // 1. Заголовок решения (RSPEC по RPNTR): смещения таблиц внутри набора
    if (!c.read(rs.rpointer_words, err)) { out.err = "RSPEC: " + err; ::close(c.fd); return out; }
    const uint8_t* hdr = c.body();
    // RSPEC: [0]=не исп., [1]=DofBlk — таблица заголовков DOF (смещение в словах, REAL*8)
    const int64_t ptr_dof_table = f64_as_word(hdr, 1);

    // 2. Таблица DOF: для каждого узлового DOF — тип, номер узла, смещение данных
    if (!c.read(ptr_dof_table, err)) { out.err = "DOF table: " + err; ::close(c.fd); return out; }
    struct DofInfo { int32_t node; int64_t data_word; };
    std::vector<DofInfo> ux, uy, uz, sxy, syz, sxz;
    int32_t max_node_seen = 0;
    for (int64_t r = 0; r + 2 < c.words; r += 3) {
        int32_t type = i32(c.body(), r + 0);
        int32_t node = i32(c.body(), r + 1);
        int64_t dw   = (int64_t)i32(c.body(), r + 2);
        auto& vec = [&]() -> std::vector<DofInfo>& {
            switch (type) {
                case 1: return ux;  case 2: return uy;  case 3: return uz;
                case 9: return sxy; case 10: return syz; case 11: return sxz;
                default: static std::vector<DofInfo> dummy; return dummy;
            }
        }();
        if (type == 1 || type == 2 || type == 3 ||
            type == 9 || type == 10 || type == 11)
            vec.push_back({node, dw});
        max_node_seen = std::max(max_node_seen, node);
    }

    // 3. Формируем список узлов вывода (фильтр компонента или все узлы)
    std::vector<int32_t> nodes_out;
    if (node_filter) nodes_out = *node_filter;
    else { nodes_out.reserve(ux.size()); for (auto& d : ux) nodes_out.push_back(d.node); }
    std::sort(nodes_out.begin(), nodes_out.end());
    nodes_out.erase(std::unique(nodes_out.begin(), nodes_out.end()), nodes_out.end());

    // Индексация: узел -> позиция в vec (упорядочены по возрастанию)
    auto build_index = [](const std::vector<DofInfo>& v) {
        std::vector<std::pair<int32_t,int64_t>> m; m.reserve(v.size());
        for (auto& d : v) m.push_back({d.node, d.data_word});
        std::sort(m.begin(), m.end());
        return m;
    };
    auto IUX = build_index(ux), IUY = build_index(uy), IUZ = build_index(uz);
    auto ISX = build_index(sxy), ISY = build_index(syz), ISZ = build_index(sxz);

    auto lookup = [&](const std::vector<std::pair<int32_t,int64_t>>& m, int32_t n) -> int64_t {
        auto it = std::lower_bound(m.begin(), m.end(), std::make_pair(n, (int64_t)0));
        return (it != m.end() && it->first == n) ? it->second : -1;
    };

    // 4. Чтение значений напряжений по одному pread на компонент (узлы отсортированы,
    //    данные в файле тоже по узлам — последовательный паттерн внутри набора)
    out.node_ids = nodes_out;
    out.stress.assign(nodes_out.size() * 6, 0.0);

    auto load_component = [&](const std::vector<std::pair<int32_t,int64_t>>& idx,
                              int comp_slot) -> bool {
        for (size_t k = 0; k < nodes_out.size(); ++k) {
            int64_t dw = lookup(idx, nodes_out[k]);
            if (dw < 0) continue;
            if (!c.read(dw, err)) return false;
            out.stress[k * 6 + comp_slot] = f64(c.body(), 0);
        }
        return true;
    };

    bool ok = load_component(IUX, 0) && load_component(IUY, 1) && load_component(IUZ, 2)
           && load_component(ISX, 3) && load_component(ISY, 4) && load_component(ISZ, 5);
    ::close(c.fd);
    out.ok = ok;
    if (!ok) out.err = err;
    return out;
}

// ---- планировщик: пул потоков, очередь наборов, упорядоченный CSV-вывод ----
int run_parallel(File& meta, Options opt) {   // копия opt — можно менять component
    int nthreads = opt.threads > 0 ? opt.threads
        : (int)std::thread::hardware_concurrency();
    if (nthreads <= 0) nthreads = 4;
    nthreads = std::min<int>(nthreads, (int)opt.sets.size());

    // фильтр узлов по компоненту
    const std::vector<int32_t>* filter = nullptr;
    if (!opt.component.empty()) {
        for (auto& c : meta.components)
            if (c.name == opt.component) { filter = &c.nodes; break; }
        if (!filter) {
            std::transform(opt.component.begin(), opt.component.end(),
                           opt.component.begin(), ::toupper);
            for (auto& c : meta.components)
                if (c.name == opt.component) { filter = &c.nodes; break; }
        }
        if (!filter) {
            std::fprintf(stderr, "Компонент '%s' не найден.\nДоступные:", opt.component.c_str());
            for (auto& c : meta.components) std::fprintf(stderr, " %s", c.name.c_str());
            std::fprintf(stderr, "\n");
            return 1;
        }
    }

    std::vector<SetData> results(opt.sets.size());
    std::atomic<size_t> next{0};

    auto worker = [&]() {
        for (;;) {
            size_t t = next.fetch_add(1);
            if (t >= opt.sets.size()) break;
            results[t] = read_set(meta, opt.sets[t], filter);
        }
    };

    std::vector<std::thread> pool;
    for (int i = 0; i < nthreads; ++i) pool.emplace_back(worker);
    for (auto& th : pool) th.join();

    // вывод
    FILE* fo = stdout;
    if (!opt.out.empty()) {
        fo = fopen(opt.out.c_str(), "w");
        if (!fo) { std::fprintf(stderr, "Не создать %s\n", opt.out.c_str()); return 1; }
    }
    std::fprintf(fo, "set,node,SX,SY,SZ,SXY,SYZ,SXZ%s\n", opt.mises ? ",MISES" : "");
    int rc = 0;
    for (auto& sd : results) {
        if (!sd.ok) {
            std::fprintf(stderr, "set %d: ошибка: %s\n", sd.set_index, sd.err.c_str());
            rc = 1; continue;
        }
        for (size_t k = 0; k < sd.node_ids.size(); ++k) {
            const double* s = &sd.stress[k * 6];
            double vm = 0;
            if (opt.mises)
                vm = std::sqrt(0.5 * ((s[0]-s[1])*(s[0]-s[1]) + (s[1]-s[2])*(s[1]-s[2]) + (s[0]-s[2])*(s[0]-s[2]))
                               + 3.0 * (s[3]*s[3] + s[4]*s[4] + s[5]*s[5]));
            std::fprintf(fo, "%d,%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e%s%.6e%s\n",
                         sd.set_index, sd.node_ids[k], s[0], s[1], s[2], s[3], s[4], s[5],
                         opt.mises ? "," : "", opt.mises ? vm : 0.0, "");
        }
    }
    if (fo != stdout) fclose(fo);
    return rc;
}

} // namespace rst
