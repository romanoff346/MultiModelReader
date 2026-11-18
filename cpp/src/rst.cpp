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
    std::fseek(fp, 0, SEEK_END);
    file_bytes_ = std::ftell(fp);
    std::fseek(fp, 0, SEEK_SET);

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

    // ---- Компоненты (name selections): узловые и элементные группы ----
    // Layout PTRCMS: одна FORTRAN-запись
    //   [NCMS,NNAM,NVAL] int32 | имена NCMS*(NNAM+1) байт | маркеры типов 'N'/'E'
    //   | (выравнивание) | для каждой группы: [cnt][id...], id = int32
    if (ptr_comp != 0) {
        if (read_record(ptr_comp, b, w, err)) {
            ncms_hdr = rd_i32(b.data() + 0);
            nnam_hdr = rd_i32(b.data() + 4);
            parse_components(b.data(), (int64_t)b.size());
        }
        err.clear();                                  // компоненты — опциональны
    }

    // ---- Смещение таблицы указателей элементов (для фильтрации по элементам) ----
    if (ptr_elem != 0 && read_record(ptr_elem, b, w, err)) {
        n_elem_ptr = rd_i64(b.data());
        ptr_elem_words = ptr_elem;
        nelem_hint_ = n_elem_ptr;
    }
    // ELMAP (обычно PTR(7)): [NELMAX][NELEM] — точное число элементов
    const int64_t ptr_elmap = ptr(7);
    if (ptr_elmap != 0 && read_record(ptr_elmap, b, w, err)) {
        int32_t nelmax = rd_i32(b.data() + 0);
        int32_t nelem  = rd_i32(b.data() + 4);
        if (nelem > 0) nelem_hint_ = nelem;
        else if (nelmax > 0) nelem_hint_ = nelmax;
    }
    err.clear();
    return true;
}

// Разбор склеенного блока компонентов: определяет тип каждой группы и
// заполняет Component::nodes / Component::elems.
void File::parse_components(const uint8_t* p, int64_t total) {
    components.clear();
    if (total < 12 || ncms_hdr <= 0) return;
    const int32_t ncms = ncms_hdr, nnam = nnam_hdr;
    const int64_t names_bytes = (int64_t)ncms * (nnam + 1);

    // Маркеры типов идут сразу после имён (по 4 байта на компонент).
    auto type_at = [&](int32_t c) -> char {
        int64_t off = 12 + names_bytes + (int64_t)c * 4;
        if (off + 4 > total) return '?';
        char ch = (char)p[off];
        if (ch == 'N' || ch == 'E') return ch;   // узловой / элементный блок
        return '?';                              // неизвестно — считаем узловым
    };

    struct Grp { std::string name; char type; };
    std::vector<Grp> grps;
    for (int32_t c = 0; c < ncms; ++c) {
        std::string nm((const char*)p + 12 + (int64_t)c * (nnam + 1), nnam);
        // Обрезаем хвостовые пробелы и \0, сохраняя внутренние пробелы имён
        auto end = std::find_if(nm.rbegin(), nm.rend(),
                                [](char ch){ return ch != ' ' && ch != '\0'; });
        if (end == nm.rend()) nm.clear();
        else nm.resize((size_t)(nm.rend() - end));
        if (nm.empty()) continue;
        grps.push_back({nm, type_at(c)});
    }

    // Списки id идут после блока типов. Стартовый сдвиг подбирается подбором
    // padding (выравнивание на слово): читаем группы последовательно и сверяем
    // счётчики с лимитами (узлы <= Glbnnod, элементы <= NELEM) и размером буфера.
    int64_t lists_base = 12 + names_bytes + (int64_t)ncms * 4;
    bool parsed = false;
    for (int64_t pad : {0LL, 1LL, 2LL, 3LL}) {         // выравнивание на слово
        int64_t try_off = lists_base + pad;
        if (try_off + 4 > total) continue;
        int64_t pos = try_off;
        bool ok = true;
        std::vector<std::vector<int32_t>> ids(grps.size());
        for (size_t g = 0; g < grps.size() && ok; ++g) {
            if (pos + 4 > total) { ok = false; break; }
            int32_t cnt = rd_i32(p + pos); pos += 4;
            int64_t lim = (grps[g].type == 'E') ? nelem_hint_ : n_nodes;
            if (lim <= 0) lim = 100000000;
            if (cnt < 0 || cnt > lim || pos + (int64_t)cnt * 4 > total) { ok = false; break; }
            ids[g].resize(cnt);
            for (int32_t k = 0; k < cnt; ++k) ids[g][k] = rd_i32(p + pos + k * 4);
            pos += (int64_t)cnt * 4;
        }
        if (ok) {
            for (size_t g = 0; g < grps.size(); ++g) {
                Component comp; comp.name = grps[g].name;
                if (grps[g].type == 'E') { comp.is_elem = true; comp.elems = std::move(ids[g]); }
                else                     { comp.is_node = true; comp.nodes = std::move(ids[g]); }
                components.push_back(std::move(comp));
            }
            parsed = true; break;
        }
    }
    if (!parsed) {                                    // fallback: без маркеров типов (старые файлы)
        int64_t nodes_off = lists_base;
        for (int64_t pad : {0LL, 1LL, 2LL, 3LL}) {
            int64_t co = 12 + names_bytes + pad;      // без маркеров типов (старые файлы)
            bool ok = true;
            for (size_t g = 0; g < grps.size() && ok; ++g) {
                int64_t o = co + (int64_t)g * 4;
                if (o + 4 > total) { ok = false; break; }
                int32_t cnt = rd_i32(p + o);
                if (cnt < 0 || (int64_t)cnt > n_nodes) ok = false;
            }
            if (ok) { nodes_off = co + (int64_t)ncms * 4; break; }
        }
        int64_t pos = nodes_off;
        for (size_t g = 0; g < grps.size(); ++g) {
            if (pos + 4 > total) break;
            int32_t cnt = rd_i32(p + pos); pos += 4;
            if (cnt < 0 || pos + (int64_t)cnt * 4 > total) break;
            Component comp; comp.name = grps[g].name; comp.is_node = true;
            comp.nodes.resize(cnt);
            for (int32_t k = 0; k < cnt; ++k) comp.nodes[k] = rd_i32(p + pos + k * 4);
            pos += (int64_t)cnt * 4;
            components.push_back(std::move(comp));
        }
    }
}

// ---- Ленивая загрузка связности элементов (ELMAP + ELCONN) ----
bool File::load_elements(std::string& err) {
    if (elements_loaded) return true;
    if (ptr_elem_words == 0) { err = "нет таблицы элементов"; return false; }
    std::vector<uint8_t> b; int64_t w = 0;
    if (!read_record(ptr_elem_words, b, w, err)) return false;
    // ELMAP: абсолютное словесное смещение таблицы элементов
    const int64_t elmmap_abs = (int64_t)llround(rd_f64(b.data()));
    if (!read_record(elmmap_abs, b, w, err)) { err = "read ELMAP fail"; return false; }
    // Запись ELMAP: строки по 5 слов? — формат варьируется; берём первый int64 как NELEM,
    // дальше id->offset таблица. Читаем defensively.
    elem_tab.max_eid = (int32_t)n_elem_ptr;
    if (elem_tab.max_eid <= 0) elem_tab.max_eid = (int32_t)rd_i64(b.data());
    // Таблица указателей соединений внутри ELMAP-блока:
    int64_t conn_abs = 0;
    for (int probe = 1; probe + 1 < w; ++probe) {
        int64_t v = (int64_t)llround(rd_f64(b.data() + probe * 8));
        if (v > elmmap_abs && v * kWord < file_bytes_) { conn_abs = v; break; }
    }
    if (conn_abs == 0) { err = "ELCONN не найден"; return false; }
    if (!read_record(conn_abs, b, w, err)) { err = "read ELCONN fail"; return false; }
    // ELCONN: [max eid][node1..node8 per element], int32, stride 9 слов? —
    // стандартный layout: 8 узлов на элемент + id. Парсим по 8 int32 на запись.
    elem_tab.e2n.assign(((size_t)elem_tab.max_eid + 1) * 8, 0);
    elem_tab.ncomp.assign((size_t)elem_tab.max_eid + 1, 0);
    size_t avail = (size_t)w;                         // слов в теле
    int32_t eid = 0;
    for (size_t i = 0; i + 8 <= avail && eid < elem_tab.max_eid; i += 8) {
        ++eid;
        for (int k = 0; k < 8; ++k)
            elem_tab.e2n[(size_t)eid * 8 + k] = rd_i32(b.data() + (i + k) * 4);
        for (int k = 7; k >= 0 && elem_tab.e2n[(size_t)eid * 8 + k] == 0; --k) {}
    }
    elements_loaded = true;
    err.clear();
    return true;
}

bool File::expand_elem_component(const Component& c, std::vector<int32_t>& out_nodes,
                                 std::string& err) {
    if (!load_elements(err)) return false;
    out_nodes.clear();
    for (int32_t e : c.elems) {
        if (e <= 0 || e > elem_tab.max_eid) continue;
        for (int k = 0; k < 8; ++k) {
            int32_t n = elem_tab.e2n[(size_t)e * 8 + k];
            if (n > 0) out_nodes.push_back(n);
        }
    }
    std::sort(out_nodes.begin(), out_nodes.end());
    out_nodes.erase(std::unique(out_nodes.begin(), out_nodes.end()), out_nodes.end());
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
