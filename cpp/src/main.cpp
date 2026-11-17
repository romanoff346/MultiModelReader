// main.cpp — CLI: параллельное чтение узловых напряжений из ANSYS .rst
//
// Примеры:
//   ./rst_reader model.rst --list                       # список наборов и CM
//   ./rst_reader model.rst --substep 5 --component FIXTURE_NODES --mises --out s5.csv
//   ./rst_reader model.rst --sets 1-8 --all-nodes --threads 8 --out all.csv
#include "rst.h"
#include "reader.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using namespace rst;

static void usage() {
    std::printf(
        "Использование: rst_reader FILE.rst [опции]\n"
        "  --list                 показать наборы результатов и name selections\n"
        "  --set N                прочитать набор с индексом N (0-based)\n"
        "  --sets A-B             прочитать наборы A..B включительно\n"
        "  --substep N            найти набор по номеру подшага\n"
        "  --time T               найти ближайший набор по времени\n"
        "  --component NAME       узлы из name selection NAME\n"
        "  --all-nodes            все узлы модели\n"
        "  --mises                вычислять эквивалентное напряжение Мизеса\n"
        "  --threads K            число потоков чтения (по умолчанию: ядра CPU)\n"
        "  --out FILE             CSV-выход (по умолчанию stdout)\n");
}

int main(int argc, char** argv) {
    if (argc < 2) { usage(); return 2; }
    Options opt;
    opt.path = argv[1];
    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&](const char* what) -> const char* {
            if (i + 1 >= argc) { std::fprintf(stderr, "Нужен аргумент после %s\n", what); std::exit(2); }
            return argv[++i];
        };
        if      (a == "--list")       opt.list = true;
        else if (a == "--set")        opt.sets.push_back(std::atoi(next("--set")));
        else if (a == "--sets") {
            int lo, hi; sscanf(next("--sets"), "%d-%d", &lo, &hi);
            for (int s = lo; s <= hi; ++s) opt.sets.push_back(s);
        }
        else if (a == "--substep")    opt.find_substep  = std::atoi(next("--substep"));
        else if (a == "--time")       opt.find_time     = std::atof(next("--time"));
        else if (a == "--component")  opt.component     = next("--component");
        else if (a == "--all-nodes")  opt.all_nodes     = true;
        else if (a == "--mises")      opt.mises         = true;
        else if (a == "--threads")    opt.threads       = std::atoi(next("--threads"));
        else if (a == "--out")        opt.out           = next("--out");
        else { std::fprintf(stderr, "Неизвестная опция: %s\n", a.c_str()); usage(); return 2; }
    }

    File f; std::string err;
    if (!f.open(opt.path, err)) {
        std::fprintf(stderr, "Ошибка открытия RST: %s\n", err.c_str());
        return 1;
    }

    if (opt.list || (!opt.has_sets() && !opt.find_substep && opt.find_time < 0)) {
        std::printf("Файл: %s\nНаборов (sets): %d, узлов: %lld\n",
                    f.path.c_str(), f.nsets, (long long)f.n_nodes);
        std::printf("%-6s %-6s %-8s %s\n", "idx", "step", "substep", "time");
        for (size_t i = 0; i < f.sets.size(); ++i)
            std::printf("%-6zu %-6d %-8d %.6g\n", i, f.sets[i].step,
                        f.sets[i].substep, f.sets[i].time);
        std::printf("\nName selections (%zu):\n", f.components.size());
        for (auto& c : f.components)
            std::printf("  %-32s %zu nodes\n", c.name.c_str(), c.nodes.size());
        if (opt.list) { f.close(); return 0; }
        std::printf("\nУкажите --set/--sets/--substep/--time.\n");
        f.close(); return 0;
    }

    // Разрешение выбора набора
    if (opt.sets.empty()) {
        int idx = opt.find_substep ? f.find_set_by_substep(opt.find_substep)
                                   : f.find_set_by_time(opt.find_time);
        if (idx < 0) { std::fprintf(stderr, "Набор не найден.\n"); f.close(); return 1; }
        opt.sets.push_back(idx);
    }

    int rc = 0;
    for (int s : opt.sets) {
        Options one = opt;
        one.sets = {s};
        rc |= run_parallel(f, one);
    }
    f.close();
    return rc;
}
