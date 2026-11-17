// reader.h — параллельное (по наборам результатов) чтение узловых напряжений
#pragma once

#include "rst.h"

#include <atomic>
#include <string>
#include <thread>
#include <vector>

namespace rst {

struct Options {
    std::string path;
    std::vector<int> sets;        // индексы наборов (0-based)
    int  find_substep = 0;        // >0 — найти набор по подшагу
    double find_time = -1.0;      // >=0 — найти набор по времени
    std::string component;        // имя name selection (пусто = все узлы)
    bool all_nodes = true;
    bool mises = false;
    bool list = false;
    int  threads = 0;             // 0 = auto (hardware_concurrency)
    std::string out;              // "" = stdout

    bool has_sets() const { return !sets.empty(); }
};

// Результат чтения одного набора
struct SetData {
    int set_index = -1;
    bool ok = false;
    std::string err;
    std::vector<int32_t> node_ids;                 // узлы в порядке вывода
    std::vector<double>  stress;                   // 6 компонентов на узел (SX..SXZ)
};

// Прочитать один набор: заголовок решения -> таблица DOF -> узловые данные.
// Потокобезопасно за счёт отдельного FILE* на поток (posix pread/fopen).
SetData read_set(File& meta, int set_index, const std::vector<int32_t>* node_filter);

// Параллельный прогон по всем наборам из opt; запись CSV в opt.out.
int run_parallel(File& meta, Options opt);

} // namespace rst
