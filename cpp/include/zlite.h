// zlite.h — минимальный DEFLATE-декодер (raw/deflate stream), без внешних зависимостей.
// Используется для чтения сжатых записей ANSYS RST (компрессия ANSCOMP = zlib deflate).
#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>

namespace zl {

// Распаковать deflate-stream из src[0..src_len) в dst (dst очищается, затем дозаполняется).
// level_hint игнорируется (zlib-заголовок уже снят вызывающим кодом).
// Возвращает true при успехе.
bool inflate_raw(const uint8_t* src, size_t src_len, std::vector<uint8_t>& dst);

// Проверить zlib-обёртку и распаковать: src указывает на [CMF FLG][deflate...][adler32].
bool inflate_zlib(const uint8_t* src, size_t src_len, std::vector<uint8_t>& dst);

} // namespace zl
