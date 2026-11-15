#!/usr/bin/env python3
"""
nodal_stress_reader.py — быстрое чтение узловых напряжений из файлов
результатов ANSYS MAPDL (*.rst / *.rnn) по имени Name Selection (Component).

Программа НЕ требует установленного ANSYS: используется быстрый бинарный
ридер ansys-mapdl-reader (Cython + numpy). Имена name selections и списки
узлов читаются напрямую из заголовков RST-файла (без построения VTK-меша),
напряжения извлекаются нативным фильтром ридера `nodes=` — он читает с диска
только элементы, содержащие узлы выбранного компонента.

Источники name selections (пробуются по порядку):
  1. Файл результатов (.rst/.rnn)          — если модель сохранена с CM
     (MAPDL: `CMSLV,PICK` + `CM,...,N` до сохранения или /POST1 с CMSflg=1)
  2. Архив модели / CDB-файл (.cdb/.arc/.txt) — блок *CM,N,...
     Передаётся ключом --archive; если не задан, автоматически ищется
     одноимённый .cdb рядом с RST.

Примеры:
    python nodal_stress_reader.py file.rst                       # интерактивно
    python nodal_stress_reader.py file.rst --list                # список CM
    python nodal_stress_reader.py file.rst -a model.cdb -c FIXTURE_NODES -t 5 --mises
    python nodal_stress_reader.py file.rst -a model.cdb --select # выбрать из списка
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import time
import warnings

import numpy as np

try:
    from ansys.mapdl.reader.rst import Result
except ImportError:  # pragma: no cover
    sys.exit(
        "Не найден пакет 'ansys-mapdl-reader'. Установите его командой:\n"
        "    pip install ansys-mapdl-core    (или pip install ansys-mapdl-reader)"
    )

warnings.filterwarnings("ignore", category=DeprecationWarning)

STRESS_LABELS = ["SX", "SY", "SZ", "SXY", "SYZ", "SXZ"]


# ---------------------------------------------------------------------------
# Чтение name selections (components)
# ---------------------------------------------------------------------------

def read_components_from_archive(cdb_path: str) -> dict[str, np.ndarray]:
    """Разбирает блоки '*CM,NAME,N,nval,offset' из архива/CDB (MapdlBinaryReader).

    Быстро: только header/table1/table2, без триangulation/VTK.
    """
    from ansys.mapdl.reader.mesh import MapdlBinaryReader

    reader = MapdlBinaryReader(cdb_path)
    comps = {}
    for name, nodes in reader._node_comps.items():
        key = name.strip().upper()
        if key:
            comps[key] = np.asarray(nodes, dtype=np.int32)
    elem_comps = {k.strip().upper(): np.asarray(v, dtype=np.int32)
                  for k, v in reader._elem_comps.items()}
    read_components_from_archive.element_components = elem_comps
    return comps


def _load_tables(rst_path: str) -> tuple[bytes, bytes]:
    """Читает table1/table2 RST-файла напрямую (без Mesh/VTK)."""
    ext = os.path.splitext(rst_path)[1].lower()
    d = os.path.dirname(os.path.abspath(rst_path))
    stem = os.path.splitext(os.path.basename(rst_path))[0]
    if ext == ".rst":
        names = ["fileh", stem + "h", "file.h"]
    else:  # .rnn
        names = ["filenn", stem + "nn", "file.nn"]
    table1 = next((os.path.join(d, n) for n in names
                   if os.path.isfile(os.path.join(d, n))), None)
    if table1 is None:
        raise FileNotFoundError(
            f"Не найден файл таблиц '{names[0]}' рядом с {rst_path}")
    table2 = table1[:-1] + "2"
    with open(table1, "rb") as f:
        data1 = f.read()
    data2 = b""
    if os.path.isfile(table2):
        with open(table2, "rb") as f:
            data2 = f.read()
    return data1, data2


_CM_RE = re.compile(rb"\*CM,(.{0,90})")


def read_components_from_rst(rst_path: str) -> dict[str, np.ndarray]:
    """Name selections прямо из заголовков RST (таблица геометрии, record 1).

    Возвращает {} если модель была сохранена без компонентов (CMSflg=0).
    """
    try:
        data1, data2 = _load_tables(rst_path)
    except FileNotFoundError:
        return {}
    if b"*CM," not in data1:
        return {}

    node_comps: dict[str, list[int]] = {}
    elem_comps: dict[str, list[int]] = {}
    for match in _CM_RE.finditer(data1):
        raw = match.group(1).split(b",")
        if len(raw) < 4:
            continue
        name = raw[0].decode("utf-8", "ignore").strip().upper()
        ctype = raw[1].decode("utf-8", "ignore").strip().upper()
        try:
            nval, offset = int(raw[2]), int(raw[3])
        except ValueError:
            continue
        if not name or nval <= 0 or not data2:
            continue
        tag = f"<{offset}N<".encode()
        start = data2.find(tag)
        if start == -1:
            continue
        start += len(tag)
        end = data2.find(b">>", start)
        if end == -1:
            end = min(start + 4 * nval, len(data2))
        values = np.frombuffer(data2[start:end], dtype="<i4")
        target = node_comps if ctype == "N" else elem_comps if ctype == "E" else None
        if target is not None:
            target.setdefault(name, []).extend(values.tolist())

    read_components_from_rst.element_components = {
        k: np.array(v, dtype=np.int32) for k, v in elem_comps.items()}
    return {k: np.array(v, dtype=np.int32) for k, v in node_comps.items()}


def gather_components(rst_path: str, archive: str | None,
                      verbose: bool = True) -> tuple[dict[str, np.ndarray], str]:
    """Собирает name selections из RST и/или архива. Возвращает (comps, source)."""
    t0 = time.perf_counter()
    comps = read_components_from_rst(rst_path)
    src = "RST"
    # архив: явный ключ или автопоиск одноимённого .cdb
    cands = [archive] if archive else []
    if not cands:
        base = os.path.splitext(rst_path)[0]
        cands = [base + ".cdb", base + ".arc"]
    for c in cands:
        if c and os.path.isfile(c):
            arc_comps = read_components_from_archive(c)
            if arc_comps:
                for k, v in arc_comps.items():
                    comps.setdefault(k, v)
                src = os.path.basename(c)
                break
    if verbose:
        print(f"[тайминг] разбор name selections: {(time.perf_counter()-t0)*1e3:.0f} мс"
              f" (источник: {src}, найдено: {len(comps)})")
    return comps, src


# ---------------------------------------------------------------------------
# Открытие RST и извлечение напряжений (максимально быстро)
# ---------------------------------------------------------------------------

def open_result(rst_path: str) -> Result:
    """Открывает RST: mesh нужен для фильтрации узлов, VTK-разбор отключён."""
    t0 = time.perf_counter()
    res = Result(rst_path, read_mesh=True, parse_vtk=False)
    print(f"[тайминг] открытие RST (mesh без VTK): "
          f"{(time.perf_counter()-t0)*1e3:.0f} мс")
    return res


def extract_stresses(res: Result, rnum: int,
                     node_numbers: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Напряжения [n,6] только для узлов component (фильтр nodes= в Cython).

    Возвращает (nnum, stress). Сортировка узлов component → ускорение кэша.
    """
    sel = np.unique(node_numbers.astype(np.int32))
    t0 = time.perf_counter()
    nnum, stress = res.nodal_stress(rnum, nodes=sel)
    print(f"[тайминг] извлечение напряжений ({sel.size} узлов): "
          f"{(time.perf_counter()-t0)*1e3:.0f} мс")
    return nnum, stress


def principal_stresses(stress: np.ndarray) -> np.ndarray:
    """Главные напряжения S1>=S2>=S3 и von Mises (пакетный eigvalsh)."""
    sx, sy, sz, sxy, syz, sxz = stress.T
    t = np.zeros((stress.shape[0], 9), dtype=np.float64)
    t[:, 0], t[:, 4], t[:, 8] = sx, sy, sz
    t[:, 1] = t[:, 3] = sxy
    t[:, 5] = t[:, 7] = syz
    t[:, 2] = t[:, 6] = sxz
    eig = np.linalg.eigvalsh(t.reshape(-1, 3, 3))[:, ::-1]
    s1, s2, s3 = eig.T
    mises = np.sqrt(0.5 * ((s1 - s2) ** 2 + (s2 - s3) ** 2 + (s1 - s3) ** 2))
    return np.column_stack([eig, mises])


# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

def choose_record(res: Result, interactive: bool, rnum: int | None) -> int:
    total = res.nsets
    if rnum is not None:
        if not (0 <= rnum < total):
            sys.exit(f"Запись должна быть в диапазоне 0..{total-1} (всего {total}).")
        return rnum
    if not interactive or total <= 1:
        return total - 1
    ls_table = res._resultheader["ls_table"]
    times = res.time_values
    print("\nДоступные записи результатов:")
    for i in range(min(total, 30)):
        step, sub = int(ls_table[i, 0]), int(ls_table[i, 1])
        print(f"  [{i:>4}] step={step:>4} substep={sub:>4}  time={times[i]:.6g}")
    sel = input(f"Номер записи (0..{total-1}, по умолчанию {total-1}): ").strip()
    return int(sel) if sel else total - 1


def select_component(interactive: bool, chosen: str | None,
                     components: dict[str, np.ndarray]) -> str:
    names = sorted(components)
    if not names:
        sys.exit("Name selections не найдены. Сохраните компоненты в RST "
                 "(CMSLV,PICK + CM,...,N) или укажите архив модели ключом --archive.")
    if chosen is not None:
        key = chosen.strip().upper()
        if key not in components:
            sys.exit(f"Компонент '{chosen}' не найден. Доступные: {', '.join(names)}")
        return key
    if not interactive:
        sys.exit("Укажите name selection ключом --component/-c (неинтерактивный режим).")
    print("\nName selections (компоненты узлов):")
    for i, n in enumerate(names, start=1):
        print(f"  [{i:>3}] {n:<32} узлов: {components[n].size}")
    while True:
        sel = input("Выберите name selection (номер или имя): ").strip()
        if sel.isdigit() and 1 <= int(sel) <= len(names):
            return names[int(sel) - 1]
        if sel.upper() in components:
            return sel.upper()
        print("Неверный ввод, попробуйте ещё раз.")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="Быстрое чтение узловых напряжений ANSYS по Name Selection.")
    p.add_argument("rst", help="Файл результатов MAPDL (.rst / .rnn)")
    p.add_argument("-a", "--archive", metavar="CDB",
                   help="архив/CDB модели с блоками *CM (если в RST их нет)")
    p.add_argument("-l", "--list", action="store_true",
                   help="только вывести список name selections")
    p.add_argument("--select", action="store_true",
                   help="только выбрать name selection (интерактивно)")
    p.add_argument("-c", "--component", help="имя name selection (без диалога)")
    p.add_argument("-t", "--time", "--rnum", dest="rnum", type=int,
                   help="номер записи результатов (0-based)")
    p.add_argument("--principal", action="store_true",
                   help="добавить главные напряжения S1..S3")
    p.add_argument("--mises", action="store_true",
                   help="добавить интенсивные напряжения (von Mises)")
    p.add_argument("--csv", metavar="FILE", help="сохранить результат в CSV")
    p.add_argument("--head", type=int, default=20,
                   help="строк в консоли (0 = все)")
    p.add_argument("--no-input", action="store_true",
                   help="пакетный режим: не задавать вопросов")
    args = p.parse_args(argv)

    if not os.path.isfile(args.rst):
        sys.exit(f"Файл не найден: {args.rst}")
    interactive = sys.stdin.isatty() and not args.no_input

    # 1. Список name selections (быстрый парсер заголовков, без VTK)
    components, _src = gather_components(args.rst, args.archive)

    if args.list:
        if not components:
            sys.exit("Name selections не найдены (см. справку --help).")
        print("\nName selections:")
        for n in sorted(components):
            print(f"  {n:<32} узлов: {components[n].size}")
        return 0

    # 2. Выбор name selection
    chosen = select_component(interactive, args.component, components)
    if args.select:
        print(f"Выбрано: {chosen}")
        return 0

    # 3. Открытие RST (без меша — максимально быстро)
    res = open_result(args.rst)
    rnum = choose_record(res, interactive, args.rnum)

    # 4. Напряжения только для узлов выбранного name selection
    nnum, stress = extract_stresses(res, rnum, components[chosen])

    missing = int(np.setdiff1d(np.unique(components[chosen]), nnum).size)

    header = ["Node"] + STRESS_LABELS
    cols = [nnum[:, None].astype(np.float64), stress]
    if args.principal or args.mises:
        pr = principal_stresses(stress)
        if args.principal:
            header += ["S1", "S2", "S3"]
            cols.append(pr[:, :3])
        if args.mises:
            header += ["SEQV"]
            cols.append(pr[:, 3:4])
    table = np.hstack(cols)

    times = res.time_values
    ls = res._resultheader["ls_table"]
    print(f"\nФайл: {os.path.basename(args.rst)} | запись #{rnum} "
          f"(step={int(ls[rnum,0])}, substep={int(ls[rnum,1])}, time={times[rnum]:.6g})"
          f" | name selection: {chosen} | узлов с данными: {nnum.size}"
          + (f" | без данных: {missing}" if missing else ""))

    limit = args.head if args.head > 0 else max(table.shape[0], 1)
    widths = [max(len(h), 12) for h in header]
    line = " ".join(h.rjust(w) for h, w in zip(header, widths))
    print(line)
    print("-" * len(line))
    for row in table[:limit]:
        cells = [f"{int(row[0]):>{widths[0]}}"]
        cells += [f"{v:>{w}.5g}" for v, w in zip(row[1:], widths[1:])]
        print(" ".join(cells))
    if table.shape[0] > limit:
        print(f"... всего строк: {table.shape[0]}")

    if args.csv:
        np.savetxt(args.csv, table, delimiter=",",
                   header=",".join(header), comments="", fmt="%.8e")
        print(f"Сохранено: {args.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
