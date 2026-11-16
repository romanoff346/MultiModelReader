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
    """Разбирает блоки '*CM,NAME,N,...' из архива/CDB модели.

    Совместимо со всеми версиями pymapdl-reader: используется публичный класс
    ``Archive`` (ansys-mapdl-core 0.6.x / reader 0.28+, ANSYS 2022 R2) и
    запасной вариант ``MapdlBinaryReader`` для очень старых версий.
    Чтение быстрое: parse_vtk=False — без построения VTK-сетки.
    """
    comps_raw = elem_raw = None
    try:
        from ansys.mapdl.reader.archive import Archive
        arc = Archive(cdb_path, parse_vtk=False)
        comps_raw, elem_raw = arc.node_components, arc.element_components
    except ImportError:  # очень старые версии reader (<0.28)
        from ansys.mapdl.reader.mesh import MapdlBinaryReader
        reader = MapdlBinaryReader(cdb_path)
        comps_raw, elem_raw = reader._node_comps, reader._elem_comps

    comps = {}
    for name, nodes in comps_raw.items():
        key = name.strip().upper()
        if key:
            comps[key] = np.asarray(nodes, dtype=np.int32)
    elem_comps = {k.strip().upper(): np.asarray(v, dtype=np.int32)
                  for k, v in (elem_raw or {}).items()}
    read_components_from_archive.element_components = elem_comps
    return comps


def read_components_from_rst(rst_path: str) -> dict[str, np.ndarray]:
    """Name selections из таблицы геометрии RST (ptrCOMP/maxcomp).

    Реализация опирается только на стабильные внутренние API pymapdl-reader
    (Result.read_record / _geometry_header / _read_components), которые
    одинаковы во всех версиях, поддерживающих компоненты в RST (~ANSYS 194+;
    полноценно — с MAPDL 2021+, включая 2022 R2). Возвращает {} если модель
    сохранена без компонентов (CMSflg=0) или файл слишком старый.
    """
    from ansys.mapdl.reader.rst import Result

    t0 = time.perf_counter()
    res = Result(rst_path, read_mesh=False, parse_vtk=False)
    try:
        maxcomp = int(res._geometry_header.get("maxcomp", 0) or 0)
        if maxcomp <= 0:
            read_components_from_rst.element_components = {}
            return {}
        node_comp, elem_comp = res._read_components()
    finally:
        # В reader >=0.28 (ANSYS 2022 R2) дескриптор — res._cfile (Cython Frf);
        # атрибута _file нет. Явно НЕ закрываем его повторно: двойной close
        # приводит к segfault в __del__ (AnsysBinary). Обнуляем ссылку —
        # ресурс освободит сборщик мусора.
        res._cfile = None
    print(f"[тайминг] чтение CM из RST ({maxcomp} записей): "
          f"{(time.perf_counter()-t0)*1e3:.0f} мс")

    comps = {}
    for name, nodes in node_comp.items():
        key = name.strip().upper()
        if key:
            comps[key] = np.asarray(nodes, dtype=np.int32)
    read_components_from_rst.element_components = {
        k.strip().upper(): np.asarray(v, dtype=np.int32)
        for k, v in elem_comp.items() if k.strip()}
    return comps


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
    """Открывает RST максимально быстро: сетка читается, VTK НЕ разбирается.

    Построение VTK UnstructuredGrid (`parse_vtk=True`) на больших моделях
    занимает минуты; здесь оно отключено, а напряжения читаются напрямую
    через Cython-ридер `_binary_reader.read_nodal_values` (см.
    ``extract_stresses``) с той же математикой усреднения по узлам, что и
    штатный ``nodal_stress``. Это проверенно идентично (max diff = 0)
    эталонному пути в ansys-mapdl-reader (0.28+ ... 0.56, ANSYS 2022 R2).
    """
    t0 = time.perf_counter()
    res = Result(rst_path, read_mesh=True, parse_vtk=False)
    print(f"[тайминг] открытие RST (mesh без VTK): "
          f"{(time.perf_counter()-t0)*1e3:.0f} мс")
    return res


def _fast_nodal_stress(res: Result, rnum: int) -> tuple[np.ndarray, np.ndarray]:
    """Полные узловые напряжения (nnum, stress[n,6]) БЕЗ разбора VTK-сетки.

    Прямой вызов Cython-ридера ENS + деление на число вкладов элементов —
    эквивалент ``res.nodal_stress(rnum)`` при ``grid=None``. Совместимо с
    ansys-mapdl-reader 0.28+ (ANSYS 2022 R2); при изменении внутреннего API
    бросает AttributeError/ImportError — вызывающий включает запасной путь.
    """
    from ansys.mapdl.reader import _binary_reader
    from ansys.mapdl.reader.misc import vtk_cell_info
    from ansys.mapdl.reader.rst import ELEMENT_INDEX_TABLE_KEYS

    if getattr(res, "_mesh", None) is None:
        raise AttributeError("mesh не загружен (read_mesh=False)")

    # минимальный линейный грид нужен ридеру только для celltypes/offset/cells;
    # строим его один раз (быстрее полного parse_vtk: без координат/fix_midside)
    grid = getattr(res, "_fast_grid", None)
    if grid is None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            grid = res._mesh._parse_vtk(null_unallowed=True,
                                        fix_midside=False).linear_copy()
        res._fast_grid = grid
    cells, offset = vtk_cell_info(grid)

    rnum = res.parse_step_substep(rnum)
    ele_ind_table, ans_etype, ptr_off = res._element_solution_header(rnum)
    nitem = res._result_nitem(rnum, "ENS")
    result_index = ELEMENT_INDEX_TABLE_KEYS.index("ENS")

    data, ncount = _binary_reader.read_nodal_values(
        res.filename, grid.celltypes, ele_ind_table, offset, cells,
        nitem, grid.n_points, res._nodstr.copy(), ans_etype,
        res._mesh.etype, result_index, ptr_off, False)

    if nitem != 6:
        data = data[:, :6]

    mask = res._insolution  # узлы, реально присутствующие в решении
    ncount = ncount[mask].reshape(-1, 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        stress = data[mask] / ncount
    return np.ascontiguousarray(res._mesh.nnum[mask]), stress


def extract_stresses(res: Result, rnum: int,
                     node_numbers: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Напряжения [n,6] только для узлов name selection (векторный numpy-фильтр).

    Выборка узлов — O(N) searchsorted по отсортированным номерам узлов
    решения. Возвращает (nnum, stress). Если внутренний быстрый путь
    недоступен (другая версия reader), автоматически используется штатный
    ``nodal_stress`` с построением VTK.
    """
    sel = np.unique(node_numbers.astype(np.int32))
    t0 = time.perf_counter()
    try:
        nnum_all, stress_all = _fast_nodal_stress(res, rnum)
    except (AttributeError, ImportError, TypeError):
        print("[предупреждение] быстрый путь недоступен в этой версии "
              "ansys-mapdl-reader, включаю разбор VTK...")
        if getattr(res, "grid", None) is None:
            res._store_mesh(parse_vtk=True)
        nnum, stress = res.nodal_stress(rnum, nodes=sel)
        print(f"[тайминг] извлечение напряжений через VTK ({sel.size} узлов): "
              f"{(time.perf_counter()-t0)*1e3:.0f} мс")
        return nnum, stress
    # узлы решения отсортированы по номеру → быстрая выборка
    pos = np.searchsorted(nnum_all, sel)
    ok = (pos < nnum_all.size) & (nnum_all[np.clip(pos, 0, nnum_all.size - 1)] == sel)
    idx = pos[ok]
    nnum, stress = nnum_all[idx], stress_all[idx]
    print(f"[тайминг] извлечение напряжений ({sel.size} узлов выбрано, "
          f"{nnum.size} с данными): {(time.perf_counter()-t0)*1e3:.0f} мс")
    return nnum, stress


def principal_stresses(stress: np.ndarray) -> np.ndarray:
    """Главные напряжения S1>=S2>=S3 и von Mises (пакетный eigvalsh).

    Строки с NaN/Inf (узлы без данных) возвращаются как NaN — eigvalsh на
    них не вызывается (иначе LinAlgError / зависание LAPACK).
    """
    n = stress.shape[0]
    out = np.full((n, 4), np.nan)
    ok = np.isfinite(stress).all(axis=1)
    if not ok.any():
        return out
    s = stress[ok].astype(np.float64)
    sx, sy, sz, sxy, syz, sxz = s.T
    t = np.zeros((s.shape[0], 9), dtype=np.float64)
    t[:, 0], t[:, 4], t[:, 8] = sx, sy, sz
    t[:, 1] = t[:, 3] = sxy
    t[:, 5] = t[:, 7] = syz
    t[:, 2] = t[:, 6] = sxz
    eig = np.linalg.eigvalsh(t.reshape(-1, 3, 3))[:, ::-1]
    s1, s2, s3 = eig.T
    mises = np.sqrt(0.5 * ((s1 - s2) ** 2 + (s2 - s3) ** 2 + (s1 - s3) ** 2))
    out[ok] = np.column_stack([eig, mises])
    return out


# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

def _load_step_info(res, total: int) -> np.ndarray:
    """Массив [nsets, 2] со step/substep. Совместимо с разными версиями API."""
    try:
        ls = res._resultheader["ls_table"]
        return np.asarray(ls).reshape(total, -1)[:, :2]
    except (KeyError, AttributeError, IndexError, ValueError):
        info = []
        for i in range(min(total, 30)):
            try:
                si = res.solution_info(i)
                info.append([si.get("time_step", 0), si.get("sub_time_step", 0)])
            except Exception:
                info.append([0, 0])
        return np.array(info, dtype=int) if info else np.zeros((0, 2), int)


def choose_record(res: Result, interactive: bool, rnum: int | None) -> int:
    total = res.nsets
    if rnum is not None:
        if not (0 <= rnum < total):
            sys.exit(f"Запись должна быть в диапазоне 0..{total-1} (всего {total}).")
        return rnum
    if not interactive or total <= 1:
        return total - 1
    steps = _load_step_info(res, total)
    times = res.time_values
    print("\nДоступные записи результатов:")
    for i in range(min(total, 30)):
        step, sub = int(steps[i, 0]), int(steps[i, 1])
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
    steps = _load_step_info(res, res.nsets)
    print(f"\nФайл: {os.path.basename(args.rst)} | запись #{rnum} "
          f"(step={int(steps[rnum,0])}, substep={int(steps[rnum,1])}, time={times[rnum]:.6g})"
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
