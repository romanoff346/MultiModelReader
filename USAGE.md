# Инструкция по использованию `nodal_stress_reader.py`

Чтение узловых напряжений из файлов результатов ANSYS MAPDL (`.rst`)
с выбором по Name Selection (*CM). Оптимизировано для максимальной скорости.
Проверено с ANSYS 2022 R2 (v222), `ansys-mapdl-reader >= 0.52.4`.

## Требования и установка

```bash
pip install ansys-mapdl-reader numpy
```

Python 3.8+. Файл результатов должен быть версии v222 или новее
(для старых версий пересоберите RST в MAPDL 2022 R2).

## Интерактивный режим (по умолчанию)

```bash
python nodal_stress_reader.py file.rst
```

Программа:
1. выводит список Name Selection из RST;
2. предлагает выбрать компоненту;
3. показывает доступные записи результатов;
4. выводит напряжения узлов выбранной компоненты.

## Основные сценарии CLI

| Задача | Команда |
|---|---|
| Только список name selections | `python nodal_stress_reader.py file.rst -l` |
| Только выбор компоненты | `python nodal_stress_reader.py file.rst --select` |
| Напряжения компоненты, запись №2 | `python nodal_stress_reader.py file.rst -c TOP_FACE -t 2` |
| Главные напряжения S1..S3 + von Mises | `python nodal_stress_reader.py file.rst -c TOP_FACE --principal --mises` |
| Сохранить в CSV | `python nodal_stress_reader.py file.rst -c TOP_FACE --csv out.csv` |
| Все строки в консоли | `... --head 0` (по умолчанию показывается 20) |
| Без вопросов (скрипты / CI) | `... --no-input -c TOP_FACE -t 0` |
| Name selection из архива/CDB модели | `python nodal_stress_reader.py file.rst -a file.cdb -l` |

## Пример вывода

```
Файл: file.rst | запись #2 (step=1, substep=2, time=0.5) | name selection: TOP_FACE | узлов с данными: 145

        Node           SX           SY           SZ           SXY          SYZ          SXZ
------------------------------------------------------------------------------------------
          12      1.2345e+07   -2.1000e+05    ...
```

Колонки напряжений: SX, SY, SZ, SXY, SYZ, SXZ (компоненты тензора);
при `--principal` добавляются S1, S2, S3; при `--mises` — SEQV (von Mises).

## Использование как библиотеки

```python
from nodal_stress_reader import gather_components, open_result, extract_stresses

comps, _ = gather_components("file.rst")          # {имя: массив номеров узлов}
res = open_result("file.rst")
nnum, stress = extract_stresses(res, rnum=0, node_array=comps["TOP_FACE"])
# stress: Nx6 — [SX, SY, SZ, SXY, SYZ, SXZ]
```

## Как обеспечена высокая скорость

- Список компонент читается парсером заголовков RST **без загрузки меша/VTK**.
- Напряжения читаются напрямую из блоков ENS через `_fast_nodal_stress()`
  (внутренний Cython-ридер) с фильтрацией узлов через `np.searchsorted` —
  без полного разбора VTK. Результат побитово совпадает с эталоном.
- Безопасное освобождение файла (`res._cfile = None`) вместо `close()`,
  который вызывал segfault в `__del__` некоторых версий reader.
- Типичное время отклика: единицы миллисекунд – десятки мс на запись.

## Ограничения

- Если в RST нет блоков *CM, укажите архив/CDB модели ключом `-a`.
- Главные напряжения устойчивы к NaN-строкам (пропущенным результатам узлов).
