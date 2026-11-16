# MultiModelReader
Чтение и отображение нескольких моделей ANSYS без использования ANSYS

## Nodal Stress Reader
Быстрое чтение узловых напряжений из `.rst` с выбором по Name Selection.
Полная инструкция: [USAGE.md](USAGE.md).

Кратко:
```bash
python nodal_stress_reader.py file.rst              # интерактивный режим
python nodal_stress_reader.py file.rst -l           # список name selections
python nodal_stress_reader.py file.rst -c TOP_FACE --principal --mises
```  
