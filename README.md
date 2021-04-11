### Анализ данных с устройств  

Данные лежат в папке `analytics/data`. Данные о ЧСС - файл `HR.csv`. Данные о уровнях глюкозы - файл `BG.txt`.  
Для их анализа нужно запустить:  
```
python .\analytics\analytics.py
```
В результате в `results/bg_real/` и `results/hr/` появятся графики ЧСС и глюкозы, а в `experiments/` будут лежать входные данные для моделирования и реальные значения глюкозы.  

### Запуск моделирования  

Для сравнения реальных данных с данными моделирования нужно запустить:
```
python .\model.py
```
В результате в `results/bg_compare/` появятся графики сравнения, в количественные оценки будут выведены в консоль и файл `results.xlsx`.