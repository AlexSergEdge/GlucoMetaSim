import csv
import numpy as np

import matplotlib.pyplot as plt

from model import simulate_enhanced, Simulation
from coefficients import *
from functions import *

from PyQt5 import QtWidgets, QtCore
import gui

# Константы
MAX_ALLOWED_GLUCOSE_MMOLL = 15  # Максимальный допустимый уровень глюкозы в начале тренировки
MIN_ALLOWED_GLUCOSE_MMOLL = 4  # Минимальный допустимый уровень глюкозы в начале тренировки
GLUCOSE_STEP_MMOLL = 1  # Шаг изменения глюкозы при создании таблицы

EXERCISE_DURATION_RANGE = range(15, 60 + 15, 15)  # Массив длительностей тренировок при создании таблицы
EXERCISE_HEART_RATE_RANGE = range(90, 150 + 30, 30)  # Массив ЧСС при создании таблицы
POST_TRAINING_TIME = 120  # Запас в 2 часа после тренировки (через 2 часа глюкоза должна быть в норме)
EXERCISE_START_TIME = 0  # Время начала упражнений
START_MEAL = 12000  # Начальное значение подбираемой пищи
START_INSULIN = 0.1  # Начальное значение подбираемого инсулина

# Перевод из ммоль/л в мг/дл
def to_mgdl(mmoll):
    return mmoll * MG_DL_TO_MMOL_L_CONVENTION_FACTOR

# Перевод ммоль/л в массу глюкозы
def to_gmass(mmoll):
    return mmoll * MG_DL_TO_MMOL_L_CONVENTION_FACTOR * VG

# Массив начальных значений глюкозы при создании таблицы
GLUCOSE_RANGE = np.arange(to_gmass(MIN_ALLOWED_GLUCOSE_MMOLL), to_gmass(MAX_ALLOWED_GLUCOSE_MMOLL + GLUCOSE_STEP_MMOLL), to_gmass(GLUCOSE_STEP_MMOLL))


# Перевод из массы глюкоы в ммоль/л
def to_mmoll(mgdl):
    return mgdl / (MG_DL_TO_MMOL_L_CONVENTION_FACTOR)

# Перевод инсулина из ед. в пмоль/л
def to_pmoll(u):
    return int(u / (kd + ka1) * 6.0)

# Перевод из пмоль/л в ед.
def to_u(pmoll):
    return pmoll * (kd + ka1) / 6.0

# Перевод хлебной единицы (12 г углеводов) в массу углеводов
def be_to_mg(be):
    return be * 12000

# Функция подбора пищи
def get_meal(simulation, Dcurr, Dprev, Gh, Gl):
    simulation.set_meal(Dcurr)
    res, t = simulate_enhanced(simulation)
    Gres = G(0, res[:,0])[-1]
    print(f'Gres = {to_mmoll(Gres)}')
    # Если G низкий
    if Gres < Gl:
        # Если до этого еды было меньше, увеличиваем на разницу текущей и предыдущей еды
        if Dcurr > Dprev:
            tmp = Dcurr
            Dcurr += (Dcurr - Dprev)
            Dprev = tmp
        # Если до этого еды было больше, увеличивем на половину разницы текущей и предыдущей еды
        else:
            tmp = Dcurr
            Dcurr += (Dprev - Dcurr) / 2
            Dprev = tmp
        # Рекурсивно вызываем функцию с обновленными параметрами
        return get_meal(simulation, Dcurr, Dprev, Gh, Gl)
    # Если G высокий
    if Gres > Gh:
        # Если до этого еды было меньше, уменьшаем пищу (делаем среднее между текущим и предыдущим)
        if Dcurr > Dprev:
            tmp = Dcurr
            Dcurr = (Dcurr + Dprev) / 2
            Dprev = tmp
        # Если до этого еды было больше, уменьшаем на половину разницы между текущим и предыдущим
        else:
            tmp = Dcurr
            Dcurr -= (Dprev - Dcurr) / 2
            Dprev = tmp
        # print(f'    Food decreased: {Dprev} -> {Dcurr}')
        # Рекурсивно вызываем функцию с обновленными параметрами
        return get_meal(simulation, Dcurr, Dprev, Gh, Gl)
    # Выходим, когда подобранная пища приводит к правильному уровню глюкозы
    return Dcurr


# Функция подбора пищи
def get_insulin(simulation, Icurr, Iprev, Gh, Gl):
    simulation.set_insulin(to_pmoll(Icurr))
    res, t = simulate_enhanced(simulation)
    Gres = G(0, res[:,0])[-1]
    print(f'Gres = {to_mmoll(Gres)}')
    # Если G высокий
    if Gres > Gh:
        # Если до этого инсулина было меньше, увеличиваем на разницу текущего и предыдущего инсулина
        if Icurr > Iprev:
            tmp = Icurr
            Icurr += (Icurr - Iprev)
            Iprev = tmp
        # Если до этого инсулина было больше, увеличивем на половину разницы текущего и предыдущего инсулина
        else:
            tmp = Icurr
            Icurr += (Iprev - Icurr) / 2
            Iprev = tmp
        # Рекурсивно вызываем функцию с обновленными параметрами
        return get_insulin(simulation, Icurr, Iprev, Gh, Gl)
    # Если G низкий
    if Gres < Gl:
        # Если до этого инсулина было меньше, уменьшаем инсулин (делаем среднее между текущим и предыдущим)
        if Icurr > Iprev:
            tmp = Icurr
            Icurr = (Icurr + Iprev) / 2
            Iprev = tmp
        # Если до этого инсулина было больше, уменьшаем на половину разницы между текущим и предыдущим
        else:
            tmp = Icurr
            Icurr -= (Iprev - Icurr) / 2
            Iprev = tmp
        # Рекурсивно вызываем функцию с обновленными параметрами
        return get_insulin(simulation, Icurr, Iprev, Gh, Gl)
    # Выходим, когда подобранный инсулин приводит к правильному уровню глюкозы    
    return Icurr


# Функция создает таблицу заранее рассчитанных тренировок
def precalculate_table(body_weight, heart_rate_basal, G_high, G_low):
    with open('table.csv', 'w', newline='') as f:
        csv_writer = csv.writer(f)
        csv_writer.writerow(['Gp0','Duration','HR', 'Meal', 'Insulin', 'G'])
        # проходимся по всем значениям
        for glucose in GLUCOSE_RANGE:
            for duration in EXERCISE_DURATION_RANGE:
                time = duration + POST_TRAINING_TIME
                ex_finish = EXERCISE_START_TIME + duration
                for hr in EXERCISE_HEART_RATE_RANGE:
                    meal, insulin, sim = run_selection(time, body_weight, glucose, EXERCISE_START_TIME, ex_finish, hr, heart_rate_basal, G_low, G_high)

                    sim_log = sim.get_sim_log()[-1]
                    g_col = sim_log[0] 
                    final_g = g_col[:,0]

                    print(f'Gp0={glucose}, dur={duration}, hr={hr}, food={meal}, insulin={insulin}, G={to_mmoll(G(0, final_g))[-1]}')
                    csv_writer.writerow([glucose, duration, hr, meal, insulin, to_mmoll(G(0, final_g))[-1]])


def get_recommendation_from_model(body_weight, heart_rate_basal, G_high, G_low, glucose, duration, hr):
    time = duration + POST_TRAINING_TIME
    ex_finish = EXERCISE_START_TIME + duration
    return run_selection(time, body_weight, glucose, EXERCISE_START_TIME, ex_finish, hr, heart_rate_basal, G_low, G_high)


def run_selection(time, bw, glucose, ex_start, ex_finish, hr, heart_rate_basal, G_low, G_high):
    simulation = Simulation(
        time=time, samples=time, BW=bw, D=0, Gpb=glucose, Gp0=glucose, Djins=0, IIRb=1.0, meal_time=0, injection_time=0,
        ex_on=True, ex_start=ex_start, ex_finish=ex_finish, ex_hr=hr, HRb=heart_rate_basal
    )
    res, t = simulate_enhanced(simulation)
    Gres = G(0, res[:,0])  # Концентрация глюкозы в плазме
    # Если глюкозв в целевом диапазоне, идем дальше
    if Gres[-1] >= G_low and Gres[-1] <= G_high:
        recommended_meal = 0
        recommended_insulin = 0
    # Если ниже, запускаем подбор пищи
    elif Gres[-1] < G_low:
        recommended_meal = get_meal(simulation, START_MEAL, 0, G_high, G_low)
        recommended_insulin = 0
    # Если выше, запускаем подбор инсулина
    elif Gres[-1] > G_high:
        recommended_meal = 0
        recommended_insulin = get_insulin(simulation, START_INSULIN, 0, G_high, G_low)
        
    return recommended_meal, recommended_insulin, simulation


# Получение рекомендации
def get_recommendations_from_table(glucose, duration, hr):
    # Если глюкозы слишком высокая или низкая - рекомендуем воздержаться от тренировки и дождаться нормализации
    if glucose >= to_gmass(MAX_ALLOWED_GLUCOSE_MMOLL):
        print('Too high glucose, training is not recommended, make sure it stays below {MAX_ALLOWED_GLUCOSE_MMOLL} mmol/l')
        return None
    if glucose <= to_gmass(MIN_ALLOWED_GLUCOSE_MMOLL):
        print('Too low glucose, training is not recommended, make sure is stays above {MIN_ALLOWED_GLUCOSE_MMOLL} mmol/l')
        return None
    # Читаем файл
    with open('table.csv', 'r', newline='') as f:
        csv_reader = csv.reader(f)
        next(csv_reader, None)  # skip the headers
        # create 2D array
        data = []
        for row in csv_reader:
            data.append(row)
    
    data = np.array(data)
    gvals = data[:, 0]
    prev_g = None
    upper_g = None
    glucose_to_check = None
    # Находим отрезок, в котором лежит значение глюкозы (либо при точном равенстве оставляем это значение)
    for g in gvals:
        if glucose == float(g):
            glucose_to_check = float(g)
            break
        if prev_g and float(g) > glucose:
            upper_g = float(g)
            glucose_to_check = (prev_g, upper_g)
            break
        prev_g = float(g)
    
    # Если значение глюкозы точно совпало со значением в таблице - сразу берем еду и инсулин
    if isinstance(glucose_to_check, float):
        for row in data:
            if float(row[0]) == glucose_to_check and int(row[1]) == duration and int(row[2]) == hr:
                recommended_meal = float(row[3])
                recommended_insulin = float(row[4])
                res_G = float(row[5])
    # Иначе интерполируем и получаем промежуточные значения инсулина/еды
    else:
        for row in data:
            if float(row[0]) == glucose_to_check[0] and int(row[1]) == duration and int(row[2]) == hr:
                meal1 = float(row[3])
                ins1 = float(row[4])
                res_G1 = float(row[5])
            if float(row[0]) == glucose_to_check[1] and int(row[1]) == duration and int(row[2]) == hr:
                meal2 = float(row[3])
                ins2 = float(row[4])
                res_G2 = float(row[5])
        
        recommended_meal = find_value_by_ratio(glucose_to_check[0], glucose, glucose_to_check[1], meal1, meal2)
        recommended_insulin = find_value_by_ratio(glucose_to_check[0], glucose, glucose_to_check[1], ins1, ins2)
        res_G = find_value_by_ratio(glucose_to_check[0], glucose, glucose_to_check[1], res_G1, res_G2)
    
    return recommended_meal, recommended_insulin, res_G

# https://brilliant.org/wiki/section-formula/           
def find_value_by_ratio(y1, y, y2, val1, val2):
    # y1, y2 - known (glucose levels from table)
    # y is known and divides [y1, y2]
    # function divides [val1,  val2] in same ratio as y divides [y1, y2]
    m = y - y1
    n = y2 - y
    return (m * val2 + n * val1) / (m + n)



def plot_selection_results(sim, Gh, Gl):
    
    sel_sim_log = sim.get_sim_log()[:-1]
    final_res = sim.get_sim_log()[-1]

    colors = plt.cm.rainbow(np.linspace(0,1,len(sim.get_sim_log())))

    results = []
    times = []
    labels = []
    linestyles = []
    name = 'Selection process'

    for simulation in sel_sim_log:
        x = simulation[0]
        t = simulation[1]
        D = simulation[2]
        Djins = to_u(simulation[3])
        bg = G(0, x[:,0])
        results.append(bg)
        times.append(t)
        labels.append(f'D={D}, Djins={Djins:.1f}')
        linestyles.append('dashed')

    
    x = final_res[0]
    t = final_res[1]
    D = final_res[2]
    Djins = to_u(final_res[3])
    bg = G(0, x[:,0])
    results.append(bg)
    times.append(t)
    labels.append(f'D={D}, Djins={Djins:.1f}')
    linestyles.append('solid')

    print_multiple_graphs(results, times, labels, colors, linestyles, name, (to_mmoll(Gl), to_mmoll(Gh)), (sim.ex_start, sim.ex_finish), is_mmoll=True)

        
def main():

    bw = 60
    hrb = 60
    Gh = to_mgdl(7.4)
    Gl = to_mgdl(4.8)


    glucose = to_gmass(5.6)
    glucose = 203.24529599999997
    glucose = to_gmass(12.6)
    print(f'Glucose mass {glucose}')
    duration = 30
    hr = 120

    # precalculate_table(
    #     body_weight=bw, heart_rate_basal=hrb, G_high=Gh, G_low=Gl
    # )
    
    meal, ins, res_G = get_recommendations_from_table(
        glucose, duration, hr
    )
    print(f'Gres = {res_G}\nMeal: {meal}, Insulin {ins}')

    meal, ins, sim = get_recommendation_from_model(
        bw, hrb, Gh, Gl, glucose, duration, hr
    )
    plot_selection_results(sim, Gh, Gl)
    print(f'Meal: {meal}, Insulin {ins}')


class CalculatorApp(QtWidgets.QMainWindow, gui.Ui_MainWindow):
    def __init__(self):
        super().__init__()
        self.setupUi(self)

if __name__ == "__main__":
    # main()
    app = QtWidgets.QApplication([])
    window = CalculatorApp()
    window.show()
    app.exec_()