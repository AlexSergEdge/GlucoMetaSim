import csv
import numpy as np

from model import simulate_enhanced, Simulation
from coefficients import *
from functions import *


# Перевод из ммоль/л в мг/дл
def to_mgdl(mmoll):
    return mmoll * MG_DL_TO_MMOL_L_CONVENTION_FACTOR

# Перевод ммоль/л в массу глюкозы
def to_gmass(mmoll):
    return mmoll * MG_DL_TO_MMOL_L_CONVENTION_FACTOR * VG

# Перевод инсулина из ед. в пмоль/л
def to_pmoll(u):
    return int(u / (kd + ka1) * 6.0)

# Перевод хлебной единицы (12 г углеводов) в массу углеводов
def be_to_mg(be):
    return be * 12000

# Функция создает таблицу заранее рассчитанных тренировок
# На вход получает (параметры пользователя):
# - вес тела
# - базальный ЧСС
# - верхняя граница глюкозы 
# - нижняя граница глюкозы
def precalculate_table(body_weight, heart_rate_basal, G_high, G_low):
    
    # Фиксированные константные параметры
    posttraining_time = 120  # время, добавляемое после тренировки
    ex_on = True  # модель физических нагрузок подключена
    ex_start = 0  # начало тренировки в 0 минут
    IIRb = 1.0  # базальный инсулин
    meal_time = 0  # прием пищи в 0 минут
    injection_time = 0  # введение инсулина в 0 минут 
    
    # Кортежи значений параметров, которые меняются от тренировки к тренировке:
    glucose_range = np.arange(to_gmass(3), to_gmass(15), to_gmass(1))  # начальный уровень глюкозы в плазме
    ex_duration_range = range(15, 60 + 5, 5)  # длительность тренировки
    ex_hr_range = range(90, 180 + 30, 30)  # средний пульс в ходе тренировки
    
    with open('table.csv', 'w', newline='') as f:
        csv_writer = csv.writer(f)
        csv_writer.writerow(['Gp0','Duration','HR', 'Meal', 'Insulin'])

        # проходимся по всем значениям
        for glucose in glucose_range:
            for duration in ex_duration_range:
                time = duration + posttraining_time
                ex_finish = ex_start + duration
                for hr in ex_hr_range:
                    # ВАЖНО: еда и инсулин сначала на 0
                    # print(f't={time}, dur={duration}, Gp0={glucose}, hr={hr}')
                    simulation = Simulation(
                        time=time, 
                        samples=time, 
                        BW=body_weight, 
                        D=0, 
                        Gpb=glucose, 
                        Gp0=glucose,
                        Djins=0,
                        IIRb=IIRb,
                        meal_time=meal_time,
                        injection_time=injection_time,
                        ex_on=ex_on,
                        ex_start=ex_start,
                        ex_finish=ex_finish,
                        ex_hr=hr,
                        HRb=heart_rate_basal
                    )

                    res, t = simulate_enhanced(simulation)
                    Gres = G(0, res[:,0])  # Концентрация глюкозы в плазме
                    # Если глюкозв в целевом диапазоне, идем дальше
                    if Gres[-1] >= G_low and Gres[-1] <= G_high:
                        # print(f'Good: {Gres[-1]}')
                        recommended_food = 0
                        recommended_insulin = 0
                    # Если ниже, запускаем подбор пищи
                    elif Gres[-1] < G_low:
                        recommended_food = get_food(simulation, 6000, 0, G_high, G_low)
                        recommended_insulin = 0
                    # Если выше, запускаем подбор инсулина
                    elif Gres[-1] > G_high:
                        recommended_food = 0
                        recommended_insulin = get_insulin(simulation, 0.1, 0, G_high, G_low)

                    print(f'Gp0={glucose}, dur={duration}, hr={hr}, food={recommended_food}, insulin={recommended_insulin}')
                    csv_writer.writerow([glucose, duration, hr, recommended_food, recommended_insulin])
                

# Функция подбора пищи
def get_food(simulation, Dcurr, Dprev, Gh, Gl):
    
    simulation.set_meal(Dcurr)
    res, t = simulate_enhanced(simulation)
    Gres = G(0, res[:,0])[-1]

    # Если G низкий
    if Gres < Gl:
        # print(f'    LOW: {Gres} < {Gl}')
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
        # print(f'    Food increased: {Dprev} -> {Dcurr}')
        # Рекурсивно вызываем функцию с обновленными параметрами
        return get_food(simulation, Dcurr, Dprev, Gh, Gl)
    # Если G высокий
    if Gres > Gh:
        # print(f'    HIGH: {Gres} > {Gh}')
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
        return get_food(simulation, Dcurr, Dprev, Gh, Gl)
    # Выходим, когда подобранная пища приводит к правильному уровню глюкозы
    # print(f'Good: {Gres} with food {Dcurr}')
    return Dcurr


# Функция подбора пищи
def get_insulin(simulation, Icurr, Iprev, Gh, Gl):
    
    simulation.set_insulin(to_pmoll(Icurr))
    res, t = simulate_enhanced(simulation)
    Gres = G(0, res[:,0])[-1]
    
    # Если G высокий
    if Gres > Gh:
        # print(f'    HIGH: {Gres} > {Gh}')
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
        # print(f'    Insulin increased: {Iprev} -> {Icurr}')
         

        # Рекурсивно вызываем функцию с обновленными параметрами
        return get_insulin(simulation, Icurr, Iprev, Gh, Gl)
    # Если G низкий
    if Gres < Gl:
        # print(f'    LOW: {Gres} < {Gl}')
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
        # print(f'    Insulin decreased: {Iprev} -> {Icurr}')
    
        # Рекурсивно вызываем функцию с обновленными параметрами
        return get_insulin(simulation, Icurr, Iprev, Gh, Gl)
    # Выходим, когда подобранный инсулин приводит к правильному уровню глюкозы
    # print(f'Good: {Gres} with insulin {Icurr}')
    return Icurr



def main():
    # precalculate_table(
    #     body_weight=62,  # вес тела
    #     heart_rate_basal=70,  # ЧСС в состоянии покоя
    #     G_high=to_mgdl(7.4),
    #     G_low=to_mgdl(4.8)
    # )


if __name__ == "__main__":
    main()