import numpy as np

from model import simulate_enhanced, Simulation
from coefficients import *
from functions import *

def to_mgdl(mmoll):
    return mmoll * MG_DL_TO_MMOL_L_CONVENTION_FACTOR

def to_gmass(mmoll):
    return mmoll * MG_DL_TO_MMOL_L_CONVENTION_FACTOR * VG

def to_pmoll(u):
    return int(u / (kd + ka1) * 6.0)

def be_to_mg(be):
    return be * 12000

# Функция создает таблицу заранее рассчитанных тренировок
def precalculate_table(body_weight, heart_rate_basal, G_high, G_low):
    
    # Фиксированные константные параметры
    posttraining_time = 120  # время, добавляемое после тренировки
    ex_on = True  # модель физических нагрузок подключена
    ex_start = 0  # начало тренировки в 0 минут
    IIRb = 1.0  # базальный инсулин
    meal_time = 0  # прием пищи в 0 минут
    injection_time = 0  # введение инсулина в 0 минут 
    
    # Кортежи значений параметров, которые меняются от тренировки к тренировке:
    glucose_range = np.arange(to_gmass(3), to_gmass(15), to_gmass(1))

    ex_duration_range = range(15, 60, 5)  # длительность тренировки
    ex_hr_range = range(90, 180, 30)  # средний пульс в ходе тренировки
    
    # Кортежи значений параметров, которые будут подбираться
    # ins_range = range(0, to_pmoll(1), to_pmoll(0.1))  # введенный инсулин
    ins_range = np.arange(0, to_pmoll(1), to_pmoll(0.1)) 
    digestion_range = range(0, be_to_mg(4), be_to_mg(1))  # принятая пища (be - хлебная единица, 1 be = 12 г)


    for glucose in glucose_range:
        for duration in ex_duration_range:
            time = duration + posttraining_time
            ex_finish = ex_start + duration
            for hr in ex_hr_range:
                # WARNING: WITH 0 FOOD AND INSULIN
                print(f't={time}, dur={duration}, Gp0={glucose}, hr={hr}')
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
                
                if Gres[-1] > G_low and Gres[-1] < G_high:
                    print(f'Good: {Gres[-1]}')
                    continue
                if Gres[-1] < G_low:
                    get_food(simulation, Gres[-1], 6000, 0, G_high, G_low)
                elif Gres[-1] > G_high:
                    print(f'HIGH: {Gres[-1]}')
                

def get_food(simulation, Gres, Dcurr, Dprev, Gh, Gl):
    
    if Gres < Gl:
        print(f'    LOW: {Gres} < {Gl}')
        
        if Dcurr > Dprev:
            tmp = Dcurr
            Dcurr += (Dcurr - Dprev)
            Dprev = tmp
        else:
            tmp = Dcurr
            Dcurr += (Dprev - Dcurr) / 2
            Dprev = tmp

        print(f'    Food increased: {Dprev} -> {Dcurr}')

        simulation.set_meal(Dcurr)
        res, t = simulate_enhanced(simulation)
        Gres = G(0, res[:,0]) 
    
        return get_food(simulation, Gres[-1], Dcurr, Dprev, Gh, Gl)

    if Gres > Gh:
        print(f'    HIGH: {Gres} > {Gh}')
        
        if Dcurr > Dprev:
            tmp = Dcurr
            Dcurr = (Dcurr + Dprev) / 2
            Dprev = tmp
        else:
            tmp = Dcurr
            Dcurr -= (Dprev - Dcurr) / 2
            Dprev = tmp

        print(f'    Food decreased: {Dprev} -> {Dcurr}')

        simulation.set_meal(Dcurr)
        res, t = simulate_enhanced(simulation)
        Gres = G(0, res[:,0])
    
        return get_food(simulation, Gres[-1], Dcurr, Dprev, Gh, Gl)

    print(f'Good: {Gres} with food {Dcurr}')
    return Dcurr


def main():
    precalculate_table(
        body_weight=62,  # вес тела
        heart_rate_basal=70,  # ЧСС в состоянии покоя
        G_high=to_mgdl(7.4),
        G_low=to_mgdl(4.8)
    )


if __name__ == "__main__":
    main()