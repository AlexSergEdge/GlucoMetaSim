#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
from math import tanh

import numpy as np
from scipy.integrate import odeint

from utils import solve_quadratic_equasion
from coefficients import *
from functions import *


def simulate( 
    time,  # длительность моделирования, мин 
    samples,  # число временных отсчетов
    BW,  # вес тела
    D,  # масса принятых углеводов, мг
    Gpb,  # базальный уровень глюкозы в крови
    Gp0,  # начальный уровень глюкозы в крови
    Djins,  # доза быстрого инсулина
    IIRb,  # базальный уровень инсулина
    meal_time,  # время приема пищи, мин с начала моделирования
    injection_time,  # время введения нисулина, мин с начала моделирования
    ex_on,  # включение/выключение физической активности
    ex_start,  # время начала упражнения
    ex_finish,  # время завершения упражнения
    ex_hr,  # ЧСС во время упражнения
    HRb  # базальное значение ЧСС
):
    '''Разные параметры'''
    # Параметры для ЖКТ модели, добавлено s чтобы не путать это с alpha и beta для инсулина
    if D != 0:
        alphas = (5.0 / 2.0) * (1 / (D * (1 - b)))
        betas = (5.0 / 2.0) * (1 / (D * d))
    else:
        alphas = 0
        betas = 0
    m3 = HEb * m1 / (1 - HEb)

    '''Базальные значения глюкозы'''
    Gb = G(0, Gpb)  # базальная концентрация инсулина
    # Рассчитаем базальный уровень глюкозы Gtb (в тканях с медленным всасыванием глюкозы)    
    '''
    Есть формула 20: 
      Gtb = (Fcns - EGPb + k1 * Gp) / k2 (20)
    EGP есть в формуле 21:
      EGPb = Fcns + (Vm0 * Gtb) / (Km0 + Gtb) (21)
    Поэтому решаем квадратное уравнение:
      k2 * Gtb ^ 2 + (k2*Km0 - k1*Gpb + Vm0) * Gtb - k1 * Km0 * Gpb = 0
    '''
    Gtb = solve_quadratic_equasion(
        a_eq=k2,
        b_eq=k2*Km0 - k1*Gpb + Vm0,
        c_eq=-k1 * Km0 * Gpb)

    '''Параметры инсулина'''
    # Начальное значение подкожного инсулина:
    if injection_time == 0:  # ВАЖНО: При введении в момент времени 0
        Isc1ss = Djins + IIRb / (kd + ka1)
    else:  # ВАЖНО: При введении в другое время
        Isc1ss = IIRb / (kd + ka1)
    Isc2ss = kd / ka2 * Isc1ss
    Ipb = IIRb / (m2 + m4 - (m1 * m2)/(m1 + m3))  # базальный инсулин в плазме здорового человека
    Ilb = Ipb * m2/(m1 + m3)
    Ib = I(0, Ipb) # базальная кнцентрация инсулина в плазме

    '''Разные параметры'''
    EGPb = Fcns + (Vm0 * Gtb) / (Km0 + Gtb) # базальный EGP (21)
    global kp1
    kp1 = EGPb + kp2 * Gpb + kp3 * Ib   # глобальный параметр переопределяется в модели подкожного инсулина

    '''Начальные значения'''
    Gt0 = Gtb
    Il0 = Ilb
    Ip0 = Ipb
    Isc10 = Isc1ss
    Isc20 = Isc2ss
    Ione0 = Ib
    Id0 = Ib
    if meal_time == 0:  # ВАЖНО: При приеме пищи в момент времени 0
        Qsto10 = D
    else:  # ВАЖНО: При приеме в другое время
        Qsto10 = 0
    Qsto20 = 0
    Qgut0 = 0
    X0 = 0
    Y0 = 0
    Z0 = 0

    '''Система дифф. уравнений'''
    def ode_system(x, t):
        '''Входные параметры'''
        Gp = x[0]
        Gt = x[1]
        Il = x[2]
        Ip = x[3]
        Isc1 = x[4]
        Isc2 = x[5]
        Ione = x[6]
        Id = x[7]
        Qsto1 = x[8]
        Qsto2 = x[9]
        Qgut = x[10]
        X = x[11]
        Y = x[12]
        Z = x[13]
        '''Подсистема глюкозы'''
        dGpdt = EGP(t, Gp, Id) + Ra(t, Qgut, BW) - Uii(t) - E(t, Gp) - k1 * Gp + k2 * Gt
        if not ex_on:  # Если упражнения не включены - используем Uid а не Uid_ex
            dGtdt = -Uid(t, X, Gt) + k1 * Gp - k2 * Gt
        else:
            dGtdt = -Uid_ex(t, X, Y, Z, Gt, Ib, HRb, ex_start, ex_finish, ex_hr) + k1 * Gp - k2 * Gt
        '''Подсистема инсулина'''
        dIldt = -(m1 + m3) * Il + m2 * Ip
        dIpdt = -(m2 + m4) * Ip + m1 * Il + R(t, Isc1, Isc2)  # У здорового нет R
        ''' Подсистема подкожного инсулина'''
        # ВАЖНО: Если введение инсулина происходит в момент времени 0, то Isc10 = Djins
        dIsc1dt = -(kd + ka1) * Isc1 + IIR(t, IIRb)
        dIsc2dt = kd * Isc1 - ka2 * Isc2
        '''Эндогенная продукция глюкозы'''
        dIonedt = -ki * (Ione - I(t, Ip))
        dIddt = -ki * (Id - Ione)
        '''Глюкоза в пищеварительном тракте'''
        # ВАЖНО: Если прием пищи мгновенный (по умолчанию), ingestion(t) всегда 0, а Qsto10 = D
        dQsto1dt = -kgri *  Qsto1 + ingestion(t)
        dQsto2dt = -kempt(t, Qsto1, Qsto2, D, meal_time, alphas, betas) * Qsto2 + kgri * Qsto1
        dQgutdt = -kabs * Qgut + kempt(t, Qsto1, Qsto2, D, meal_time, alphas, betas) * Qsto2        
        '''Инсулин в межклеточной жидкости'''
        dXdt = -p2U * X + p2U * (I(t, Ip) - Ib)  # was deleted - Ib
        '''Физическая активность'''
        dYdt = -(1/Thr) * (Y - (HR(t, ex_start, ex_finish, ex_hr, HRb) - HRb))
        dZdt = -(fex(t, Y, HRb)/ Tin + 1 / Tex) * Z + fex(t, Y, HRb)

        return [dGpdt, dGtdt, dIldt, dIpdt, dIsc1dt, dIsc2dt, dIonedt, dIddt, dQsto1dt, dQsto2dt, dQgutdt, dXdt, dYdt, dZdt]
    
    if meal_time > 0 and meal_time == injection_time:
        t = np.linspace(0, time, samples)  # начало, завершение моделирования, число отсчетов
        x10 = [Gp0, Gt0, Il0, Ip0, Isc10, Isc20, Ione0, Id0, Qsto10, Qsto20, Qgut0, X0, Y0, Z0]  # начальные условия
        t1 = np.linspace(0, meal_time, meal_time)
        x1 = odeint(ode_system, x10, t1)
        # Когда происходит прием пищи и ввод инсулина - "разрываем" моделирование
        x20 = get_updated_init_conditions(x1, D, Djins)
        t2 = np.linspace(meal_time, time + 1, time - meal_time + 1)
        x2 = odeint(ode_system, x20, t2, hmax=1)
        
        x2 = np.delete(x2, (0), axis=0)  # delete first row

        for row in (x2):
            x1 = np.append(x1, [row], axis=0)
        x = x1 
    else:
        t = np.linspace(0, time, samples)  # начало, завершение моделирования, число отсчетов
        x0 = [Gp0, Gt0, Il0, Ip0, Isc10, Isc20, Ione0, Id0, Qsto10, Qsto20, Qgut0, X0, Y0, Z0] # начальные значения
        x = odeint(ode_system, x0, t, hmax=1)
    return x, t


# Простой пример для проверки
def test_example():
    test_Gpb = (4.8 * MG_DL_TO_MMOL_L_CONVENTION_FACTOR) * VG
    test_Gp0 = (10.2 * MG_DL_TO_MMOL_L_CONVENTION_FACTOR) * VG
    food_and_insulin_time = 25  # Пока что жестко прошито, в планах разделение

    res, t = simulate(
        time=420,  # длительность моделирования, мин 
        samples=420,  # число временных отсчетов
        BW=77,  # вес тела
        D=0,  # масса принятых углеводов, мг
        Gpb=test_Gp0,  # базальный уровень глюкозы в крови
        Gp0=test_Gp0,  # начальный уровень глюкозы в крови
        Djins=0,  # доза быстрого инсулина
        IIRb=2,  # базальный уровень инсулина
        meal_time=food_and_insulin_time,  # время приема пищи, мин с начала моделирования
        injection_time=food_and_insulin_time,  # время введения нисулина, мин с начала моделирования
        ex_on=True,  # включение/выключение физической активности
        ex_start=20,  # время начала упражнения
        ex_finish=80,  # время завершения упражнения
        ex_hr=120,  # ЧСС во время упражнения
        HRb=60  # базальное значение ЧСС
    )
    print_graphs(res, t, 'test')


def str2bool(v):
    return v.lower() in ("yes", "true", "t", "1")


# Функция открывает файлы с данными (входными для модели и реальными)
# и запускает моделирование для каждой из тренировок
# Получив результаты моделирования, функция вызывает print_graphs из functions.py
# print_graphs сохраняет необходимые графики и выводит количественные оценки
def test_on_text_data():
    with open(os.path.join('experiments','data.txt'), 'r') as file:
        text = file.read()
    with open(os.path.join('experiments','real_bg.txt'), 'r') as file:
        real_data = file.read()
    counter = 0
    for line, real in zip(text.split('\n')[:-1], real_data.split('\n')[:-1]):
        counter += 1
        line = line.split(',')
        res, t = simulate(
            time=int(line[0]),  # длительность моделирования, мин 
            samples=int(line[1]),  # число временных отсчетов
            BW=int(line[2]),  # вес тела
            D=int(line[3]),  # масса принятых углеводов, мг
            Gpb=float(line[4]) * MG_DL_TO_MMOL_L_CONVENTION_FACTOR * VG,  # базальный уровень глюкозы в крови
            Gp0=float(line[5]) * MG_DL_TO_MMOL_L_CONVENTION_FACTOR * VG,  # начальный уровень глюкозы в крови
            Djins=float(line[6]),  # доза быстрого инсулина
            IIRb=float(line[7]),  # базальный уровень инсулина
            meal_time=int(line[8]),  # время приема пищи, мин с начала моделирования
            injection_time=int(line[9]),  # время введения нисулина, мин с начала моделирования
            ex_on=str2bool(line[10]),  # включение/выключение физической активности
            ex_start=int(line[11]),  # время начала упражнения
            ex_finish=int(line[12]),  # время завершения упражнения
            ex_hr=float(line[13]),  # ЧСС во время упражнения
            HRb=float(line[14])  # базальное значение ЧСС
        )
        real = real.split(',')
        food = None
        insulin = None
        if int(line[3]):
            food = (int(line[8]), int(line[3]))
        if int(line[6]):
            insulin = (int(line[9]), int(line[6]))

        print_graphs(res, t, counter, real, food, insulin)


if __name__ == "__main__":
    test_example()
    test_on_text_data()
    