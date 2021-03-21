#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from utils import solve_quadratic_equasion
from coefficients import *


# Эндогенное производство глюкозы
def EGP(t, Gp, Id):
    return kp1 - kp2 * Gp - kp3 * Id





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

    # Начальные значения
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


'''
Система дифф. уравнений
'''
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
    dGpdt = EGP(t, Gp, Id) + Ra(t, Qgut) - Uii(t) - E(t, Gp) - k1 * Gp + k2 * Gt
    if not ex_on:  # Если упражнения не включены - используем Uid а не Uid_ex
        dGtdt = -Uid(t, X, Gt) + k1 * Gp - k2 * Gt
    else:
        dGtdt = -Uid_ex(t, X, Y, Z, Gt) + k1 * Gp - k2 * Gt
    '''Подсистема инсулина'''
    dIldt = -(m1 + m3) * Il + m2 * Ip
    dIpdt = -(m2 + m4) * Ip + m1 * Il + R(t, Isc1, Isc2)  # У здорового нет R
    ''' Подсистема подкожного инсулина'''
    # ВАЖНО: Если введение инсулина происходит в момент времени 0, то Isc10 = Djins
    dIsc1dt = -(kd + ka1) * Isc1 + IIR(t)
    dIsc2dt = kd * Isc1 - ka2 * Isc2
    '''Эндогенная продукция глюкозы'''
    dIonedt = -ki * (Ione - I(t, Ip))
    dIddt = -ki * (Id - Ione)
    '''Глюкоза в пищеварительном тракте'''
    # ВАЖНО: Если прием пищи мгновенный (по умолчанию), ingestion(t) всегда 0, а Qsto10 = D
    dQsto1dt = -kgri *  Qsto1 + ingestion(t)
    dQsto2dt = -kempt(t, Qsto1, Qsto2) * Qsto2 + kgri * Qsto1
    dQgutdt = -kabs * Qgut + kempt(t, Qsto1, Qsto2) * Qsto2        
    '''Инсулин в межклеточной жидкости'''
    dXdt = -p2U * X + p2U * (I(t, Ip) - Ib)  # was deleted - Ib
    '''Физическая активность'''
    dYdt = -(1/Thr) * (Y - (HR(t) - HRb))
    dZdt = -(fex(t, Y)/ Tin + 1 / Tex) * Z + fex(t, Y)

    return [dGpdt, dGtdt, dIldt, dIpdt, dIsc1dt, dIsc2dt, dIonedt, dIddt, dQsto1dt, dQsto2dt, dQgutdt, dXdt, dYdt, dZdt]




if __name__ == "__main__":
    simulate()