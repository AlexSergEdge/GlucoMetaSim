#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
from math import tanh

import numpy as np
from scipy.integrate import odeint

from utils import solve_quadratic_equasion
from coefficients import *
from functions import *

class Simulation():
    def __init__(
        self,
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
        self.time = time
        self.samples = samples
        self.BW = BW
        self.D = D
        self.Gpb = Gpb
        self.Gp0 = Gp0
        self.Djins = Djins
        self.IIRb = IIRb
        self.meal_time = meal_time
        self.injection_time = injection_time
        self.ex_on = ex_on
        self.ex_start = ex_start
        self.ex_finish = ex_finish
        self.ex_hr = ex_hr
        self.HRb = HRb
    
    def set_meal(self, meal):
        self.D = meal

    def set_insulin(self, ins):
        self.Djins = ins


def simulate_enhanced(simulation): 
    
    time = simulation.time
    samples = simulation.samples
    BW = simulation.BW
    D = simulation.D
    Gpb = simulation.Gpb
    Gp0 = simulation.Gp0
    Djins = simulation.Djins
    IIRb = simulation.IIRb
    meal_time = simulation.meal_time
    injection_time = simulation.injection_time
    ex_on = simulation.ex_on
    ex_start = simulation.ex_start
    ex_finish = simulation.ex_finish
    ex_hr = simulation.ex_hr
    HRb = simulation.HRb

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