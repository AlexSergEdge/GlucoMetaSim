#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from coefficients import *
from scipy.integrate import quad
import matplotlib.pyplot as plt

'''Функции'''
# Эндогенное производство глюкозы
def EGP(t, Gp, Id):
    return kp1 - kp2 * Gp - kp3 * Id

# Скорость увеличения концентрации глюкозы в крови
def Ra(t, Qgut, BW):
    return (f * kabs * Qgut) / BW

# Инсулинонезависимая утилизация
def Uii(t):
    return Fcns

# Почечная экскреция
def E(t, Gp): 
    if Gp > ke2:
        return ke1 * (Gp - ke2)
    else:
        return 0

# Инсулинозависимая утилизация
def Uid(t, X, Gt):
    Vm = Vm0 + Vmx * X
    Km = Km0 # + Kmx * X
    return (Vm * Gt) / (Km + Gt)

'''Функции физ. активности'''
# Инсулинозависимая утилизация при упражнениях
def Uid_ex(t, X, Y, Z, Gt, Ib, HRb, ex_start, ex_finish, ex_hr):
    dividend = Vm0 * (1 + betaex * Y) + Vmx * (1 + alphaex * Z * W(t, Z, HRb, ex_start, ex_finish, ex_hr)) * (X + Ib) - Vmx * Ib
    divisor = Km0 * (1 - upsilonex * Z * W(t, Z, HRb, ex_start, ex_finish, ex_hr) * (X + Ib)) + Gt
    return (dividend / divisor) * Gt 

# Вычисляем площадь под графиком ЧСС
def W(tt, Z, HRb, ex_start, ex_finish, ex_hr):
    if Z > 0.0:
        y = lambda x: HR(x, ex_start, ex_finish, ex_hr, HRb) - HRb
        y, err = quad(y, 0, tt, limit=100)
        return y
    return 0

# Функция ЧСС (ступенька)
def HR(tt, ex_start, ex_finish, ex_hr, HRb):
    if tt > ex_start and tt < ex_finish:
        return ex_hr
    return HRb

# Вспомогательная функция
def fex(t, Y, HRb):
    dividend =  (Y / (a * HRb)) ** n
    divisor = 1 + (Y / (a * HRb)) ** n
    return dividend / divisor

'''Другие функции'''
# Утилизация глюкозы суммарная
def U(t, Uid, Uii):
    return Uid + Uii

# Концентрация инсулина в плазме крови
def I(t, Ip):
    return Ip / VI

# Концентрация глюкозы в плазме крови
def G(t, Gp):
    return Gp / VG

# Общий объем глюкозы (твердой + мягкой) в желудке
def Qsto(t, Qsto1, Qsto2):
    return Qsto1 + Qsto2

# Динамический параметр опустошения желудка
def kempt(t, Qsto1, Qsto2, D, meal_time):
    if D != 0 and t > meal_time:
        return kmin + (kmax - kmin) / 2.0 * (tanh(alphas * (Qsto(t, Qsto1, Qsto2) - b * D)) - tanh(betas * (Qsto(t, Qsto1, Qsto2) - d * D)) + 2.0)
    else:
        return kmin

# Скорость появления инсулина в плазме
def R(t, Isc1, Isc2):
    return ka1 * Isc1 + ka2 * Isc2

# При необходимости ставить временные рамки
def ingestion(time):
    return 0
# При необходимости ставить временные рамки
def infusion(time):
    return 0
# Доза инсулина
def IIR(t, IIRb):
    return infusion(t) + IIRb



def print_graphs(x, t):
    Gp = x[:,0] # Масса глюкозы в плазме и быстро-наполняюющихся тканях, mg/kg
    Gt = x[:,1] # Масса глюкозы в медленно-наполняющихся тканях, mg/kg
    Gres = G(0, Gp) # Концентрация глюкозы в плазме
    plt.plot(t, Gres / MG_DL_TO_MMOL_L_CONVENTION_FACTOR, label='Plasma glucose', color='r')
    plt.title('Plasma glucose, mmol/l')
    plt.ylabel('G, mmol/l')
    plt.xlabel('time, min')
    plt.legend()
    plt.show()