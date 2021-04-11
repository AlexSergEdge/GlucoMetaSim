#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import numpy as np
import math

import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.stats import pearsonr, spearmanr, median_absolute_deviation
from sklearn.metrics import mean_absolute_error, mean_absolute_percentage_error, mean_squared_error

from coefficients import *

import xlwings

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
def kempt(t, Qsto1, Qsto2, D, meal_time, alphas, betas):
    if D != 0 and t > meal_time:
        return kmin + (kmax - kmin) / 2.0 * (math.tanh(alphas * (Qsto(t, Qsto1, Qsto2) - b * D)) - math.tanh(betas * (Qsto(t, Qsto1, Qsto2) - d * D)) + 2.0)
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


'''
Функции для моделирования и отображения
'''
# Функция обновляет входные данные для моделирования следующей части (при введении инсулина или приеме пищи)
def get_updated_init_conditions(x, food, insulin):
    return [
        x[-1,0], x[-1,1], x[-1,2], x[-1,3], 
        x[-1,4] + insulin, 
        x[-1,5] + kd / ka2 * insulin,
        x[-1,6], x[-1,7], 
        x[-1,8] + food,
        x[-1,9], x[-1,10], x[-1,11], x[-1,12], x[-1,13]
    ]


# Функция сравнивает реальные данные с данными моделирования и сохраняет графики и количественные параметры
def print_graphs(x, t, name, real=None, food=None, insulin=None):
    Gp = x[:,0] # Масса глюкозы в плазме и быстро-наполняюющихся тканях, mg/kg
    Gt = x[:,1] # Масса глюкозы в медленно-наполняющихся тканях, mg/kg
    Gres = G(0, Gp) # Концентрация глюкозы в плазме
    if not real:
        plt.plot(t, Gres / MG_DL_TO_MMOL_L_CONVENTION_FACTOR, label='Глюкоза в плазме', color='r')
        plt.title('Глюкоза в плазме крови, ммоль/л')
        plt.ylabel('Концентрация, ммоль/л')
        plt.xlabel('время, мин')
        plt.legend()
        plt.savefig(os.path.join('results','bg_model',f'bg_model_{name}.png'))
        plt.close()
    else:
        real_data = np.array([])
        for val in real:
            real_data = np.append(real_data, float(val))

        print(f'Training {name}')
        mae = mean_absolute_error(real_data * MG_DL_TO_MMOL_L_CONVENTION_FACTOR, Gres)
        print(f'MAE (mmol/l): {mae}')
        mape = mean_absolute_percentage_error(real_data * MG_DL_TO_MMOL_L_CONVENTION_FACTOR, Gres)
        print(f'MAPE (%): {mape}')
        rmse = mean_squared_error(real_data * MG_DL_TO_MMOL_L_CONVENTION_FACTOR, Gres, squared=False)
        print(f'RMSE (mg/dl): {rmse}')
        rhop, pvalp = pearsonr(real_data * MG_DL_TO_MMOL_L_CONVENTION_FACTOR, Gres)
        print(f'Pearson: rho = {rhop}, pval = {pvalp}')
        rhos, pvals = spearmanr(real_data * MG_DL_TO_MMOL_L_CONVENTION_FACTOR, Gres)
        print(f'Spearman: rho = {rhos}, pval = {pvals}')

        # Количественные параметры сохраним в файл
        app = xlwings.App(visible=False)
        wb = xlwings.Book('results.xlsx')
        ws = wb.sheets[0]
        ws.range(int(name) + 1, 2).value = mae
        ws.range(int(name) + 1, 3).value = mape
        ws.range(int(name) + 1, 4).value = rmse
        ws.range(int(name) + 1, 5).value = rhop
        ws.range(int(name) + 1, 6).value = rhos
        wb.save()
        wb.close()
        app.quit()

        plt.plot(t, Gres, label='Глюкоза в плазме (модель)', color='r')
        plt.plot(t, real_data * MG_DL_TO_MMOL_L_CONVENTION_FACTOR, label='Глюкоза в плазме (реальная)', color='b')
        plt.title('Глюкоза в плазме крови, мг/дл')
        plt.ylabel('Концентрация, мг/дл')
        plt.xlabel('время, мин')

        bottom, top = plt.ylim()
        food_y_text = (bottom + top) / 2
        ins_y_text = bottom + (top - bottom) / 4

        # Перевод из пмоль/л в единицы инсулина 1 μIU/mL = 6.00 pmol/L
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6501531/
        if food:
            plt.axvline(x=food[0], c='orange', ymax=1.0, lw=1)
            food_info = "Прием пищи\n({} мг.)".format(food[1])
            plt.text(food[0] + 7, food_y_text, food_info, rotation=90, fontsize=9, ha='center')
        if insulin:
            plt.axvline(x=insulin[0], c='cyan', ymax=1.0, lw=1)
            ins_info = 'Введение инсулина\n({:.2f} ед.)'.format(insulin[1] * (kd + ka1) / 6.00)
            plt.text(insulin[0] + 4, ins_y_text, ins_info, rotation=90, fontsize=9, ha='center')

        plt.legend()
        plt.savefig(os.path.join('results','bg_compare',f'bg_compare_{name}.png'))
        plt.close()