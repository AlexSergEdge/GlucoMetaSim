#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# Есть формула 20: 
#   Gtb = (Fcns - EGPb + k1 * Gp) / k2 (20)
# EGP есть в формуле 21:
#   EGPb = Fcns + (Vm0 * Gtb) / (Km0 + Gtb) (21)
# Поэтому решаем квадратное уравнение:
#   k2 * Gtb ^ 2 + (k2*Km0 - k1*Gpb + Vm0) * Gtb - k1 * Km0 * Gpb = 0
def solve_quadratic_equasion(a_eq, b_eq, c_eq):
    discr = (b_eq**2) - (4*a_eq*c_eq)

    if discr == 0:
        eq_res1 = (-b_eq + sqrt(discr))
        eq_res2 = (-b_eq + sqrt(discr))
    elif discr > 0:
        eq_res1 = (-b_eq + sqrt(discr)) / (2*a_eq)
        eq_res2 = (-b_eq - sqrt(discr)) / (2*a_eq)
    else:
        print("No real solutions for Gtb, exiting")
        exit(1)
    if eq_res1 >= 0:
        return eq_res1
    elif eq_res2 >= 0:
        return eq_res2
    else:
        print("No positive solution for Gtb, exiting")
        exit(1)