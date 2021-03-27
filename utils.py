#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from math import sqrt

def solve_quadratic_equasion(a_eq, b_eq, c_eq):
    discr = (b_eq**2) - (4*a_eq*c_eq)

    if discr == 0:
        eq_res1 = (-b_eq + sqrt(discr))
        eq_res2 = (-b_eq + sqrt(discr))
    elif discr > 0:
        eq_res1 = (-b_eq + sqrt(discr)) / (2*a_eq)
        eq_res2 = (-b_eq - sqrt(discr)) / (2*a_eq)
    else:
        print("No real solutions exiting")
        exit(1)
    if eq_res1 >= 0:
        return eq_res1
    elif eq_res2 >= 0:
        return eq_res2
    else:
        print("No positive solution exiting")
        exit(1)