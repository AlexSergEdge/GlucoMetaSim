#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import datetime
import numpy as np
import matplotlib.pyplot as plt
import os

def to_minutes(timestr):
    minutes = 0
    for part in timestr.split(':'):
        minutes = minutes*60 + int(part)
    return minutes


def plot_glucose(glucose_info, show_plot=True):

    bg_values = []
    time_values = []


    for info in glucose_info:

        start_time = info['glucose_list'][0][1]
        end_time = info['glucose_list'][-1][1]

        final_glucose_list = []
        final_time_list = []
        g_prev = info['glucose_list'][0][0]
        t_prev = info['glucose_list'][0][1]

        final_glucose_list.append(g_prev)
        final_time_list.append(t_prev)

        for g, t in info['glucose_list'][1:]:
            curr_t = t_prev + 1
            while(curr_t <= t):
                # линейная интерполяция
                glucose_new = g_prev + (curr_t - t_prev) * (g - g_prev)/(t - t_prev)
                final_glucose_list.append(glucose_new)
                final_time_list.append(curr_t)
                curr_t += 1
            g_prev = g
            t_prev = t

        if show_plot:
            plt.plot(final_time_list, final_glucose_list, label='Glucose', color='y')
            plt.title('Plasma glucose')
            plt.ylabel('G, mmol/l')
            plt.xlabel('time, min')
            plt.legend()
            plt.show()

        # print(start_time)
        # print(end_time)
        # print(final_time_list)
        # exit(0)

        bg_values.append(final_glucose_list)
        time_values.append(final_time_list)

    return bg_values, time_values

def plot_hr(heart_rate_info, show_plot=True):
    
    hr_values = []
    time_values = []
    
    for info in heart_rate_info:
        hr_arr = []
        t_arr = []
        for hr, t in info['hr_list']:
            hr_arr.append(hr)
            t_arr.append(t)
        
        if show_plot:
            plt.plot(t_arr, hr_arr, label='HR', color='r')
            plt.title('HR')
            plt.ylabel('bpm')
            plt.xlabel('time, min')
            plt.legend()
            plt.show()

        hr_values.append(hr_arr)
        time_values.append(t_arr)
    return hr_values, time_values

TRAININGS = {
    '1': {
        'date': '2020-12-11',
        'start': '13:29',
        'finish': '14:05'
    },
    '2': {
        'date': '2020-12-11',
        'start': '18:10',
        'finish': '19:09'
    },
    '3': {
        'date': '2020-12-12',
        'start': '14:18',
        'finish': '14:48'
    },
    '4': {
        'date': '2020-12-13',
        'start': '17:19',
        'finish': '18:18'
    },
    '5': {
        'date': '2020-12-23',
        'start': '13:04',
        'finish': '13:54' 
    },
    '6': {
        'date': '2020-12-25',
        'start': '19:39',
        'finish': '21:35'
    },
    '7': {
        'date': '2020-12-27',
        'start': '17:33',
        'finish': '18:47'
    }
}

# [mass, time]
CARB_INTAKE = [
    [0, 0], 
    [0, 0], 
    [0, 0], 
    [0, 0], 
    [0, 0], 
    [0, 0], 
    [0, 0]
]
# [dose, time]
INSULIN_INTAKE = [
    [0, 0], 
    [0, 0], 
    [0, 0], 
    [0, 0], 
    [0, 0], 
    [0, 0], 
    [0, 0]
]



def get_glucose_info():
    glucose_info = []
    with open(os.path.join('analytics','data','BG.txt'), 'r') as bg_file:
        for _ in range(3):
            next(bg_file)
        for training in TRAININGS.items():
            required_date = training[1]['date']
            required_start = to_minutes(training[1]['start'])
            required_finish = to_minutes(training[1]['finish'])
            
            glucose_list = []
            
            for line in bg_file:
                line = line.split()
                date = line[1]
                time = to_minutes(line[2])
                glucose = line[4]
                if date != required_date.replace('-','/'):
                    continue
                if (time <= required_start - 30):
                    continue
                if (time >= required_finish + 30):
                    continue
                glucose_list.append((float(glucose), int(time)))
        
            glucose_list = sorted(glucose_list, key=lambda x: x[1])                       
            return_dict = {
                'name': training[0],
                'date': training[1]['date'],
                'glucose_list': glucose_list 
            }
            glucose_info.append(return_dict)
            
            bg_file.seek(0)
            for _ in range(3):
                next(bg_file)
    return glucose_info


def get_heart_rate_info():
    average_hrs = []      
    heart_rate_info = []   
    for training in TRAININGS.items():
        required_date = training[1]['date']
        required_start = to_minutes(training[1]['start'])
        required_finish = to_minutes(training[1]['finish'])
        
        heart_rate_list = []
        
        with open(os.path.join('analytics','data','HR.csv'), 'r') as hr_file:
            csv_reader = csv.reader(hr_file, delimiter=',')
            count = 0
            hr_accum = 0
            hr_samples = 0
            for row in csv_reader:
                if count == 0:
                    count += 1
                    continue
                count += 1
                date = row[0]
                time = to_minutes(row[1])
                heart_rate = row[2]
                if date != required_date:
                    continue
                if (time < required_start):
                    continue
                if (time > required_finish):
                    continue
                hr_accum += float(heart_rate)
                hr_samples += 1
                heart_rate_list.append((float(heart_rate), int(time)))


            #heart_rate_list = sorted(heart_rate_list, key=lambda x: x[1])

            average_hrs.append(hr_accum / hr_samples)

            return_dict = {
                'name': training[0],
                'date': training[1]['date'],
                'hr_list': heart_rate_list 
            }
            heart_rate_info.append(return_dict)
    return heart_rate_info, average_hrs


def get_timespans(time_values):
    start_times = [i[0] for i in time_values]
    finish_times = [i[-1] for i in time_values]
    timespans = []
    for st, ft in zip(start_times, finish_times):
        timespans.append(ft - st + 1)
    return timespans, start_times, finish_times


if __name__ == "__main__":
    glucose_info = get_glucose_info()
    
    bg_values, time_values = plot_glucose(glucose_info, False)
    print(len(time_values[0]))

    start_glucoses = [i[0] for i in bg_values]
    print('Start glucoses: ', start_glucoses)

    timespans_bg, start_times_bg, finish_times_bg = get_timespans(time_values)
    print(timespans_bg[0])

    heart_rate_info, average_hrs = get_heart_rate_info()
    print('Average HRs:', average_hrs)
    hr_values, time_values = plot_hr(heart_rate_info, False)

    timespans_hr, start_times_hr, finish_times_hr = get_timespans(time_values)

    modelling_start_times = []
    ex_start_times = []
    ex_real_start_times = []
    
    for sbg, shr in zip(start_times_bg, start_times_hr): 
        modelling_start_times.append(0)
        if sbg < shr:
            ex_start_times.append(shr - sbg)
        else:
            ex_start_times.append(0)

    modelling_finish_times = []
    ex_finish_times = []

    for fbg, fhr, timespan_bg in zip(finish_times_bg, finish_times_hr, timespans_bg):
        modelling_finish_times.append(timespan_bg)
        if fbg > fhr:
            ex_finish_times.append(timespan_bg - (fbg - fhr))
        else:
            ex_finish_times.append(timespan_bg)

    print('modelling_start_times: ', modelling_start_times)
    print('modelling_finish_times: ', modelling_finish_times)
    print('ex_start_times: ', ex_start_times)
    print('ex_finish_times: ', ex_finish_times)

    with open(os.path.join('experiments', 'real_bg.txt'), 'w+') as file:
        for i in range(len(TRAININGS.keys())):
            line = ','.join(str(v) for v in bg_values[i])
            file.write(line + '\n')


    with open(os.path.join('experiments', 'data.txt'), 'w+') as file:
        for i in range(len(TRAININGS.keys())):
            time = modelling_finish_times[i]
            samples = modelling_finish_times[i]
            BW = 62
            D = CARB_INTAKE[i][0]
            Gpb = start_glucoses[i]
            Gp0 = start_glucoses[i]
            Djins = INSULIN_INTAKE[i][0]
            IIRb = 2
            meal_time = CARB_INTAKE[i][1]
            injection_time = INSULIN_INTAKE[i][1]
            ex_on = True
            ex_start = ex_start_times[i]
            ex_finish = ex_finish_times[i]
            ex_hr = average_hrs[i]
            HRb = 60
            file.write(f'{time},{samples},{BW},{D},{Gpb},{Gp0},{Djins},{IIRb},{meal_time},{injection_time},{ex_on},{ex_start},{ex_finish},{ex_hr},{HRb}\n')
         


