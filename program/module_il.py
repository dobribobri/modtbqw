#coding=utf-8
from __future__ import division
from __future__ import print_function
from matplotlib import  pyplot as plt
import easygui
import math
from rec_itu_p676_3 import Model
from parser_ import Parser as Parser2
from proc import Proc
import numpy as np


def do_():
    path = easygui.fileopenbox()
    p = Parser2(path)
    p.parse()
    print("Данные загружены. ")
    # p.d_print()
    print("--------")

    THETA = int(path[len(path) - 11 : len(path) - 9])
    theta = (math.pi/180) * THETA

    '''
    model = Model(Temperature=15, Pressure=1013, Rho=7.5)
    print("Q и W для условий стандартной атмосферы (зенит), по встроенной модели:")
    q, w = model.get_QW_model_zenith(freq1=21, freq2=27, T_obl=-2)
    print("Q = ", q, "\tW = ", w, "\n--------\n")
    print("Q и W для условий стандартной атмосферы (theta = ", theta, "рад.), по встроенной модели:")
    q, w = model.get_QW_model_(freq1=21, freq2=27, T_obl=-2, theta=theta)
    print("Q = ", q, "\tW = ", w, "\n--------\n")
    '''

    model_ = Model(Temperature=15, Pressure=768*1.33322, Rho=7.5)
    pc = Proc(p.data, model_, T_avg=2.3, T_obl=-2, v_wind=None, time_interval=None)

    print("Температура ", model_.T, "цельс.,\tДавление ", model_.P, "гПа,\tАбс. влажность ", model_.rho, "г/м^3")
    print("T_avg ", pc.T_avg, "цельс.,\tТемпература облака ", pc.T_obl, "цельс.")
    print("Скорость ветра ", pc.wind, "км/ч,\tВременной интервал (период) ", pc.time, "ч")
    print("Начальная координата ", pc.start, "км")
    print("\nTHETA = ", THETA, "град.,\t", theta, "рад.")

    
    print("cos(Theta) = ", math.cos(theta))
    print("int_rho = ", model_.get_Q())

    # pc.set_T_t0_model_Tb(theta)
    pc.set_T_t0_Tb_avg()

    plt.figure(1)
    plt.title(u"Поглощение в кислороде")
    Freqs, Tbs = [], []
    for f in np.arange(18.0, 27.2, 0.2):
        Freqs.append(f)
        Tbs.append(model_.tauO_theory(f, theta))
    plt.plot(Freqs, Tbs, "r+")

    plt.figure(2)
    plt.title(u"Поглощение в водяном паре")
    Freqs, Tbs = [], []
    for f in np.arange(18.0, 27.2, 0.2):
        Freqs.append(f)
        Tbs.append(model_.tauRho_theory(f, theta))
    plt.plot(Freqs, Tbs, "r+")

    plt.figure(3)
    plt.title(u"$T_{b}$ - модельные значения (безоблачная атм.)\n Температура " + str(model_.T) +
              u" °C, Давление " + str(int(model_.P / 1.33322)) + u" мм.рт.ст., Абс. влажность " +
              str(model_.rho) + u" $g/m^{3}$")
    Freqs, Tbs = [], []
    for f in np.arange(18.0, 27.2, 0.2):
        Freqs.append(f)
        Tbs.append(model_.get_Tb(f, pc.T_avg, theta=theta))
    plt.plot(Freqs, Tbs, "r+")
    plt.show()
    plt.close()
    

    print("\n==== Средние за период яркостные температуры: ====")
    print("Здесь и далее точки {T_ярк. = 0} пропускаются")
    A = pc.get_Tb_avg()
    B = pc.get_calibr_Tb_avg(theta)
    print("18.7 ГГц:", A[18.7])
    print("Calibr.  :", B[18.7])
    print("19.1 ГГц:", A[19.1])
    print("Calibr.  :", B[19.1])
    print("21.1 ГГц:", A[21.1])
    print("Calib.  :", B[21.1])
    print("21.7 ГГц:", A[21.7])
    print("Calibr.  :", B[21.7])
    print("23.9 ГГц:", A[23.9])
    print("Calibr.  :", B[23.9])
    print("24.1 ГГц:", A[24.1])
    print("Calibr.  :", B[24.1])
    print("26.9 ГГц:", A[26.9])
    print("Calibr.  :", B[26.9])
    Freqs = [18.7, 19.1, 21.1, 21.7, 23.9, 24.1, 26.9]
    Tbs_model, Tbs_given, Tbs_calibr = [], [], []
    for freq in Freqs:
        Tbs_model.append(model_.get_Tb(freq, pc.T_avg, theta))
        Tbs_given.append(A[freq])
        Tbs_calibr.append(B[freq])
    plt.plot(Freqs, Tbs_model, "b+", label='Model')
    plt.plot(Freqs, Tbs_given, "g+", label='Data')
    plt.plot(Freqs, Tbs_calibr, "r+", label='Calibr')
    plt.legend()
    plt.show()

    print("\n----Двухчастотный метод определения Q и W----")
    print("(По средним за период ярк. температурам)")
    print("freq1 = 18.7 ГГц, \tfreq2 = 21.7 ГГц")
    q, w = pc.get_QW4Tb_avg(freq1=18.7, freq2=21.7, theta=theta)
    print("Q = ", q, "\t W = ", w)

    print("freq1 = 18.7 ГГц, \tfreq2 = 23.9 ГГц")
    q, w = pc.get_QW4Tb_avg(freq1=18.7, freq2=23.9, theta=theta)
    print("Q = ", q, "\t W = ", w)

    print("freq1 = 21.7 ГГц, \tfreq2 = 26.9 ГГц")
    q, w = pc.get_QW4Tb_avg(freq1=21.7, freq2=26.9, theta=theta)
    print("Q = ", q, "\t W = ", w)

    print("freq1 = 23.9 ГГц, \tfreq2 = 26.9 ГГц")
    q, w = pc.get_QW4Tb_avg(freq1=23.9, freq2=26.9, theta=theta)
    print("Q = ", q, "\t W = ", w)


    print("\n----Оптимизационный метод----")
    print("(По средним за период ярк. температурам)")
    print("4 частоты: 18.7, 21.7, 23.9, 26.9 ГГц")
    frequencies = [18.7, 21.7, 23.9, 26.9]
    q, w = pc.get_opt_QW4Tb_avg(frequencies, theta=theta)
    print("\nQ = ", q, "\t W = ", w)



    print("\n==== Средние за период Q и W ====")
    print("----Двухчастотный метод----")
    print("(По всем ярк. температурам)")
    print("freq1 = 18.7 ГГц, \tfreq2 = 21.7 ГГц")
    q, w = pc.get_QW_avg(freq1=18.7, freq2=21.7, theta=theta)
    print("Q = ", q, "\t W = ", w)

    print("freq1 = 18.7 ГГц, \tfreq2 = 23.9 ГГц")
    q, w = pc.get_QW_avg(freq1=18.7, freq2=23.9, theta=theta)
    print("Q = ", q, "\t W = ", w)

    print("freq1 = 21.7 ГГц, \tfreq2 = 26.9 ГГц")
    q, w = pc.get_QW_avg(freq1=21.7, freq2=26.9, theta=theta)
    print("Q = ", q, "\t W = ", w)

    print("freq1 = 23.9 ГГц, \tfreq2 = 26.9 ГГц")
    q, w = pc.get_QW_avg(freq1=23.9, freq2=26.9, theta=theta)
    print("Q = ", q, "\t W = ", w)

    '''
    pc = Proc(defaultdict(list), model)
    q, w = pc._get_QW(18, 7.57, 21, 15.9, (math.pi/180) * 51)
    print("Q = ", q, "\t W = ", w)
    q, w = pc._get_QW(18, 7.57, 22, 19.39, (math.pi/180) * 51)
    print("Q = ", q, "\t W = ", w)
    q, w = pc._get_QW(21, 15.9, 27, 13.45, (math.pi/180) * 51)
    print("Q = ", q, "\t W = ", w)
    q, w = pc._get_QW(22, 19.39, 27, 13.45, (math.pi/180) * 51)
    print("Q = ", q, "\t W = ", w)
    '''