#coding=utf-8
from __future__ import division
from __future__ import print_function
from matplotlib import  pyplot as plt
import easygui
import math
from rec_itu_p676_3 import Model
from parser import Parser as Parser1
from proc import Proc
import numpy as np
from drawer import Drawer
from collections import defaultdict
from matplotlib import rc
from configer import Configer


def do_():
    rc('font', **{'family': 'serif'})
    rc('text', usetex=True)
    rc('text.latex', unicode=True)
    rc('text.latex', preamble=r"\usepackage[T2A]{fontenc}")
    rc('text.latex', preamble=r"\usepackage[utf8]{inputenc}")
    rc('text.latex', preamble=r"\usepackage[russian]{babel}")

    path = easygui.fileopenbox()
    p = Parser1("tb_timeP", path)
    DATA = p.parse()

    p.shift_tb(DATA)

    QW_freqs = [18, 21, 22, 27]
    freqs2show = [18, 19.2, 20.4, 21.6, 22.2, 22.4, 23.2, 24.4, 25.6, 26.8]

    Drawer().draw(DATA, freqs2show, title=u"$T_{b}$ - без калибровки (экспериментальные данные)",
                  xlabel=u"время (сек.)", ylabel="")


    config_path = easygui.fileopenbox()
    info = Configer(config_path)
    info.get_info()


    # ==================================================================
    theta = info.theta
    model_ = Model(Temperature=info.T, Pressure=info.P, Rho=info.rho)
    proc = Proc(DATA, model_, T_avg=info.Tavg, T_obl=info.Tobl,
                start_coordinate=info.start, stop_coordinate=info.stop,
                is_data_calibred=False)
    # ==================================================================


    plt.title("")
    Freqs, TauO, TauH2O = [], [], []
    for f in np.arange(5, 350, 0.5):
        Freqs.append(f)
        TauO.append(model_.tauO_theory(f, theta=51) / model_.dB2np)
        TauH2O.append(model_.tauRho_theory(f, theta=51) / model_.dB2np)
    plt.plot(Freqs, TauO, label=u"Кислород", color="red")
    plt.plot(Freqs, TauH2O, label=u"Водяной пар", color="blue")
    ax = plt.axes()
    ax.set_xscale("log", nonposx='clip')
    ax.set_yscale("log", nonposy='clip')
    plt.xlabel(u"ГГц")
    plt.ylabel(u"Дб")
    plt.legend(loc="best")
    plt.show()
    plt.close()


    plt.title(u"$T_{b}$ - модельные значения (безоблачная атм.)\n Температура " + str(model_.T) +
              u" °C, Давление " + str(int(model_.P/1.33322)) +  u" мм.рт.ст., Абс. влажность " +
              str(model_.rho) + u" $g/m^{3}$")
    FREQS, TBR, TBRnQ = [], [], []
    for freq in np.arange(18, 27.2, 0.2):
        T_br = model_.get_Tb_Q(freq, proc.T_avg, theta)
        Tbrnq = model_.get_Tb(freq, proc.T_avg, theta)
        FREQS.append(freq)
        TBR.append(T_br)
        TBRnQ.append(Tbrnq)
    Q = model_.get_Q()
    plt.plot(FREQS, TBRnQ, "b+", label=u'$T_{b}$ : $T_{avg}(1-e^{-tau})$')
    plt.plot(FREQS, TBR, "r+", label=u'$T_{b}$ :  Q = ' + str(Q) + u' $g/m^{2}$')
    print("Q_model = ", Q)
    plt.xlabel(u'ГГц')
    plt.ylabel(u'K')
    plt.legend()
    plt.show()
    plt.close()

    plt.title(u"Поглощение в кислороде - модель (безобл. атм.)")
    Freqs, Tbs = [], []
    for f in np.arange(18.0, 27.2, 0.2):
        Freqs.append(f)
        Tbs.append(model_.tauO_theory(f, theta))
    plt.plot(Freqs, Tbs, "r+")
    plt.xlabel(u"ГГц")
    plt.ylabel(u"непер")
    plt.show()
    plt.close()

    plt.title(u"Поглощение в водяном паре - модель (безобл. атм.)")
    Freqs, Tbs = [], []
    for f in np.arange(18.0, 27.2, 0.2):
        Freqs.append(f)
        Tbs.append(model_.tauRho_theory(f, theta))
    plt.plot(Freqs, Tbs, "r+")
    plt.xlabel(u"ГГц")
    plt.ylabel(u"непер")
    plt.show()
    plt.close()

    plt.title("$k_{w}$")
    Freqs, K = [], []
    for f in np.arange(18.0, 27.2, 0.01):
        Freqs.append(f)
        K.append(model_.kw(f, -2))
    plt.plot(Freqs, K)
    plt.xlabel(u"ГГц")
    plt.ylabel(u"Значение")
    plt.show()
    plt.close()


    print(info.nclbeg, info.nclend)
    proc.set_T_t0_nocluds_interval(time_begin=info.nclbeg, time_end=info.nclend)
    proc.calibr_Tb(theta=theta)

    Drawer().draw(proc.data, freqs2show, title=u"$T_{b}$ - яркостные температуры",
                  xlabel=u"время", ylabel=u"температура [K]")

    # ==================================================================== #
    # ======================= Generate DataSet =========================== #
    # ==================================================================== #
    proc.generate_dataset("QW__f1_18GHz__f2_21GHz.txt", freq1=18, freq2=21, theta=theta)
    proc.generate_dataset("QW__f1_18GHz__f2_22GHz.txt", freq1=18, freq2=22, theta=theta)
    proc.generate_dataset("QW__f1_21GHz__f2_27GHz.txt", freq1=21, freq2=27, theta=theta)
    proc.generate_dataset("QW__f1_22GHz__f2_27GHz.txt", freq1=22, freq2=27, theta=theta)
    # ==================================================================== #


    QW = defaultdict(list)
    QW["18, 21"] = proc.get_resolved_QW(freq1=18, freq2=21, theta=theta)
    QW["18, 22"] = proc.get_resolved_QW(freq1=18, freq2=22, theta=theta)
    QW["21, 27"] = proc.get_resolved_QW(freq1=21, freq2=27, theta=theta)
    QW["22, 27"] = proc.get_resolved_QW(freq1=22, freq2=27, theta=theta)
    # QW["18, 27"] = proc.get_resolved_QW(freq1=18, freq2=27, theta=theta)
    # QW["18, 27"] = proc.get_QW(freq1=18, freq2=27, theta=theta)

    QW_err = defaultdict(list)
    QW_err["18, 21"] = proc.get_QW_errors(freq1=18, freq2=21, theta=theta, Tbr_error=1, Tavg_error=2)
    QW_err["18, 22"] = proc.get_QW_errors(freq1=18, freq2=22, theta=theta, Tbr_error=1, Tavg_error=2)
    QW_err["21, 27"] = proc.get_QW_errors(freq1=21, freq2=27, theta=theta, Tbr_error=1, Tavg_error=2)
    QW_err["22, 27"] = proc.get_QW_errors(freq1=22, freq2=27, theta=theta, Tbr_error=1, Tavg_error=2)
    # QW_err["18, 27"] = proc.get_QW_errors(freq1=18, freq2=27, theta=theta, Tbr_error=1, Tavg_error=2)


    plt.figure(1)
    plt.title(u"Q - полная масса водяного пара")
    k = 0
    for key in QW.keys():
        TIME, Q = [], []
        dTIME, dQ, dQ_err = [], [], []
        i = 0
        for time, q, _ in QW[key]:
            TIME.append(time)
            Q.append(q)
            if (i - k) % 100 == 0:
                dTIME.append(time)
                _, dq, _ = QW_err[key][i]
                dQ.append(q)
                dq = dq / 10 #!!!
                dQ_err.append(dq)
            i += 1
        plt.plot(TIME, Q, label=key+u" ГГц")
        ecolor = ""
        if k == 0: ecolor = "blue"
        if k == 25: ecolor = "orange"
        if k == 50: ecolor = "green"
        if k == 75: ecolor = "red"
        if k == 100: ecolor = "purple"
        # plt.errorbar(dTIME, dQ, yerr=dQ_err, fmt='o', ecolor="black", mfc=ecolor, mec=ecolor, ms=2, mew=3, capsize=2, capthick=3, elinewidth=3)
        k += 25
    plt.xlabel(u"время")
    plt.ylabel("$g/cm^{2}$")
    axes = plt.gca()
    # axes.set_ylim([1., 2.5])
    plt.legend()

    plt.figure(2)
    plt.title(u"W - водозапас облаков")
    k = 0
    for key in QW.keys():
        TIME, W = [], []
        dTIME, dW, dW_err = [], [], []
        i = 0
        for time, _, w in QW[key]:
            TIME.append(time)
            w = w * 10
            W.append(w)
            if (i - k) % 100 == 0:
                dTIME.append(time)
                _, _, dw = QW_err[key][i]
                dW.append(w)
                dW_err.append(dw)
            i += 1
        plt.plot(TIME, W, label=key+u" ГГц")
        ecolor = ""
        if k == 0: ecolor = "blue"
        if k == 25: ecolor = "orange"
        if k == 50: ecolor = "green"
        if k == 75: ecolor = "red"
        if k == 100: ecolor = "purple"
        # plt.errorbar(dTIME, dW, yerr=dW_err, fmt='o', ecolor="black", mfc=ecolor, mec=ecolor, ms=2, mew=3, capsize=2, capthick=3, elinewidth=3)
        k += 25
    axes = plt.gca()
    # axes.set_ylim([-0.5, 1])
    plt.xlabel(u"время")
    plt.ylabel("$kg/m^{2}$")
    plt.legend()

    plt.show()


    QW_opt = proc.get_opt_QW(QW_freqs, theta)

    plt.figure(1)
    TIME, Q, W = [], [], []
    for time, q, w in QW_opt:
        TIME.append(time)
        Q.append(q)
        W.append(w * 10)
    plt.plot(TIME, Q, color="red", label="Q")
    plt.plot(TIME, W, color="blue", label="W")
    plt.xlabel(u"время")
    plt.ylabel("Q [$g/cm^{2}$], W [$kg/m^{2}$]")
    plt.legend()

    '''
    plt.figure(3)
    plt.title("dQ")
    for key in QW_err:
        TIME, dQ = [], []
        for time, dq, _ in QW_err[key]:
            TIME.append(time)
            dQ.append(dq)
        plt.plot(TIME, dQ, label=key + " GHz")
    plt.xlabel("time")
    plt.ylabel("$value$")
    plt.legend()

    plt.figure(4)
    plt.title("dW")
    for key in QW_err:
        TIME, dW = [], []
        for time, _, dw in QW_err[key]:
            TIME.append(time)
            dW.append(dw)
        plt.plot(TIME, dW, label=key + " GHz")
    plt.xlabel("time")
    plt.ylabel("$value$")
    plt.legend()
    '''

    plt.show()

    '''
    plt.title("")
    TIME, TAU1821, TAU1822, TAU2721, TAU2722, TAU1827 = [], [], [], [], [], []
    min_len = len(proc.data[18])
    for f in [18, 21, 22, 27]:
        if len(proc.data[f]) < min_len: min_len = len(proc.data[f])
    for i in range(min_len):
        t18, Tb18 = proc.data[18][i]
        t21, Tb21 = proc.data[21][i]
        t22, Tb22 = proc.data[22][i]
        t27, Tb27 = proc.data[27][i]
        
        tau18 = model_.tau_experiment(Tb18, proc.T_avg + 273, theta)
        tau21 = model_.tau_experiment(Tb21, proc.T_avg + 273, theta)
        tau22 = model_.tau_experiment(Tb22, proc.T_avg + 273, theta)
        tau27 = model_.tau_experiment(Tb27, proc.T_avg + 273, theta)
        
        tau18 = model_.tauRho_theory(18, theta)
        tau21 = model_.tauRho_theory(21, theta)
        tau22 = model_.tauRho_theory(22, theta)
        tau27 = model_.tauRho_theory(27, theta)
        
        
        theta = 0
        tau1821 = math.fabs(model_.krho(18, theta) * model_.kw(21, -2) - model_.krho(21, theta) * model_.kw(18, -2))
        tau1822 = math.fabs(model_.krho(18, theta) * model_.kw(22, -2) - model_.krho(22, theta) * model_.kw(18, -2))
        tau2721 = math.fabs(model_.krho(21, theta) * model_.kw(27, -2) - model_.krho(27, theta) * model_.kw(21, -2))
        tau2722 = math.fabs(model_.krho(22, theta) * model_.kw(27, -2) - model_.krho(27, theta) * model_.kw(22, -2))
        tau1827 = math.fabs(model_.krho(18, theta) * model_.kw(27, -2) - model_.krho(27, theta) * model_.kw(18, -2))

        t = (t18 + t21 + t22 +t27) / 4
        TIME.append(t)
        TAU1821.append(tau1821)
        TAU1822.append(tau1822)
        TAU1827.append(tau1827)
        TAU2721.append(tau2721)
        TAU2722.append(tau2722)
    plt.plot(TIME, TAU1821, label="18 - 21")
    plt.plot(TIME, TAU1822, label="18 - 22")
    plt.plot(TIME, TAU2721, label="21 - 27")
    plt.plot(TIME, TAU2722, label="22 - 27")
    plt.plot(TIME, TAU1827, label="18 - 27")
    plt.legend()
    plt.show()

    print("18-21", TAU1821[0])
    print("18-22", TAU1822[0])
    print("21-27", TAU2721[0])
    print("22-27", TAU2722[0])
    print("18-27", TAU1827[0])
    '''

    print("QW avg's - no clouds interval")
    avgQ, avgW = 0, 0
    qavg, wavg = proc.get_QW_err_avg_noclouds(18, 21, theta, 1, 2)
    avgQ += math.fabs(qavg / 10)
    avgW += math.fabs(wavg)
    print("18-21:\tQ = ", qavg / 10, "\tW = ", wavg)
    qavg, wavg = proc.get_QW_err_avg_noclouds(18, 22, theta, 1, 2)
    avgQ += math.fabs(qavg / 10)
    avgW += math.fabs(wavg)
    print("18-22:\tQ = ", qavg / 10, "\tW = ", wavg)
    qavg, wavg = proc.get_QW_err_avg_noclouds(21, 27, theta, 1, 2)
    avgQ += math.fabs(qavg / 10)
    avgW += math.fabs(wavg)
    print("21-27:\tQ = ", qavg / 10, "\tW = ", wavg)
    qavg, wavg = proc.get_QW_err_avg_noclouds(22, 27, theta, 1, 2)
    avgQ += math.fabs(qavg / 10)
    avgW += math.fabs(wavg)
    print("22-27:\tQ = ", qavg / 10, "\tW = ", wavg)
    avgQ /= 4
    avgW /= 4
    print("AVG:\tQ = ", avgQ, "\tW = ", avgW)

    print("QW avg's - clouds interval")
    avgQ, avgW = 0, 0
    qavg, wavg = proc.get_QW_err_avg_clouds(18, 21, theta, 1, 2)
    avgQ += math.fabs(qavg / 10)
    avgW += math.fabs(wavg)
    print("18-21:\tQ = ", qavg / 10, "\tW = ", wavg)
    qavg, wavg = proc.get_QW_err_avg_clouds(18, 22, theta, 1, 2)
    avgQ += math.fabs(qavg / 10)
    avgW += math.fabs(wavg)
    print("18-22:\tQ = ", qavg / 10, "\tW = ", wavg)
    qavg, wavg = proc.get_QW_err_avg_clouds(21, 27, theta, 1, 2)
    avgQ += math.fabs(qavg / 10)
    avgW += math.fabs(wavg)
    print("21-27:\tQ = ", qavg / 10, "\tW = ", wavg)
    qavg, wavg = proc.get_QW_err_avg_clouds(22, 27, theta, 1, 2)
    avgQ += math.fabs(qavg / 10)
    avgW += math.fabs(wavg)
    print("22-27:\tQ = ", qavg / 10, "\tW = ", wavg)
    avgQ /= 4
    avgW /= 4
    print("AVG:\tQ = ", avgQ, "\tW = ", avgW)