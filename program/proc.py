#coding=utf-8
from __future__ import division
from __future__ import print_function
from collections import defaultdict
import math
import numpy as np
import sys
from rec_itu_p676_3 import Model


class Proc:
    # v_wind - скорость ветра [км/ч]
    # time_interval - временной интервал [ч]
    # T_avg = T0 - dT - [цельс.] - средняя абсолютная температура атмосферы, равная температуре изотермической атмосферы,
    #                   которая даёт то же излучение, что и реальная атмосфера, T0 - абсолютная температура
    #                   вблизи поверхности Земли, dT - поправка на неизотермичность атмосферы (от 5 до 32 градусов)
    # T_obl - предполагаемая температура облака [цельс]
    # theta - зенитный угол [рад.]

    # инициализация
    def __init__(self, data, MODEL = Model(),
                 T_avg = 15,
                 T_obl = -2,
                 v_wind = None,
                 time_interval = None,
                 start_coordinate = 0,
                 stop_coordinate = None,
                 is_data_calibred = False):
        self.data = data
        self.MODEL = MODEL
        self.T_avg = T_avg
        self.T_obl = T_obl
        self.start = start_coordinate
        self.upp_x = None
        if v_wind and time_interval:
            self.upp_x = start_coordinate + v_wind * time_interval
        self.stop = None
        if stop_coordinate:
            self.stop = stop_coordinate
        else:
            self.stop = self.upp_x
        self.wind = v_wind
        self.time = time_interval
        self.__cut_data()
        self.Tbr_t0 = {}
        self.set_T_t0_min()
        self.A = {}
        self.isDataCalibred = is_data_calibred
        self.noclouds_begin = None
        self.noclouds_end = None

    # обрезать данные по времени
    def __cut_data(self):
        if self.stop is None: return
        data = defaultdict(list)
        for key in self.data.keys():
            for x, temp in self.data[key]:
                if self.start <= x <= self.stop:
                    data[key].append((x, temp))
        self.data = data

    # запомнить первое значение яркостной температуры с целью дальнейшей калибровки по нему
    def set_T_t0_first(self):
        T_t0 = {}
        for key in sorted(self.data.keys()):
            for _, temp in self.data[key]:
                if temp == 0: continue
                else:
                    T_t0[key] = temp
                    break
        self.Tbr_t0 = T_t0

    # запомнить минимальные значения яркостных температур с целью дальнейшей калибровки по ним
    def set_T_t0_min(self):
        T_t0 = {}
        for key in sorted(self.data.keys()):
            min = sys.maxint
            for _, temp in self.data[key]:
                if temp == 0: continue
                if temp < min:
                    min = temp
            T_t0[key] = min
        self.Tbr_t0 = T_t0

    # запомнить средние значения яркостных температур на интервале безоблачной погоды с целью дальнейшей калибровки
    def set_T_t0_nocluds_interval(self, time_begin = 0, time_end = None):
        T_t0 = {}
        for key in sorted(self.data.keys()):
            avg_t, k = 0, 0
            for time, temp in self.data[key]:
                if (time_end and (time_begin <= time <= time_end)) or (not(time_end) and (time_begin <= time)):
                    avg_t += temp
                    k += 1
            if k != 0: avg_t /= k
            T_t0[key] = avg_t
        self.Tbr_t0 = T_t0
        self.noclouds_begin = time_begin
        self.noclouds_end = time_end

    # запомнить модельные яркостные температуры для данных метеоусловий с целью дальнейшей калибровки по ним
    def set_T_t0_model_Tb(self, theta):
        T_t0 = {}
        for key in sorted(self.data.keys()):
            T_t0[key] = self.MODEL.get_Tb(key, self.T_avg, theta)
        self.Tbr_t0 = T_t0

    # запомнить средние по полученным данным яркостные температуры с целью дальнейшей калибровки по ним
    def set_T_t0_Tb_avg(self):
        A = self.get_Tb_avg()
        self.Tbr_t0 = A

    # получить средние яркостные температуры за период
    def get_Tb_avg(self):
        A = {}
        for key in sorted(self.data.keys()):
            t_avg, k = 0, 0
            for _, temp in self.data[key]:
                if temp == 0: continue
                t_avg += temp
                k += 1
            t_avg /= k
            A[key] = t_avg
        return A

    # калибровка средних значений яркостных температур за период
    def calibr_Tb_avg(self, theta):
        if self.A == {}:
            self.A = self.get_Tb_avg()
            for key in self.A.keys():
                Tb_k = (self.T_avg + 273) * (1 - math.exp(-self.MODEL.tau_theory(key, theta)))
                # print(Tb_k)
                self.A[key] = Tb_k + (290 - Tb_k) / (290 - self.Tbr_t0[key]) * (self.A[key] - self.Tbr_t0[key]) #!!!
                # self.A[key] = Tb_k + (290 - Tb_k) / (290 - self.Tbr_t0[key]) * (self.A[key])
                # self.A[key] = Tb_k / self.Tbr_t0[key] * self.A[key]
        return

    def get_calibr_Tb_avg(self, theta):
        self.calibr_Tb_avg(theta)
        return self.A

    # двухчастотный метод определения Q и W
    def _get_QW(self, freq1, T_br1, freq2, T_br2, theta = 0.):
        M = np.array([[self.MODEL.krho(freq1, theta), self.MODEL.kw(freq1, self.T_obl)],
                      [self.MODEL.krho(freq2, theta), self.MODEL.kw(freq2, self.T_obl)]])
        v = np.array([self.MODEL.tau_experiment(T_br1, self.T_avg + 273, theta) - self.MODEL.tauO_theory(freq1, theta),
                      self.MODEL.tau_experiment(T_br2, self.T_avg + 273, theta) - self.MODEL.tauO_theory(freq2, theta)])
        s = np.linalg.solve(M, v)
        return s.tolist()[0], s.tolist()[1]

    # определение Q и W по откалиброванным средним значениям яркостных температур за период
    def get_QW4Tb_avg(self, freq1, freq2, theta = 0.):
        self.calibr_Tb_avg(theta)
        # A = self.get_Tb_avg()
        if not(freq1 in self.A.keys() and freq2 in self.A.keys()): return None
        a = self._get_QW(freq1, self.A[freq1], freq2, self.A[freq2], theta)
        return a

    # определение Q и W по откалиброванным средним значениям яркостных температур за период (оптимизационный метод)
    def get_opt_QW4Tb_avg(self, frequencies, theta):
        self.calibr_Tb_avg(theta)
        for freq in frequencies:
            if not(freq in self.A.keys()): return None
        m = self.MODEL
        minimum, sQ, sW = sys.maxint, 0, 0
        i = 0
        for Q in np.arange(0.2, 6, 0.01):
            if i % 1000 == 0: print("#", sep='..', end='')
            for W in np.arange(0.00, 1, 0.01):
                i += 1
                J = 0
                for freq in frequencies:
                    J += pow(
                        (m.tau_experiment(self.A[freq], self.T_avg + 273, theta) -
                            (m.tauO_theory(freq, theta) + m.krho(freq, theta) * Q + m.kw(freq, self.T_obl) * W)),
                        2)
                # print "Q = ", Q, "\tW = ", W, "\tJ = ", J
                if J < minimum:
                    minimum = J
                    sQ = Q
                    sW = W
        return sQ, sW

    # получение временного хода Q и W по ходу откалиброванных яркостных температур оптимизационным методом
    def get_opt_QW(self, frequencies, theta):
        self.calibr_Tb(theta)
        for freq in frequencies:
            if not(freq in self.data.keys()): return None
        m = self.MODEL
        min_len = sys.maxint
        for key in self.data.keys():
            if len(self.data[key]) < min_len: min_len = len(self.data[key])
        a = []
        infQ, supQ, stepQ, BQ = 0.8, 3.5, 0.01, 0.3
        infW, supW, stepW, BW = 0., 1.2, 0.001, 0.1
        print("Total: ", min_len)
        f = open("qw_opt.txt", "w")
        for i in range(min_len):
            print(i, " processed. ", end="")
            minimum, sQ, sW, st = sys.maxint, 0, 0, 0
            TB = {}
            for freq in frequencies:
                t, Tb = self.data[freq][i]
                st += t
                TB[freq] = Tb
            st /= len(frequencies)
            for Q in np.arange(infQ, supQ, stepQ):
                for W in np.arange(infW, supW, stepW):
                    J = 0
                    for freq in frequencies:
                        J += pow(
                            (m.tau_experiment(TB[freq], self.T_avg + 273, theta) -
                             (m.tauO_theory(freq, theta) + m.krho(freq, theta) * Q + m.kw(freq, self.T_obl) * W)),
                            2)
                    if J < minimum:
                        minimum = J
                        sQ = Q
                        sW = W
            a.append((st, sQ, sW))
            print("\tt = ", st, "\tQ = ", sQ, "\tW = ", sW)
            f.write(str(st) + "\t" + str(sQ) + "\t" + str(sW) + "\n")
            infQ = sQ - BQ
            supQ = sQ + BQ
            infW = sW - BW
            supW = sW + BW
        f.close()
        return a

    # калибровка яркостных температур (получение временного хода откалиброванных яркостных температур)
    def calibr_Tb(self, theta):
        if not self.isDataCalibred:
            for key in self.data.keys():
                for i in range(len(self.data[key])):
                    t, Tb = self.data[key][i]
                    if Tb == 0: continue
                    Tbk = (self.T_avg + 273) * (1 - math.exp(-self.MODEL.tau_theory(key, theta=theta)))
                    Tb = Tbk + (290 - Tbk) / (290 - self.Tbr_t0[key]) * (Tb - self.Tbr_t0[key]) #!!!
                    # Tb = Tbk + (290 - Tbk) / (290 - self.Tbr_t0[key]) * (Tb)
                    # Tb = Tbk / self.Tbr_t0[key] * Tb
                    self.data[key][i] = (t, Tb)
            self.isDataCalibred = True
        return

    # получение временного хода Q и W по ходу откалиброванных яркостных температур двухчастотным методом
    def get_QW(self, freq1, freq2, theta = 0.):
        self.calibr_Tb(theta)
        a = []
        l = len(self.data[freq1]) if len(self.data[freq1]) < len(self.data[freq2]) else len(self.data[freq2])
        for i in range(l):
            t1, Tb1 = self.data[freq1][i]
            t2, Tb2 = self.data[freq2][i]
            if Tb1 == 0 or Tb2 == 0: continue
            q, w = self._get_QW(freq1, Tb1, freq2, Tb2, theta)
            t = (t1 + t2) / 2
            a.append((t, q, w))
        return a

    # средние по временному ходу Q и W
    def get_QW_avg(self, freq1, freq2, theta = 0.):
        a = self.get_QW(freq1, freq2, theta)
        Q, W, k = 0, 0, 0
        for _, q, w in a:
            Q += q
            W += w
            k += 1
        Q /= k
        W /= k
        return Q, W

    # корректировка значений Q и W (временной ход)
    def get_W0corrected_QW(self, freq1, freq2, theta = 0.):
        b, e = self.noclouds_begin, self.noclouds_end
        if not b: return None
        a = self.get_QW(freq1, freq2, theta)
        w_avg, k = 0, 0
        for t, _, w in a:
            if (e and b <= t <= e) or (not(e) and b <= t):
                w_avg += w
                k += 1
        w_avg /= k
        for i in range(len(a)):
            t, q, w = a[i]
            w -= w_avg
            a[i] = (t, q, w)
        return a

    # resolve-корректировка значений Q и W (временной ход)
    def get_resolved_QW(self, freq1, freq2, theta = 0.):
        b, e = self.noclouds_begin, self.noclouds_end
        if not b: return None
        a = self.get_QW(freq1, freq2, theta)
        w_avg, k = 0, 0
        for t, _, w in a:
            if (e and b <= t <= e) or (not (e) and b <= t):
                w_avg += w
                k += 1
        w_avg /= k
        a = []
        l = len(self.data[freq1]) if len(self.data[freq1]) < len(self.data[freq2]) else len(self.data[freq2])
        for i in range(l):
            t1, T_br1 = self.data[freq1][i]
            t2, T_br2 = self.data[freq2][i]
            if T_br1 == 0 or T_br2 == 0: continue

            krho_f1, krho_f2 = self.MODEL.krho(freq1, theta), self.MODEL.krho(freq2, theta)
            kw_f1, kw_f2 = self.MODEL.kw(freq1, self.T_obl), self.MODEL.kw(freq2, self.T_obl)
            tau_ex_t1, tau_ex_t2 = self.MODEL.tau_experiment(T_br1, self.T_avg + 273, theta), \
                                   self.MODEL.tau_experiment(T_br2, self.T_avg + 273, theta)
            tauO_th_f1, tauO_th_f2 = self.MODEL.tauO_theory(freq1, theta), self.MODEL.tauO_theory(freq2, theta)

            M = np.array([[krho_f1, kw_f1],
                          [krho_f2, kw_f2]])
            v = np.array(
                [tau_ex_t1 - tauO_th_f1 - kw_f1 * w_avg,
                 tau_ex_t2 - tauO_th_f2 - kw_f2 * w_avg])
            s = np.linalg.solve(M, v)

            t = (t1 + t2) / 2
            a.append((t, s.tolist()[0], s.tolist()[1]))
        return a

    def generate_dataset(self, fname, freq1, freq2, theta = 0.):
        b, e = self.noclouds_begin, self.noclouds_end
        if not b: return None
        a = self.get_QW(freq1, freq2, theta)
        w_avg, k = 0, 0
        for t, _, w in a:
            if (e and b <= t <= e) or (not (e) and b <= t):
                w_avg += w
                k += 1
        w_avg /= k
        a = []
        fout = open(fname, "w")
        l = len(self.data[freq1]) if len(self.data[freq1]) < len(self.data[freq2]) else len(self.data[freq2])
        for i in range(l):
            t1, T_br1 = self.data[freq1][i]
            t2, T_br2 = self.data[freq2][i]
            if T_br1 == 0 or T_br2 == 0: continue

            krho_f1, krho_f2 = self.MODEL.krho(freq1, theta), self.MODEL.krho(freq2, theta)
            kw_f1, kw_f2 = self.MODEL.kw(freq1, self.T_obl), self.MODEL.kw(freq2, self.T_obl)
            tau_ex_t1, tau_ex_t2 = self.MODEL.tau_experiment(T_br1, self.T_avg + 273, theta), \
                                   self.MODEL.tau_experiment(T_br2, self.T_avg + 273, theta)
            tauO_th_f1, tauO_th_f2 = self.MODEL.tauO_theory(freq1, theta), self.MODEL.tauO_theory(freq2, theta)

            M = np.array([[krho_f1, kw_f1],
                          [krho_f2, kw_f2]])
            v = np.array(
                [tau_ex_t1 - tauO_th_f1 - kw_f1 * w_avg,
                 tau_ex_t2 - tauO_th_f2 - kw_f2 * w_avg])
            s = np.linalg.solve(M, v)

            t = (t1 + t2) / 2
            Q, W = s.tolist()[0], s.tolist()[1]
            a.append((t, Q, W))
            # fout.write(str(self.MODEL.T) + " " + str(self.MODEL.P) + " " + str(self.MODEL.rho) + " " +
            #            str(freq1) + " " + str(freq2) + " " + str(T_br1) + " " + str(T_br2) + " " +
            #            str(Q) + " " + str(W) + "\n")
            fout.write(str(t) + " " + str(T_br1) + " " + str(T_br2) + " " + str(Q) + " " + str(W) + "\n")
        fout.close()
        return a

    def get_QW_errors(self, freq1, freq2, theta = 0., Tbr_error = 1., Tavg_error = 2.):
        self.calibr_Tb(theta)
        a = []
        l = len(self.data[freq1]) if len(self.data[freq1]) < len(self.data[freq2]) else len(self.data[freq2])
        for i in range(l):
            t1, T_br1 = self.data[freq1][i]
            t2, T_br2 = self.data[freq2][i]
            if T_br1 == 0 or T_br2 == 0: continue

            krho_f1, krho_f2 = self.MODEL.krho(freq1, theta), self.MODEL.krho(freq2, theta)
            kw_f1, kw_f2 = self.MODEL.kw(freq1, self.T_obl), self.MODEL.kw(freq2, self.T_obl)
            tau_ex_t1, tau_ex_t2 = self.MODEL.tau_experiment(T_br1, self.T_avg + 273, theta), \
                                   self.MODEL.tau_experiment(T_br2, self.T_avg + 273, theta)
            tauO_th_f1, tauO_th_f2 = self.MODEL.tauO_theory(freq1, theta), self.MODEL.tauO_theory(freq2, theta)

            d_tau_t1 = self.MODEL.tau_experiment(T_br1 + Tbr_error, self.T_avg + 273 + Tavg_error, theta) - \
                                 tau_ex_t1
            d_tau_t2 = self.MODEL.tau_experiment(T_br2 + Tbr_error, self.T_avg + 273 + Tavg_error, theta) - \
                                 tau_ex_t2
            delta_tau_f1 = self.MODEL.tau_theory(freq1, theta) - tauO_th_f1
            delta_tau_f2 = self.MODEL.tau_theory(freq2, theta) - tauO_th_f2

            q, w = self._get_QW(freq1, T_br1, freq2, T_br2, theta)

            dQ = math.sqrt(kw_f2**2 * d_tau_t1**2 + kw_f1**2 * d_tau_t2**2) / \
                 (delta_tau_f1 * kw_f2**2 - delta_tau_f2 * kw_f1**2) * \
                 q

            dW = math.sqrt(krho_f2**2 * d_tau_t1**2 + krho_f1**2 * d_tau_t2**2) / \
                 (delta_tau_f1 * krho_f2**2 - delta_tau_f2 * krho_f1**2) * \
                 w

            t = (t1 + t2) / 2
            a.append((t, dQ, dW))
        return a

    def get_QW_avg_resolved_noclouds(self, freq1, freq2, theta):
        a = self.get_resolved_QW(freq1, freq2, theta)
        time_begin = self.noclouds_begin
        time_end = self.noclouds_end
        avgQ, avgW, k = 0, 0, 0
        for t, q, w in a:
            if (time_end and (time_begin <= t <= time_end)) or (not(time_end) and (time_begin <= t)):
                avgQ += q
                avgW += w
                k += 1
        avgQ /= k
        avgW /= k
        return avgQ, avgW

    def get_QW_err_avg_noclouds(self, freq1, freq2, theta, Tb_err, Tavg_err):
        a = self.get_QW_errors(freq1, freq2, theta, Tb_err, Tavg_err)
        time_begin = self.noclouds_begin
        time_end = self.noclouds_end
        avgQ, avgW, k = 0, 0, 0
        for t, q, w in a:
            if (time_end and (time_begin <= t <= time_end)) or (not(time_end) and (time_begin <= t)):
                avgQ += q
                avgW += w
                k += 1
        avgQ /= k
        avgW /= k
        return avgQ, avgW

    def get_QW_err_avg_clouds(self, freq1, freq2, theta, Tb_err, Tavg_err):
        a = self.get_QW_errors(freq1, freq2, theta, Tb_err, Tavg_err)
        time_begin = self.noclouds_begin
        time_end = self.noclouds_end
        avgQ, avgW, k = 0, 0, 0
        for t, q, w in a:
            if not((time_end and (time_begin <= t <= time_end)) or (not(time_end) and (time_begin <= t))):
                avgQ += q
                avgW += w
                k += 1
        avgQ /= k
        avgW /= k
        return avgQ, avgW