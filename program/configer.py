#coding=utf-8
from __future__ import division
from __future__ import print_function
import re as regexp
import math


class Configer:
    def __init__(self, path):
        self.path = path

        # стандартные условия
        self.theta = (math.pi / 180) * 51  # rad
        self.T = 15                        # cels
        self.P = 1013                      # hPa
        self.rho = 7.5                     # g/m^3
        self.Tavg = 15                     # cels
        self.Tobl = -2                     # cels
        self.start = 0                     # sec
        self.stop = None                   # sec
        self.nclbeg = 0                    # sec
        self.nclend = None                 # sec
        self.Vwind = None                  # km/h
        self.interv = None                 # km

        # Пример config-файла:
        # ==========================================
        # Theta   51        # deg
        # T       17.3      # cels
        # P       748       # мм.рт.ст
        # Rho     13        # г/м^3
        # Tavg    2.3       # cels
        # Tobl    -2        # cels
        # start   49600     # sec
        # stop    -1        # sec
        # nclbeg  49600     # sec
        # nclend  50600     # sec
        # Vwind   -1        # km/h
        # interv  -1        # km
        # ==========================================

    def get_info(self):
        file = open(self.path)
        for line in file:
            # print(line)
            d = regexp.split("[\t ]", regexp.sub("[\r\n,]", '', line))
            if len(d) != 2: continue
            x, val = str(d[0]), float(d[1])
            if x in ["Theta", "theta"]:       self.theta = (math.pi / 180) * val
            if x in ["T", "Temperature"]:     self.T = val
            if x in ["P", "Pressure"]:        self.P = 1.33322 * val
            if x in ["Rho", "rho"]:           self.rho = val
            if x in ["Tavg", "T_avg"]:        self.Tavg = val
            if x in ["Tobl", "T_obl"]:        self.Tobl = val
            if x in ["start", "Start"]:       self.start = val
            if x in ["stop", "Stop"]:
                if val == -1:                 self.stop = None
                else:                         self.stop = val
            if x in ["nclbeg", "nclBeg"]: self.nclbeg = val
            if x in ["nclend", "nclEnd"]:
                if val == -1:                 self.nclend = None
                else:                         self.nclend = val
            if x in ["V", "Vwind", "vwind", "v_wind", "V_wind"]:
                if val == -1:                 self.Vwind = None
                else:                         self.Vwind = val
            if x in ["interv", "int", "interval"]:
                if val == -1:                 self.interv = None
                else:                         self.interv = val

        file.close()
        return
