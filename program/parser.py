# -*- coding: utf-8 -*-
from __future__ import division
from __future__ import print_function
import re
import os
import glob
from switch import Switch
from collections import defaultdict
from converter import Converter
import numpy as np


class Parser:
    def __init__(self, mode, path):
        self.mode = mode
        self.path = path

    @staticmethod
    def parse_line(text):
        l_data = re.split("[\t ]", re.sub("[\r\n]", '', text))
        l_data = [elem for elem in l_data if elem]
        return l_data

    def __file_not_empty(self, filePath):
        return os.path.getsize(filePath) > 0

    def parse(self):
        DATA = defaultdict(list)
        with Switch(self.mode) as case:
            if case("tb"):
                print(".parse() ~tb(file)")
                t_file = open(self.path, "r")
                t_file.readline()
                for line in t_file:
                    l_data = Parser.parse_line(line)
                    hh, mm, ss = l_data[1].split(':')
                    ms = l_data[2]
                    timestamp = Converter.timeToX(hh, mm, ss, ms)
                    DATA[int(l_data[5])].append((timestamp, float(l_data[6])))
                    DATA[int(l_data[7])].append((timestamp, float(l_data[8])))
                t_file.close()
                print("done!")
            if case("tb_timeP"):
                print(".parse() ~tb(file) --- TIME_P")
                t_file = open(self.path, "r")
                t_file.readline()
                for line in t_file:
                    l_data = Parser.parse_line(line)
                    timestamp = float(l_data[3])
                    DATA[float(l_data[5])/1000].append((timestamp, float(l_data[6])))
                    DATA[float(l_data[7])/1000].append((timestamp, float(l_data[8])))
                t_file.close()
                print("done!")
            if case("tb_g_1"):
                print(".parse() ~tb(file) - гетеродин")
                t_file = open(self.path, "r")
                t_file.readline()
                for line in t_file:
                    l_data = Parser.parse_line(line)
                    hh, mm, ss = l_data[1].split(':')
                    ms = l_data[2]
                    timestamp = Converter.timeToX(hh, mm, ss, ms)
                    #print(l_data[1], " ", l_data[2], " ", l_data[4], " ", l_data[6])
                    DATA[int(l_data[4])].append((timestamp, float(l_data[6])))
                t_file.close()
                print("done!")
            if case("tb_g_2"):
                print(".parse() ~tb(file) - гетеродин")
                t_file = open(self.path, "r")
                t_file.readline()
                for line in t_file:
                    l_data = Parser.parse_line(line)
                    hh, mm, ss = l_data[1].split(':')
                    ms = l_data[2]
                    timestamp = Converter.timeToX(hh, mm, ss, ms)
                    #print(l_data[1], " ", l_data[2], " ", l_data[4], " ", l_data[6])
                    DATA[int(l_data[4])].append((timestamp, float(l_data[8])))
                t_file.close()
                print("done!")
            if case("th"):
                print(".parse() ~th(dir)")
                for d_filename in glob.glob(self.path + "*.dat"):
                    print(d_filename)
                    d_file = open(d_filename, "r")
                    d_file.readline()
                    for line in d_file:
                        l_data = Parser.parse_line(line.replace('-', ''))
                        hh, mm, ss, ms = l_data[2:6]
                        timestamp = Converter.timeToX(hh, mm, ss, ms)
                        DATA["v1"].append((timestamp, float(l_data[9])))
                        # DATA["v2"].append((timestamp, float(l_data[10])))
                        DATA["az"].append((timestamp, float(l_data[11])))
                    d_file.close()
                for key in DATA.keys():
                    DATA[key] = sorted(DATA[key], key=lambda (time, value): time)
                print("done!")
        return DATA

    def shift_tb(self, DATA):
        ssize = len(DATA[18])
        msize = len(DATA[21.2])
        s_arr = [21.2, 21.4, 21.6, 21.8, 22.0, 22.2, 22.4, 22.6, 22.8, 23.0, 23.2, 23.4, 23.6, 23.8, 24.0]
        for f in s_arr:
            if not DATA[f]:
                print(f)
                continue
            l = len(DATA[f])
            if l > 1.25 * ssize:
                k = 0
                av = []
                for i in range(0, l-1, 2):
                    t_i, x_i = DATA[f][i]
                    t_i1, x_i1 = DATA[f][i+1]
                    t = (t_i + t_i1) / 2
                    x = (x_i + x_i1) / 2
                    av.append((t, x))
                    k += 1
                DATA[f] = av
                if k < msize: msize = k
        return

    def prepare(self, DATA):
        if 0 in DATA.keys(): del DATA[0]
        return