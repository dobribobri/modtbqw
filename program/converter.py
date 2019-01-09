# -*- coding: utf-8 -*-
from __future__ import division
from __future__ import print_function


class Converter:

    @staticmethod
    def timeToX(hh, mm, ss, ms):
        return float(ms) + float(ss) * 1000 + \
               float(mm) * 60 * 1000 + float(hh) * 3600 * 1000