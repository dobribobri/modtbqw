from collections import defaultdict
import re as regexp


class Parser:
    def __init__(self, path):
        self.path = path
        self.data = defaultdict(list)

    def parse(self):
        t = []
        file = open(self.path)
        for k, line in enumerate(file):
            d = regexp.split("[\t ]", regexp.sub("[\r\n]", '', line))
            d = [e for e in d if e]
            if not k:
                for e in d[1:]:
                    t.append(float(e))
                continue
            for i in range(len(d) - 1):
                self.data[float(d[0])].append((t[i], float(d[i + 1])))

    def d_print(self):
        for key in sorted(self.data.keys()):
            print(key)
            print(self.data[key])
            print(' ')