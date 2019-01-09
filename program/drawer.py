from matplotlib import pyplot as plt


class Drawer:
    def __init__(self):
        self.color = {
            18 : "navy",
            19.2 : "dodgerblue",
            20.4 : "darkgreen",
            21.6 : "limegreen",
            22.2 : "darkorange",
            22.4 : "orange",
            23.2 : "m",
            24.4 : "darkviolet",
            25.6 : "brown",
            26.8 : "red"
        }

    def draw(self, data, frequencies = None, title = "", xlabel = "", ylabel = ""):
        plt.title(title)
        for key in sorted(data.keys()):
            if frequencies:
                if not(key in frequencies): continue
            TIME, TEMP = [], []
            for time, temp in data[key]:
                TIME.append(time)
                TEMP.append(temp)
            plt.plot(TIME, TEMP, label=str(key)+" GHz", color=self.color[key])
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.legend(loc='right')
        plt.show()