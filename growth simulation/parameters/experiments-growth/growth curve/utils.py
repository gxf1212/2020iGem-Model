import re
import matplotlib.pyplot as plt

def t2s(t):
    # d,h,m,s = re.split(':| ',t.strip())  # file modifying time
    y, mo, d, h, m, s = re.split(':| |/', t.strip())
    return int(d) * 86400 + int(h) * 3600 + int(m) * 60 + int(s)


def plot1(value, time):
    for i in range(value.shape[1]):
        plt.plot(time, value[:, i])
        plt.legend(["1", "2", "3", "4", "5", "6", "7", "8"])
        plt.show()


