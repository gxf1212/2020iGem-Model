import xlrd
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt


def read_data(filename):
    sheet = xlrd.open_workbook(filename).sheet_by_index(0)
    x = sheet.col_values(0)
    y = sheet.col_values(1)
    return np.array(x[2:]), np.array(y[2:])


def fit_mm(x, y):
    # x1 = (1/x).reshape(-1, 1)
    # y1 = (1/y).reshape(-1, 1)
    x1 = 1 / x
    y1 = 1 / y
    slope, intercept, r_value, p_value, std_err = linregress(x1, y1)
    print('coefficient of determination(𝑅²) :', r_value**2)
    plt.figure()
    plt.scatter(x1, y1)
    t = 1
    t = np.linspace(0, t, t*10+1)
    plt.plot(t, slope * t + intercept)

    V = 1 / intercept
    K = V * slope
    return V, K


#%% main
file = "west.xls"
# file = "central.xls"
# file = "east.xls"
y, x = read_data(file)
v, k = fit_mm(x, y)

#%% plot
t = np.linspace(0,40,401)
plt.figure()
plt.title(file)
plt.scatter(x, y)
plt.plot(t, v*t/(t+k), '--', linewidth=1)
plt.show()

