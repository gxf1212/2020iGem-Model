# try to get analytical solution
# impossible
import numpy as np
import sympy as sy
import sys
sys.path.append(r'./main')
from utils import *

# %% preparation
# parameters
# 增长率，decay rate
r1 = np.log(2) / 30
r2 = np.log(2) / 20
m1 = 0.001
m2 = 0.001

# 半饱和值
Kc1 = 1
Kn1 = 9.31
Kp1 = 1.30  # 没确定
Kn2 = 9.31
Kp2 = 1.30

# 转化速率，暂时认为和时间无关
# 芽孢呼吸，maintenance energy coefficient of 0.67 mmol glucose/g CDW/h
Re1 = 0.67e-3 * 12 / 60
# 蓝藻净光合，4.41+-0.16μmol(mg·min)，1000/54为光强
NPh = 4.41e-3 * 12.01 * np.tanh((1000 / 54) * 0.016 / 4.41)
# 蓝藻固氮,1.72+-0.25 fmol N/cell/h，细胞质量是估计的
Nf = 1.72e-15 * 14.01 / 60 / (1e-9)
# 芽孢溶磷，无机+有机
Pf = 8.87e-4 + 0

# Q值，假设为常数
# 芽孢
Q1c = 0.303  # 待定，能用
Q1n = 0.0734
Q1p = 0.00607
# 蓝藻
Q2c = 0.2523
Q2n = 0.0846
Q2p = 0.00248

## 模拟
# 迭代变量，设定为初值
n1 = 1
n2 = 1
# 沙丘
rc = 7.2e-3
rn = 2.6e-4
rp = 8.88e-5
# 植被
# rc = 2.35e-2
# rn = 2.8e-4
# rp = 3.676e-4

# variables
t = sy.symbols("t")
N1, N2, Rc, Rn, Rp = sy.symbols("N1 N2 Rc Rn Rp")
G1, G2, D1, D2 = sy.symbols("G1 G2 D1 D2")
f1, f2 = sy.symbols("f1 f2")


# %% define functions
def equations():
    return [G1 == f1 * r1 * N1,
            D1 == m1 * N1,
            G2 == f2 * r2 * N2,
            D2 == m2 * N2,
            f1 == np.min([MM(Rc, Kc1, N1), MM(Rn, Kn1, N1), MM(Rp, Kp1, N1)]),
            f2 == np.min([MM(Rn, Kn2, N2), MM(Rp, Kp2, N2)]),
            sy.diff(N1) == G1 - D1,
            sy.diff(N2) == G2 - D2,
            sy.diff(Rc) == -Re1 * N1 + NPh * N2 - Q1c * G1 - Q2c * sy.diff(N2),
            sy.diff(Rn) == Nf * N2 - Q1n * G1 - Q2n * sy.diff(N2),
            sy.diff(Rp) == Pf * N1 - Q1p * G1 - Q2p * sy.diff(N2)]


# %%
sy.dsolve(equations())

# plt.plot(x, y[:, 0])  # y数组（矩阵）的第一列，（因为维度相同，plt.plot(x, y)效果相同）
# plt.grid()
# plt.show()
