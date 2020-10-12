# try to obtain numerical solution

import sys

import pylab as pl
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, subplot, figure, legend, xlabel, ylabel, title

sys.path.append(r'E:\undergraduate_study\study\research&competition\igem及合成生物学\土壤改造\建模\2020iGem-Model\growth simulation\main')
from utils import *
# the project is defined in the folder 'growth stimulation',
# we have to append the path 'main', and just ignore the error
# we have to use date type 'float64'


# %% preparation
# parameters
# 增长率，decay rate
r1 = np.log(2) / 30
r2 = np.log(2) / 20
# m1=0.196/60 # 14.812 min^{-1}
m1 = 0.001
m2 = 0.001

# 分泌毒素的抑制，gamma，单位和N^{-2}相同
ga1 = 0.1  # 2对1的限制
ga2 = 0.1

# 半饱和值，无量纲
Kc1 = 2.865 * 0.4
Kn1 = 0.638
Kp1 = 0.057
Kn2 = 0.050993  # 蓝丝菌
Kp2 = 0.018739  # 铜绿的

# Q值，假设为常数
# 芽孢
# Q1c = 0.303 # 估算的质量分数
Q1c = 0.34 * 8 / 0.3 * 1e-3 * 72  # reciprocal of growth yield
Q1n = 0.0734
Q1p = 0.00607
# 蓝藻
# Q2c=0.2523 # 估算的质量分数
Q2c = 0.488 / 0.65
Q2n = 0.0846
Q2p = 0.00248
# Q2p = Q2p*20;  # consider P luxury uptake by Nostoc
Q1n=Q1n*20

# 念珠藻两种细胞的比例，为了更真实
f_Nf = 0.1
f_Ph = 1 - f_Nf

# 转化速率，暂时认为和时间无关
# 芽孢呼吸，maintenance energy coefficient of 0.67 mmol glucose/g CDW/h
Re1 = 0.67e-3 * 12 / 60
# 蓝藻（用了铜绿）净光合，4.41+-0.16μmol(mg·min)，1000/54为光强
NPh = 4.41e-3 * 12.01 * np.tanh((1000 / 54) * 0.016 / 4.41)
# NPh=NPh*f_Ph
# bacterial N fixation
Nf1 = 1e-4
# 蓝藻固氮,1.72+-0.25 fmol N/cell/h，细胞质量是估计的
Nf2 = 1.72e-15 * 14.01 / 60 * (7722e6 / 0.42)
# Nf=Nf*f_Nf
Nf2 = Nf2 * 7
# 芽孢溶磷，无机+有机
Pf = 8.87e-4 / 1440 + 0
Pf = Pf * 10

# pH值
pH0 = 9.1

# Variables
# biomass: 1 for 芽孢，2 for 蓝藻
# c,n,p表示各自元素

# 迭代变量，设定为初值
# 初始菌浓度，通常藻比菌多
n1 = 0.01
n2 = 0.2
# 初始营养浓度，单位g/L

# 土壤容重，自然状态的体积，单位g/L
# 沙丘
# rc=7.2e-3
# rn=2.6e-4
# rp=8.88e-7
# rho=1.426e3
# 植被
rc = 2.35e-2
rn = 2.8e-4
rp = 3.676e-6
rho = 1.446e3

rc = rc / 1.724  # 源数据为有机质，非碳
rn = rn * 0.05; # 氮要用可利用氮（无机，而部分有机氮可转化）
# rp = rp / 0.06  # 磷别用总磷

# 自己添加资源
# rc=0.02
# rn=0.001
# rp=0.001
# rho=1.33e3

rc = rc * rho
rn = rn * rho
rp = rp * rho  # 此时单位为g/L


#%%useful functions
def MM(r, k, n):
    return r / (r + k * n)


def toxin(ga, n1, n2):
    return 1 - ga * n1 * n2


def numerical_simulation(n1, n2, rc, rn, rp, t_max):
    # 变量
    # 生物量：1表示芽孢，2表示蓝藻
    N1 = []
    N2 = []
    # c, n, p表示各自元素
    Rc = []
    Rn = []
    Rp = []
    # 记录生长速率等，作为监控
    G1 = []
    G2 = []

    # 输入参数：几个变量的初值

    step = 0.01
    time = [(t+1)*step for t in range(int(t_max/step))]
    for t in time:
        # 更新的顺序可能无所谓，因为step取得很小
        # 更新资源限制因素
        f1 = min([MM(rc, Kc1, n1), MM(rn, Kn1, n1), MM(rp, Kp1, n1)])
        # f1 = min([MM_ori(rc, Kc1), MM_ori(rn, Kn1), MM_ori(rp, Kp1)]) # 不考虑N的资源限制
        f2 = min([MM(rn, Kn2, n2), MM(rp, Kp2, n2)])

        # 更新资源量，并定义Grow
        Grow1 = f1 * r1 * n1
        Grow2 = f2 * r2 * n2
        # considering toxin
        # Grow1 = f1 * r1 * n1 * toxin(ga1, n1, n2)
        # Grow2 = f2 * r2 * n2 * toxin(ga2, n1, n2)
        # resources
        rc = rc + (-Re1 * n1 + NPh * n2 - Q1c * Grow1 - Q2c * (Grow2 - m2 * n2)) * step
        rn = rn + (Nf1 * n1 + Nf2 * n2 - Q1n * Grow1 - Q2n * (Grow2 - m2 * n2)) * step
        rp = rp + (Pf * n1 - Q1p * Grow1 - Q2p * (Grow2 - m2 * n2)) * step
        # 模拟芽孢物质也归还
        # rc = rc + (-Re1 * n1 + NPh * n2 - Q1c * (Grow1 - m1 * n1) - Q2c * (Grow2 - m2 * n2)) * step
        # rn = rn + (Nf1 * n1 + Nf2 * n2 - Q1n * (Grow1 - m1 * n1) - Q2n * (Grow2 - m2 * n2)) * step
        # rp = rp + (Pf * n1 - Q1p * (Grow1 - m1 * n1) - Q2p * (Grow2 - m2 * n2)) * step

        # 更新生物量
        n1 = n1 + (Grow1 - m1 * n1) * step
        n2 = n2 + (Grow2 - m2 * n2) * step

        # 储存需要的变量
        N1.append(n1)
        N2.append(n2)
        Rc.append(rc)
        Rn.append(rn)
        Rp.append(rp)
        # 监控变量
        G1.append(Grow1)
        G2.append(Grow2)

    return np.array(N1), np.array(N2), np.array(Rc), np.array(Rn), \
           np.array(Rp), np.array(time), np.array(G1), np.array(G2)


def plot_result(see_toxin=False):
    # growth comparation between the two
    fig1 = figure(1, figsize=(8, 6))
    plt.clf()
    # gs = gridspec.GridSpec(1,2,width_ratios=[2,1])
    # ax0=plt.subplot(gs[0])
    plt_growth = plot(time, N1, 'b', time, N2, 'r')
    title('growth curve')
    xlabel('time/min')
    ylabel('biomass concentration/($g\cdot L^{-1}$)')
    legend(handles=plt_growth, labels=['B.S', 'Nostoc'], loc='best')

    # ax1 = plt.subplot(gs[1])
    # gs1 = gridspec.GridSpec(2, 1, width_ratios=[1, 1])
    # ax2 = plt.subplot(gs1[0])
    # ax3 = plt.subplot(gs1[1])

    # growth and decay
    fig2 = figure(2, figsize=(10, 4))
    plt.clf()
    pl.gca().yaxis.get_major_formatter().set_powerlimits((0, 1))  # 将坐标轴的base number设置为一位。1是指科学计数法时的位数
    subplot(1, 2, 1)
    plot(time, G1, time, N1 * m1)
    title('B.S')
    legend(['growth', 'decay'])

    subplot(1, 2, 2)
    title('Nostoc')
    plot(time, G2, time, N2 * m2)
    legend(['growth', 'decay'])
    plt.show()


    # % 资源
    fig3 = figure(3, figsize=(14, 10.5))
    plt.clf()
    # plot(time, Rc, time, Rn, time, Rp)
    # legend(['C','N','P'])
    subplot(2, 1, 1)
    plot(time, Rc)
    xlabel('time/min')
    ylabel('carbon concentration/($g\cdot L^{-1}$)')
    legend(labels='C')

    subplot(2, 2, 3)
    plot(time, Rn)
    legend(labels='N')
    xlabel('time/min')
    ylabel('nitrogen concentration/($g\cdot L^{-1}$)')

    subplot(2, 2, 4)
    plot(time, Rp)
    legend(labels='P')
    xlabel('time/min')
    ylabel('phosphorus concentration/($g\cdot L^{-1}$)')

    plt.show()


    # limiting factor analysis
    figure(4, figsize=(12, 9))
    plt.clf()
    subplot(2,2,1)
    title('f of B.S')
    # ha1=scatter(time, G1/N1/r1/toxin(ga1, N1, N2),color='k',marker='o',s=5)
    ha1=plot(time, G1/N1/r1,linewidth=2.5)
    ha2=plot(time, MM(Rc, Kc1, N1), time, MM(Rn, Kn1, N1), time, MM(Rp, Kp1, N1), linewidth=1)
    legend(handles=[ha1,ha2[0],ha2[1],ha2[2]],labels=['f', 'C', 'N', 'P'])

    subplot(2,2,2)
    title('f of Nostoc')
    # ha1=scatter(time, G2/N2/r2/toxin(ga2, N1, N2), color='k',marker='o',s=5)
    ha1=plot(time, G2/N2/r2, linewidth=2.5)
    # # ha2=plot(time, MM(Rn, Kn2, N2), time, MM(Rp, Kp2, N2))
    # # legend(handles=[ha1,ha2[0],ha2[1]],labels=['f','N','P'])

    if see_toxin:
        subplot(2, 2, 3)
        title('toxin factor')
        plot(time, toxin(ga1, N1, N2), time, toxin(ga2, N1, N2))
        legend(['B.S', 'Nostoc'])



# fig2：重复运行图会变色和图例不同
# fig1：不能控大小
# fig3：控制坐标轴科学计数法


#%% 调用模拟函数
days = 4
t_max = 1440 * days  # 1天是1440min
N1, N2, Rc, Rn, Rp, time, G1, G2 = numerical_simulation(n1, n2, rc, rn, rp, t_max)

# store variables to visualize in MATLAB
# save_files(N1, N2, Rc, Rn, Rp, time, G1, G2)
# %% 展示结果
# plot_settings()
plot_result()
# plot_result(see_toxin=True)

# Python用来比较正式的模拟，MATLAB做定性、单点的

#%% 调参









