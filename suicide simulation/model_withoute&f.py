import scipy
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def arabinoseInsideFunction(A_in, V_e, A_ex, E_mRNA, K_e, V_fgh, FGH_mRNA, K_fgh, Mu, Gamma_ain):
    return V_e * A_ex * E_mRNA / (K_e + A_ex) + V_fgh * A_ex * FGH_mRNA / (K_fgh + A_ex) - (Mu + Gamma_ain) * A_in


def araEmRNADynamicsFunction(E_mRNA, Alpha_e, V_me, A_in, K_me, Mu, Gamma_e):
    return Alpha_e + V_me * (A_in ** 3) / (K_me ** 3 + A_in ** 3) - (Mu + Gamma_e) * E_mRNA


def araFGHmRNADynamicsFunction(FGH_mRNA, Alpha_fgh, V_mfgh, A_in, K_mfgh, Mu, Gamma_fgh):
    return Alpha_fgh + V_mfgh * (A_in ** 3) / (K_mfgh ** 3 + A_in ** 3) - (Mu + Gamma_fgh) * FGH_mRNA


def mRNAFunction(M, D, k_m, A_in, K_1, K, d_m, d_large, MazF, K_t):
    return D * k_m * (1 + K_1 * (A_in)) / (K + K_1 * (A_in)) - d_m * M

def mazfFunction(b_2, M, MazF, d_c):
    return b_2 * M - d_c * MazF


def mazeFunction(b_1, M, a_T, MazF, MazE, d_a):
    return b_1 * M - a_T * MazF * MazE - d_a * MazE


def mazefFunction(a_T, MazE, MazF, d_c, MazEF, d_a2):
    return a_T * MazE * MazF - d_c * MazEF - d_a2 * MazEF


def totalModel(y, t, parameters):
    A_in, E_mRNA, FGH_mRNA, M, MazF = y
    V_e, A_ex, K_e, V_fgh, K_fgh, Mu, Gamma_ain, \
    Alpha_e, V_me, K_me, Gamma_e, \
    Alpha_fgh, V_mfgh, K_mfgh, Gamma_fgh, \
    D, k_m, K_1, K, d_m, d_large, K_t, \
    b_2, a_T, d_c, d_a2, b_1, d_a \
        = parameters
    _ain = arabinoseInsideFunction(A_in, V_e, A_ex, E_mRNA, K_e, V_fgh, FGH_mRNA, K_fgh, Mu, Gamma_ain)
    _emRNA = araEmRNADynamicsFunction(E_mRNA, Alpha_e, V_me, A_in, K_me, Mu, Gamma_e)
    _fghmRNA = araFGHmRNADynamicsFunction(FGH_mRNA, Alpha_fgh, V_mfgh, A_in, K_mfgh, Mu, Gamma_fgh)
    _mRNA = mRNAFunction(M, D, k_m, A_in, K_1, K, d_m, d_large, MazF, K_t)
    # _maze = mazeFunction(b_1, M, a_T, MazF, MazE, d_a)
    _mazf = mazfFunction(b_2, M, MazF, d_c)
    # _mazef = mazefFunction(a_T, MazE, MazF, d_c, MazEF, d_a2)
    return np.array([_ain, _emRNA, _fghmRNA, _mRNA, _mazf])








if '__main__':
    V_e = 48
    A_ex = 1E-7

    K_e = 1.4E-4
    V_fgh = 4.8
    K_fgh = 6E-6
    Mu = 5.78E-4
    Gamma_ain = 0.0155
    Alpha_e = 2.83E-11
    V_me = 4.28E-9
    K_me = 1.06E-4
    Gamma_e = 0.011552
    Alpha_fgh = 7.55E-11
    V_mfgh = 1.1266E-8
    K_mfgh = 1E-9
    Gamma_fgh = 0.0115
    D = 2
    k_m = 0.03  # it is 0.18/min=0.03/s
    K_1 = 2.52E10
    K = 7200
    d_m = 0.002
    d_large = 0.198
    K_t = 15
    b_2 = 0.009
    a_T = 0.0365
    d_c = 0.00028
    d_a2 = 0.000231
    b_1 = 0.122
    d_a = 0.00231
    A_in = 0
    E_mRNA = 0
    FGH_mRNA = 0
    M = 0.00416634
    MazF = 0.1335502
    parameters = [
        V_e, A_ex, K_e, V_fgh, K_fgh, Mu, Gamma_ain,
        Alpha_e, V_me, K_me, Gamma_e,
        Alpha_fgh, V_mfgh, K_mfgh, Gamma_fgh,
        D, k_m, K_1, K, d_m, d_large, K_t,
        b_2, a_T, d_c, d_a2, b_1, d_a]
    time = np.arange(0, 5000, 1)


    solve = odeint(totalModel, [A_in, E_mRNA, FGH_mRNA, M, MazF], time, (parameters,), atol=1e-7, rtol=1e-11, full_output=1)
    A_in, E_mRNA, FGH_mRNA, M, MazF = solve[0][:, 0], solve[0][:, 1], solve[0][:, 2], solve[0][:, 3], \
                                                   solve[0][:, 4]
    # M_1 = M.reshape(3000,1)
    print(M)
    '''
    print(time)
    print(A_in)
    A_in_1=A_in.reshape(1000,1)
    print(A_in_1)
    '''
    # print(MazF)
    # print(M)

    with open("withoutef_A_in.txt","w") as f:
        for i in A_in:
            f.write(str(i))
            f.write('\n')
    with open("withoutef_M.txt","w") as f:
        for i in M:
            f.write(str(i))
            f.write('\n')
    with open("withoutef.txt","w") as f:
        for i in MazF:
            f.write(str(i))
            f.write('\n')


    plt.subplot(3, 2, 1)
    plt.plot(time, A_in)
    plt.subplot(3, 2, 2)
    plt.plot(time, E_mRNA)
    plt.subplot(3, 2, 3)
    plt.plot(time, FGH_mRNA)
    plt.subplot(3, 2, 4)
    plt.plot(time, M)
    plt.subplot(3, 2, 5)
    plt.plot(time, MazF)
    plt.subplot(3, 2, 6)
    plt.plot(time, MazF, M)

    plt.show()
    print(2)
