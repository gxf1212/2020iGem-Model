import scipy
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def arabinoseInsideFunction(A_in,V_e,A_ex,E_mRNA,K_e,V_fgh,FGH_mRNA,K_fgh,Mu,Gamma_ain):
    return V_e*A_ex*E_mRNA/(K_e+A_ex)+V_fgh*A_ex*FGH_mRNA/(K_fgh+A_ex)-(Mu+Gamma_ain)*A_in
def araEmRNADynamicsFunction(E_mRNA,Alpha_e,V_me,A_in,K_me,Mu,Gamma_e):
    return Alpha_e+V_me*(A_in**3)/(K_me**3+A_in**3)-(Mu+Gamma_e)*E_mRNA
def araFGHmRNADynamicsFunction(FGH_mRNA,Alpha_fgh,V_mfgh,A_in,K_mfgh,Mu,Gamma_fgh):
    return Alpha_fgh+V_mfgh*(A_in**3)/(K_mfgh**3+A_in**3)-(Mu+Gamma_fgh)*FGH_mRNA
def mRNAFunction(M,D,k_m,A_in,K_1,K,d_m,d_large,MazF,K_t):
    return D*k_m*(1+K_1*(A_in**2))/(K+K_1*(A_in**2))-d_m*M-d_large*(MazF**2)/(MazF**2+K_t**2)*M
def mazfFunction(b_2,M,a_T,MazF,MazE,d_c,d_a2,MazEF):
    return b_2*M-a_T*MazF*MazE-d_c*MazF+d_a2*MazEF
def mazeFunction(b_1,M,a_T,MazF,MazE,d_a):
    return b_1*M-a_T*MazF*MazE-d_a*MazE
def mazefFunction(a_T,MazE,MazF,d_c,MazEF,d_a2):
    return a_T*MazE*MazF-d_c*MazEF-d_a2*MazEF
def totalModel(y,t,parameters):
    A_in,E_mRNA,FGH_mRNA,M,MazF,MazE,MazEF=y
    V_e, A_ex, K_e, V_fgh, K_fgh, Mu, Gamma_ain,\
    Alpha_e,V_me,K_me,Gamma_e, \
    Alpha_fgh, V_mfgh, K_mfgh, Gamma_fgh, \
    D, k_m, K_1, K, d_m, d_large,  K_t, \
    b_2, a_T,  d_c, d_a2,  b_1,d_a   \
        =parameters
    _ain=arabinoseInsideFunction(A_in,V_e,A_ex,E_mRNA,K_e,V_fgh,FGH_mRNA,K_fgh,Mu,Gamma_ain)
    _emRNA=araEmRNADynamicsFunction(E_mRNA,Alpha_e,V_me,A_in,K_me,Mu,Gamma_e)
    _fghmRNA=araFGHmRNADynamicsFunction(FGH_mRNA,Alpha_fgh,V_mfgh,A_in,K_mfgh,Mu,Gamma_fgh)
    _mRNA=mRNAFunction(M,D,k_m,A_in,K_1,K,d_m,d_large,MazF,K_t)
    _maze=mazeFunction(b_1,M,a_T,MazF,MazE,d_a)
    _mazf=mazfFunction(b_2,M,a_T,MazF,MazE,d_c,d_a2,MazEF)
    _mazef=mazefFunction(a_T,MazE,MazF,d_c,MazEF,d_a2)
    return np.array([_ain,_emRNA,_fghmRNA,_mRNA,_maze,_mazef,_mazf])
if '__main__':
    V_e=48
    A_ex=1E-6
    K_e=1.4E-4
    V_fgh =4.8
    K_fgh=6E-6
    Mu=5.78E-4
    Gamma_ain=0.0155
    Alpha_e=2.83E-11
    V_me=4.28E-9
    K_me=1.06E-4
    Gamma_e=0.011552
    Alpha_fgh=7.55E-11
    V_mfgh=1.1266E-8
    K_mfgh=1E-9
    Gamma_fgh=0.0115
    D=2
    k_m=0.3
    K_1=2.52E10
    K=7200
    d_m=0.002
    d_large=0.198
    K_t=15
    b_2=0.009
    a_T=0.000365
    d_c=0.00028
    d_a2=0.000231
    b_1=0.122
    d_a=0.00231
    parameters=[
        V_e, A_ex, K_e, V_fgh, K_fgh, Mu, Gamma_ain,\
        Alpha_e,V_me,K_me,Gamma_e, \
        Alpha_fgh, V_mfgh, K_mfgh, Gamma_fgh, \
        D, k_m, K_1, K, d_m, d_large,  K_t, \
        b_2, a_T,  d_c, d_a2,  b_1,d_a]
    time = np.arange(0, 7000, 1)
    solve=odeint(totalModel,[0,0,0,0,0,0,0],time,(parameters,),atol=1e-7, rtol=1e-11,full_output=1)
    A_in, E_mRNA, FGH_mRNA, M, MazF, MazE, MazEF=solve[0][:,0],solve[0][:,1],solve[0][:,2],solve[0][:,3],solve[0][:,4],solve[0][:,5],solve[0][:,6]
    plt.subplot(4,2,1)
    plt.plot(time,A_in)
    plt.subplot(4, 2, 2)
    plt.plot(time, E_mRNA)
    plt.subplot(4, 2, 3)
    plt.plot(time, FGH_mRNA)
    plt.subplot(4, 2, 4)
    plt.plot(time, M)
    plt.subplot(4, 2, 5)
    plt.plot(time, MazE)
    plt.subplot(4, 2, 6)
    plt.plot(time, MazF)
    plt.subplot(4, 2, 7)
    plt.plot(time, MazEF)
    plt.show()
    print(1)