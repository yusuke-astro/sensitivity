import  numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

GHz = 1e9
p = 1e-12
W = 0.13106*p

F = 50*GHz
def P(f):
    global GHz
    global F
    WL = 299792458/F############################################
    WL
    u = 1e-6
    m = 1e-3
    e_0 = 8.85418782*10**(-12)
    h = 6.62607004*10**(-34)
    k_B = 1.38064852*10**(-23)
    T = 2.725

    HWP = 0.8952#0.8962

    Ape = 1.0 - np.exp(-((np.pi)**2 / 2.0) * ((24*m) / (2.75 * 3.0 * WL))**2)
    Ape

    e_Mir = 4.0 * np.sqrt(np.pi * F * e_0 * 1.39*10**(-8))###############################
    e_Mir
    r_Mir = 1.0 - np.exp(-((4.0 * np.pi * 2.0*u) / WL)**2)
    r_Mir

    Pri = 1.0 - e_Mir - r_Mir
    Pri
    Sec = Pri

    #K_20 =1.0# 0.0130

    e_2KF = 1.0 - np.exp((-2.0 * np.pi * 5.0*m * 1.5 * 2.3*10**(-4)) / WL)
    e_2KF
    KF_2 = 1.0 - 0.05 - e_2KF
    KF_2

    e_Lens = 1.0 - np.exp((-2.0 * np.pi * 9.0*m * 3.4 * 5.0*10**(-5)) / WL)
    e_Lens
    Lens = 1.0 - 0.05 - e_Lens
    Lens

    det = 1.0 - 0.32
    det
    fermi = (h * f) / (np.exp((h * f) / (k_B * T)) - 1.0)

    return HWP * Ape * Pri * Sec * KF_2 * Lens * det * fermi

print("center=",P(F))###################################



width = (F * 0.3)/2.0
width
F-width
F+width

freq = np.arange( F - width, F + width, 0.001*GHz)

P(F)
#plt.plot(freq, P(freq),"x")
value = integrate.quad(P, F - width, F + width)
print(value)
print(W)
