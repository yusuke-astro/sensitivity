import  numpy as np
from scipy import integrate

GHz = 1e9
p = 1e-12
W = 0.0883*p

def P(f):
    global GHz
    WL = 299792458/(50*GHz)
    WL
    u = 1e-6
    m = 1e-3
    e_0 = 1.25663706212*10**(-6)
    h = 6.62607004*10**(-34)
    k_B = 1.38064852*10**(-23)
    T = 2.725

    HWP = 0.8952   ########0.8962

    Ape = 1.0 - np.exp(-((np.pi)**2 / 2.0) * (2.75 / (3.0 * WL))**2)

    e_Mir = 4.0 * np.sqrt(np.pi * f * e_0 * 1.39*10**(-8))
    r_Mir = 1.0 - np.exp(-((4.0 * np.pi * 2.0*u) / WL)**2)
    Pri = 1.0 - e_Mir - r_Mir

    Sec = Pri

    K_20 = 0.0084  ####0.0130

    e_2KF = 1.0 - np.exp((-2.0 * np.pi * 5.0*m * 1.5 * 2.3*10**(-4)) / WL)
    KF_2 = 1.0 - 0.05 - e_2KF

    e_Lens = 1.0 - np.exp((-2.0 * np.pi * 9.0*m * 3.4 * 5.0*10**(-4)) / WL)
    Lens = 1.0 - 0.05 - e_Lens

    det = 1.0 - 0.32

    fermi = (h * f) / (np.exp((h * f) / (k_B * T)) - 1.0)

    return HWP * Ape * Pri * Sec * K_20 * KF_2 * Lens * det * fermi

value = integrate.quad(P, 34*GHz, 46*GHz)
print(value)
print(W)
