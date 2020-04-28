import  numpy as np
from scipy import integrate

GHz = 1e9
p = 1e-12
W = 0.08883*p
u = 1e-6
m = 1e-3
e_0 = 8.85418782*10**(-12)
h = 6.62607004*10**(-34)
k_B = 1.38064852*10**(-23)

f = 40*GHz
WL = 299792458/f
F = f
T = 2.725

def width(f):
    global GHz
    global p
    global W
    global u
    global m
    global e_0
    global h
    global k_B

    if f == 40*GHz:
        frac = 0.30
    if f == 50*GHz:
        frac = 0.30
    if f == 60*GHz:
        frac = 0.23
    if f == 68*GHz:
        frac = 0.23
    if f == 78*GHz:
        frac = 0.23
    if f == 89*GHz:
        frac = 0.23
    if f == 100*GHz:
        frac = 0.23
    if f == 119*GHz:
        frac = 0.30
    if f == 140*GHz:
        frac = 0.30
    BW = (f * frac)/2.0
    return BW

def r_HWP_LFT(f, BW):
    global GHz
    global p
    global W
    global u
    global m
    global e_0
    global h
    global k_B
    if f == 40*GHz:
        r_HWP = 0.1000
    if f == 50*GHz:
        r_HWP = 0.1000
    if f == 60*GHz:
        r_HWP = 0.0682
    if f == 68*GHz:
        r_HWP = 0.0514
    if f == 78*GHz:
        r_HWP = 0.0415
    if f == 89*GHz:
        r_HWP = 0.0368
    if f == 100*GHz:
        r_HWP = 0.0275
    if f == 119*GHz:
        r_HWP = 0.013
    if f == 140*GHz:
        r_HWP = 0.0014
    return r_HWP
r_HWP_LFT(f, BW)
def epsi_HWP_LFT(f, BW):
    global GHz
    global p
    global W
    global u
    global m
    global e_0
    global h
    global k_B
    if f == 40*GHz:
        epsi_HWP = 0.0038
    if f == 50*GHz:
        epsi_HWP = 0.0048
    if f == 60*GHz:
        epsi_HWP = 0.058
    if f == 68*GHz:
        epsi_HWP = 0.066
    if f == 78*GHz:
        epsi_HWP = 0.075
    if f == 89*GHz:
        epsi_HWP = 0.097
    if f == 100*GHz:
        epsi_HWP = 0.0115
    if f == 119*GHz:
        epsi_HWP = 0.0135
    if f == 140*GHz:
        epsi_HWP = 0.0161
    return epsi_HWP
epsi_HWP_LFT(f, BW)

def ita_HWP_LFT(f):
    global GHz
    global p
    global W
    global u
    global m
    global e_0
    global h
    global k_B
    if f == 40*GHz:
        ita_HWP = 0.8962
    if f == 50*GHz:
        ita_HWP = 0.8952
    if f == 60*GHz:
        ita_HWP = 0.9260
    if f == 68*GHz:
        ita_HWP = 0.9420
    if f == 78*GHz:
        ita_HWP = 0.9510
    if f == 89*GHz:
        ita_HWP = 0.9535
    if f == 100*GHz:
        ita_HWP = 0.9610
    if f == 119*GHz:
        ita_HWP = 0.9735
    if f == 140*GHz:
        ita_HWP = 0.9825
    return ita_HWP

def ita_Apt_LFT(f, D):
    global GHz
    global p
    global W
    global u
    global m
    global e_0
    global h
    global k_B
    global WL
    ita_Apt = 1.0 - np.exp(-((np.pi)**2 / 2.0) * ((24*m) / (2.75 * 3.0 * WL))**2)
    return ita_Apt

def epsi_Mir(f):
    global GHz
    global p
    global W
    global u
    global m
    global e_0
    global h
    global k_B
    global WL
    e_Mir = 4.0 * np.sqrt(np.pi * f * e_0 * 1.39*10**(-8))
    return e_Mir
epsi_Mir(f)

def r_Mir(f):
    global GHz
    global p
    global W
    global u
    global m
    global e_0
    global h
    global k_B
    global WL
    r_Mir = 1.0 - np.exp(-((4.0 * np.pi * 2.0*u) / WL)**2)
    return r_Mir
r_Mir(f)

def ita_Mir(f):
    global GHz
    global p
    global W
    global u
    global m
    global e_0
    global h
    global k_B
    global WL
    ita_Mir = 1.0 - epsi_Mir(f) - r_Mir(f)
    return ita_Mir
ita_Mir(f)

def epsi_2KF(f):
    global GHz
    global p
    global W
    global u
    global m
    global e_0
    global h
    global k_B
    global WL
    e_2K = 1.0 - np.exp((-2.0 * np.pi * 5.0*m * 1.5 * 2.3*10**(-4)) / WL)
    return e_2K
epsi_2KF(f)

def ita_2K(f):
    global GHz
    global p
    global W
    global u
    global m
    global e_0
    global h
    global k_B
    global WL
    ita_2KF = 1.0 - 0.05 - epsi_2KF(f)
    return ita_2KF
ita_2K(f)

def ita_5K(f, D):
    return 1.0 - ita_Apt_LFT(f, D)

def epsi_Lenslet(f):
    global GHz
    global p
    global W
    global u
    global m
    global e_0
    global h
    global k_B
    global WL
    epsi_Lens = 1.0 - np.exp((-2.0 * np.pi * 9.0*m * 3.4 * 5.0*10**(-5)) / WL)
    return epsi_Lens
epsi_Lenslet(f)

def ita_Lenslet(f):
    global GHz
    global p
    global W
    global u
    global m
    global e_0
    global h
    global k_B
    global WL
    ita_Lens = 1.0 - 0.05 - epsi_Lenslet(f)
    return ita_Lens
ita_Lenslet(f)

def r_det_LFT():
    global GHz
    global p
    global W
    global u
    global m
    global e_0
    global h
    global k_B
    r_det = 0.32
    return r_det

def ita_det_LFT():
    global GHz
    global p
    global W
    global u
    global m
    global e_0
    global h
    global k_B
    ita_det = 1.0 - r_det_LFT()
    return ita_det

def Bose(f,T):
    global GHz
    global p
    global W
    global u
    global m
    global e_0
    global h
    global k_B
    global F
    bose = (h * f) / (np.exp((h * f) / (k_B * T)) - 1.0)
    return bose

def P_CMB(f):
    global GHz
    global p
    global W
    global u
    global m
    global e_0
    global h
    global k_B
    global F
    global T
    global WL
    CMB = ita_HWP_LFT(F)*ita_Apt_LFT(f,24.0)*ita_Mir(f)*ita_Mir(f)*ita_2K(f)*ita_Lenslet(f)*ita_det_LFT()*Bose(f,T) * 1.000
    return CMB

BW = width(F)
BW
f
P_CMB(F)
ita_HWP_LFT(F)
ita_Apt_LFT(f,24.0)
ita_Mir(f)
ita_2K(f)
ita_Lenslet(f)
ita_det_LFT()
"""
| p_cmb  7.44698507836e-24 ,
|
| eff_hwp 0.896 ,
|
| eff_apt 0.525 ,
|
| eff_mirror1 0.999 ,
|
| eff_mirror2 0.999 ,
|
| eff_filter 0.949 ,
|
| eff_lenslet 0.949 ,
|
| eff_detector 0.68

"""
F
BW
value = integrate.quad(P_CMB, F-BW, F+BW)
print(value)
