import  numpy as np
from scipy import integrate
from my_module import my_sens_func as st
import matplotlib.pyplot as plt

GHz = 1e9
p = 1e-12
Tera=1e12
W = 0.08883*p
u = 1e-6
m = 1e-3
e_0 = 8.85418782*10**(-12)
h = 6.62607004*10**(-34)
k_B = 1.38064852*10**(-23)
c = 299792458
############################################

f =40*GHz
WL = c/f
F = f

#freqに対してDが2つあるところがある
if f == 40*GHz:
    D = 24.0
if f == 60*GHz:
    D = 24.0
if f == 78*GHz:
    #D = 24.0
    D = 16.0
if f == 50*GHz:
    D = 24.0
if f == 68*GHz:
    #D = 24.0
    D = 16.0
if f == 89*GHz:
    D = 24.0
    D = 16.0
if f == 100*GHz:
    D = 16.0
if f == 119*GHz:
    D = 16.0
if f == 140*GHz:
    D = 16.0

st.D = D
#T = 2.725

def P_CMB(f):
    global F
    global D
    T = 2.725
    return st.ita_HWP_LFT(F)*st.ita_Apt_LFT(f,D)*st.ita_Mir(f)*st.ita_Mir(f)*st.ita_2K(f)*st.ita_Lenslet(f)*st.ita_det_LFT()*st.Bose(f,T) * 1.000


def P_HWP(f):
    global F
    global D
    T = 20.00
    T_r = 5.000
    HWP_e = st.ita_Apt_LFT(f,D)*st.ita_Mir(f)*st.ita_Mir(f)*st.ita_2K(f)*st.ita_Lenslet(f)*st.ita_det_LFT()*st.Bose(f,T) * st.epsi_HWP_LFT(F)
    HWP_r = st.ita_Apt_LFT(f,D)*st.ita_Mir(f)*st.ita_Mir(f)*st.ita_2K(f)*st.ita_Lenslet(f)*st.ita_det_LFT()*st.Bose(f,T_r) * st.r_HWP_LFT(F)
    return HWP_e + HWP_r

def P_Ape(f):
    global F
    global D
    T = 5.000
    return st.ita_5K(f,D)*st.ita_Mir(f)*st.ita_Mir(f)*st.ita_2K(f)*st.ita_Lenslet(f)*st.ita_det_LFT()*st.Bose(f,T)

def P_m1(f):
    global F
    T = 5.000
    Pri = st.ita_Mir(f)*st.ita_2K(f)*st.ita_Lenslet(f)*st.ita_det_LFT()*st.Bose(f,T)
    Pri_e = Pri * st.epsi_Mir(f)
    Pri_r = Pri * st.r_Mir(f)
    return Pri_e + Pri_r

def P_m2(f):
    global F
    T = 5.000
    Sec = st.ita_2K(f)*st.ita_Lenslet(f)*st.ita_det_LFT()*st.Bose(f,T)
    Sec_e = Sec * st.epsi_Mir(f)
    Sec_r = Sec * st.r_Mir(f)
    return Sec_e + Sec_r

def P_20K(f):
    global F
    global D
    T = 20.0
    return st.ita_20K(F)*st.ita_2K(f)*st.ita_Lenslet(f)*st.ita_det_LFT()*st.Bose(f,T) * 1.000

def P_2KF(f):
    global F
    T = 2.000
    T_r = 0.100
    Flt2K_e = st.ita_Lenslet(f)*st.ita_det_LFT()*st.Bose(f,T) * st.epsi_2KF(f)
    Flt2K_r = st.ita_Lenslet(f)*st.ita_det_LFT()*st.Bose(f,T_r) * st.r_2KF()
    return Flt2K_e + Flt2K_r

def P_Lens(f):
    global F
    T = 0.100
    T_r = 5.000
    Lens_e = st.ita_det_LFT()*st.Bose(f,T) * st.epsi_Lenslet(f)
    Lens_r = st.ita_det_LFT()*st.Bose(f,T_r) * st.r_Lenslet()
    return Lens_e + Lens_r

def P_Det():
    T = 0.100
    Det_e = st.Bose(f,T) * st.r_det_LFT()
    Det_r = st.Bose(f,T) * st.ita_det_LFT()
    return Det_e + Det_r

def P_Opt(f):
    return P_CMB(f)+P_HWP(f)+P_Ape(f)+P_m1(f)+P_m2(f)+P_20K(f)+P_2KF(f)+P_Lens(f)+P_Det()


BW = st.width(F)

value_pcmb = integrate.quad(P_CMB, F-BW, F+BW)
value = integrate.quad(P_Opt, F-BW, F+BW)
print("\nf=",f/GHz)
#print("BW=", BW/GHz)
print("P_CMB({:.0})={}".format(f, P_CMB(f)))
print("P_HWP({:.0})={}".format(f, P_HWP(f)))
print("P_APT({:.0})={}".format(f, P_Ape(f)))
print("P_20K({:.0})={}".format(f, P_20K(f)))
print("P_m1({:.0})={}".format(f, P_m1(f)))
print("P_m2({:.0})={}".format(f, P_m2(f)))
print("P_2KF({:.0})={}".format(f, P_2KF(f)))
print("P_Lens({:.0})={}".format(f, P_Lens(f)))
print("P_Det()={:.0}".format(P_Det()))
print("P_Opt({:.0})={}".format(f, P_Opt(f)))
"""
freq 40.0 , 
p_cmb  7.44698507836e-24 ,
p_hwp 2.13635191675e-24 , 
p_apt 1.6458882942e-23 , 
p_20K 2.09293683272e-24 , 
p_m1 1.76179775095e-26 , 
p_m2 1.76269440977e-26 , 
p_filter 1.53348853609e-26 , 
p_lens 1.92529436274e-24 ,
p_opt 3.01110309396e-23

"""
print("int_P_cmb={0:.6}".format(value_pcmb[0]/p))
print("int_P_Opt={0:.6}".format(value[0]/p))


freq = np.array([40,60,78,50,68,89,68,89,119,78,100,140])*GHz
#freq = np.array([40,50,60,68,78,89,100,119,140,166,195,235,280,337,402])*GHz
print("\nIn the LFT freqency...")
INTEGRATE = []
for i in range(len(freq)):
    #print(freq[i])
    WL = c/freq[i]
    F = freq[i]
    BW = st.width(F)
    #print(BW)

    INTEGRATE.append(integrate.quad(P_Opt, F-BW, F+BW)[0])
    print("Freq={}, P_Opt={:.2}".format(freq[i]/GHz,INTEGRATE[i]/p))
INTEGEATE = np.array(INTEGRATE)

hasebe = np.array([0.36,0.29,0.28,0.38,0.28,0.28,0.36,0.34,0.40,0.36,0.33,0.38])*p
delta = hasebe - INTEGEATE
plt.figure()
plt.plot(freq/GHz,hasebe,"o",label="Hasebe-san")
plt.plot(freq/GHz,INTEGEATE,"o",label="Takase")
plt.grid()
plt.legend()

plt.figure()
plt.plot(freq/GHz,delta,"o", label="Hase-Taka")
plt.grid()
plt.legend()

plt.show()
