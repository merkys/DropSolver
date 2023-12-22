#!/usr/bin/python3

from scipy.interpolate import interp1d
from scipy.optimize import Bounds, minimize
import numpy
import sympy

Kd = 0.001      # water viscosity, [Pa*s] *)
Kvisc = 0.0014  # oil consistency index, [Pa*s^n] *)
EtaZero = 0.0014
EtaInf = 0.0014
wn = 7 * 10 ** -5   # width of the focussing nozzle, [m] *)(* 20/20/H11; 28x15.5/H17; 33/5/H18; 37.5/4.5/20; 10x17; 33.5x5//17; 37x4//H20*)
Ln = 7 * 10 ** -5   # length of the focussing nozzle, [m] *)
H = 8 * 10 ** -5    # Height of the chip, [m] *)
sigmaEQ = (52) * 0.001  # surface tension of water/oil interface with no surfactant *)
wcont = 6 * 10 ** -5    # Width of the Continuous-phase channel *)
wdisp = 7 * 10 ** -5    # Width of the Dispersed-phase channel *)
wout = 11 * 10 ** -5    # Width of the outlet channel *)
RhoO = 1614     # Oil density [kg/m^3] *)
RhoW = 1000     # Water density, [kg/m^3] *)
Ms = 7.5    # Molar mass of surfactant [kg/mol] *)
omega = 0.006   # surfactant fraction, USED IN OIL FOR EXPERIMENT, w/w *)
NA = 6.023 * 10 ** 23 # Avogadro Number *)
R = 8.314   # Universal gas-constant *)
T = 293     # Absolute temperature, [K] *)
GAMMAinf = 8 * 10 ** -6 # Surface excess of surfactant at complete-coverage of interface [mol/m^2] *)
Kads = 1400
Kdes = 0.0002
CMC = 0.128     # Critical Micelle Concentration, [mol/m^3] *)
n = 1.0     # NONLINEARITY INDEX for continuous phase*)
p = 1.0     # NONLINEARITY INDEX for dispersed phase*)
EpsilonHFE = 5.8
F = 96400
z = 1
r0 = 0.2 * 10 ** -9 * 45 ** 0.6

Tau0c = 0.0026
def GammaDOTc(Qoil, wjet0sol):
    return Qoil/H ** 2/wn - wjet0sol
def etaoNN(Qoil, wjet0sol):
    return EtaInf + (EtaZero - EtaInf)/(1 + (GammaDOTc(Qoil, wjet0sol)/77) ** n)
# 0.03267+(0.02946-0.03267)/(1+(GammaDOTc(Qoil, wjet0sol)/99.7) ** n); M8410 Min.Oil Carreau-model
# etaoNN=Tau0c/(GammaDOTc(Qoil,wjet0sol)+Kvisc*GammaDOTc(Qoil, wjet0sol) ** (n-1); M5904 Min.Oil HB-model

def GammaDOTd(Qoil, wjet0sol):
    return Qw / H ** 2 / wjet0sol
# def etaw(wjet0sol, Qw, H, wn):
#     return etaINF+(Kd-etaINF)/(1+(B1*GammaDOTd(wjet0sol, Qw, H, wn)) ** p)/.etaINF→0.00532/.Kd→1.276/.B1→4.578;
def etaw(Qoil, wjet0sol):
    etaINF1 = 0.001
    B1 = 4.691
    return etaINF1 + (Kd - etaINF1)/(1 + (B1 * GammaDOTd(Qoil, wjet0sol)) ** p)

def lhs(Qoil, wjet0sol):
    return wjet0sol * (1 + etaw(Qoil, wjet0sol) * ((Qw/Qoil)/etaoNN(Qoil, wjet0sol)))
def rhs(Qoil, wjet0sol):
    return wn * (etaw(Qoil, wjet0sol) * Qw / Qoil / etaoNN(Qoil, wjet0sol))

Qw        = 220 * 2.78 * 10 ** -13;
QoilStart =  90 * 2.78 * 10 ** -13;
QoilEnd   = 600 * 2.78 * 10 ** -13;
QoilStep  = 100 * 2.78 * 10 ** -13;
def lhs_rhs_diff(wjet0sol, Qoil):
    return (lhs(Qoil, wjet0sol) - rhs(Qoil, wjet0sol))**2
data_x = numpy.arange(QoilStart, QoilEnd, QoilStep)
data_y = []
for Qoil in data_x:
    data_y.append(*list(minimize(lhs_rhs_diff, [10 ** -5], args=(Qoil,), bounds=Bounds(10 ** -5, 10 ** -4)).x))
print(data_x)
print(data_y)

wjet0solF = interp1d(data_x, data_y, kind='quadratic')

# Second step: with known wjet0solF[Qoil] function solve ....
# In[69]:=
def GammaDOTc(Qoil):
    return (Qoil/H**2)/(wn - wjet0solF(Qoil))
def etaoNN(Qoil):
    return EtaInf + (EtaZero - EtaInf)/(1+(GammaDOTc(Qoil)/77)**n)
def GammaDOTd(Qoil):
    return (Qw/H**2)/wjet0solF(Qoil)
# (*etaw[wjet0sol_,Qw_,H_,wn_]:=etaINF+(Kd-etaINF)/(1+(B1*GammaDOTd[wjet0sol,Qw,H,wn])^p)/.etaINF\[Rule]0.00532/.Kd\[Rule]1.276/.B1\[Rule]4.578;*)
def etaw(Qoil):
    etaINF = 0.001
    B1 = 4.691
    return etaINF+(Kd-etaINF)/(1+(B1*GammaDOTd(Qoil))**p)

# In[73]:=
def dPO(Qoil):
    return abs(Qoil/H**3/wn*12*etaoNN(Qoil)*Ln/(1-0.63*Ln/wn))
def dPW(Qoil):
    return abs(Qw/H**3/wn*12*etaw(Qoil)*Ln/(1-0.63*Ln/wn))
def HF(Qoil):
    Y = 1*10**6
    alfa = 0.7
    return 1+(3/2*alfa*(dPO(Qoil)+dPW(Qoil))*wn/Y/H)**0.25

# In[76]:=
print([dPO(QoilStart), dPW(QoilStart), HF(QoilStart)])
# Out[76]= {2.21883,3.87416,1.04864}

# Tdrop is simply global var, we don't need to pass it as function argument
Tdrop = sympy.Symbol('Tdrop')
Pi = sympy.pi

# In[77]:=
Ldrop = H+(Tdrop*Qw-Pi/6*H**3)/wout/H
def L(Qoil):
    return 2*((wdisp - wjet0solF(Qoil))**2+4 * wcont**2)**0.5+2*Ln+2*Ldrop+Pi*H

# In[79]:= L[QoilStart]
# Out[79]= 0.000636683 +2 (1/12500+1250000000/11 (-(\[Pi]/11718750000000)+6.116*10^-11 Tdrop))

# In[80]:=
Cbulk=RhoO*omega/Ms
Nmon0 = (1+2*((Cbulk/CMC)*(Cbulk-CMC))**0.5).real # CHECK: What is Re()?
Dmon = (R*T/NA)/(3*Pi*Kvisc*r0)
Dmic = (R*T/NA)/(3*Pi*Kvisc*r0*Nmon0**0.3333)
Deff = ((Cbulk-CMC)*Dmic+CMC*Dmon)/Cbulk

print(Cbulk/CMC)
print(Dmon.evalf())
print(Dmic.evalf())
print(Deff.evalf())

## The bottom part I do not understand yet (A.M.)

# In[86]:= (*Micelle kinetics and surfactant adsorption*)
CBmic = (Cbulk-CMC) # (*Qoil*Tdrop*NA/.wjet0sol->sol1;*)
LambdaAgg = 0.5*(r0+r0*Nmon0**0.333)*(Dmon+Dmic)
# (* Cmic is a function of variable Tdrop!*)
CDECmic=Cmic(Tdrop)*(LambdaF+LambdaS)(*Tdrop*)/.LambdaF->10**6*(2/Nmon0)/.LambdaS->10**2*(Nmon0-2)/Nmon0 # (*/.LambdaS\[Rule]0.01*LambdaF*); # (*F\[Rule] Fast, S\[Rule]Slow *)
CAGGmic = (LambdaAgg/Nmon0)*(Cmon(Tdrop)) # (*Nmic/Nmon0*LambdaAgg*Tdrop*)

CBmon=CMC(*Qoil*Tdrop*NA/.wjet0sol->sol1*);
CRELmon=Module[{LambdaF=10^6*(2/Nmon0),LambdaS=10^2*(Nmon0-2)/Nmon0},Cmic[Tdrop]*(LambdaF+LambdaS)]
CAGGmon=LambdaAgg*(Cmon(Tdrop)) # (*Nmon*(1+LambdaAgg/Nmon0*Tdrop)*)

def Pe(Qoil):
    return (Qoil/H/wcont)*wout/Dmon
