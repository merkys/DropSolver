#!/usr/bin/python3

from scipy.optimize import Bounds, minimize
import numpy

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
    return lhs(Qoil, wjet0sol) - rhs(Qoil, wjet0sol)
print(minimize(lhs_rhs_diff, [10 ** -5], args=(QoilStart,), bounds=Bounds(10 ** -5, 10 ** -4)))
