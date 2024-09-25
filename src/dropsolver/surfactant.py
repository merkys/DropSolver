from scipy.interpolate import interp1d
from scipy.optimize import newton
from dropsolver.util import inclusive_range
import numpy
from sympy import Eq, Function, Symbol, dsolve, exp, log, pi

def calculate(Kd=0.001, etaINF1=0.001, B1=4.691, p=1.0, Kvisc=0.0014, EtaZero=0.0014, EtaInf=0.0014, B2=77.0, n=1.0, wn=7e-05, Ln=6e-05, H=8e-05, wcont=7e-05, wdisp=6e-05, wout=0.00011, omega=0.006, sigmaEQ=0.052, Qw=6.116e-11, QoilStart=2.502e-11, QoilEnd=1.668e-10, QoilStep=2.78e-11, Ms=7.5, GAMMAinf=8e-06, Kads=600, Kdes=8.57e-05, Enth=5.629, m=0.06, CMC=0.128, is_ionic=True, debug=False, reporter=None):
    """
    Calculate the droplet volume for given microfluid chip parameters.
    """
    # PARAMETERS SET BY THE CHOICE OF OIL/SURFACTANT TYPE; further not changed by user
    # SURFACTANT AND ADSORPTION PARAMETERS

    if not is_ionic:
        Ms = 12.5
        GAMMAinf = 3.4 * 10 ** -6
        Kads = 18
        Kdes = 0.006
        m = 0
        Enth = 0
        CMC = 0.15

    # CONSTANTS
    RhoO = 1614             # Oil density [kg/m^3]
    RhoW = 1000             # Water density, [kg/m^3]
    NA = 6.023 * 10 ** 23   # Avogadro Number
    R = 8.314               # Universal gas-constant
    T = 293                 # Absolute temperature, [K]
    EpsilonHFE = 5.8        # Dielectric constant of oil
    F = 96400               # Faraday's number *)
    z = 1
    r0 = 0.2 * 10 ** -9 * 45 ** 0.6

    def GammaDOTc(Qoil, wjet0sol):
        return Qoil/H ** 2/(wn - wjet0sol)
    def etaoNN(Qoil, wjet0sol):
        return EtaInf + (EtaZero - EtaInf)/(1 + (B2 * GammaDOTc(Qoil, wjet0sol)) ** n)

    # In[55]:=
    def GammaDOTd(Qoil, wjet0sol):
        return Qw / H ** 2 / wjet0sol
    # def etaw(wjet0sol, Qw, H, wn):
    #     return etaINF+(Kd-etaINF)/(1+(B1*GammaDOTd(wjet0sol, Qw, H, wn)) ** p)/.etaINF→0.00532/.Kd→1.276/.B1→4.578;
    def etaw(Qoil, wjet0sol):
        return etaINF1 + (Kd - etaINF1)/(1 + (B1 * GammaDOTd(Qoil, wjet0sol)) ** p)

    def lhs(Qoil, wjet0sol):
        return wjet0sol * (1 + etaw(Qoil, wjet0sol) * ((Qw/Qoil)/etaoNN(Qoil, wjet0sol)))
    def rhs(Qoil, wjet0sol):
        return wn * (etaw(Qoil, wjet0sol) * Qw / Qoil / etaoNN(Qoil, wjet0sol))

    # In[62]:=
    def lhs_rhs_diff(wjet0sol, Qoil):
        return lhs(Qoil, wjet0sol) - rhs(Qoil, wjet0sol)
    data_x = inclusive_range(QoilStart, QoilEnd, QoilStep)
    data_y = []
    for Qoil in data_x:
        data_y.append(newton(lhs_rhs_diff, x0=(Kd*Qw)/(EtaZero*QoilEnd)*wn, x1=wn, args=(Qoil,)))
    if debug:
        print(data_x)
        print(data_y)

    wjet0solF = interp1d(data_x, data_y, kind='quadratic')

    # Second step: with known wjet0solF[Qoil] function solve ....
    # In[69]:=
    def GammaDOTc(Qoil):
        return (Qoil/H**2)/(wn - wjet0solF(Qoil))
    def etaoNN(Qoil):
        return EtaInf + (EtaZero - EtaInf)/(1 + (B2 * GammaDOTc(Qoil))**n)
    def GammaDOTd(Qoil):
        return Qw / H ** 2 / wjet0solF(Qoil)
    def etaw(Qoil):
        return etaINF1 + (Kd - etaINF1)/(1 + (B1 * GammaDOTd(Qoil)) ** p)

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
    if debug:
        print([dPO(QoilStart), dPW(QoilStart), HF(QoilStart)])
    # Out[76]= {2.21883,3.87416,1.04864}

    # Tdrop is simply global var, we don't need to pass it as function argument
    Tdrop = Symbol('Tdrop')

    # In[77]:=
    def Ldrop(Tdrop=Tdrop):
        return H+(Tdrop*Qw-pi/6*H**3)/wout/H
    def L(Qoil, Tdrop=Tdrop):
        return 2*((wdisp - wjet0solF(Qoil))**2+4 * wcont**2)**0.5+2*Ln+2*Ldrop(Tdrop)+pi*H

    # In[79]:=
    if debug:
        print(L(QoilStart))
    # Out[79]= 0.000636683 +2 (1/12500+1250000000/11 (-(\[Pi]/11718750000000)+6.116*10^-11 Tdrop))

    # In[80]:=
    Cbulk=RhoO*omega/Ms
    Nmon0 = (1+2*((Cbulk/CMC)*(Cbulk-CMC))**0.5).real
    Dmon = (R*T/NA)/(3*pi*Kvisc*r0)
    Dmic = (R*T/NA)/(3*pi*Kvisc*r0*Nmon0**(1/3))
    Deff = ((Cbulk-CMC)*Dmic+CMC*Dmon)/Cbulk

    if debug:
        print(Cbulk/CMC)
        print(Dmon.evalf())
        print(Dmic.evalf())
        print(Deff.evalf())
    # Out[81]= 10.0875
    # Out[83]= 1.56138*10^-10
    # Out[84]= 7.85654*10^-11
    # Out[85]= 8.62553*10^-11

    # In[86]:= (*Micelle kinetics and surfactant adsorption*)
    CBmic = (Cbulk-CMC) # (*Qoil*Tdrop*NA/.wjet0sol->sol1;*)
    if is_ionic:
        LambdaAgg = 0.5*(r0+r0*Nmon0**(1/3))*(Dmon+Dmic)
    else:
        LambdaAgg = Cbulk*(1+Nmon0**(1/3))**2/Nmon0**(1/3)
    # (* Cmic is a function of variable Tdrop!*)
    Cmic = Function('Cmic')
    # (*LambdaS->0.01*LambdaF*); # (*F->Fast, S->Slow *)
    LambdaF = 10 ** 6 * (2/Nmon0)
    if is_ionic:
        LambdaS = 10 ** 2 * (Nmon0-2)/Nmon0
    else:
        LambdaS = 10 ** 3 * (Nmon0-2)/Nmon0
    CDECmic = Cmic(Tdrop)*(LambdaF+LambdaS)
    Cmon = Function('Cmon')
    if is_ionic:
        CAGGmic = (LambdaAgg/Nmon0)*(Cmon(Tdrop)) # (*Nmic/Nmon0*LambdaAgg*Tdrop*)
    else:
        CAGGmic = LambdaAgg/Nmon0*(Cmon(Tdrop)) # CHECKME: Same as for non-ionic?

    CBmon = CMC # (*Qoil*Tdrop*NA/.wjet0sol->sol1*);
    CRELmon=Cmic(Tdrop) * (LambdaF + LambdaS)
    CAGGmon=LambdaAgg * Cmon(Tdrop) # (*Nmon*(1+LambdaAgg/Nmon0*Tdrop)*)

    def Pe(Qoil):
        return (Qoil/H/wcont)*wout/Dmon

    # In[95]:=
    Qoil = QoilStart
    eq = [Eq(CBmic+CAGGmic-CDECmic-(1+Pe(Qoil)) * Cmic(Tdrop).diff(Tdrop), 0),
          Eq(CBmon+CRELmon-CAGGmon-(1+Pe(Qoil)) * Cmon(Tdrop).diff(Tdrop), 0)]
    ABC = dsolve(eq, [Cmon(Tdrop), Cmic(Tdrop)], ics={Cmon(0): CMC, Cmic(0): Cbulk-CMC})
    ABC = [ABC[0].rhs.evalf(), ABC[1].rhs.evalf()]

    # In[96]:=
    def CmonINT(Qoil, Tdrop=Tdrop): # Return values do not seem to depend on Qoil
        if isinstance(Tdrop, Symbol):
            return ABC[0]
        else:
            return ABC[0].evalf(subs={'Tdrop': Tdrop})
    def CmicINT(Qoil, Tdrop=Tdrop): # Return values do not seem to depend on Qoil
        if isinstance(Tdrop, Symbol):
            return ABC[1]
        else:
            return ABC[1].evalf(subs={'Tdrop': Tdrop})

    # In[98]:=
    def Psi(Tdrop=Tdrop):
        return (9 * 10 ** 9) / EpsilonHFE * (1.6 * 10 ** -19) / (pi * Dmon * Tdrop)**0.5
    def Cfactor(Tdrop=Tdrop):
        return exp(-z * Psi(Tdrop) * (F / R / T))

    if is_ionic:
        TauI=(Kads*Cbulk**2+Kdes) ** -1
    else:
        TauI=(Kads*Cbulk+Kdes) ** -1
    def zz1(Qoil, Tdrop=Tdrop):
        if is_ionic:
            return Kads/Kdes*(Cfactor(Tdrop)*CmonINT(Qoil, Tdrop))**2
        else:
            return Kads/Kdes*CmonINT(Qoil, Tdrop)
    def zz2(Qoil, Tdrop=Tdrop):
        return zz1(Qoil, Tdrop)/(1+zz1(Qoil, Tdrop))

    def DGamma(Qoil, Tdrop=Tdrop):
        return 1-exp(-(Tdrop/TauI)/(1+Pe(Qoil)))*zz1(Qoil, Tdrop) * exp(-Enth*zz2(Qoil, Tdrop)**m)/(1+zz1(Qoil, Tdrop)*exp(-Enth*zz2(Qoil, Tdrop)**m))
    # DGammaEl[Qoil_]:=27.72*0.138/1.138*(DGamma[Qoil])^1.138;
    # SIGMAel[Qoil_]:=4*R*T/F*(2*EpsilonHFE*(8.85*10^-11)*R*T*CmonINT[Qoil])^0.5*(1-Cosh[z*Psi*(F/R/T)])
    def SIGMAio(Qoil, Tdrop=Tdrop):
        if is_ionic:
            return sigmaEQ+0.5*(1-exp(-Tdrop/TauI/(1+Pe(Qoil))))*GAMMAinf*R*T*log(1.005-DGamma(Qoil, Tdrop))
        else:
            return sigmaEQ+(1-exp(-Tdrop/TauI/(1+Pe(Qoil))))*GAMMAinf*R*T*log(1.003-DGamma(Qoil, Tdrop))

    # In[109]:= (* Capillary force at the interface *)
    def Pup(Qoil, Tdrop=Tdrop):
        return SIGMAio(Qoil, Tdrop)*(1/wdisp+1/H) # (* Acting upwards the outlet *)
    def Pdown(Qoil, Tdrop=Tdrop):
        return -SIGMAio(Qoil, Tdrop)*(2/wdisp+1/H) # (* Acting upwards the dispersed channel *)
    def FgammaIO(Qoil, Tdrop=Tdrop):
        return -(SIGMAio(Qoil, Tdrop)/wdisp)*H*wdisp

    # (* The sum of two *)

    # In[112]:=
    def Uc(Qoil):
        return ((Qoil/H)/HF(Qoil))/(wn-wjet0solF(Qoil))
    def Ud(Qoil):
        return (Qw/H)/HF(Qoil)/wjet0solF(Qoil)
    def CaNC(Qoil, Tdrop=Tdrop):
        return (etaoNN(Qoil)*Uc(Qoil))/SIGMAio(Qoil, Tdrop) # (*(Qw/Qoil)^0.333*)
    def CaND(Qoil, Tdrop=Tdrop):
        return (etaw(Qoil)*Ud(Qoil))/SIGMAio(Qoil, Tdrop)*(Qoil/Qw) ** (1/3)
    def Ljet(Qoil, Tdrop=Tdrop):
        return (etaw(Qoil)/SIGMAio(Qoil, Tdrop))*(8/pi/H/(CaNC(Qoil, Tdrop)+CaND(Qoil, Tdrop)))*(Qw*Qoil/2) ** 0.5

    def We(Qoil, Tdrop=Tdrop):
        return RhoO*Uc(Qoil) ** 2 * (1+1/CaND(Qoil) ** (1/3)) ** 2 * (wout*CaNC(Qoil) ** (2/3)) ** 2/SIGMAio(Qoil, Tdrop) # (*Weber number*)

    # In[122]:= (* LOSS OF ENERGY/FORCE AS A RESULT OF VISCOSITY\[Equal]> DRAG COEFFICIENT; Dv=8*pi/(Re*Log[7.4/Re])?????  *)

    # (*Drag coefficient & DISSIPATION *)
    # (*Reyn[Qoil_]:=RhoW*(Qoil/H/wcont)*wn/etaoNN[Qoil]*)
    def Reyn(Qoil):
        return RhoW*(Qoil/H/wcont)*wn/etaoNN(Qoil)
    def xi(Qoil):
        z1 = 2.345 * n / (2.423 * n + 0.918)
        return 24/Reyn(Qoil)*(1+0.15*(Reyn(Qoil)) ** z1) # (*+Del(Qoil)+Drel*)

    def CAF(Qoil, Tdrop=Tdrop):
        return 1.34 * CaNC(Qoil, Tdrop) ** (2/3) / CaND(Qoil, Tdrop) ** (2/3) * (1 + 3.35 * CaND(Qoil, Tdrop) ** (1/3)) ** 2

    def VcontSq(Qoil, Tdrop=Tdrop):
        return (Qoil/H/HF(Qoil)/(wn-wjet0solF(Qoil))) ** 2 * (wcont+Ln)/(Ln+Ljet(Qoil, Tdrop))+(Qoil/H/HF(Qoil)/(wn-wjet0solF(Qoil))) ** 2 * (wn/wout) ** 2 * (Ljet(Qoil, Tdrop)+Ldrop(Tdrop))/(Ln+Ljet(Qoil, Tdrop))

    def Ffrict(Qoil, Tdrop=Tdrop):
        return CAF(Qoil, Tdrop) * xi(Qoil) * etaoNN(Qoil) * Tdrop * (VcontSq(Qoil, Tdrop)) * (Ln/wn) ** 2.25

    # In[128]:= (* SHEAR FORCE-RELATED TERM *)
    def Fshear(Qoil, Tdrop=Tdrop):
        return etaoNN(Qoil) * (Qoil/(wn-wjet0solF(Qoil)) ** 2) * (2*Ljet(Qoil, Tdrop)+L(Qoil, Tdrop))
    # (*Plot[Fshear[Qw,Kvisc,Qoil, etaw, wn, Ln,H,wcont,wdisp,n,Tdrop],{Tdrop,10^-6,0.001}]*)

    # In[129]:= (* RESISTANCE TO OIL FLOW *)
    def Fresist(Qoil, Tdrop=Tdrop):
        return etaoNN(Qoil) * (Qoil/(wn-wjet0solF(Qoil)) ** 3) * (2*Ljet(Qoil, Tdrop)+L(Qoil, Tdrop)) ** 2
    # (*Plot[Fresist[Qw, etaw,Kvisc, Qoil, wn, Ln,H, wcont,wdisp,n,Tdrop],{Tdrop,10^-6,0.001}]*)

    # In[134]:=
    def LHS(Qoil, Tdrop=Tdrop):
        return Fshear(Qoil, Tdrop) + Fresist(Qoil, Tdrop) - Ffrict(Qoil, Tdrop)
    def RHS(Qoil, Tdrop=Tdrop):
        return FgammaIO(Qoil, Tdrop)

    # (* TOLIAU JAU SEKA SPRENDIMAS Tdrop VERTĖMS RASTI SU GAUTOMIS wjet0sol (ir Qoil) VERTĖMIS; FindRoot[LHS-RHS\[Equal]0,{Tdrop,10^-4,0.04}] *)

    # In[407]:= (* "FindRoot" ribos- apatinė ir viršutinė, kad nereiktų įvedinėti atsitiktinių skaičių *)
    LIM1 = (QoilStart / QoilEnd) * (etaw(QoilStart) * Ud(QoilStart)/sigmaEQ) ** (1/3)
    LIM2 = (etaw(QoilStart) * Ud(QoilStart)/sigmaEQ) ** (1/3)

    # Out[407]= 0.0102476
    # Out[408]= 0.0683174
    if debug:
        print(LIM1)
        print(LIM2)

    # In[409]:=
    def LHS_RHS_diff(Tdrop_sol, Qoil):
        return float(LHS(Qoil, Tdrop_sol).evalf() - RHS(Qoil, Tdrop_sol).evalf())
    data22_x = inclusive_range(QoilStart, QoilEnd, QoilStep)
    data22_y = []
    for Qoil in data22_x:
        try:
            data22_y.append(newton(LHS_RHS_diff, x0=LIM1, x1=LIM2, args=(Qoil,)))
        except RuntimeError:
            data22_y.append(numpy.nan)
        if reporter:
            reporter.ping()
    data22_y = numpy.array(data22_y)
    if debug:
        print(data22_x)
        print(data22_y)

    # In[137]:=
    # Out[137]= 0.034013
    if debug:
        print(SIGMAio(QoilStart, 0.0627).evalf())

    # In[153]:=
    data2_x = data22_x / (2.78 * 10 ** -13)
    data2_y = Qw/(10 ** -15) * data22_y
    if debug:
        print(data2_x)
        print(data2_y)
    # Out[153]= {{90.,4003.19},{190.,3704.71},{290.,3413.34},{390.,3156.99},{490.,2940.5},{590.,2758.15}}

    # In[139]:=
    # numpy.savetxt( '/dev/stdout', numpy.transpose(numpy.vstack((data2_x, data2_y))), delimiter="\t" )

    # In[142]:=
    def Vd(Qoil):
        return wn**3 / 0.7*(wn/Ln)**2.5*(etaw(Qoil,wjet0sol)*Qw/wdisp/etaoNN(Qoil,wjet0sol)/Qoil*wcont)**(2/3)/(1+(etaw(Qoil,wjet0sol)*Qw/wdisp/H/sigmaEQ)**(1/3))**2

    return numpy.transpose(numpy.vstack((data2_x, data2_y)))
