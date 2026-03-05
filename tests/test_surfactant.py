from dropsolver import calculate
import numpy
import sys
import unittest

Q_SI_multiplier = 2.78 * 10 ** -13
micrometre = 1e-6

class TestSurfactant(unittest.TestCase):

    def test_newtonian_ionic_0_6_perc(self):
        target = numpy.array([[50, 2062.1939], [100, 1646.4591], [150, 1337.4537], [200, 1123.3634], [250, 968.5950], [300, 852.1015]])
        result = calculate(is_ionic=True, B2=1, Qw=50*Q_SI_multiplier, QoilStart=50*Q_SI_multiplier, QoilEnd=300*Q_SI_multiplier, QoilStep=50*Q_SI_multiplier, wn=37*micrometre, Ln=13*micrometre, H=18*micrometre, wcont=120*micrometre, wdisp=80*micrometre, wout=70*micrometre)
        self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)

    def test_newtonian_nonionic_0_6_perc(self):
        target = numpy.array([[50, 2147.7976], [100, 1674.5007], [150, 1348.7122], [200, 1127.4445], [250, 969.8363], [300, 852.4417]])
        result = calculate(is_ionic=False, B2=1, Qw=50*Q_SI_multiplier, QoilStart=50*Q_SI_multiplier, QoilEnd=300*Q_SI_multiplier, QoilStep=50*Q_SI_multiplier, wn=37*micrometre, Ln=13*micrometre, H=18*micrometre, wcont=120*micrometre, wdisp=80*micrometre, wout=70*micrometre)
        self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)

    def test_non_newtonian_nonionic_0_6_perc(self):
        target = numpy.array([[50, 2706.1604], [100, 1669.6661], [150, 1285.9275], [200, 1075.5646], [250, 937.1200], [300, 836.1838]])
        result = calculate(is_ionic=False, Kd=0.381, etaINF1=0.00375, B1=4.691, B2=1, p=0.529, wn=37*micrometre, Ln=13*micrometre, H=18*micrometre, wcont=120*micrometre, wdisp=80*micrometre, wout=70*micrometre, Qw=50*Q_SI_multiplier, QoilStart=50*Q_SI_multiplier, QoilEnd=300*Q_SI_multiplier, QoilStep=50*Q_SI_multiplier)
        self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)

    def test_newtonian_ionic_5_perc(self):
        target = numpy.array([[90, 6540.0668], [190, 9227.8656], [290, 9741.0197], [390, 9337.6029], [490, 8616.6547], [590, 7838.7299]])
        result = calculate(is_ionic=True, omega=0.05)
        self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)

    def test_newtonian_ionic_6_perc(self):
        target = numpy.array([[90, 5454.1951], [190, 8074.2408], [290, 8867.4801], [390, 8759.1669], [490, 8256.4023], [590, 7619.1919]])
        result = calculate(is_ionic=True, omega=0.06)
        self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)

    def test_newtonian_nonionic_6_perc(self):
        target = numpy.array([[90, 11290.8562], [190, 12483.4485], [290, 11529.4405], [390, 10278.0363], [490, 9115.0400], [590, 8110.3990]])
        result = calculate(is_ionic=False, omega=0.06)
        self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)

    # ~ def test_non_newtonian_nonionic(self):
        # ~ target = numpy.array([[100, 72.5572], [200, 52.1084], [300, 36.1085], [400, 25.5143], [500, 18.6754], [600, 14.186], [700, 11.1508], [800, 9.032], [900, 7.5064], [1000, 6.3763], [1100, 5.5175], [1200, 4.8496], [1300, 4.3197], [1400, 3.8917], [1500, 3.5403]])
        # ~ result = calculate(is_ionic=False, Kd=1.276, etaINF1=0.00532, B1=4.578, p=0.672, EtaZero=0.006, EtaInf=0.016, B2=0.0001, n=1.48, wn=40*micrometre, Ln=10*micrometre, H=28*micrometre, wcont=30*micrometre, wdisp=40*micrometre, wout=70*micrometre, sigmaEQ=0.0375, omega=0.026, Qw=100*Q_SI_multiplier, QoilStart=100*Q_SI_multiplier, QoilEnd=1500*Q_SI_multiplier, QoilStep=100*Q_SI_multiplier)
        # ~ print(result)
        # ~ self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)
