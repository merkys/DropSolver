from dropsolver import calculate
import numpy
import unittest

Q_SI_multiplier = 2.78 * 10 ** -13
micrometre = 1e-6

class TestSurfactant(unittest.TestCase):

    def test_newtonian_ionic_0_6_perc(self):
        target = numpy.array([[50, 74.8777], [100, 56.3615], [150, 44.3066], [200, 36.4674], [250, 31.0657], [300, 27.146]])
        result = calculate(is_ionic=True, B2=1, Qw=50*Q_SI_multiplier, QoilStart=50*Q_SI_multiplier, QoilEnd=300*Q_SI_multiplier, QoilStep=50*Q_SI_multiplier, wn=37*micrometre, Ln=13*micrometre, H=18*micrometre, wcont=120*micrometre, wdisp=80*micrometre, wout=70*micrometre)
        self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)

    def test_newtonian_nonionic_0_6_perc(self):
        target = numpy.array([[50, 74.8777], [100, 56.3615], [150, 44.3066], [200, 36.4674], [250, 31.0657], [300, 27.146]])
        result = calculate(is_ionic=False, B2=1, Qw=50*Q_SI_multiplier, QoilStart=50*Q_SI_multiplier, QoilEnd=300*Q_SI_multiplier, QoilStep=50*Q_SI_multiplier, wn=37*micrometre, Ln=13*micrometre, H=18*micrometre, wcont=120*micrometre, wdisp=80*micrometre, wout=70*micrometre)
        self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)

    def test_non_newtonian_nonionic_0_6_perc(self):
        target = numpy.array([[50, 134.08], [100, 95.4445], [150, 77.2938], [200, 65.8209], [250, 57.5501], [300, 51.1632]])
        result = calculate(is_ionic=False, Kd=0.381, etaINF1=0.00375, B1=4.691, B2=1, p=0.529, wn=37*micrometre, Ln=13*micrometre, H=18*micrometre, wcont=120*micrometre, wdisp=80*micrometre, wout=70*micrometre, Qw=50*Q_SI_multiplier, QoilStart=50*Q_SI_multiplier, QoilEnd=300*Q_SI_multiplier, QoilStep=50*Q_SI_multiplier)
        self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)

    def test_newtonian_ionic_5_perc(self):
        target = numpy.array([[90, 2279.93], [190, 2546.71], [290, 2402.18], [390, 2190.66], [490, 1985.67], [590, 1802.76]])
        result = calculate(is_ionic=True, omega=0.05)
        self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)

    def test_newtonian_ionic_6_perc(self):
        target = numpy.array([[90, 2168.96], [190, 2495.85], [290, 2380.01], [390, 2180.28], [490, 1980.41], [590, 1799.91]])
        result = calculate(is_ionic=True, omega=0.06)
        self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)

    def test_newtonian_nonionic_6_perc(self):
        target = numpy.array([[90, 2446.99], [190, 2606.96], [290, 2425.79], [390, 2201.11], [490, 1990.78], [590, 1805.47]])
        result = calculate(is_ionic=False, omega=0.06)
        self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)

    # ~ def test_non_newtonian_nonionic(self):
        # ~ target = numpy.array([[100, 72.5572], [200, 52.1084], [300, 36.1085], [400, 25.5143], [500, 18.6754], [600, 14.186], [700, 11.1508], [800, 9.032], [900, 7.5064], [1000, 6.3763], [1100, 5.5175], [1200, 4.8496], [1300, 4.3197], [1400, 3.8917], [1500, 3.5403]])
        # ~ result = calculate(is_ionic=False, Kd=1.276, etaINF1=0.00532, B1=4.578, p=0.672, EtaZero=0.006, EtaInf=0.016, B2=0.0001, n=1.48, wn=40*micrometre, Ln=10*micrometre, H=28*micrometre, wcont=30*micrometre, wdisp=40*micrometre, wout=70*micrometre, sigmaEQ=0.0375, omega=0.026, Qw=100*Q_SI_multiplier, QoilStart=100*Q_SI_multiplier, QoilEnd=1500*Q_SI_multiplier, QoilStep=100*Q_SI_multiplier)
        # ~ print(result)
        # ~ self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)
