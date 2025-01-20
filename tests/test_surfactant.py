from dropsolver import calculate
import numpy
import unittest

Q_SI_multiplier = 2.78 * 10 ** -13
micrometre = 1e-6

class TestSurfactant(unittest.TestCase):

    # From 'notebooks/2024-03-22/HFE7500-IonicSurfactant (0.6perc) [2].txt'
    def test_newtonian_ionic_0_6_perc(self):
        target = numpy.array([[50, 74.8777], [100, 56.3615], [150, 44.3066], [200, 36.4674], [250, 31.0657], [300, 27.146]])
        result = calculate(is_ionic=True, B2=1, Qw=50*Q_SI_multiplier, QoilStart=50*Q_SI_multiplier, QoilEnd=300*Q_SI_multiplier, QoilStep=50*Q_SI_multiplier, wn=37*micrometre, Ln=13*micrometre, H=18*micrometre, wcont=120*micrometre, wdisp=80*micrometre, wout=70*micrometre)
        self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)

    # From 'notebooks/2024-03-22/HFE7500-NonionicSurfactant (0.6perc) [2].txt'
    def test_newtonian_nonionic_0_6_perc(self):
        target = numpy.array([[50, 74.8777], [100, 56.3615], [150, 44.3066], [200, 36.4674], [250, 31.0657], [300, 27.146]])
        result = calculate(is_ionic=False, B2=1, Qw=50*Q_SI_multiplier, QoilStart=50*Q_SI_multiplier, QoilEnd=300*Q_SI_multiplier, QoilStep=50*Q_SI_multiplier, wn=37*micrometre, Ln=13*micrometre, H=18*micrometre, wcont=120*micrometre, wdisp=80*micrometre, wout=70*micrometre)
        self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)

    # From 'notebooks/2024-03-22/HFE7500-NonionicSurfactant (0.6perc) non-Newtonian disp.phase [2].txt'
    def test_non_newtonian_nonionic_0_6_perc(self):
        target = numpy.array([[50, 134.08], [100, 95.4445], [150, 77.2938], [200, 65.8209], [250, 57.5501], [300, 51.1632]])
        result = calculate(is_ionic=False, Kd=0.381, etaINF1=0.00375, B1=4.691, B2=1, p=0.529, wn=37*micrometre, Ln=13*micrometre, H=18*micrometre, wcont=120*micrometre, wdisp=80*micrometre, wout=70*micrometre, Qw=50*Q_SI_multiplier, QoilStart=50*Q_SI_multiplier, QoilEnd=300*Q_SI_multiplier, QoilStep=50*Q_SI_multiplier)
        self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)

    def test_newtonian_ionic_5_perc(self):
        target = numpy.array([[90, 2279.93], [190, 2546.71], [290, 2402.18], [390, 2190.66], [490, 1985.67], [590, 1802.76]])
        result = calculate(is_ionic=True, omega=0.05)
        self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)

    # From 'notebooks/2024-11-27/Ionic Surfactant 6perc. DefaultParameters.txt'
    def test_newtonian_ionic_6_perc(self):
        target = numpy.array([[90, 2168.96], [190, 2495.85], [290, 2380.01], [390, 2180.28], [490, 1980.41], [590, 1799.91]])
        result = calculate(is_ionic=True, omega=0.06)
        self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)

    # From 'notebooks/2024-11-27/NonIonic Surfactant 6perc. DefaultParameters.txt'
    def test_newtonian_nonionic_6_perc(self):
        target = numpy.array([[90, 2446.99], [190, 2606.96], [290, 2425.79], [390, 2201.11], [490, 1990.78], [590, 1805.47]])
        result = calculate(is_ionic=False, omega=0.06)
        self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)
