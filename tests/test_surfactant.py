from dropsolver import calculate
import numpy
import unittest

Q_SI_multiplier = 2.78 * 10 ** -13
micrometre = 1e-6

class TestSurfactant(unittest.TestCase):

    # From 'notebooks/2024-03-22/HFE7500-IonicSurfactant (0.6perc) [2].txt'
    def test_newtonian_ionic_0_6_perc(self):
        target = numpy.array([[50, 107.429], [100, 76.5061], [150, 59.4823], [200, 48.9298], [250, 41.8043], [300, 36.6827]])
        result = calculate(is_ionic=True, B2=1, Qw=50*Q_SI_multiplier, QoilStart=50*Q_SI_multiplier, QoilEnd=300*Q_SI_multiplier, QoilStep=50*Q_SI_multiplier, wn=37*micrometre, Ln=13*micrometre, H=18*micrometre, wcont=120*micrometre, wdisp=80*micrometre, wout=70*micrometre)
        self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)

    # From 'notebooks/2024-03-22/HFE7500-NonionicSurfactant (0.6perc) [2].txt'
    def test_newtonian_nonionic_0_6_perc(self):
        target = numpy.array([[50, 107.429], [100, 76.506], [150, 59.4823], [200, 48.9298], [250, 41.8043], [300, 36.6827]])
        result = calculate(is_ionic=False, B2=1, Qw=50*Q_SI_multiplier, QoilStart=50*Q_SI_multiplier, QoilEnd=300*Q_SI_multiplier, QoilStep=50*Q_SI_multiplier, wn=37*micrometre, Ln=13*micrometre, H=18*micrometre, wcont=120*micrometre, wdisp=80*micrometre, wout=70*micrometre)
        self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)

    # From 'notebooks/2024-03-22/HFE7500-NonionicSurfactant (0.6perc) non-Newtonian disp.phase [2].txt'
    def test_non_test_newtonian_nonionic_0_6_perc(self):
        target = numpy.array([[50, 1146.36], [100, 290.37], [150, 182.068], [200, 136.282], [250, 110.076], [300, 92.742]])
        result = calculate(is_ionic=False, Kd=0.381, etaINF1=0.00375, B1=4.691, B2=1, p=0.529, wn=37*micrometre, Ln=13*micrometre, H=18*micrometre, wcont=120*micrometre, wdisp=80*micrometre, wout=70*micrometre, Qw=50*Q_SI_multiplier, QoilStart=50*Q_SI_multiplier, QoilEnd=300*Q_SI_multiplier, QoilStep=50*Q_SI_multiplier)
        self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)
