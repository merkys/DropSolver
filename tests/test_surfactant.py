from surfactant.surfactant import calculate
import numpy
import unittest

Q_SI_multiplier = 2.78 * 10 ** -13
micrometre = 1e-6

class TestStringMethods(unittest.TestCase):

    # This test used to pass in 37fdcd25ae0582d80540af974dbdefa56c7810e6 and earlier
    def test_ionic(self):
        target = numpy.array([[90, 4003.19], [190, 3704.71], [290, 3413.34], [390, 3156.99], [490, 2940.5], [590, 2758.15]])
        result = calculate(is_ionic=True)
        self.assertLess(abs((target - result).min()), 20)

    # From 'notebooks/2024-03-21/Newtonian HFE7500- IonicSurfactant (0.6perc).txt'
    # This test uses m = 5.6 and Enth = 0.06
    def test_newtonian_ionic_0_6_perc(self):
        target = numpy.array([[50, 385.954], [100, 269.174], [150, 213.79], [200, 179.768], [250, 156.414], [300, 139.246]])
        result = calculate(is_ionic=True, Qw=50*Q_SI_multiplier, QoilStart=50*Q_SI_multiplier, QoilEnd=300*Q_SI_multiplier, QoilStep=50*Q_SI_multiplier, B1=1, B2=1, wn=37*micrometre, Ln=13*micrometre, H=18*micrometre, wcont=120*micrometre, wdisp=80*micrometre, wout=70*micrometre)
        self.assertLess(abs((target - result).min()), 20)

    # From 'notebooks/2024-03-21/Newtonian HFE7500- NonionicSurfactant (0.6perc).txt'
    def test_newtonian_nonionic_0_6_perc(self):
        target = numpy.array([[50, 385.955], [100, 269.174], [150, 213.79], [200, 179.768], [250, 156.414], [300.,139.246]])
        result = calculate(is_ionic=False, Qw=50*Q_SI_multiplier, QoilStart=50*Q_SI_multiplier, QoilEnd=300*Q_SI_multiplier, QoilStep=50*Q_SI_multiplier, B1=1, B2=1, wn=37*micrometre, Ln=13*micrometre, H=18*micrometre, wcont=120*micrometre, wdisp=80*micrometre, wout=70*micrometre)
        self.assertLess(abs((target - result).min()), 20)

    # From 'notebooks/2024-03-22/HFE7500-NonionicSurfactant (0.6perc) non-Newtonian disp.phase [2].txt'
    def test_non_test_newtonian_nonionic_0_6_perc(self):
        target = numpy.array([[50, 1146.36}, [100, 290.37], [150, 182.068], [200, 136.282], [250, 110.076], [300, 92.742]])
        result = calculate(is_ionic=False, Kd=0.381, etaINF1=0.00375, B1=4.691, B2=1, p=0.529, wn=37*micrometre, Ln=13*micrometre, H=18*micrometre, wcont=120*micrometre, wdisp=80*micrometre, wout=70*micrometre, Qw=50*Q_SI_multiplier, QoilStart=50*Q_SI_multiplier, QoilEnd=301*Q_SI_multiplier, QoilStep=50*Q_SI_multiplier)
        self.assertLess(abs((target - result).min()), 20)
