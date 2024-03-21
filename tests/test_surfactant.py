from surfactant.surfactant import calculate
import numpy
import unittest

Q_SI_multiplier = 2.78 * 10 ** -13
micrometre = 10e-6

class TestStringMethods(unittest.TestCase):

    # This test used to pass in 37fdcd25ae0582d80540af974dbdefa56c7810e6 and earlier
    def test_ionic(self):
        target = numpy.array([[90, 4003.19], [190, 3704.71], [290, 3413.34], [390, 3156.99], [490, 2940.5], [590, 2758.15]])
        result = calculate(is_ionic=True)
        self.assertLess(abs((target - result).min()), 20)

    def test_newtonian_ionic_0_6_perc(self):
        target = numpy.array([[50, 385.954], [100, 269.174], [150, 213.79], [200, 179.768], [250, 156.414], [300, 139.246]])
        result = calculate(is_ionic=True, Qw=50*Q_SI_multiplier, QoilStart=50*Q_SI_multiplier, QoilEnd=300*Q_SI_multiplier, QoilStep=50*Q_SI_multiplier, B1=1, B2=1, wn=37*micrometre, Ln=13*micrometre, H=18*micrometre, wcont=120*micrometre, wdisp=80*micrometre, wout=70*micrometre)
        self.assertLess(abs((target - result).min()), 20)
