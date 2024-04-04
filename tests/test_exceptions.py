from surfactant.surfactant import calculate
import numpy
import unittest

Q_SI_multiplier = 2.78 * 10 ** -13
micrometre = 1e-6

class TestExceptions(unittest.TestCase):

    # From 'notebooks/2024-03-29/Išimtiniai atvejai.txt', first example
    def test_1(self):
        result = calculate(Kd=0.0236, etaINF1=0.0236, B1=1, B2=0.01, p=0, omega=0.01, Qw=50*Q_SI_multiplier, QoilStart=50*Q_SI_multiplier, QoilEnd=351*Q_SI_multiplier, QoilStep=50*Q_SI_multiplier, wn=25*micrometre, Ln=15*micrometre, H=20*micrometre, wcont=30*micrometre, wdisp=40*micrometre, wout=70*micrometre, is_ionic=False)
        self.assertEqual(len(list(filter(numpy.isnan, result[:,1]))), 6)
