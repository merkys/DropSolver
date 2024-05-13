from dropsolver import calculate
import numpy
import unittest

Q_SI_multiplier = 2.78 * 10 ** -13
micrometre = 1e-6

class TestExceptions(unittest.TestCase):

    # From 'notebooks/2024-03-29/Išimtiniai atvejai.txt', example [1]
    def test_1(self):
        result = calculate(Kd=0.0236, etaINF1=0.0236, B1=1, B2=0.01, p=0, omega=0.01, Qw=50*Q_SI_multiplier, QoilStart=50*Q_SI_multiplier, QoilEnd=350*Q_SI_multiplier, QoilStep=50*Q_SI_multiplier, wn=25*micrometre, Ln=15*micrometre, H=20*micrometre, wcont=30*micrometre, wdisp=40*micrometre, wout=70*micrometre, is_ionic=False)
        self.assertEqual(len(list(filter(numpy.isnan, result[:,1]))), 6)

    # From 'notebooks/2024-03-29/Išimtiniai atvejai.txt', example [2]
    def test_2(self):
        result = calculate(B1=1, B2=0.01, p=0, omega=0.01, Qw=300*Q_SI_multiplier, QoilStart=50*Q_SI_multiplier, QoilEnd=350*Q_SI_multiplier, QoilStep=50*Q_SI_multiplier, wn=25*micrometre, Ln=15*micrometre, H=20*micrometre, wcont=30*micrometre, wdisp=40*micrometre, wout=70*micrometre, is_ionic=False)
        self.assertEqual(len(list(filter(numpy.isnan, result[:,1]))), 7)

    # From 'notebooks/2024-03-29/Išimtiniai atvejai.txt', example [3]
    def test_3(self):
        target = numpy.array([[50, 0.152369], [100, 0.126003], [150, 0.113378], [200, 0.105204], [250,0.0992039], [300, 0.0944858]])
        result = calculate(B1=1, B2=0.78, p=1, Kvisc=0.0294, EtaZero=0.0294, EtaInf=0.0326, n=2, omega=0.01, Qw=50*Q_SI_multiplier, QoilStart=50*Q_SI_multiplier, QoilEnd=300*Q_SI_multiplier, QoilStep=50*Q_SI_multiplier, wn=20*micrometre, Ln=20*micrometre, H=20*micrometre, wcont=30*micrometre, wdisp=40*micrometre, wout=70*micrometre, is_ionic=False)
        self.assertLess(max(abs(target[:,1] - result[:,1])/target[:,1]), 0.1)
