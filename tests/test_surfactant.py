from surfactant.surfactant import calculate
import numpy
import unittest

class TestStringMethods(unittest.TestCase):

    # This test used to pass in 37fdcd25ae0582d80540af974dbdefa56c7810e6 and earlier
    def test_ionic(self):
        target = numpy.array([[90, 4003.19], [190, 3704.71], [290, 3413.34], [390, 3156.99], [490, 2940.5], [590, 2758.15]])
        result = calculate(is_ionic=True)
        self.assertLess(abs((target - result).min()), 20)
