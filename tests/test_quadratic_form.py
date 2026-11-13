"""Independent probability checks and the experimental binding contract."""
import math
import unittest
from concurrent.futures import ThreadPoolExecutor

from fastlmmrslib import qf


class TestQuadraticForm(unittest.TestCase):
    def probability(self, statistic, coefficients, **kwargs):
        value, fault, trace = qf(statistic, coefficients, acc=1e-7, **kwargs)
        self.assertEqual(fault, 0, (value, fault, trace))
        self.assertIsNotNone(value)
        self.assertTrue(0 <= value <= 1)
        self.assertEqual(len(trace), 7)
        self.assertTrue(all(math.isfinite(x) for x in trace))
        return value

    def test_scaled_exponential(self):
        for scale in [0.5, 2.0]:
            for statistic in [0.0, 0.5, 3.0, 10.0]:
                with self.subTest(scale=scale, statistic=statistic):
                    value = self.probability(statistic, [scale], dof=[2])
                    self.assertAlmostEqual(value, math.exp(-statistic/(2*scale)), delta=1e-7)

    def test_unequal_positive_weights(self):
        for a in [1e-3, 1e-8]:
            with self.subTest(a=a):
                value = self.probability(20, [a, a, 1.0, 1.0])
                expected = math.exp(-10)/(1-a)
                self.assertAlmostEqual(value, expected, delta=1e-7)

    def test_signed_exponentials(self):
        for statistic in [-10.0, 0.0, 1.0, 10.0]:
            with self.subTest(statistic=statistic):
                value = self.probability(statistic, [1.0, -2.0], dof=[2,2])
                expected = (math.exp(-statistic/2)/3 if statistic >= 0
                            else 1-2*math.exp(statistic/4)/3)
                self.assertAlmostEqual(value, expected, delta=1e-7)

    def test_normal_component(self):
        for statistic in [-2.0, 0.0, 2.0, 12.0]:
            with self.subTest(statistic=statistic):
                value = self.probability(statistic, [0.0], sigma=2.0)
                expected = math.erfc(statistic/(2*math.sqrt(2)))/2
                self.assertAlmostEqual(value, expected, delta=1e-7)

    def test_noncentral_squared_normal(self):
        # a*(Z+sqrt(lambda))**2: sum the two elementary normal tails.
        for a in [0.5, 2.0]:
            for noncentrality in [0.5, 5.0]:
                for statistic in [1.0, 9.0]:
                    with self.subTest(a=a, nc=noncentrality, statistic=statistic):
                        root = math.sqrt(statistic/a)
                        shift = math.sqrt(noncentrality)
                        expected = (math.erfc((root-shift)/math.sqrt(2))
                                    +math.erfc((root+shift)/math.sqrt(2)))/2
                        value = self.probability(statistic, [a], noncentrality=[noncentrality])
                        self.assertAlmostEqual(value, expected, delta=1e-7)

    def test_failures_have_no_probability(self):
        cases = [
            (1.0, [1.0], dict(dof=[]), 3),
            (1.0, [1.0], dict(noncentrality=[-1.0]), 3),
            (math.nan, [1.0], {}, 3),
            (1.0, [math.inf], {}, 3),
            (1.0, [1.0], dict(acc=0.0), 3),
            (1.0, [1.0], dict(lim=1), 4),
            (11.498286033063403, [1.4256991860196153],
             dict(dof=[0], noncentrality=[1.2791211646129712]), 1),
        ]
        for statistic, coefficients, kwargs, expected_fault in cases:
            with self.subTest(kwargs=kwargs):
                probability, fault, trace = qf(statistic, coefficients, **kwargs)
                self.assertIsNone(probability)
                self.assertEqual(fault, expected_fault)
                self.assertEqual(len(trace), 7)

    def test_degenerate_boundary_convention(self):
        self.assertEqual(self.probability(0.0, [0.0]), 1.0)
        self.assertEqual(self.probability(1.0, [0.0]), 0.0)

    def test_threaded_binding_calls(self):
        def calculate(index):
            scale = 0.5 + index % 7
            return qf(2+scale, [0.1*scale, -0.3, 0.5, 0.8],
                      dof=[1,2,3,4], noncentrality=[0.0,0.2*scale,0.0,1.0], sigma=0.5)
        expected = list(map(calculate, range(48)))
        with ThreadPoolExecutor(max_workers=4) as pool:
            self.assertEqual(list(pool.map(calculate,range(48))), expected)


if __name__ == '__main__':
    unittest.main()
