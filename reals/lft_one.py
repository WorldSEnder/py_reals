import fractions
import math
from .defs import POWER_2, EXPONENT_2


class LFTOne():
    MODE_INCREASING = 0
    MODE_DECREASING = 1

    @staticmethod
    def is_plusminus_same_sign(a, b):
        """Checks that a*(-1) + b has the same sign as a*(1) + b"""
        # same sign <=> (b - a) * (b + a) > 0
        # (b - a) * (b + a) = b**2 - a**2
        # same sign <=> b**2 > a**2 <=> abs(b) > abs(a)
        return abs(b) > abs(a)

    @staticmethod
    def digit_from_lower_bound(a, b):
        """returns a digit suitable for extraction (given the interval is small enough),
        assuming the lower bound is given as the fraction a / b."""
        # the digit num is equivalent to the interval
        # [num - 1, num + 1] / POWER_2
        # we can extract a digit num if both lower and upper bound
        # fall into the interval [num - 1, num + 1] / POWER_2
        # so we just need to pick the biggest num such that
        # (num - 1) // POWER_2 <= lowerBound
        # EXCEPT that would result in num = POWER_2, in which case we return num - 1
        # TODO: this is biased against negative numbers
        # we might as well pick the smallest num such that
        # upperBound <= num // POWER_2
        # but alas, we would have to use a ceiling division there, so we dont
        num = 1 + (a << EXPONENT_2) // b
        if num == POWER_2:
            return num - 1
        return num

    @staticmethod
    def is_small_enough(a, b):
        """returns if the interval length given by the fraction (2 * a) / b (must be positive)
        is small enough to successfully extract a digit, assuming the bounds are contracting"""
        # small enough means that the interval_length <= 1 / POWER_2
        return a <= b >> (EXPONENT_2 + 1)

    @classmethod
    def identity(cls):
        return cls(1, 0, 0, 1)

    @classmethod
    def digit(cls, num):
        assert -POWER_2 < num < POWER_2
        return cls(1, 0, num, POWER_2)

    @classmethod
    def from_fraction(cls, frac):
        return cls(frac.numerator, 0, 0, frac.denominator)

    def __init__(self, a, b, c, d):
        self._matrix = [a, b, c, d]
        self._calculateCharacteristics()

    def clone(self):
        [a, b, c, d] = self._matrix
        return LFTOne(a, b, c, d)

    def __str__(self):
        [a, b, c, d] = self._matrix
        return "[{a}\t{c}\n{b}\t{d}]".format(a=a, b=b, c=c, d=d)

    def _calculateCharacteristics(self):
        [a, b, c, d] = self._matrix
        self._lft_type = mode = LFTOne.MODE_INCREASING if self._determineMP else LFTOne.MODE_DECREASING
        self._interval_length_num, self._interval_length_denom = {
            LFTOne.MODE_INCREASING: lambda: (a * d - c * b, d ** 2 - b ** 2),
            LFTOne.MODE_DECREASING: lambda: (c * b - a * d, d ** 2 - b ** 2),
        }[mode]()
        self._lowest_bound_num, self._lowest_bound_denom = {
            LFTOne.MODE_INCREASING: lambda: (c - a, d - b),
            LFTOne.MODE_DECREASING: lambda: (c + a, d + b),
        }[mode]()

    def times(self, other):
        # calculates self * other
        [a, b, c, d] = self._matrix
        [u, v, w, x] = other._matrix
        self._matrix[0] = a * u + c * v
        self._matrix[1] = b * u + d * v
        self._matrix[2] = a * w + c * x
        self._matrix[3] = b * w + d * x
        self._calculateCharacteristics()

    def timesdigit(self, digit):
        # special cases times(LFTOne.digit(digit))
        assert -POWER_2 < digit < POWER_2
        [a, b, c, d] = self._matrix
        w, exp = digit, EXPONENT_2
        self._matrix[0] = a
        self._matrix[1] = b
        self._matrix[2] = a * w + (c << exp)
        self._matrix[3] = b * w + (d << exp)
        # TODO: inline
        self._calculateCharacteristics()

    def invtimes(self, other):
        # calculates inv(other) * self
        [u, v, w, x] = self._matrix
        [d, b, c, a] = other._matrix
        b = -b
        c = -c
        self._matrix[0] = a * u + c * v
        self._matrix[1] = b * u + d * v
        self._matrix[2] = a * w + c * x
        self._matrix[3] = b * w + d * x
        self._calculateCharacteristics()

    def invtimesdigit(self, digit):
        # special cases invtimes(LFTOne.digit(digit))
        assert -POWER_2 < digit < POWER_2
        [u, v, w, x] = self._matrix
        c, exp = -digit, EXPONENT_2
        self._matrix[0] = (u << exp) + c * v
        self._matrix[1] = v
        self._matrix[2] = (w << exp) + c * x
        self._matrix[3] = x
        # TODO: inline
        self._calculateCharacteristics()

    @property
    def _determineMP(self):
        [a, b, c, d] = self._matrix
        # L(x) = (ax + c) / (bx + d)
        # L(-1) < L(1) ~~
        # (c - a) * (d + b) < (c + a) * (d - b)
        # 2 * c * b - 2 * a * d < 0
        # c * b - a * d
        return c * b < a * d

    @property
    def lft_type(self):
        """determines if the LFT is increasing or decreasing"""
        return self._lft_type

    def normalize(self):
        [a, b, c, d] = self._matrix
        ab = math.gcd(a, b)
        cd = math.gcd(c, d)
        abcd = math.gcd(ab, cd)
        gcd = max(1, abcd)
        self._matrix[:] = [a // gcd, b // gcd, c // gcd, d // gcd]
        # TODO: inline
        self._calculateCharacteristics()

    @property
    def next_index_to_pull(self):
        # assert self.is_contracting
        # L(1) - L(-1) = (c + a) / (d + b) - (c - a) / (d - b)
        #              = [(c + a) * (d - b) - (c - a) * (d + b)] / (d - b) * (d + b)
        #              = 2 * (a * d - c * b) / (d * d - b * b)
        is_small_enough = LFTOne.is_small_enough(self._interval_length_num, self._interval_length_denom)
        return None if is_small_enough else 0

    def extract(self):
        assert self.next_index_to_pull is None
        extracted_digit = LFTOne.digit_from_lower_bound(self._lowest_bound_num, self._lowest_bound_denom)
        assert -POWER_2 < extracted_digit < POWER_2
        self.invtimesdigit(extracted_digit)
        # assert self.is_contracting
        return extracted_digit

    @property
    def is_bounded(self):
        [_a, b, _c, d] = self._matrix
        return LFTOne.is_plusminus_same_sign(b, d)

    @property
    def bounds(self):
        assert self.is_bounded
        [a, b, c, d] = self._matrix
        at_m1 = fractions.Fraction(c - a, d - b)
        at_p1 = fractions.Fraction(c + a, d + b)
        return at_m1, at_p1

    @property
    def _determinant(self):
        [a, b, c, d] = self._matrix
        return a * d - b * c

    @property
    def _signature(self):
        [_a, b, _c, d] = self._matrix
        return d ** 2 - b ** 2

    @property
    def is_contracting(self):
        if not self.is_bounded:
            raise Exception("Not bounded, will never be contracting")
        bm1, bp1 = self.bounds
        if abs(bm1) <= 1 and abs(bp1) <= 1:
            return True
        if abs(bm1) > 1 and abs(bp1) > 1 and bm1 * bp1 > 0:
            raise Exception("interval outside [-1, 1], will never be contracting")
        return False

    @property
    def interval_length(self):
        return 2 * fractions.Fraction(abs(self._determinant), self._signature)

__all__ = ["LFTOne"]
