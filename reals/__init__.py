import fractions
import decimal
import math

EXPONENT_2 = 64
POWER_2 = 2 ** EXPONENT_2


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
        num = 1 + a * POWER_2 // b
        if num == POWER_2:
            return num - 1
        return num

    @staticmethod
    def is_small_enough(a, b):
        """returns if the interval length given by the fraction a / b (must be positive)
        is small enough to successfully extract a digit, assuming the bounds are contracting"""
        # small enough means that the interval_length <= 2 / POWER_2
        return 2 * a * POWER_2 <= b

    def __init__(self, a, b, c, d):
        self._matrix = [a, b, c, d]

    def clone(self):
        [a, b, c, d] = self._matrix
        return LFTOne(a, b, c, d)

    def __str__(self):
        [a, b, c, d] = self._matrix
        return "[{a}\t{c}\n{b}\t{d}]".format(a=a, b=b, c=c, d=d)

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

    def times(self, other):
        # calculates self * other
        [a, b, c, d] = self._matrix
        [u, v, w, x] = other._matrix
        self._matrix[0] = a * u + c * v
        self._matrix[1] = b * u + d * v
        self._matrix[2] = a * w + c * x
        self._matrix[3] = b * w + d * x

    def timesdigit(self, digit):
        # special cases times(LFTOne.digit(digit))
        assert -POWER_2 < digit < POWER_2
        [a, b, c, d] = self._matrix
        w, x = digit, POWER_2
        self._matrix[0] = a
        self._matrix[1] = b
        self._matrix[2] = a * w + c * x
        self._matrix[3] = b * w + d * x

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

    def invtimesdigit(self, digit):
        # special cases invtimes(LFTOne.digit(digit))
        assert -POWER_2 < digit < POWER_2
        [u, v, w, x] = self._matrix
        c, a = -digit, POWER_2
        self._matrix[0] = a * u + c * v
        self._matrix[1] = v
        self._matrix[2] = a * w + c * x
        self._matrix[3] = x

    def normalize(self):
        [a, b, c, d] = self._matrix
        ab = math.gcd(a, b)
        cd = math.gcd(c, d)
        abcd = math.gcd(ab, cd)
        gcd = max(1, abcd)
        self._matrix[:] = [a // gcd, b // gcd, c // gcd, d // gcd]

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
        isIncreasing = self._determineMP
        return LFTOne.MODE_INCREASING if isIncreasing else LFTOne.MODE_DECREASING

    @property
    def next_index_to_pull(self):
        [a, b, c, d] = self._matrix
        # TODO: duplicated work here, when we also calculate this for lft_type
        assert self.is_contracting
        mode = self.lft_type
        # L(1) - L(-1) = (c + a) / (d + b) - (c - a) / (d - b)
        #              = [(c + a) * (d - b) - (c - a) * (d + b)] / (d - b) * (d + b)
        #              = 2 * (a * d - c * b) / (d * d - b * b)
        is_small_enough = {
            LFTOne.MODE_INCREASING: lambda: LFTOne.is_small_enough(a * d - c * b, self._signature),
            LFTOne.MODE_DECREASING: lambda: LFTOne.is_small_enough(c * b - a * d, self._signature),
        }[mode]()
        return None if is_small_enough else 0

    def extract(self):
        assert self.next_index_to_pull is None
        # after normalization, this guarantees d - b > 0
        [a, b, c, d] = self._matrix
        mode = self.lft_type
        extracted_digit = {
            LFTOne.MODE_INCREASING: lambda: LFTOne.digit_from_lower_bound(c - a, d - b),
            LFTOne.MODE_DECREASING: lambda: LFTOne.digit_from_lower_bound(c + a, d + b),
        }[mode]()
        assert -POWER_2 < extracted_digit < POWER_2
        self.invtimesdigit(extracted_digit)
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
        return self.is_bounded and all(abs(bound) <= 1 for bound in self.bounds)

    @property
    def interval_length(self):
        return 2 * fractions.Fraction(abs(self._determinant), self._signature)


class LFTTwo():
    # the mode determines which endpoint is the min and max of the output interval
    MODE_MM_PP = 0x03
    MODE_MP_PP = 0x13
    MODE_PM_PP = 0x23
    MODE_MM_PM = 0x02
    MODE_MP_PM = 0x12
    MODE_PP_PM = 0x32
    MODE_MM_MP = 0x01
    MODE_PM_MP = 0x21
    MODE_PP_MP = 0x31
    MODE_MP_MM = 0x10
    MODE_PM_MM = 0x20
    MODE_PP_MM = 0x30

    def __init__(self, a, b, c, d, e, f, g, h):
        self._matrix = [a, b, c, d, e, f, g, h]

    def clone(self):
        [a, b, c, d, e, f, g, h] = self._matrix
        return LFTTwo(a, b, c, d, e, f, g, h)

    def __str__(self):
        [a, b, c, d, e, f, g, h] = self._matrix
        return "[{a}\t{c}\t| {e}\t{g}\n{b}\t{d}\t| {f}\t{h}]".format(
            a=a, b=b, c=c, d=d, e=e, f=f, g=g, h=h)

    def timesX(self, other):
        [a, b, c, d, e, f, g, h] = self._matrix
        [u, v, w, x] = other._matrix
        self._matrix[0] = a * u + c * v
        self._matrix[1] = b * u + d * v
        self._matrix[2] = a * w + c * x
        self._matrix[3] = b * w + d * x
        self._matrix[4] = e * u + g * v
        self._matrix[5] = f * u + h * v
        self._matrix[6] = e * w + g * x
        self._matrix[7] = f * w + h * x

    def timesY(self, other):
        # suppose we have a flip operation that swaps X and Y
        # (i.e. it swaps [c, d] with [e, f])
        # then this is swap . timeX other . swap
        [a, b, c, d, e, f, g, h] = self._matrix
        [u, v, w, x] = other._matrix
        self._matrix[0] = a * u + e * v
        self._matrix[1] = b * u + f * v
        self._matrix[2] = c * u + g * v
        self._matrix[3] = d * u + h * v
        self._matrix[4] = a * w + e * x
        self._matrix[5] = b * w + f * x
        self._matrix[6] = c * w + g * x
        self._matrix[7] = d * w + h * x

    def timesDigitX(self, digit):
        [a, b, c, d, e, f, g, h] = self._matrix
        w, x = digit, POWER_2
        self._matrix[0] = a
        self._matrix[1] = b
        self._matrix[2] = a * w + c * x
        self._matrix[3] = b * w + d * x
        self._matrix[4] = e
        self._matrix[5] = f
        self._matrix[6] = e * w + g * x
        self._matrix[7] = f * w + h * x

    def timesDigitY(self, digit):
        [a, b, c, d, e, f, g, h] = self._matrix
        w, x = digit, POWER_2
        self._matrix[0] = a
        self._matrix[1] = b
        self._matrix[2] = c
        self._matrix[3] = d
        self._matrix[4] = a * w + e * x
        self._matrix[5] = b * w + f * x
        self._matrix[6] = c * w + g * x
        self._matrix[7] = d * w + h * x

    def invtimes(self, other):
        # calculates inv(other) * self
        [a, b, c, d, e, f, g, h] = self._matrix
        [w, u, v, x] = other._matrix
        u = -u
        v = -v
        self._matrix[0] = x * a + v * b
        self._matrix[1] = u * a + w * b
        self._matrix[2] = x * c + v * d
        self._matrix[3] = u * c + w * d
        self._matrix[4] = x * e + v * f
        self._matrix[5] = u * e + w * f
        self._matrix[6] = x * g + v * h
        self._matrix[7] = u * g + w * h

    def invtimesdigit(self, digit):
        [a, b, c, d, e, f, g, h] = self._matrix
        v, x = -digit, POWER_2
        self._matrix[0] = a * x + v * b
        self._matrix[1] = b
        self._matrix[2] = c * x + v * d
        self._matrix[3] = d
        self._matrix[4] = e * x + v * f
        self._matrix[5] = f
        self._matrix[6] = g * x + v * h
        self._matrix[7] = h

    @property
    def _determineXM(self):
        """determines MM < MP"""
        [a, b, c, d, e, f, g, h] = self._matrix
        # WHEN COMPARING FRACTIONS, KEEP IN MIND THAT is_bounded WILL guarantee THAT THE
        # SIGN OF THE DOMINATOR IS THE SAME FOR ALL POINTS!
        # L(x, y) = (a xy + c x + e y + g) / (b xy + d x + f y + h)
        # Lmx(t) = L(-1, t) = ((e - a) t + (g - c)) / ((f - b) t + (h - d))
        # Lmx(-1) < Lmx(1) ~~
        # ((g - c) - (e - a)) * ((h - d) + (f - b)) < ((g - c) + (e - a)) * ((h - d) - (f - b)) ~~
        # ((g - c) - (e - a)) * ((h - d) + (f - b)) < ((g - c) + (e - a)) * ((h - d) - (f - b)) ~~
        # 2 * (g - c) * (f - b) - 2 * (e - a) * (h - d) < 0 ~~
        # (g - c) * (f - b) - (e - a) * (h - d) < 0
        return (g - c) * (f - b) < (e - a) * (h - d)
        # Similarly (for later):
        # L(-1, 1) - L(-1, -1) = Lmx(1) < Lmx(-1) =
        # ((g - c) + (e - a)) / ((h - d) + (f - b)) - ((g - c) - (e - a)) / ((h - d) - (f - b)) =
        # 2 * [(e - a) * (h - d) - (g - c) * (f - b)] / ((h - d) ** 2 - (f - b) ** 2)

    @property
    def _determineXP(self):
        """determines PM < PP"""
        [a, b, c, d, e, f, g, h] = self._matrix
        # L(x, y) = (a xy + c x + e y + g) / (b xy + d x + f y + h) at x = 1
        # Lpx(t) = ((e + a) t + (g + c)) / ((f + b) t + (h + d))
        # Lpx(-1) < Lpx(1) ~~
        return (g + c) * (f + b) < (e + a) * (h + d)
        # Similarly (for later):
        # L(1, 1) - L(1, -1) = Lpx(1) < Lpx(-1) =
        # 2 * [(e + a) * (h + d) - (g + c) * (f + b)] / ((h + d) ** 2 - (f + b) ** 2)

    @property
    def _determineYM(self):
        """determines MM < PM"""
        [a, b, c, d, e, f, g, h] = self._matrix
        # L(x, y) = (a xy + c x + e y + g) / (b xy + d x + f y + h) at y = -1
        # Lmy(t) = ((c - a) t + (g - e)) / ((d - b) t + (h - f))
        # Lmy(-1) < Lmy(1) ~~
        return (g - e) * (d - b) < (c - a) * (h - f)
        # Similarly (for later):
        # L(1, -1) - L(-1, -1) = Lmy(1) < Lmy(-1) =
        # 2 * [(c - a) * (h - f) - (g - e) * (d - b)] / ((h - f) ** 2 - (d - b) ** 2)

    @property
    def _determineYP(self):
        """determines MP < PP"""
        [a, b, c, d, e, f, g, h] = self._matrix
        # L(x, y) = (a xy + c x + e y + g) / (b xy + d x + f y + h) at y = 1
        # Lpy(t) = ((c + a) t + (g + e)) / ((d + b) t + (h + f))
        # Lpy(-1) < Lpy(1) ~~
        return (g + e) * (d + b) < (c + a) * (h + f)
        # Similarly (for later):
        # L(1, 1) - L(-1, 1) = Lpy(1) < Lpy(-1) =
        # 2 * [(c + a) * (h + f) - (g + e) * (d + b)] / ((h + f) ** 2 - (d + b) ** 2)

    @property
    def _determineCrossMMPP(self):
        """determines MM < PP"""
        [a, b, c, d, e, f, g, h] = self._matrix
        # L(x, y) = (a xy + c x + e y + g) / (b xy + d x + f y + h)
        # Lmmpp(t) = L(t, t) = ((c + e) t + (g + a tt)) / ((d + f) t + (h + b tt))
        # L(-1, -1) < L(1, 1) ~~
        # Lmmpp(-1) < Lmmpp(1) ~~
        # ((g + a) - (c + e)) * ((h + b) + (d + f)) < ((g + a) + (c + e)) * ((h + b) - (d + f)) ~~
        # 2 * (g + a) * (d + f) - 2 * (c + e) * (h + b) < 0
        # (g + a) * (d + f) - (c + e) * (h + b) < 0
        return (g + a) * (d + f) < (c + e) * (h + b)
        # Similarly (for later):
        # L(1, 1) - L(-1, -1) =
        # 2 * [(c + e) * (h + b) - (g + a) * (d + f)] / ((h + b) ** 2 - (d + f) ** 2)

    @property
    def _determineCrossMPPM(self):
        """determines MP < PM"""
        [a, b, c, d, e, f, g, h] = self._matrix
        # L(x, y) = (a xy + c x + e y + g) / (b xy + d x + f y + h)
        # Lmppm(t) = L(t, -t) = ((c - e) t + (g - a tt)) / ((d - f) t + (h - b tt))
        # L(-1, 1) < L(1, -1) ~~
        # Lmppm(-1) < Lmppm(1) ~~
        # ((g - a) - (c - e)) * ((h - b) + (d - f)) < ((g - a) + (c - e)) * ((h - b) - (d - f)) ~~
        # 2 * (g - a) * (d - f) - 2 * (c - e) * (h - b) < 0 ~~
        # (g - a) * (d - f) - (c - e) * (h - b) < 0
        return (g - a) * (d - f) < (c - e) * (h - b)
        # Similarly (for later):
        # L(1, -1) - L(-1, 1) =
        # 2 * [(c - e) * (h - b) - (g - a) * (d - f)] / ((h - b) ** 2 - (d - f) ** 2)

    @property
    def lft_type(self):
        """determines if the LFT is increasing or decreasing"""
        isIncrAtXM = self._determineXM
        isIncrAtXP = self._determineXP
        if isIncrAtXP and isIncrAtXM:
            # L(-1, -1) <= L(-1, 1) and L(1, -1) <= L(1, 1)
            isIncrAtYM = self._determineYM
            isIncrAtYP = self._determineYP
            if isIncrAtYP and isIncrAtYM:
                # L(-1, -1) <= L(1, -1) and L(-1, 1) <= L(1, 1)
                return LFTTwo.MODE_MM_PP
            elif isIncrAtYP:  # and not isIncrAtYM
                # L(1, -1) <= L(-1, -1) and L(-1, 1) <= L(1, 1)
                return LFTTwo.MODE_PM_PP
            elif isIncrAtYM:  # and not isIncrAtYP
                # L(-1, -1) <= L(1, -1) and L(1, 1) <= L(-1, 1)
                return LFTTwo.MODE_MM_MP
            else:  # not isIncrAtYM and not isIncrAtYP
                # L(1, -1) <= L(-1, -1) and L(1, 1) <= L(-1, 1)
                return LFTTwo.MODE_PM_MP
        elif isIncrAtXP:  # and not isIncrAtXM
            # L(-1, 1) <= L(-1, -1) and L(1, -1) <= L(1, 1)
            isIncrCrossMMPP = self._determineCrossMMPP
            isIncrCrossMPPM = self._determineCrossMPPM
            if isIncrCrossMMPP and isIncrCrossMPPM:
                # L(-1, -1) <= L(1, 1) and L(-1, 1) <= L(1, -1)
                return LFTTwo.MODE_MP_PP
            elif isIncrCrossMMPP:  # and not isIncrCrossMPPM
                # L(-1, -1) <= L(1, 1) and L(1, -1) <= L(-1, 1)
                return LFTTwo.MODE_PM_PP
            elif isIncrCrossMPPM:  # and not isIncrCrossMMPP
                # L(1, 1) <= L(-1, -1) and L(-1, 1) <= L(1, -1)
                return LFTTwo.MODE_MP_MM
            else:  # not isIncrCrossMPPM and not isIncrCrossMMPP
                # L(1, 1) <= L(-1, -1) and L(1, -1) <= L(-1, 1)
                return LFTTwo.MODE_PM_MM
        elif isIncrAtXM:  # and not isIncrAtXP
            # L(-1, -1) <= L(-1, 1) and L(1, 1) <= L(1, -1)
            isIncrCrossMMPP = self._determineCrossMMPP
            isIncrCrossMPPM = self._determineCrossMPPM
            if isIncrCrossMMPP and isIncrCrossMPPM:
                # L(-1, -1) <= L(1, 1) and L(-1, 1) <= L(1, -1)
                return LFTTwo.MODE_MM_PM
            elif isIncrCrossMMPP:  # and not isIncrCrossMPPM
                # L(-1, -1) <= L(1, 1) and L(1, -1) <= L(-1, 1)
                return LFTTwo.MODE_MM_MP
            elif isIncrCrossMPPM:  # and not isIncrCrossMMPP
                # L(1, 1) <= L(-1, -1) and L(-1, 1) <= L(1, -1)
                return LFTTwo.MODE_PP_PM
            else:  # not isIncrCrossMPPM and not isIncrCrossMMPP
                # L(1, 1) <= L(-1, -1) and L(1, -1) <= L(-1, 1)
                return LFTTwo.MODE_PP_MP
        else:  # not isIncrAtXM and not isIncrAtXP
            # L(-1, 1) <= L(-1, -1) and L(1, 1) <= L(1, -1)
            isIncrAtYM = self._determineYM
            isIncrAtYP = self._determineYP
            if isIncrAtYP and isIncrAtYM:
                # L(-1, -1) <= L(1, -1) and L(-1, 1) <= L(1, 1)
                return LFTTwo.MODE_MP_PM
            elif isIncrAtYP:  # and not isIncrAtYM
                # L(1, -1) <= L(-1, -1) and L(-1, 1) <= L(1, 1)
                return LFTTwo.MODE_MP_MM
            elif isIncrAtYM:  # and not isIncrAtYP
                # L(-1, -1) <= L(1, -1) and L(1, 1) <= L(-1, 1)
                return LFTTwo.MODE_PP_PM
            else:  # not isIncrAtYM and not isIncrAtYP
                # L(1, -1) <= L(-1, -1) and L(1, 1) <= L(-1, 1)
                return LFTTwo.MODE_PP_MM

    def normalize(self):
        [a, b, c, d, e, f, g, h] = self._matrix
        ab = math.gcd(a, b)
        cd = math.gcd(c, d)
        abcd = math.gcd(ab, cd)
        ef = math.gcd(e, f)
        gh = math.gcd(g, h)
        efgh = math.gcd(ef, gh)
        abcdefgh = math.gcd(abcd, efgh)
        gcd = max(1, abcdefgh)
        self._matrix[:] = [
                a // gcd, b // gcd, c // gcd, d // gcd,
                e // gcd, f // gcd, g // gcd, h // gcd
            ]

    @property
    def is_bounded(self):
        [_a, b, _c, d, _e, f, _g, h] = self._matrix
        # D(x, y) = b (xy) + d x + f y + h
        # must not be 0 for x, y in [-1, 1]. By mean value theorem
        # D(+-1, +-1) must all have the same sign and != 0
        # so we test that for
        # D(1, +-1) = (f + b) y + (h + d)
        bounded_at_xp1 = LFTOne.is_plusminus_same_sign(f + b, h + d)
        # D(-1, +-1) = (f - b) y + (h - d)
        bounded_at_xm1 = LFTOne.is_plusminus_same_sign(f - b, h - d)
        # D(+-1, 1) = (d + b) x + (h + f)
        bounded_at_yp1 = LFTOne.is_plusminus_same_sign(d + b, h + f)
        is_bounded = bounded_at_xp1 and bounded_at_xm1 and bounded_at_yp1
        if is_bounded:
            # D(+-1, -1) = (d - b) x + (h - f)
            assert LFTOne.is_plusminus_same_sign(d - b, h - f)
        return is_bounded

    @property
    def bounds(self):
        [a, b, c, d, e, f, g, h] = self._matrix
        at_xm1ym1 = fractions.Fraction(g - e - c + a, h - f - d + b)
        at_xp1ym1 = fractions.Fraction(g - e + c - a, h - f + d - b)
        at_xm1yp1 = fractions.Fraction(g + e - c - a, h + f - d - b)
        at_xp1yp1 = fractions.Fraction(g + e + c + a, h + f + d + b)
        return at_xm1ym1, at_xp1ym1, at_xm1yp1, at_xp1yp1

    @property
    def next_index_to_pull(self):
        [a, b, c, d, e, f, g, h] = self._matrix
        # assert self.is_contracting
        mode = self.lft_type
        # TODO: duplicated work here, when we also calculate this for lft_type
        small_enough = {
            LFTTwo.MODE_MM_PP: lambda: LFTOne.is_small_enough(
                (c + e) * (h + b) - (g + a) * (d + f), (h + b) ** 2 - (d + f) ** 2),
            LFTTwo.MODE_MP_PP: lambda: LFTOne.is_small_enough(
                (c + a) * (h + f) - (g + e) * (d + b), (h + f) ** 2 - (d + b) ** 2),
            LFTTwo.MODE_PM_PP: lambda: LFTOne.is_small_enough(
                (e + a) * (h + d) - (g + c) * (f + b), (h + d) ** 2 - (f + b) ** 2),
            LFTTwo.MODE_MM_PM: lambda: LFTOne.is_small_enough(
                (c - a) * (h - f) - (g - e) * (d - b), (h - f) ** 2 - (d - b) ** 2),
            LFTTwo.MODE_MP_PM: lambda: LFTOne.is_small_enough(
                (c - e) * (h - b) - (g - a) * (d - f), (h - b) ** 2 - (d - f) ** 2),
            LFTTwo.MODE_PP_PM: lambda: LFTOne.is_small_enough(
                (g + c) * (f + b) - (e + a) * (h + d), (h + d) ** 2 - (f + b) ** 2),
            LFTTwo.MODE_MM_MP: lambda: LFTOne.is_small_enough(
                (e - a) * (h - d) - (g - c) * (f - b), (h - d) ** 2 - (f - b) ** 2),
            LFTTwo.MODE_PM_MP: lambda: LFTOne.is_small_enough(
                (g - a) * (d - f) - (c - e) * (h - b), (h - b) ** 2 - (d - f) ** 2),
            LFTTwo.MODE_PP_MP: lambda: LFTOne.is_small_enough(
                (g + e) * (d + b) - (c + a) * (h + f), (h + f) ** 2 - (d + b) ** 2),
            LFTTwo.MODE_MP_MM: lambda: LFTOne.is_small_enough(
                (g - c) * (f - b) - (e - a) * (h - d), (h - d) ** 2 - (f - b) ** 2),
            LFTTwo.MODE_PM_MM: lambda: LFTOne.is_small_enough(
                (g - e) * (d - b) - (c - a) * (h - f), (h - f) ** 2 - (d - b) ** 2),
            LFTTwo.MODE_PP_MM: lambda: LFTOne.is_small_enough(
                (g + a) * (d + f) - (c + e) * (h + b), (h + b) ** 2 - (d + f) ** 2),
        }[mode]()
        if small_enough:
            return None
        # TODO: find a solid way to determine which stream to pull next, instead of basically by chance
        pre_hash = hash((a, b, c, d, e, f, g, h)) % 2**32
        pre_hash = (pre_hash + 0x479ab41d) + (pre_hash << 8)
        pre_hash = (pre_hash ^ 0x5aedd67d) ^ (pre_hash >> 3)
        pre_hash = (pre_hash + 0x17bea992) + (pre_hash << 7)
        return pre_hash % 2

    def extract(self):
        [a, b, c, d, e, f, g, h] = self._matrix
        # assert self.is_contracting
        mode = self.lft_type
        # take the minimum point TODO: biased against smaller negative digits
        extracted_digit = {
            LFTTwo.MODE_MM_PP: lambda: LFTOne.digit_from_lower_bound(+ a - c - e + g,   b - d - f + h),
            LFTTwo.MODE_MP_PP: lambda: LFTOne.digit_from_lower_bound(- a - c + e + g, - b - d + f + h),
            LFTTwo.MODE_PM_PP: lambda: LFTOne.digit_from_lower_bound(- a + c - e + g, - b + d - f + h),
            LFTTwo.MODE_MM_PM: lambda: LFTOne.digit_from_lower_bound(+ a - c - e + g,   b - d - f + h),
            LFTTwo.MODE_MP_PM: lambda: LFTOne.digit_from_lower_bound(- a - c + e + g, - b - d + f + h),
            LFTTwo.MODE_PP_PM: lambda: LFTOne.digit_from_lower_bound(+ a + c + e + g,   b + d + f + h),
            LFTTwo.MODE_MM_MP: lambda: LFTOne.digit_from_lower_bound(+ a - c - e + g,   b - d - f + h),
            LFTTwo.MODE_PM_MP: lambda: LFTOne.digit_from_lower_bound(- a + c - e + g, - b + d - f + h),
            LFTTwo.MODE_PP_MP: lambda: LFTOne.digit_from_lower_bound(+ a + c + e + g,   b + d + f + h),
            LFTTwo.MODE_MP_MM: lambda: LFTOne.digit_from_lower_bound(- a - c + e + g, - b - d + f + h),
            LFTTwo.MODE_PM_MM: lambda: LFTOne.digit_from_lower_bound(- a + c - e + g, - b + d - f + h),
            LFTTwo.MODE_PP_MM: lambda: LFTOne.digit_from_lower_bound(+ a + c + e + g,   b + d + f + h),
        }[mode]()
        assert -POWER_2 < extracted_digit < POWER_2
        self.invtimesdigit(extracted_digit)
        return extracted_digit

    @property
    def is_contracting(self):
        return self.is_bounded and all(abs(bound) <= 1 for bound in self.bounds)


def zero_stream():
    while True:
        yield 0


def one_stream():
    while True:
        yield POWER_2 - 1


def bbp_formula_base_2_32():
    """this calculates pi - 3 in base 16 via the BBP formula.
    We calculate pi - 3 instead of pi, so that the result is
    in the range of [-1, 1]"""
    SHIFT = 4 * 14
    M = 1 << SHIFT
    MASK = M - 1

    def S(j, n):
        # Left sum
        s = 0
        k = 0
        while k <= n:
            r = 8 * k + j
            s = (s + (pow(16, n - k, r) << SHIFT) // r) & MASK
            k += 1
        # fractional part
        t = 0
        k = -1
        while 1:
            # int(16**(n-k) * M)
            xp = int(16 ** k * M)
            newt = t + xp // (8 * (n - k) + j)
            # Iterate until t no longer changes
            if t == newt:
                break
            else:
                t = newt
            k -= 1
        return s + t
    n = 0
    # not all of the last 6 digits are reliable, but the leading 8 are
    EXT_SHIFT = 4 * 6
    EXT_MASK = (1 << 4 * 8) - 1
    while True:
        x = ((4*S(1, n) - 2*S(4, n) - S(5, n) - S(6, n)) >> EXT_SHIFT) & EXT_MASK
        yield x
        n += 8


def adapted_bpp_arbitrary_base():
    bbp_generator = bbp_formula_base_2_32()

    if EXPONENT_2 % 32 == 0:
        digits_per_power2 = EXPONENT_2 // 32
        while True:
            n = digits_per_power2
            out = 0
            while n:
                out = (out << 32) + next(bbp_generator)
                n -= 1
            yield out
    elif 32 % EXPONENT_2 == 0:
        digits_per_gen = 32 // EXPONENT_2
        MASK = POWER_2 - 1
        n = 0
        digit = 0
        while True:
            if n == 0:
                n = digits_per_gen
                digit = next(bbp_generator)
            n -= 1
            yield (digit >> n * EXPONENT_2) & MASK
    else:
        print("Chosen exponent is not compatible with 32. This is not advised.")
        bits_remaining = 0
        digit = 0
        MASK = POWER_2 - 1
        while True:
            while bits_remaining < EXPONENT_2:
                digit = (digit << 32) + next(bbp_generator)
                bits_remaining += 32
            yield (digit >> (bits_remaining - EXPONENT_2)) & MASK
            digit = digit & ~(MASK << (bits_remaining - EXPONENT_2))
            bits_remaining -= EXPONENT_2


def transform_unary(lft, digitstream):
    assert lft.is_contracting

    def transformed():
        local_lft = lft.clone()
        for digit in digitstream():
            while local_lft.next_index_to_pull is None:
                yield local_lft.extract()
                local_lft.normalize()
            local_lft.timesdigit(digit)
    return transformed


def transform_binary(lft, xstream, ystream):
    assert lft.is_contracting

    def transformed():
        local_lft = lft.clone()
        xgen = xstream()
        ygen = ystream()
        while True:
            next_pull = local_lft.next_index_to_pull
            while next_pull is None:
                yield local_lft.extract()
                local_lft.normalize()
                next_pull = local_lft.next_index_to_pull
            if next_pull == 0:
                pulled_x = next(xgen)
                local_lft.timesDigitX(pulled_x)
            else:
                pulled_y = next(ygen)
                local_lft.timesDigitY(pulled_y)
    return transformed


def from_fraction(frac):
    if not abs(frac) <= 1:
        raise ValueError("fraction must be in the interval [-1, 1]")
    return transform_unary(LFTOne.from_fraction(frac), one_stream)


def dec_from_frac(frac):
    return frac.numerator / decimal.Decimal(frac.denominator)


def format_num(digitstream, sift_zeroes=20, sift_digits=4):
    matrix = LFTOne.identity()
    for digit in digitstream():
        if digit == 0 and sift_zeroes > 0:
            sift_zeroes -= 1
        else:
            sift_digits -= 1
        matrix.times(LFTOne.digit(digit))
        if sift_digits <= 0:
            break
    matrix.normalize()
    lower, upper = matrix.bounds
    return "[{l}, {u}]".format(l=dec_from_frac(lower), u=dec_from_frac(upper))


class RealNumber():
    def __init__(self, generator):
        self._generator = generator

    def __str__(self):
        return format_num(self._generator)


class UnaryOperation():
    def __init__(self, lft):
        self._matrix = lft

    def __call__(self, number):
        new_gen = transform_unary(self._matrix, number._generator)
        return RealNumber(new_gen)


class BinaryOperation():
    def __init__(self, lft):
        self._matrix = lft

    def __call__(self, x, y):
        new_gen = transform_binary(self._matrix, x._generator, y._generator)
        return RealNumber(new_gen)
