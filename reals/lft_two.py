import fractions
import math
from .defs import POWER_2, EXPONENT_2
from .lft_one import LFTOne


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
        self._calculateCharacteristics()

    def clone(self):
        [a, b, c, d, e, f, g, h] = self._matrix
        return LFTTwo(a, b, c, d, e, f, g, h)

    def __str__(self):
        [a, b, c, d, e, f, g, h] = self._matrix
        return "[{a}\t{c}\t| {e}\t{g}\n{b}\t{d}\t| {f}\t{h}]".format(
            a=a, b=b, c=c, d=d, e=e, f=f, g=g, h=h)

    def _calculateCharacteristics(self):
        [a, b, c, d, e, f, g, h] = self._matrix
        self._lft_type = mode = self._calc_lft_type()
        self._interval_length_num, self._interval_length_denom = {
            LFTTwo.MODE_MM_PP: lambda: ((c + e) * (h + b) - (g + a) * (d + f), (h + b) ** 2 - (d + f) ** 2),
            LFTTwo.MODE_MP_PP: lambda: ((c + a) * (h + f) - (g + e) * (d + b), (h + f) ** 2 - (d + b) ** 2),
            LFTTwo.MODE_PM_PP: lambda: ((e + a) * (h + d) - (g + c) * (f + b), (h + d) ** 2 - (f + b) ** 2),
            LFTTwo.MODE_MM_PM: lambda: ((c - a) * (h - f) - (g - e) * (d - b), (h - f) ** 2 - (d - b) ** 2),
            LFTTwo.MODE_MP_PM: lambda: ((c - e) * (h - b) - (g - a) * (d - f), (h - b) ** 2 - (d - f) ** 2),
            LFTTwo.MODE_PP_PM: lambda: ((g + c) * (f + b) - (e + a) * (h + d), (h + d) ** 2 - (f + b) ** 2),
            LFTTwo.MODE_MM_MP: lambda: ((e - a) * (h - d) - (g - c) * (f - b), (h - d) ** 2 - (f - b) ** 2),
            LFTTwo.MODE_PM_MP: lambda: ((g - a) * (d - f) - (c - e) * (h - b), (h - b) ** 2 - (d - f) ** 2),
            LFTTwo.MODE_PP_MP: lambda: ((g + e) * (d + b) - (c + a) * (h + f), (h + f) ** 2 - (d + b) ** 2),
            LFTTwo.MODE_MP_MM: lambda: ((g - c) * (f - b) - (e - a) * (h - d), (h - d) ** 2 - (f - b) ** 2),
            LFTTwo.MODE_PM_MM: lambda: ((g - e) * (d - b) - (c - a) * (h - f), (h - f) ** 2 - (d - b) ** 2),
            LFTTwo.MODE_PP_MM: lambda: ((g + a) * (d + f) - (c + e) * (h + b), (h + b) ** 2 - (d + f) ** 2),
        }[mode]()
        self._lowest_bound_num, self._lowest_bound_denom = {
            LFTTwo.MODE_MM_PP: lambda: (+ a - c - e + g,   b - d - f + h),
            LFTTwo.MODE_MP_PP: lambda: (- a - c + e + g, - b - d + f + h),
            LFTTwo.MODE_PM_PP: lambda: (- a + c - e + g, - b + d - f + h),
            LFTTwo.MODE_MM_PM: lambda: (+ a - c - e + g,   b - d - f + h),
            LFTTwo.MODE_MP_PM: lambda: (- a - c + e + g, - b - d + f + h),
            LFTTwo.MODE_PP_PM: lambda: (+ a + c + e + g,   b + d + f + h),
            LFTTwo.MODE_MM_MP: lambda: (+ a - c - e + g,   b - d - f + h),
            LFTTwo.MODE_PM_MP: lambda: (- a + c - e + g, - b + d - f + h),
            LFTTwo.MODE_PP_MP: lambda: (+ a + c + e + g,   b + d + f + h),
            LFTTwo.MODE_MP_MM: lambda: (- a - c + e + g, - b - d + f + h),
            LFTTwo.MODE_PM_MM: lambda: (- a + c - e + g, - b + d - f + h),
            LFTTwo.MODE_PP_MM: lambda: (+ a + c + e + g,   b + d + f + h),
        }[mode]()

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
        self._calculateCharacteristics()

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
        self._calculateCharacteristics()

    def timesDigitX(self, digit):
        [a, b, c, d, e, f, g, h] = self._matrix
        w, exp = digit, EXPONENT_2
        self._matrix[0] = a
        self._matrix[1] = b
        self._matrix[2] = a * w + (c << exp)
        self._matrix[3] = b * w + (d << exp)
        self._matrix[4] = e
        self._matrix[5] = f
        self._matrix[6] = e * w + (g << exp)
        self._matrix[7] = f * w + (h << exp)
        # TODO: inline
        self._calculateCharacteristics()

    def timesDigitY(self, digit):
        [a, b, c, d, e, f, g, h] = self._matrix
        w, exp = digit, EXPONENT_2
        self._matrix[0] = a
        self._matrix[1] = b
        self._matrix[2] = c
        self._matrix[3] = d
        self._matrix[4] = a * w + (e << exp)
        self._matrix[5] = b * w + (f << exp)
        self._matrix[6] = c * w + (g << exp)
        self._matrix[7] = d * w + (h << exp)
        # TODO: inline
        self._calculateCharacteristics()

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
        self._calculateCharacteristics()

    def invtimesdigit(self, digit):
        [a, b, c, d, e, f, g, h] = self._matrix
        v, exp = -digit, EXPONENT_2
        self._matrix[0] = (a << exp) + v * b
        self._matrix[1] = b
        self._matrix[2] = (c << exp) + v * d
        self._matrix[3] = d
        self._matrix[4] = (e << exp) + v * f
        self._matrix[5] = f
        self._matrix[6] = (g << exp) + v * h
        self._matrix[7] = h
        # TODO: inline
        self._calculateCharacteristics()

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

    def _calc_lft_type(self):
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

    @property
    def lft_type(self):
        return self._lft_type

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
        # TODO: inline
        self._calculateCharacteristics()

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
        # TODO: duplicated work here, when we also calculate this for lft_type
        small_enough = LFTOne.is_small_enough(self._interval_length_num, self._interval_length_denom)
        if small_enough:
            return None
        # TODO: find a solid way to determine which stream to pull next, instead of basically by chance
        pre_hash = hash((a, b, c, d, e, f, g, h)) % 2**32
        return (pre_hash >> 31) % 2

    def extract(self):
        # assert self.is_contracting
        # take the minimum point TODO: biased against smaller negative digits
        extracted_digit = LFTOne.digit_from_lower_bound(self._lowest_bound_num, self._lowest_bound_denom)
        assert -POWER_2 < extracted_digit < POWER_2
        self.invtimesdigit(extracted_digit)
        return extracted_digit

    @property
    def is_contracting(self):
        return self.is_bounded and all(abs(bound) <= 1 for bound in self.bounds)

__all__ = ["LFTTwo"]
