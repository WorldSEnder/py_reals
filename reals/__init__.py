import fractions
import decimal
import math
from .defs import EXPONENT_2, POWER_2
from .lft_one import LFTOne
from .lft_two import LFTTwo


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


def prim_from_fraction(frac):
    if not abs(frac) <= 1:
        raise ValueError("fraction must be in the interval [-1, 1]")
    return transform_unary(LFTOne.from_fraction(frac), one_stream)


def format_num(digitstream, integer_digits, precision=128):
    def dec_from_frac(frac):
        return frac.numerator / decimal.Decimal(frac.denominator)

    precision //= EXPONENT_2
    digit_gen = digitstream()
    integer_part = 0
    for _ in range(integer_digits):
        integer_part *= POWER_2
        integer_part += next(digit_gen)

    matrix = LFTOne(1, 0, integer_part, 1)
    for _ in range(precision):
        digit = next(digit_gen)
        matrix.timesdigit(digit)
        matrix.normalize()
    lower, upper = matrix.bounds
    return "[{l}, {u}]".format(l=dec_from_frac(lower), u=dec_from_frac(upper))


def dec_from_frac(frac):
    return frac.numerator / decimal.Decimal(frac.denominator)


def from_matrix_prod(lft_start, matrix_gen):
    def generator():
        lft = lft_start.clone()
        matrices = matrix_gen()
        while True:
            while not lft.is_contracting or lft.next_index_to_pull is not None:
                lft.times(next(matrices))
            while lft.is_contracting and lft.next_index_to_pull is None:
                yield lft.extract()
            lft.normalize()
    return generator


def log2_matrix_gen():
    n = 1
    while True:
        yield LFTOne(- n, 2 * n + 1, -4 * n, 7 * n + 3)
        n += 1

log2_gen = from_matrix_prod(LFTOne(1, 2, 4, 6), log2_matrix_gen)


class PrimRealNumber():
    def __init__(self, generator):
        self._generator = generator

    def __str__(self):
        return format_num(self._generator, 0)


class PrimUnaryOperation():
    def __init__(self, lft):
        self._matrix = lft

    def __call__(self, number):
        new_gen = transform_unary(self._matrix, number._generator)
        return PrimRealNumber(new_gen)


class PrimBinaryOperation():
    def __init__(self, lft):
        self._matrix = lft

    def __call__(self, x, y):
        new_gen = transform_binary(self._matrix, x._generator, y._generator)
        return PrimRealNumber(new_gen)
