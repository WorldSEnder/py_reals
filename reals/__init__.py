import fractions
import decimal
import math
from .defs import EXPONENT_2, POWER_2, PRINT_HEX
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


def convert_base(digitstream, orig_base, target_base):
    def is_power2(n):
        return not (n & (n - 1))

    def exact_log2(n):
        return n.bit_length() - 1

    def is_exactly_convertible(bf, bt):
        while bt > 1:
            if bt % bf != 0:
                return False
            bt //= bf
        return True

    def discrete_log(bf, bt):
        n = 0
        while bt > 1:
            bt //= bf
            n += 1
        return n

    def largest_shared_power(bf, bt):
        while bt > 1:
            while bf % bt == 0:
                bf //= bt
            bf, bt = bt, bf
        return bf

    def get_split_strat():
        # this is complicated by the fact that
        # each digit can either be positive or
        # negative! This leads us to the realization
        # that there is no one unique sequence
        # realizing the split but rather multiple
        if is_power2(target_base):
            mask = target_base - 1
            shift = exact_log2(target_base)

            def split_trgt_base(p, n):
                split = p >> (n * shift)
                rest = p - (split << (n * shift))
                return rest, split
        else:
            def split_trgt_base(p, n):
                base_pow = pow(target_base, n)
                split = p // base_pow
                if split < 0:
                    split += 1
                rest = p - split * base_pow
                return rest, split
        return split_trgt_base

    digit_gen = digitstream()

    if is_exactly_convertible(orig_base, target_base):
        # original base is smaller, but fits exactly
        in_digits_per_out = discrete_log(orig_base, target_base)
        if is_power2(orig_base):
            shift = exact_log2(orig_base)

            def mult_orig_base(p):
                return p << shift
        else:
            def mult_orig_base(p):
                return p * orig_base
        while True:
            out = 0
            for _ in range(in_digits_per_out):
                out = mult_orig_base(out) + next(digit_gen)
            yield out

    elif is_exactly_convertible(target_base, orig_base):
        # target base is smaller, but fits exactly
        split_trgt_base = get_split_strat()
        out_digits_per_in = discrete_log(target_base, orig_base)
        while True:
            digit = next(digit_gen)
            for n in range(out_digits_per_in - 1, -1, -1):
                digit, part = split_trgt_base(digit, n)
                yield part
    shared_power = largest_shared_power(target_base, orig_base)
    if shared_power == 1:
        raise NotImplementedError()
    yield from convert_base(lambda: convert_base(digitstream, orig_base, shared_power), shared_power, target_base)


def adapted_bpp_arbitrary_base():
    yield from convert_base(bbp_formula_base_2_32, 2**32, POWER_2)


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

    base = POWER_2
    precision //= EXPONENT_2
    digit_gen = digitstream()
    integer_part = 0
    for _ in range(integer_digits):
        integer_part *= base
        integer_part += next(digit_gen)

    matrix = LFTOne(1, 0, integer_part, 1)
    for _ in range(precision):
        digit = next(digit_gen)
        matrix.timesdigitbase(digit, base)
        matrix.normalize()
    lower, upper = matrix.bounds
    return "[{l}, {u}]".format(l=dec_from_frac(lower), u=dec_from_frac(upper))


def format_hex(digitstream, integer_digits, precision=512 // 4):
    def to_hex(digit):
        assert 0 <= digit < 16
        return "%x" % digit

    def inv_digit(digit):
        assert 0 <= digit < 16
        return 16 - digit

    digit_gen = convert_base(digitstream, POWER_2, 16)
    zeroes = 0
    digit = next(digit_gen)
    while digit == 0 and precision > 0:
        zeroes += 1
        digit = next(digit_gen)
        precision -= 1
    outstr = ("-" if digit < 0 else " ") + "." + ("0" * zeroes)
    zeroes = 0
    sign = -1 if digit < 0 else 1
    saved = abs(digit)
    while precision > 0:
        zeroes = 0
        precision -= 1
        digit = next(digit_gen)
        while digit == 0 and precision > 0:
            zeroes += 1
            digit = next(digit_gen)
            precision -= 1
        digit *= sign
        if digit < 0:
            outstr += to_hex(saved - 1) + ("f" * zeroes)
            saved = inv_digit(-digit)
        elif digit > 0:
            outstr += to_hex(saved) + ("0" * zeroes)
            saved = digit
    if zeroes > 0:
        rounding = next(digit_gen) * sign
        if rounding < 0:
            outstr += to_hex(saved - 1) + ("f" * zeroes)
        elif rounding >= 0:
            outstr += to_hex(saved) + ("0" * zeroes)
    return outstr


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
        if PRINT_HEX:
            return format_hex(self._generator, 0)
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
