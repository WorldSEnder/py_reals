from . import *

zero = PrimRealNumber(zero_stream)
three_forth = PrimRealNumber(prim_from_fraction(fractions.Fraction(3, 4)))
one = PrimRealNumber(one_stream)
xplus3Over4 = PrimUnaryOperation(LFTOne(1, 0, 3, 4))
one_over_xplus2 = PrimUnaryOperation(LFTOne(0, 1, 1, 2))
piMinusThree = PrimRealNumber(adapted_bpp_arbitrary_base)
piForth = xplus3Over4(piMinusThree)
third_of = PrimUnaryOperation(LFTOne(1, 0, 0, 3))
times = PrimBinaryOperation(LFTTwo(1, 0, 3, 0, 3, 0, 0, 10))
timesMinus = PrimBinaryOperation(LFTTwo(-1, 0, -3, 0, -3, 0, 0, 10))

midpoint = PrimBinaryOperation(LFTTwo(0, 0, 1, 0, 1, 0, 0, 2))
mult = PrimBinaryOperation(LFTTwo(1, 0, 0, 0, 0, 0, 0, 1))

seven_eightth = midpoint(three_forth, one)

log2 = PrimRealNumber(log2_gen)

piSquareMinus9Tenth = times(piMinusThree, piMinusThree)
print(piSquareMinus9Tenth)
piSquareMinus9Tenth.stream_to_stdout()
# print(one)
# print(third_of(one))
# print(one_over_xplus2(third_of(one)))
# print(piMinusThree)
# print(piForth)
# print(seven_eightth)
# print(mult(piMinusThree, zero))
# print(timesMinus(piMinusThree, piMinusThree))
# print(log2)
