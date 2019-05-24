EXPONENT_2 = 64
POWER_2 = 2 ** EXPONENT_2

if (EXPONENT_2 % 32 != 0) and (32 % EXPONENT_2 != 0):
    print("Chosen exponent is not compatible with 32. This is not advised.",
          f"Some generation algs for irrational numbers work with 32 bits at " +
          "a time and are less efficient with exponent {EXPONENT_2}")
