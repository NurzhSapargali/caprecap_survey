# Root EST data, N ~ 27000 - very bad, long tails?
# Wang, J.-P. Z., & Lindsay, B. G. (2005). A Penalized Nonparametric Maximum Likelihood Approach to Species Richness Estimation. Journal of the American Statistical Association, 100(471), 942–959. http://www.jstor.org/stable/27590625
f = Dict(
    1 => 2187,
    2 => 490,
    3 => 133,
    4 => 121,
    5 => 37,
    6 => 51,
    7 => 22,
    8 => 19,
    9 => 7,
    10 => 8,
    11 => 6,
    12 => 7,
    13 => 6,
    14 => 4,
    15 => 5,
    16 => 5,
    17 => 18,
    18 => 4,
    19 => 2,
    21 => 2,
    23 => 2,
    24 => 1,
    25 => 6
)

# Golf tees data, N = 250 - Meh
# Borchers, D. L., Buckland, S. T., and Zucchini, W. (2002). Estimating Animal Abundance. Closed Populations. London: Springer.
f = Dict(
    1 => 46,
    2 => 28,
    3 => 21,
    4 => 13,
    5 => 23,
    6 => 14,
    7 => 6,
    8 => 11
)

# Meadow voles data, N = ? - Meh
f = Dict(
    1 => 29,
    2 => 15,
    3 => 15,
    4 => 16,
    5 => 27
)

# McKendrick cholera data, N = ? - Same as alternatives, but not so interesting
f = Dict(
    1 => 32,
    2 => 16,
    3 => 6,
    4 => 1
)

# Heroin usage data in Thailand, N = ?, similar to alternatives, interesting!
f = Dict(
    1 => 2176,
    2 => 1600,
    3 => 1278,
    4 => 976,
    5 => 748,
    6 => 570,
    7 => 455,
    8 => 368,
    9 => 281,
    10 => 254,
    11 => 188,
    12 => 138,
    13 => 99,
    14 => 67,
    15 => 44,
    16 => 34,
    17 => 17,
    18 => 3,
    19 => 3,
    20 => 2,
    21 => 1    
)

# Meth usage data in Thailand, N = ? - awful, long tails
f = Dict(
    1 => 3114,
    2 => 163,
    3 => 23,
    4 => 20,
    5 => 9,
    6 => 3,
    7 => 3,
    8 => 3,
    9 => 4,
    10 => 3
)


# Colorectal cancer data (low-fiber), N = 584 - bad, overestimation
f = Dict(
    1 => 145,
    2 => 66,
    3 => 39,
    4 => 17,
    5 => 8,
    6 => 8,
    7 => 7,
    8 => 3,
    9 => 1,
    11 => 3,
    22 => 1,
    28 => 1
)

# Colorectal cancer data (high-fiber), N = 722 - awful, long tails
f = Dict(
    1 => 144,
    2 => 61,
    3 => 55,
    4 => 37,
    5 => 17,
    6 => 5,
    7 => 4,
    8 => 6,
    9 => 5,
    10 => 1,
    11 => 1,
    31 => 1,
    44 => 1,
    57 => 1,
    70 => 1,
    77 => 1
)

# Scrapie data (Great Britain), N = ?, - awful results, non-monotonic decrease?
f = Dict(
    1 => 84,
    2 => 15,
    3 => 7,
    4 => 5,
    5 => 2,
    6 => 1,
    7 => 2,
    8 => 2
)

# Bowel cancer data, N = 122 - bad, too homogeneous
f = Dict(
    1 => 8,
    2 => 12,
    3 => 16,
    4 => 21,
    5 => 12,
    6 => 31
)

# Brain vessel disease data, N = 366 - not bad, but not interesting
f = Dict(
    1 => 4,
    2 => 15,
    3 => 31,
    4 => 39,
    5 => 55,
    6 => 54,
    7 => 49,
    8 => 47,
    9 => 31,
    10 => 16,
    11 => 9,
    12 => 8,
    13 => 4,
    14 => 3
)

# Product purchase data, N = 456 - underestimation, similar to alternatives
f = Dict(
    1 => 54,
    2 => 49,
    3 => 62,
    4 => 44,
    5 => 25,
    6 => 26,
    7 => 15,
    8 => 15,
    9 => 10,
    10 => 10,
    11 => 10,
    12 => 10,
    13 => 3,
    14 => 3,
    15 => 5,
    16 => 5,
    17 => 4,
    18 => 1,
    19 => 2,
    20 => 1,
)

# Dicentric chromosome data (0.405 Gy), N = 520 - underestimation, worse than alternatives 
f = Dict(
    1 => 66,
    2 => 15,
    3 => 1,
    4 => 1
)

# Dicentric chromosome data (0.6 Gy), N = 631 - underestimation, even worse
f = Dict(
    1 => 119,
    2 => 34,
    3 => 3,
    4 => 2
)
# Hepatitis A data, N = 545 - good, similar to altenratives, interesting but homogeneous
f = Dict(
    1 => 187,
    2 => 56,
    3 => 28
)
# Cottontail rabbit data, N = 135 - bad, overestimation
f = Dict(
    1 => 43,
    2 => 16,
    3 => 8,
    4 => 6,
    6 => 2,
    7 => 1
)
# Taxicab data, N = 420 - okay may be covered by CIs but still homogeneous
f = Dict(
    1 => 142,
    2 => 81,
    3 => 49,
    4 => 7,
    5 => 3,
    6 => 1
)

# Homeless population in Utrecht, N = ? - unknown, need to find more studies with this data
f = Dict(
    1 => 36,
    2 => 11,
    3 => 6,
    4 => 11,
    5 => 5,
    6 => 7,
    7 => 6,
    8 => 11,
    9 => 3,
    10 => 8,
    11 => 7,
    12 => 12,
    13 => 22,
    14 => 77
)

# Dystrophin density data, N = ? - two maxima, one (local) is identical to previous studies and the other (global) is higher
f = Dict(
    1 => 122,
    2 => 50,
    3 => 18,
    4 => 4,
    5 => 4
)
# Gotland Deep data N = ? - awful, long tails
f = Dict(
    1 => 48,
    2 => 9,
    3 => 6,
    4 => 2,
    6 => 2,
    8 => 2,
    9 => 1,
    10 => 1,
    12 => 1,
    13 => 1,
    16 => 1,    
    17 => 2,
    18 => 1,
    20 => 1,
    29 => 1,
    42 => 1,
    53 => 1
)
# Shakespeare data, N = ? - absolutely awful, long tails
f = Dict(
    5 => 1043,
    56 => 19,
    35 => 53,
    55 => 31,
    60 => 14,
    30 => 63,
    32 => 47,
    6 => 837,
    67 => 15,
    45 => 37,
    73 => 10,
    64 => 18,
    90 => 8,
    4 => 1463,
    13 => 242,
    54 => 27,
    63 => 21,
    86 => 11,
    91 => 4,
    62 => 19,
    58 => 22,
    52 => 19,
    12 => 259,
    28 => 76,
    75 => 18,
    23 => 99,
    92 => 7,
    41 => 49,
    43 => 30,
    11 => 305,
    36 => 45,
    68 => 14,
    69 => 11,
    98 => 7,
    82 => 12,
    85 => 10,
    39 => 45,
    84 => 8,
    77 => 8,
    7 => 638,
    25 => 93,
    95 => 10,
    71 => 13,
    66 => 10,
    76 => 11,
    34 => 59,
    50 => 19,
    59 => 23,
    93 => 6,
    2 => 4343,
    10 => 364,
    18 => 130,
    26 => 74,
    27 => 83,
    42 => 41,
    87 => 7,
    100 => 5,
    79 => 12,
    16 => 181,
    20 => 128,
    81 => 13,
    19 => 127,
    49 => 28,
    44 => 35,
    9 => 430,
    31 => 73,
    74 => 16,
    61 => 30,
    29 => 72,
    94 => 7,
    46 => 21,
    57 => 19,
    70 => 16,
    21 => 104,
    38 => 49,
    88 => 12,
    78 => 15,
    72 => 12,
    24 => 112,
    8 => 519,
    17 => 179,
    37 => 34,
    1 => 14376,
    53 => 28,
    22 => 105,
    47 => 41,
    83 => 11,
    99 => 7,
    89 => 9,
    14 => 223,
    3 => 2292,
    80 => 7,
    96 => 10,
    51 => 25,
    33 => 56,
    40 => 52,
    48 => 30,
    15 => 187,
    65 => 15,
    97 => 15
)
# Illegal immigrants in the Netherlands, N = ? - really bad, long tails
f = Dict(
    1 => 1645,
    2 => 183,
    3 => 37,
    4 => 13,
    5 => 1,
    6 => 1
)
# Grizzly bear data (1986)
f = Dict(
    1 => 7,
    2 => 5,
    3 => 6,
    4 => 1,
    5 => 1,
    7 => 1,
    8 => 2,
    15 => 1
)
# Grizzly bear data (1987)
f = Dict(
    1 => 7,
    2 => 3,
    3 => 1,
    4 => 1
)
# Grizzly bear data (1988)
f = Dict(
    1 => 7,
    2 => 4,
    3 => 4,
    4 => 1,
    5 => 1
)
# California drug abuse data, N = ? - severe overestimation
f = Dict(
    1 => 11982,
    2 => 3893,
    3 => 1959,
    4 => 1002,
    5 => 575,
    6 => 340,
    7 => 214,
    8 => 90,
    9 => 72,
    10 => 36,
    11 => 21,
    12 => 14,
)
# AT&T switch error data, N = ? - possible overestimation
f = Dict(5 => 1, 2 => 11, 3 => 1, 1 => 30)
# French scrapie data, N = ? - very bad, long tails
f = Dict(1 => 121, 2 => 13, 3 => 5, 4 => 2) 
# Domestic violence data, N = ? - very bad, long tails
f = Dict(1 => 15169, 2 => 1957, 3 => 393, 4 => 99, 5 => 29, 6 => 16)
# Salmonella data, N = ? - possible underestimation, may be covered by CIs
f = Dict(1 => 17, 2 => 9, 3 => 5, 4 => 6, 5 => 5, 6 => 5, 7 => 6)
# Salmonella data (validation), N = 21 - no changes from No
f = Dict(1 => 1, 2 => 3, 3 => 2, 4 => 3, 5 => 3, 6 => 4, 7 => 2)

# Hong Kong bird data, N = ? - very bad, long tails
f = Dict(
    1 => 21,
    2 => 16,
    3 => 13,
    4 => 10,
    5 => 4,
    6 => 13,
    7 => 6,
    8 => 4,
    9 => 11,
    10 => 1,
    11 => 6,
    12 => 5,
    13 => 8,
    14 => 3,
    15 => 4,
    16 => 6,
    17 => 11,
    18 => 15,
    19 => 8,
    20 => 55
)

# Traffic accident data, N = 9461 - very bad, long tails
f = Dict(
    1 => 1317,
    2 => 239,
    3 => 42,
    4 => 14,
    5 => 4,
    6 => 4,
    7 => 1
)

# Tomato flower EST data, N = ? - very bad, long tails
f = Dict(
    1 => 1434,
    2 => 253,
    3 => 71,
    4 => 33,
    5 => 11,
    6 => 6,
    7 => 2,
    8 => 3,
    9 => 1,
    10 => 2,
    11 => 2,
    12 => 1,
    13 => 1,
    14 => 1,
    16 => 1,
    23 => 1,
    27 => 1
)

# Spina bifida data, N =  ? - okay, possible overestimation
f = Dict(
    1 => 49 + 247 + 60,
    2 => 4 + 112 + 142,
    3 => 12
)

# Diabetes data, N = ? - okay, possible overestimation
f = Dict(
    1 => 182 + 74 + 709 + 10,
    2 => 8 + 7 + 20 + 12 + 650 + 104,
    3 => 14 + 46 + 18 + 157,
    4 => 58
)

# Congenital anomaly data, N = ? - okay, possible overestimation
f = Dict(
    1 => 83 + 97 + 4 + 37 + 27,
    2 => 30 + 3 + 2 + 34 + 37 + 1 + 36 + 22 + 4 + 19,
    3 => 23 + 3 + 5 + 5 + 1 + 18 + 25 + 1,
    4 => 3 + 8 + 2 + 5,
    5 => 2
)

# Dice snakes
f = Dict(
    1 => 59,
    2 => 8,
    3 => 1,
    4 => 1,
    5 => 1
)

# Forgotten books data
f = Dict(
    1 => 42,
    2 => 8,
    3 => 5,
    4 => 5,
    5 => 2,
    6 => 1,
    7 => 6
)

# Forced labour data
f = Dict(
    1 => 4069,
    2 => 1186,
    3 => 167,
    4 => 46,
    5 => 10,
    6 => 7,
    7 => 3,
    8 => 1,
    10 => 1,
    11 => 1
)