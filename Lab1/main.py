from math import sqrt
import matplotlib.pyplot as plt
import random

K = 256
N = 1000
COEF_A = 65643
COEF_B = 65643
MOD = 2147483647
INTERVALS = 30
KOLMOGOROV_QUANTILE = 1.36
PIRSON_QUANTILE = 43.8


def empiric_dist(x, values):
    if x < min(values):
        return 0.0
    if x > max(values):
        return 1.0

    less = 0
    for val in values:
        if val <= x:
            less += 1

    return less / len(values)


def rect_dist(x, a, b):
    if x < a:
        return 0.0
    if x > b:
        return 1.0
    return (x - a) / (b - a)


def kolmogorov_stat(values):
    stat = 0
    x = min(values)
    step = (max(values) - min(values)) / INTERVALS
    for i in range(INTERVALS):
        stat = max(abs(rect_dist(x, 0, 1) - empiric_dist(x, values)), stat)
        x += step

    return stat


def pirson_stat(values):
    stat = 0
    step = (max(values) - min(values)) / INTERVALS

    x = min(values)
    y = x + step
    for i in range(1,INTERVALS):
        observed = 0
        for val in values:
            if x < val <= y:
                observed += 1

        expected = len(values) * (rect_dist(y, 0, 1) - rect_dist(x, 0, 1))
        stat += ((observed - expected)**2) / expected
        x = y
        y += step

    return stat


def drawing():
    plt.figure(figsize=(9, 4))

    plt.subplot(2, 2, 1)
    plt.hist(MCG, bins=10)
    plt.title('MCG')

    plt.subplot(2, 2, 2)
    plt.title('MCG')
    x_list = [0.0]
    x_list.extend(MCG[:-1])
    plt.scatter(x_list, MCG)

    plt.subplot(2, 2, 3)
    plt.hist(MLMG, bins=10)
    plt.title('MLMG')

    plt.subplot(2, 2, 4)
    plt.title('MLMG')
    x_list2 = [0.0]
    x_list2.extend(MLMG[:-1])
    plt.scatter(x_list2, MLMG)

    plt.show()
    # end of drawing


if __name__ == '__main__':
    # Multiplicative Congruential Generator
    a = [COEF_A]
    MCG = [COEF_A / MOD]

    for i in range(1, N):
        a.append(COEF_B * a[i - 1] % MOD)
        MCG.append(a[i] / MOD)

    # MacLaren Marsaglia (MCM and random)
    MLMG = []
    index = 0
    rand_list = [random.random() for i in range(0, N)]
    v_table = [MCG[i] for i in range(0, K)]

    for i in range(0, N - K):
        index = int(K * rand_list[i])
        MLMG.append(v_table[index])
        v_table[index] = MCG[i + K]

    stat_MCG = kolmogorov_stat(MCG)
    stat_MLMG = kolmogorov_stat(MLMG)
    print("Kolmogorov statistics for MCG: " + str(stat_MCG))
    print(sqrt(len(MCG)) * stat_MCG < KOLMOGOROV_QUANTILE)
    print("Kolmogorov statistics for MLMG: " + str(stat_MLMG))
    print(sqrt(len(MLMG)) * stat_MLMG < KOLMOGOROV_QUANTILE)

    stat_MCG = pirson_stat(MCG)
    stat_MLMG = pirson_stat(MLMG)
    print("Pirson statistics for MCG: " + str(stat_MCG))
    print(stat_MCG < PIRSON_QUANTILE)
    print("Pirson statistics for MLMG: " + str(stat_MLMG))
    print(stat_MLMG < PIRSON_QUANTILE)

    drawing()
