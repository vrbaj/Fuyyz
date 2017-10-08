import numpy as np
import matplotlib.pyplot as plt


# TODO: Pocet fuzzy mnozin
# TODO: Pocet pravidel
# TODO: Umistit stredy etc.
# TODO: Zjistit platnost pravidel
# TODO: Sestavit fuzzy system
# TODO: overeni


# DATA SERIES
beta = 0.2
gamma = 0.1
n = 10
tau = 30
series_length = 5000
x_init = 0.1

x = np.zeros(series_length)
x[0:tau] = x_init

for i in range(tau-1, series_length - 1):
    x[i+1] = x[i] + beta * x[i - tau] / (1 + x[i - tau] ** n) - gamma * x[i]

plt.plot(x)
plt.ylabel('x[k]')
plt.show()

# UMISTENI STREDU, POCET FUZZY MNOZIN
left_bound = min(x)
right_bound = max(x)
fuzzy_sets_number = 5
sets_centers = np.linspace(left_bound, right_bound, fuzzy_sets_number)
print(sets_centers)