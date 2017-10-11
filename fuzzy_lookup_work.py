import numpy as np
import matplotlib.pyplot as plt
import membership_degree as md

# DATA SERIES
beta = 0.2
gamma = 0.1
n = 10
tau = 30
series_length = 2000
x_init = 0.1
history_length = 4
training_data_size = 600

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
print("STREDY:")
print(sets_centers)
distance_from_center_to_bounds = abs(sets_centers[0] - sets_centers[1]) / 2
# print(distance_from_center_to_bounds)

# SESTAVIT VEKTOR IO dat

io_data = np.zeros((training_data_size, history_length + 1))

for i in range(0, training_data_size):
    io_data[i, :] = x[i + tau:i + history_length + 1 + tau]

# print(x[tau:tau+10])
# print(io_data)
params1 = np.zeros((2,1))
params1[0] = sets_centers[0]
params1[1] = sets_centers[1]
# TODO: zjistit hodnotu MF ke vsem mnozinam
mf = md.get_l_mf_degree(io_data[1, 1], params1)
print(mf)
for i in range(1, fuzzy_sets_number - 1):
    print(sets_centers[i])


params1[0] = sets_centers[-2]
params1[1] = sets_centers[-1]
mf = md.get_r_mf_degree(io_data[1, 1], params1)
# TODO: zjistit mnozinu s nejvetsi MF

# TODO: sestavit pravidlo a spocitat degree...

# TODO: check konzistence rule base a vybrat ty s nejvetsim degree

# TODO: sestaveni fuzzy systemu

# TODO: overeni na zbytku dat

# TODO: vykresleni vysledku