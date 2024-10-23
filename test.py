import numpy as np


def cap_termica(temp):
    return 1608 - 0.7481*temp + 3.929e-4* (temp ** 2)

potenza = 38.700
peak_factor = np.array([.572, .737, .868, .958, 1, .983, .912, .802, .658, .498])

pesi = np.sum(peak_factor)

lin_power = pesi*potenza*0.85
print(lin_power)