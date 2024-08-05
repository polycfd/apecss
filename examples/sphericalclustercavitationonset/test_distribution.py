import random
import numpy as np
import matplotlib.pyplot as plt
from math import pi, sqrt, cos, log, exp

# Small file to test if the distribution used for the polydispersed cluster follows a log-normal distribution

def normal_distribution (mu, sigma) :
    u1 = random.random()
    while u1 == 0.0 :
        u1 = random.random()
    u2 = random.random()

    mag = sigma * sqrt(-2 * log(u1))
    z1 = mag * cos(2 * pi * u2) + mu
    return z1

mu = 0.0
sigma = 0.7
nBubbles = 2500
r_ref = 10.0e-6

distrib = []
for i in range(nBubbles) :
    r = r_ref * exp(normal_distribution(mu, sigma))
    while r > 20 * r_ref and r < 0.5 * r_ref :
        r = r_ref * exp(normal_distribution(mu, sigma))
    distrib.append(r)
distrib = np.array(distrib) / r_ref

fig, ax = plt.subplots(1, 1, figsize=((15, 12.5)))
ax.set_title("Test for the proposed radii distribution")

ax.set_xlabel(r"$R/R_{ref}$")
ax.grid()

count, bins, ignored = ax.hist(distrib, 100, density=True, align='mid')
x = np.linspace(min(bins), max(bins), 10000)

pdf = (np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2)) / (x * sigma * np.sqrt(2 * np.pi)))

ax.plot(x, pdf, color="red", linewidth=2)

plt.show()