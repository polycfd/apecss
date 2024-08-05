import os
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset']='stix'
plt.rcParams['font.size']=10

cm = 1/2.54

count = 0
for file in os.listdir() :
    if "Bubble_" in file :
        count += 1

Bubbles = []
for i in range(count) :
    Bubble_file = os.listdir(os.path.join(os.getcwd(), "Bubble_{}".format(i)))[0]
    Bubble = np.genfromtxt("Bubble_{}/".format(i) + Bubble_file, delimiter=" ")
    Bubbles.append(Bubble)

Bubbles_A = [[]]
for i in range(count) :
    Bubbles_A.append([])

file_results = open("Ida2009_results.txt", "r")
lines = file_results.readlines()
file_results.close()

for line in lines[3:] :
    data = line.split(" ")
    t = float(data[0])
    Bubbles_A[0].append(t)
    for i in range(count) :
        Bubbles_A[1+i].append(float(data[1 + 2 * count + i]))

fig1 = plt.figure(figsize=(21*cm,5*cm))
ax1 = plt.subplot2grid((1,4),(0,0),colspan=1)
ax2 = plt.subplot2grid((1,4),(0,1),colspan=1)
ax3 = plt.subplot2grid((1,4),(0,2),colspan=1)
ax4 = plt.subplot2grid((1,4),(0,3),colspan=1)
plt.subplots_adjust(wspace=1.2*cm,hspace=1.2*cm)

ax1.set_xlabel(r"t ($\mu$s)")
ax1.set_xlim(xmin=10.0, xmax=60.0)
ax2.set_xlabel(r"t ($\mu$s)")
ax2.set_xlim(xmin=10.0, xmax=60.0)
ax3.set_xlabel(r"t ($\mu$s)")
ax3.set_xlim(xmin=10.0, xmax=60.0)
ax4.set_xlabel(r"t ($\mu$s)")
ax4.set_xlim(xmin=10.0, xmax=60.0)

ax1.set_ylabel(r"R ($\mu$m)")
ax2.set_ylabel(r"U (m.$s^{-1}$)")
ax3.set_ylabel(r"$p_{\infty}$ / $p_{0}$ (-)")
ax3.set_ylim(ymin=-0.25, ymax=0.0)
ax4.set_ylabel(r"A (m.$s^{-2}$)")

ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()

for i in range(count) :
    ax1.plot(Bubbles[i][:, 1]*1e6, Bubbles[i][:, 3]*1e6)
    ax2.plot(Bubbles[i][:, 1]*1e6, Bubbles[i][:, 4])
    ax3.plot(Bubbles[i][:, 1]*1e6, Bubbles[i][:, 7] * 1 / (101.3e03))
    ax4.plot(np.array(Bubbles_A[0])*1e6, np.array(Bubbles_A[1+i]))

plt.show()