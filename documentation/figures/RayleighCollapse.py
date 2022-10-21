import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset']='stix'
plt.rcParams['font.size']=10

Bubble = np.genfromtxt("examples/rayleighcollapse/RP_R1.000e+00.txt", delimiter=" ")

cm = 1/2.54
fig1, (ax2,ax3,ax1) = plt.subplots(nrows=1,ncols=3,figsize=(19*cm,4*cm))

ax1.set_yscale('log')
ax1.set(xlabel=r'$t$ [s]',ylabel=r'$p_\mathrm{G}(t)$ [MPa]')
ax1.set_xlim(xmin=0,xmax=1)
ax1.set_xticks([0,0.2,0.4,0.6,0.8,1])
ax1.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax1.plot(Bubble[:, 1], (Bubble[:,5])*1e-6, linestyle='solid', linewidth=1.0,color='steelblue')

ax2.set_yscale('log')
ax2.set(xlabel=r'$t$ [s]',ylabel=r'$R(t)$ [m]')
ax2.set_xlim(xmin=0,xmax=1)
ax2.set_xticks([0,0.2,0.4,0.6,0.8,1])
ax2.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax2.plot(Bubble[:, 1], Bubble[:, 3], linestyle='solid', linewidth=1,color='steelblue')

ax3.set(xlabel=r'$t$ [s]',ylabel=r'$\dot{R}(t)$ [m/s]')
ax3.set_xlim(xmin=0,xmax=1)
ax3.set_xticks([0,0.2,0.4,0.6,0.8,1])
ax3.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax3.plot(Bubble[:, 1], Bubble[:, 4], linestyle='solid', linewidth=1,color='steelblue')

plt.subplots_adjust(wspace=1.3*cm)
ax2.minorticks_off()
ax1.xaxis.set_label_coords(0.5,-0.27)
ax2.xaxis.set_label_coords(0.5,-0.27)
ax3.xaxis.set_label_coords(0.5,-0.27)
ax1.yaxis.set_label_coords(-0.32, 0.5)
ax2.yaxis.set_label_coords(-0.32, 0.5)
ax3.yaxis.set_label_coords(-0.3, 0.5)
fig1.savefig('documentation/figures/RayleighCollapse.pdf', bbox_inches='tight',pad_inches=0.005)