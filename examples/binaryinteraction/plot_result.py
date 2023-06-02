import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset']='stix'
plt.rcParams['font.size']=10

cm = 1/2.54

Bubble0 = np.genfromtxt("Bubble_0/KellerMiksis_R5.000e-06_fa1.570e+04_pa-1.200e+05.txt", delimiter=" ")
Bubble1 = np.genfromtxt("Bubble_1/KellerMiksis_R1.000e-05_fa1.570e+04_pa-1.200e+05.txt", delimiter=" ")

fig1 = plt.figure(figsize=(17*cm,5*cm))
ax1 = plt.subplot2grid((1,3),(0,0),colspan=1)
ax2 = plt.subplot2grid((1,3),(0,1),colspan=1)
ax3 = plt.subplot2grid((1,3),(0,2),colspan=1)
plt.subplots_adjust(wspace=1.2*cm,hspace=1.2*cm)

ax1.set(xlabel=r'$t$ [$\mu$s]',ylabel=r'$R(t)$ [$\mu$m]')
ax1.set_xlim(xmin=600,xmax=750)
ax1.set_ylim(ymin=0,ymax=60)
ax1.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax1.plot(Bubble0[:, 1]*1e6, Bubble0[:, 3]*1e6, linestyle='solid', linewidth=1,color='steelblue', label=r'$R_0 = 5 \ \mu \mathrm{m}$')
ax1.plot(Bubble1[:, 1]*1e6, Bubble1[:, 3]*1e6, linestyle='solid', linewidth=1,color='goldenrod', label=r'$R_0 = 10 \ \mu \mathrm{m}$')

ax1.legend(ncol=1,labelspacing=0.2,markerfirst=True,loc='upper right',fontsize='x-small',facecolor='None',edgecolor='None',framealpha=1,frameon=True,bbox_to_anchor=(1, 1))

ax2.set(xlabel=r'$t$ [$\mu$s]',ylabel=r'$\dot{R}(t)$[m/s]')
ax2.set_xlim(xmin=600,xmax=750)
ax2.set_ylim(ymin=-400,ymax=300)
ax2.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax2.plot(Bubble0[:, 1]*1e6, Bubble0[:, 4], linestyle='solid', linewidth=1,color='steelblue')
ax2.plot(Bubble1[:, 1]*1e6, Bubble1[:, 4], linestyle='solid', linewidth=1,color='goldenrod')

ax3.set_yscale('log')
ax3.set(xlabel=r'$t$ [$\mu$s]',ylabel=r'$p_\mathrm{G}(t)$ [Pa]')
ax3.set_xlim(xmin=600,xmax=750)
ax3.set_ylim(ymin=1e1,ymax=1e10)
ax3.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax3.plot(Bubble0[:, 1]*1e6, Bubble0[:, 5], linestyle='solid', linewidth=1.0,color='steelblue')
ax3.plot(Bubble1[:, 1]*1e6, Bubble1[:, 5], linestyle='solid', linewidth=1.0,color='goldenrod')

ax1.xaxis.set_label_coords(0.5,-0.24)
ax2.xaxis.set_label_coords(0.5,-0.24)
ax3.xaxis.set_label_coords(0.5,-0.24)

ax1.yaxis.set_label_coords(-0.25, 0.5)
ax2.yaxis.set_label_coords(-0.25, 0.5)
ax3.yaxis.set_label_coords(-0.25, 0.5)

fig1.savefig('binaryinteraction.pdf', bbox_inches='tight',pad_inches=0.035)
