import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset']='stix'
plt.rcParams['font.size']=10

cm = 1/2.54

Bubble = np.genfromtxt("RP_R5.000e-06_fa6.366e+05_pa2.533e+04.txt", delimiter=" ")

fig1 = plt.figure(figsize=(17*cm,5*cm))
ax1 = plt.subplot2grid((1,3),(0,0),colspan=1)
ax2 = plt.subplot2grid((1,3),(0,1),colspan=1)
ax3 = plt.subplot2grid((1,3),(0,2),colspan=1)
plt.subplots_adjust(wspace=1.2*cm,hspace=1.2*cm)

ax1.set(xlabel=r'$t\, f$',ylabel=r'$R(t)/R_0$')
ax1.set_xlim(xmin=0,xmax=15)
ax1.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax1.plot(Bubble[:, 1]/1.570796e-6, Bubble[:, 3]/5e-6, linestyle='solid', linewidth=1,color='steelblue')

ax2.set(xlabel=r'$t\, f$',ylabel=r'$\dot{R}(t)$[m/s]')
ax2.set_xlim(xmin=0,xmax=15)
ax2.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax2.plot(Bubble[:, 1]/1.570796e-6, Bubble[:, 4], linestyle='solid', linewidth=1,color='steelblue')

ax3.set_yscale('log')
ax3.set(xlabel=r'$t\, f$',ylabel=r'$p_\mathrm{G}(t)$ [Pa]')
ax3.set_xlim(xmin=0,xmax=15)
ax3.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax3.plot(Bubble[:, 1]/1.570796e-6, Bubble[:, 5]/1e3, linestyle='solid', linewidth=1.0,color='steelblue')

ax1.xaxis.set_label_coords(0.5,-0.24)
ax2.xaxis.set_label_coords(0.5,-0.24)
ax3.xaxis.set_label_coords(0.5,-0.24)

ax1.yaxis.set_label_coords(-0.28, 0.5)
ax2.yaxis.set_label_coords(-0.25, 0.5)
ax3.yaxis.set_label_coords(-0.25, 0.5)

fig1.savefig('ultrasound_powerlaw.pdf', bbox_inches='tight',pad_inches=0.035)
