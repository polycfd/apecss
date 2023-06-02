import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset']='stix'
plt.rcParams['font.size']=10

cm = 1/2.54

Bubble = np.genfromtxt("Gilmore_R5.000e-06_fa2.350e+04_pa1.450e+05.txt", delimiter=" ")
pLmax = np.genfromtxt("EmissionsNode_pLmax.txt", delimiter=" ")
Holzfuss_p = np.genfromtxt("reference-data/Holzfuss2010_p.csv", delimiter=",")
Holzfuss_u = np.genfromtxt("reference-data/Holzfuss2010_u.csv", delimiter=",")

fig1 = plt.figure(figsize=(17*cm,10*cm))
ax1 = plt.subplot2grid((2,9),(0,0),colspan=3)
ax2 = plt.subplot2grid((2,9),(0,3),colspan=3)
ax3 = plt.subplot2grid((2,9),(0,6),colspan=3)
ax4 = plt.subplot2grid((2,21),(1,0),colspan=8)
ax5 = plt.subplot2grid((2,21),(1,10),colspan=8)
plt.subplots_adjust(wspace=30*cm,hspace=1.2*cm)

ax1.set(xlabel=r'$t$ [$\mu$s]',ylabel=r'$R(t)$ [$\mu$m]')
ax1.set_xlim(xmin=0,xmax=40)
ax1.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax1.plot(Bubble[:, 1]*1e6, Bubble[:, 3]*1e6, linestyle='solid', linewidth=0.75,color='navy')

ax2.set(xlabel=r'$t$ [$\mu$s]',ylabel=r'$\dot{R}(t)$ [m/s]')
ax2.set_xlim(xmin=0,xmax=40)
ax2.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax2.plot(Bubble[:, 1]*1e6, Bubble[:, 4], linestyle='solid', linewidth=0.75,color='navy')

ax3.set_yscale('log')
ax3.set(xlabel=r'$t$ [$\mu$s]',ylabel=r'$p_\mathrm{G}(t)$ [Pa]')
ax3.set_xlim(xmin=0,xmax=40)
ax3.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax3.plot(Bubble[:, 1]*1e6, (Bubble[:, 5]), linestyle='solid', linewidth=0.75,color='navy')

ax4.set_xscale('log')
ax4.set_yscale('log')
ax4.set(xlabel=r'$r$ [m]',ylabel=r'$\Delta p (r,t)$ [MPa]')
ax4.set_xlim(xmin=0.5e-6,xmax=1e-3)
ax4.set_xticks([1e-6,1e-5,1e-4,1e-3])
ax4.set_ylim(ymin=0.5,ymax=3e4)
ax4.set_yticks([1,1e1,1e2,1e3,1e4])
ax4.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax4.plot([2e-6,5e-4], [250,1], color='gray', linestyle='dashed', lw=0.7)
ax4.annotate('$-1$', fontsize = 8, color='gray', xy=(0.455,0.18), xycoords='axes fraction')
ax4.plot(Holzfuss_p[:, 0], (Holzfuss_p[:, 1]-1)*1e-6, linestyle='solid',linewidth=2.5,  color='mediumseagreen')
ax4.plot(pLmax[:, 1], (pLmax[:, 3]-pLmax[:, 6])*1e-6,  linestyle='solid', linewidth=0.75,color='navy')

ax5.set_xscale('log')
ax5.set_yscale('log')
ax5.set(xlabel=r'$r$ [m]',ylabel=r'$u (r,t)$ [m/s]')
ax5.set_xlim(xmin=0.5e-6,xmax=1e-3)
ax5.set_xticks([1e-6,1e-5,1e-4,1e-3])
ax5.set_ylim(ymin=2e-1,ymax=7e2)
ax5.set_yticks([1e0,1e1,1e2])
ax5.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax5.plot([4e-6,7e-4], [105,0.6], color='gray', linestyle='dashed', lw=0.7)
ax5.annotate('$-1$', fontsize = 8, color='gray', xy=(0.5,0.32), xycoords='axes fraction')
ax5.plot(Holzfuss_u[:, 0], Holzfuss_u[:, 1],  linestyle='solid',linewidth=2.5,  color='mediumseagreen',label=r'Holzfuss (2010)')
ax5.plot(pLmax[:, 1], pLmax[:, 4], linestyle='solid', linewidth=0.75,color='navy',label=r'APECSS')

ax5.legend(ncol=1,labelspacing=0.2,markerfirst=False,loc=(1.1,0.7),fontsize='x-small',facecolor='None',edgecolor='None',framealpha=1,frameon=True)

ax1.xaxis.set_label_coords(0.5,-0.24)
ax2.xaxis.set_label_coords(0.5,-0.24)
ax3.xaxis.set_label_coords(0.5,-0.24)
ax4.xaxis.set_label_coords(0.5,-0.18)
ax5.xaxis.set_label_coords(0.5,-0.18)

ax1.yaxis.set_label_coords(-0.25, 0.5)
ax2.yaxis.set_label_coords(-0.25, 0.5)
ax3.yaxis.set_label_coords(-0.25, 0.5)
ax4.yaxis.set_label_coords(-0.18, 0.5)
ax5.yaxis.set_label_coords(-0.23, 0.5)

fig1.savefig('ultrasound_sonolum.pdf', bbox_inches='tight',pad_inches=0.035)
