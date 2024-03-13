import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset']='stix'
plt.rcParams['font.size']=10

cm = 1/2.54
rsh = 997 * (1478**3) / (2*np.pi*4.075*1e3*1e6)  

E1 = np.genfromtxt("EmissionsTime_9.0000000000e-02.txt", delimiter=" ")
E3 = np.genfromtxt("EmissionsTime_2.6000000000e-01.txt", delimiter=" ")
E5 = np.genfromtxt("EmissionsTime_4.3000000000e-01.txt", delimiter=" ")

fig1 = plt.figure(figsize=(17*cm,10*cm))
ax1 = plt.subplot2grid((2,9),(1,0),colspan=3)
ax2 = plt.subplot2grid((2,9),(1,3),colspan=3)
ax3 = plt.subplot2grid((2,9),(1,6),colspan=3)
plt.subplots_adjust(wspace=30*cm,hspace=1.2*cm)

ax1.set(xlabel=r'$(r-R_0)/r_\mathrm{sh}$',ylabel=r'$\Delta p(r,t)$ [MPa]')
ax1.set_xlim(xmin=0.98,xmax=1.02)
ax1.set_ylim(ymin=-1,ymax=1)
ax1.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax1.plot((E1[:, 1]-1)/rsh, (E1[:, 2]-1e5)/1e6, linestyle='solid', linewidth=0.75,color='navy')


ax2.set(xlabel=r'$(r-R_0)/r_\mathrm{sh}$',ylabel=r'$\Delta p(r,t)$ [MPa]')
ax2.set_xlim(xmin=2.98,xmax=3.02)
ax2.set_ylim(ymin=-1,ymax=1)
ax2.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax2.plot((E3[:, 1]-1)/rsh, (rsh/(rsh + (E3[:, 1]-1))*np.pi), linestyle='solid', linewidth=2,color='goldenrod')
ax2.plot((E3[:, 1]-1)/rsh, (-rsh/(rsh + (E3[:, 1]-1))*np.pi), linestyle='solid', linewidth=2,color='goldenrod')
ax2.plot((E3[:, 1]-1)/rsh, (E3[:, 2]-1e5)/1e6, linestyle='solid', linewidth=0.75,color='navy')

ax3.set(xlabel=r'$(r-R_0)/r_\mathrm{sh}$',ylabel=r'$\Delta p(r,t)$ [MPa]')
ax3.set_xlim(xmin=4.98,xmax=5.02)
ax3.set_ylim(ymin=-1,ymax=1)
ax3.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax3.plot((E5[:, 1]-1)/rsh, (rsh/(rsh + (E5[:, 1]-1))*np.pi), linestyle='solid', linewidth=2,color='goldenrod', label='Fay solution')
ax3.plot((E5[:, 1]-1)/rsh, (-rsh/(rsh + (E5[:, 1]-1))*np.pi), linestyle='solid', linewidth=2,color='goldenrod')
ax3.plot((E5[:, 1]-1)/rsh, (E5[:, 2]-1e5)/1e6, linestyle='solid', linewidth=0.75,color='navy', label='APECSS')


ax3.legend(ncol=2,labelspacing=0.2,markerfirst=False,loc=(-0.45,1.04),fontsize='small',facecolor='None',edgecolor='None',framealpha=1,frameon=True)

ax1.yaxis.set_label_coords(-0.27, 0.5)
ax2.yaxis.set_label_coords(-0.27, 0.5)
ax3.yaxis.set_label_coords(-0.27, 0.5)

fig1.savefig('planaremitter.pdf', bbox_inches='tight',pad_inches=0.035)
