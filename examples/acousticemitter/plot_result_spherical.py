import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset']='stix'
plt.rcParams['font.size']=10

cm = 1/2.54
L = 1478/1500
k = 2*np.pi/L
R0 = 1/(2*np.pi)

kbres = np.genfromtxt("EmissionsTime_4.0000000000e-03.txt", delimiter=" ")


fig1 = plt.figure(figsize=(17*cm,10*cm))
ax1 = plt.subplot2grid((2,2),(0,0),colspan=1)
ax2 = plt.subplot2grid((2,2),(0,1),colspan=1)
plt.subplots_adjust(wspace=0.95*cm,hspace=1*cm)


ax1.plot((kbres[:,1] - R0)/L, ((R0/kbres[:,1])*np.cos(12*np.pi-k*(kbres[:,1] - R0)-0.5*np.pi)),linestyle='solid', linewidth=2, color='goldenrod')
ax1.plot((kbres[:,1] - R0)/L, (kbres[:,2]-1e5),linestyle='solid', linewidth=1.0, color='black')
ax1.set(xlabel=r'${(r-R)}/{\lambda_\mathrm{a}}$',ylabel=r'$p_1/\Delta p_\mathrm{a}$')
ax1.set_xlim(xmin=0,xmax=4)
ax1.set_ylim(ymin=-0.6,ymax=0.3)
ax1.set_yticks([-0.6,-0.3,0,0.3])
ax1.grid(color='gainsboro', linestyle='-', linewidth=0.5)

ax2.plot((kbres[:,1] - R0)/L, ((R0*R0/(kbres[:,1]*kbres[:,1]))*(np.cos(12*np.pi-k*(kbres[:,1]-R0)-0.5*np.pi)-np.cos(12*np.pi-k*(kbres[:,1]-R0)-0.5*np.pi-(0.5*np.pi-np.arctan(k*kbres[:,1])))/(np.cos(0.5*np.pi-np.arctan(k*kbres[:,1])))))+((R0/kbres[:,1])*np.cos(12*np.pi-k*(kbres[:,1]-R0)-0.5*np.pi-(0.5*np.pi-np.arctan(k*kbres[:,1])))/(np.cos(0.5*np.pi-np.arctan(k*kbres[:,1])))), linewidth=2, color='goldenrod',label='Analytical')
ax2.plot((kbres[:,1] - R0)/L, (kbres[:,3]*997*1478),linestyle='solid', linewidth=1, color='black',label='Kirkwood-Bethe')
ax2.set(xlabel=r'${(r-R)}/{\lambda_\mathrm{a}}$',ylabel=r'$u \rho c/\Delta p_\mathrm{a}$')
ax2.set_xlim(xmin=0,xmax=4)
ax2.set_ylim(ymin=-0.6,ymax=0.3)
ax2.set_yticks([-0.6,-0.3,0,0.3])
ax2.grid(color='gainsboro', linestyle='-', linewidth=0.5)

ax2.legend(fontsize='small', loc=(0,1.04), ncol=3, frameon=False)

ax1.xaxis.set_label_coords(0.5,-0.19)
ax2.xaxis.set_label_coords(0.5,-0.19)

ax1.yaxis.set_label_coords(-0.23,0.5)
ax2.yaxis.set_label_coords(-0.23,0.5)

fig1.savefig('sphericalemitter.pdf', bbox_inches='tight', pad_inches=0.035)