import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl

plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset']='stix'
plt.rcParams['font.size']=10

cm = 1/2.54

Bubble = np.genfromtxt("Gilmore_R1.330e-06.txt", delimiter=" ")
Resm27 = np.genfromtxt("EmissionsTime_6.4953804088e-06.txt", delimiter=" ")
Resm07 = np.genfromtxt("EmissionsTime_6.4955804088e-06.txt", delimiter=" ")
Res000 = np.genfromtxt("EmissionsTime_6.4956504088e-06.txt", delimiter=" ")
Res009 = np.genfromtxt("EmissionsTime_6.4957404088e-06.txt", delimiter=" ")
Res030 = np.genfromtxt("EmissionsTime_6.4959504088e-06.txt", delimiter=" ")
Res076 = np.genfromtxt("EmissionsTime_6.4964104088e-06.txt", delimiter=" ")
Res320 = np.genfromtxt("EmissionsTime_6.4988504088e-06.txt", delimiter=" ")

Liang_R = np.genfromtxt("reference-data/Liang_R.csv", delimiter=",")
Liang_dotR = np.genfromtxt("reference-data/Liang_dotR.csv", delimiter=",")
Liang_pG = np.genfromtxt("reference-data/Liang_pG.csv", delimiter=",")
Lp_m27 = np.genfromtxt("reference-data/LiangP_m27.csv", delimiter=",")
Lp_m07 = np.genfromtxt("reference-data/LiangP_m07.csv", delimiter=",")
Lp_000 = np.genfromtxt("reference-data/LiangP_000.csv", delimiter=",")
Lp_009 = np.genfromtxt("reference-data/LiangP_009.csv", delimiter=",")
Lp_030 = np.genfromtxt("reference-data/LiangP_030.csv", delimiter=",")
Lp_076 = np.genfromtxt("reference-data/LiangP_076.csv", delimiter=",")
Lp_320 = np.genfromtxt("reference-data/LiangP_320.csv", delimiter=",")
Lu_m27 = np.genfromtxt("reference-data/LiangU_m27.csv", delimiter=",")
Lu_m07 = np.genfromtxt("reference-data/LiangU_m07.csv", delimiter=",")
Lu_000 = np.genfromtxt("reference-data/LiangU_000.csv", delimiter=",")
Lu_009 = np.genfromtxt("reference-data/LiangU_009.csv", delimiter=",")
Lu_030 = np.genfromtxt("reference-data/LiangU_030.csv", delimiter=",")
Lu_076 = np.genfromtxt("reference-data/LiangU_076.csv", delimiter=",")
Lu_320 = np.genfromtxt("reference-data/LiangU_320.csv", delimiter=",")

fig1 = plt.figure(figsize=(17*cm,10*cm))
ax1 = plt.subplot2grid((2,9),(0,0),colspan=3)
ax2 = plt.subplot2grid((2,9),(0,3),colspan=3)
ax3 = plt.subplot2grid((2,9),(0,6),colspan=3)
ax4 = plt.subplot2grid((2,21),(1,0),colspan=8)
ax5 = plt.subplot2grid((2,21),(1,10),colspan=8)
plt.subplots_adjust(wspace=30*cm,hspace=1.2*cm)


ax1.set(xlabel=r'$t$ [$\mu$s]',ylabel=r'$R(t)$ [$\mu$m]')
ax1.set_xlim(xmin=0,xmax=9)
ax1.set_xticks([0,3,6,9])
ax1.set_ylim(ymin=-1,ymax=40)
ax1.set_yticks([0,10,20,30,40])
ax1.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax1.plot(Liang_R[:, 0], (Liang_R[:, 1]), linestyle='solid', linewidth=2.0,color='mediumseagreen')
ax1.plot(Bubble[:, 1]*1e6, Bubble[:, 3]*1e6, linestyle='solid', linewidth=0.75,color='navy')

ax2.set(xlabel=r'$t$ [$\mu$s]',ylabel=r'$\dot{R}(t)$ [m/s]')
ax2.set_xlim(xmin=0,xmax=9)
ax2.set_xticks([0,3,6,9])
ax2.set_ylim(ymin=-2000,ymax=1000)
ax2.set_yticks([-2000,0,1000])
ax2.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax2.plot(Liang_dotR[:, 0], (Liang_dotR[:, 1]), linestyle='solid', linewidth=2.0,color='mediumseagreen', label=r'Liang et al. (2022)')
ax2.plot(Bubble[:, 1]*1e6, Bubble[:, 4], linestyle='solid', linewidth=0.75,color='navy', label=r'APECSS')

ax2.legend(ncol=2,labelspacing=0.2,markerfirst=True,loc=(0.5,1),fontsize='x-small',facecolor='None',edgecolor='None',framealpha=1,frameon=True)

ax3.set_yscale('log')
ax3.set(xlabel=r'$t$ [$\mu$s]',ylabel=r'$p_\mathrm{G}(t)$ [Pa]')
ax3.set_xlim(xmin=0,xmax=9)
ax3.set_xticks([0,3,6,9])
ax3.set_ylim(ymin=1,ymax=1e11)
ax3.set_yticks([1e0,1e2,1e4,1e6,1e8,1e10])
ax3.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax3.plot(Liang_pG[:, 0], (Liang_pG[:, 1]*1e6), linestyle='solid', linewidth=2.0,color='mediumseagreen')
ax3.plot(Bubble[:, 1]*1e6, (Bubble[:, 5]), linestyle='solid', linewidth=0.75,color='navy')

n = 7
colors = pl.cm.cividis(np.linspace(0,1,n))

ax4.set_xscale('log')
ax4.set_yscale('log')
ax4.set(xlabel=r'$r$ [m]',ylabel=r'$\Delta p (r)$ [MPa]')
ax4.set_xlim(xmin=4e-7,xmax=1e-5)
ax4.set_ylim(ymin=10,ymax=2e4)
ax4.set_yticks([1e1,1e2,1e3,1e4])
ax4.grid(color='gainsboro', which='major', linestyle='-', linewidth=0.5,alpha=1)
ax4.grid(color='gainsboro', which='minor', linestyle='-', linewidth=0.3,alpha=1)
ax4.plot(Resm27[:, 1], (Resm27[:, 2]-Resm27[:, 5])*1e-6, linestyle='solid', linewidth=0.75, color=colors[0])
ax4.plot(Lp_m27[:, 0], Lp_m27[:, 1],  linestyle='None',markeredgewidth=0.5,marker='o',markeredgecolor=colors[0],markersize=1.75, color=colors[0])
ax4.plot(Resm07[:, 1], (Resm07[:, 2]-Resm07[:, 5])*1e-6, linestyle='solid', linewidth=0.75, color=colors[1])
ax4.plot(Lp_m07[:, 0], Lp_m07[:, 1],  linestyle='None',markeredgewidth=0.5,marker='o',markeredgecolor=colors[1],markersize=1.75, color=colors[1])
ax4.plot(Res000[:, 1], (Res000[:, 2]-Res000[:, 5])*1e-6, linestyle='solid', linewidth=0.75, color=colors[2])
ax4.plot(Lp_000[:, 0], Lp_000[:, 1],  linestyle='None',markeredgewidth=0.5,marker='o',markeredgecolor=colors[2],markersize=1.75, color=colors[2])
ax4.plot(Res009[:, 1], (Res009[:, 2]-Res009[:, 5])*1e-6, linestyle='solid', linewidth=0.75, color=colors[3])
ax4.plot(Lp_009[:, 0], Lp_009[:, 1],  linestyle='None',markeredgewidth=0.5,marker='o',markeredgecolor=colors[3],markersize=1.75, color=colors[3])
ax4.plot(Res030[:, 1], (Res030[:, 2]-Res030[:, 5])*1e-6, linestyle='solid', linewidth=0.75, color=colors[4])
ax4.plot(Lp_030[:, 0], Lp_030[:, 1],  linestyle='None',markeredgewidth=0.5,marker='o',markeredgecolor=colors[4],markersize=1.75, color=colors[4])
ax4.plot(Res076[:, 1], (Res076[:, 2]-Res076[:, 5])*1e-6, linestyle='solid', linewidth=0.75, color=colors[5])
ax4.plot(Lp_076[:, 0], Lp_076[:, 1],  linestyle='None',markeredgewidth=0.5,marker='o',markeredgecolor=colors[5],markersize=1.75, color=colors[5])
ax4.plot(Res320[:, 1], (Res320[:, 2]-Res320[:, 5])*1e-6, linestyle='solid', linewidth=0.75, color=colors[6])
ax4.plot(Lp_320[:, 0], Lp_320[:, 1],  linestyle='None',markeredgewidth=0.5,marker='o',markeredgecolor=colors[6],markersize=1.75, color=colors[6])

ax5.set_xscale('log')
ax5.set(xlabel=r'$r$ [m]',ylabel=r'$u (r)$ [m/s]')
ax5.set_xlim(xmin=4e-7,xmax=1e-5)
ax5.set_ylim(ymin=-2000,ymax=500)
ax5.grid(color='gainsboro', which='major', linestyle='-', linewidth=0.5)
ax5.grid(color='gainsboro', which='minor', axis='x',linestyle='-', linewidth=0.3)
ax5.plot(Resm27[:, 1], (Resm27[:, 3]), linestyle='solid', linewidth=0.75, color=colors[0], label=r'$t_0 - 0.27$ ns')
ax5.plot(Lu_m27[:, 0], Lu_m27[:, 1],  linestyle='None',markeredgewidth=0.5,marker='o',markeredgecolor=colors[0],markersize=1.75, color=colors[0])
ax5.plot(Resm07[:, 1], (Resm07[:, 3]), linestyle='solid', linewidth=0.75, color=colors[1], label=r'$t_0 - 0.07$ ns')
ax5.plot(Lu_m07[:, 0], Lu_m07[:, 1],  linestyle='None',markeredgewidth=0.5,marker='o',markeredgecolor=colors[1],markersize=1.75, color=colors[1])
ax5.plot(Res000[:, 1], (Res000[:, 3]), linestyle='solid', linewidth=0.75, color=colors[2], label=r'$t_0$')
ax5.plot(Lu_000[:, 0], Lu_000[:, 1],  linestyle='None',markeredgewidth=0.5,marker='o',markeredgecolor=colors[2],markersize=1.75, color=colors[2])
ax5.plot(Res009[:, 1], (Res009[:, 3]), linestyle='solid', linewidth=0.75, color=colors[3], label=r'$t_0 + 0.09$ ns')
ax5.plot(Lu_009[:, 0], Lu_009[:, 1],  linestyle='None',markeredgewidth=0.5,marker='o',markeredgecolor=colors[3],markersize=1.75, color=colors[3])
ax5.plot(Res030[:, 1], (Res030[:, 3]), linestyle='solid', linewidth=0.75, color=colors[4], label=r'$t_0 + 0.30$ ns')
ax5.plot(Lu_030[:, 0], Lu_030[:, 1],  linestyle='None',markeredgewidth=0.5,marker='o',markeredgecolor=colors[4],markersize=1.75, color=colors[4])
ax5.plot(Res076[:, 1], (Res076[:, 3]), linestyle='solid', linewidth=0.75, color=colors[5], label=r'$t_0 + 0.76$ ns')
ax5.plot(Lu_076[:, 0], Lu_076[:, 1],  linestyle='None',markeredgewidth=0.5,marker='o',markeredgecolor=colors[5],markersize=1.75, color=colors[5])
ax5.plot(Res320[:, 1], (Res320[:, 3]), linestyle='solid', linewidth=0.75, color=colors[6], label=r'$t_0 + 3.20$ ns')
ax5.plot(Lu_320[:, 0], Lu_320[:, 1],  linestyle='None',markeredgewidth=0.5,marker='o',markeredgecolor=colors[6],markersize=1.75, color=colors[6])

ax5.legend(ncol=1,labelspacing=0.2,markerfirst=True,loc='center left',fontsize='x-small',facecolor='None',edgecolor='None',framealpha=1,frameon=True,bbox_to_anchor=(1, 0.475))

ax1.xaxis.set_label_coords(0.5,-0.24)
ax2.xaxis.set_label_coords(0.5,-0.24)
ax3.xaxis.set_label_coords(0.5,-0.24)
ax4.xaxis.set_label_coords(0.5,-0.18)
ax5.xaxis.set_label_coords(0.5,-0.18)

ax1.yaxis.set_label_coords(-0.25, 0.5)
ax2.yaxis.set_label_coords(-0.25, 0.5)
ax3.yaxis.set_label_coords(-0.25, 0.5)
ax4.yaxis.set_label_coords(-0.18, 0.5)
ax5.yaxis.set_label_coords(-0.27, 0.5)

fig1.savefig('laserinducedcavitation.pdf', bbox_inches='tight',pad_inches=0.035)
