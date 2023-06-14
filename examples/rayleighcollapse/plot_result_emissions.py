import sys
import numpy as np
import matplotlib.pyplot as plt

# Reading an optional argument is required to run this script as a part of a GitHub action
if len(sys.argv) > 1:
    path = str(sys.argv[1])
else:
    path = "./"

plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset']='stix'
plt.rcParams['font.size']=10

cm = 1/2.54
tc = 0.09136264691 # Rayleigh collapse time

# Load the computed results.
Bubble = np.genfromtxt(path + "Gilmore_R1.000e+00.txt", delimiter=" ")
Ep2 = np.genfromtxt(path + "EmissionsSpace_2.000e-01.txt", delimiter=" ")
Ep5 = np.genfromtxt(path + "EmissionsSpace_5.000e-01.txt", delimiter=" ")
EpX = np.genfromtxt(path + "EmissionsSpace_1.000e+00.txt", delimiter=" ")

# Check if the computed results match the reference results of APECSS.
rtol = 1.0e-6
atol = 3.0e-16 # Based (approximately) on the machine precision of double precision.
res = 0
res += np.allclose(np.genfromtxt(path + "reference-data_emissions/APECSS_Gilmore_R1.000e+00.txt", delimiter=" "), Bubble, rtol, atol, False)
res += np.allclose(np.genfromtxt(path + "reference-data_emissions/APECSS_EmissionsSpace_2.000e-01.txt", delimiter=" "), Ep2, rtol, atol, False)
res += np.allclose(np.genfromtxt(path + "reference-data_emissions/APECSS_EmissionsSpace_5.000e-01.txt", delimiter=" "), Ep2, rtol, atol, False)
res += np.allclose(np.genfromtxt(path + "reference-data_emissions/APECSS_EmissionsSpace_1.000e+00.txt", delimiter=" "), EpX, rtol, atol, False)
if res == 4:
    print('All results of the Rayleigh collapse with emissions match the corresponding reference!')
else:
    print('-!!!- WARNING -!!!-', 4-res, 'result(s) of the Rayleigh collapse with emissions DO(ES) NOT match the corresponding reference!')
    sys.exit(1)

# Load the Navier-Stokes reference results
Mp = np.genfromtxt(path + "reference-data_emissions/NavierStokes_pressure_0p05-0p1-0p2-0p5-1-2", delimiter=" ")
Mu = np.genfromtxt(path + "reference-data_emissions/NavierStokes_velocity_0p05-0p1-0p2-0p5-1-2", delimiter=" ")
MbR = np.genfromtxt(path + "reference-data_emissions/NavierStokes_bubbleVolumeRadius", delimiter=" ")
Mbpg = np.genfromtxt(path + "reference-data_emissions/NavierStokes_bubbleAvgPres", delimiter=" ")

fig1 = plt.figure(figsize=(17*cm,15*cm))
ax1 = plt.subplot2grid((3,3),(0,0),colspan=1)
ax3 = plt.subplot2grid((3,3),(0,1),colspan=1)
ax4 = plt.subplot2grid((3,3),(1,0),colspan=1)
ax5 = plt.subplot2grid((3,3),(1,1),colspan=1)
ax6 = plt.subplot2grid((3,3),(1,2),colspan=1)
ax7 = plt.subplot2grid((3,3),(2,0),colspan=1)
ax8 = plt.subplot2grid((3,3),(2,1),colspan=1)
ax9 = plt.subplot2grid((3,3),(2,2),colspan=1)
plt.subplots_adjust(wspace=1.2*cm,hspace=1.2*cm)

ax1.set(xlabel=r'$t/t_\mathrm{c}$',ylabel=r'$R(t)/R_0$')
ax1.set_xlim(xmin=0,xmax=2)
ax1.set_xticks([0,1,2])
ax1.set_ylim(ymin=0,ymax=1)
ax1.set_yticks([0,0.5,1])
ax1.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax1.plot(MbR[:, 0]/tc, MbR[:, 2], linestyle='solid',linewidth=2,  color='mediumseagreen')
ax1.plot(Bubble[:, 1]/tc, Bubble[:, 3], linestyle='solid',linewidth=0.75,  color='navy')

ax3.set_yscale('log')
ax3.set(xlabel=r'$t/t_\mathrm{c}$',ylabel=r'$p_\mathrm{G}(t)/p_\mathrm{G,0}$')
ax3.set_xlim(xmin=0,xmax=2)
ax3.set_xticks([0,1,2])
ax3.set_ylim(ymin=1,ymax=1e6)
ax3.set_yticks([1,1e2,1e4,1e6])
ax3.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax3.plot(Mbpg[:, 0]/tc, Mbpg[:, 1]*1e-3, linestyle='solid',linewidth=2,  color='mediumseagreen', label=r'Navier-Stokes')
ax3.plot(Bubble[:, 1]/tc, Bubble[:, 5]*1e-3, linestyle='solid',linewidth=0.75,  color='navy', label=r'APECSS')

ax3.legend(ncol=1,labelspacing=0.2,markerfirst=True,loc=(1.1,0.7),fontsize='small',facecolor='None',edgecolor='None',framealpha=1,frameon=True)

ax4.set_yscale('log')
ax4.set(xlabel=r'$t/t_\mathrm{c}$',ylabel=r'$\Delta p(t)/p_\infty$')
ax4.set_xlim(xmin=1.004,xmax=1.0325)
ax4.set_ylim(ymin=5e-1,ymax=2e3)
ax4.set_yticks([1,1e1,1e2,1e3])
ax4.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax4.plot(Mp[:, 0]/tc, (Mp[:, 3]-1e5)*1e-5, linestyle='solid',linewidth=3,  color='mediumseagreen')
ax4.plot(Ep2[:, 0]/tc, (Ep2[:, 1]-1e5)*1e-5, linestyle='solid',linewidth=0.75,  color='navy')

ax5.set_yscale('log')
ax5.set(xlabel=r'$t/t_\mathrm{c}$',ylabel=r'$\Delta p(t)/p_\infty$')
ax5.set_xlim(xmin=0.92,xmax=1.12)
ax5.set_xticks([0.92,1.02,1.12])
ax5.set_ylim(ymin=5e-1,ymax=2e3)
ax5.set_yticks([1,1e1,1e2,1e3])
ax5.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax5.plot(Mp[:, 0]/tc, (Mp[:, 4]-1e5)*1e-5, linestyle='solid',linewidth=3,  color='mediumseagreen')
ax5.plot(Ep5[:, 0]/tc, (Ep5[:, 1]-1e5)*1e-5, linestyle='solid',linewidth=0.75,  color='navy')

ax6.set_yscale('log')
ax6.set(xlabel=r'$t/t_\mathrm{c}$',ylabel=r'$\Delta p(t)/p_\infty$')
ax6.set_xlim(xmin=0.92,xmax=1.12)
ax6.set_xticks([0.92,1.02,1.12])
ax6.set_ylim(ymin=5e-1,ymax=2e3)
ax6.set_yticks([1,1e1,1e2,1e3])
ax6.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax6.plot(Mp[:, 0]/tc, (Mp[:, 5]-1e5)*1e-5, linestyle='solid',linewidth=3,  color='mediumseagreen', label=r'Navier-Stokes')
ax6.plot(EpX[:, 0]/tc, (EpX[:, 1]-1e5)*1e-5,linestyle='solid',linewidth=0.75,  color='navy', label=r'APECSS')

ax7.set(xlabel=r'$t/t_\mathrm{c}$',ylabel=r'$u(t)/u_\mathrm{c}$')
ax7.set_xlim(xmin=1.004,xmax=1.0325)
ax7.set_ylim(ymin=-8,ymax=8)
ax7.set_yticks([-8,0,8])
ax7.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax7.plot(Mu[:, 0]/tc, (Mu[:, 3]*tc), linestyle='solid',linewidth=3,  color='mediumseagreen')
ax7.plot(Ep2[:, 0]/tc, (Ep2[:, 2]*tc), linestyle='solid',linewidth=0.75,  color='navy')

ax8.set(xlabel=r'$t/t_\mathrm{c}$',ylabel=r'$u(t)/u_\mathrm{c}$')
ax8.set_xlim(xmin=0.92,xmax=1.12)
ax8.set_xticks([0.92,1.02,1.12])
ax8.set_ylim(ymin=-2,ymax=2)
ax8.set_yticks([-2,0,2])
ax8.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax8.plot(Mu[:, 0]/tc, (Mu[:, 4]*tc), linestyle='solid',linewidth=3,  color='mediumseagreen')
ax8.plot(Ep5[:, 0]/tc, (Ep5[:, 2]*tc), linestyle='solid',linewidth=0.75,  color='navy')

ax9.set(xlabel=r'$t/t_\mathrm{c}$',ylabel=r'$u(t)/u_\mathrm{c}$')
ax9.set_xlim(xmin=0.92,xmax=1.12)
ax9.set_xticks([0.92,1.02,1.12])
ax9.set_ylim(ymin=-0.6,ymax=0.6)
ax9.set_yticks([-0.6,0,0.6])
ax9.grid(color='gainsboro', linestyle='-', linewidth=0.5)
ax9.plot(Mu[:, 0]/tc, (Mu[:, 5]*tc), linestyle='solid',linewidth=3,  color='mediumseagreen', label=r'Navier-Stokes')
ax9.plot(EpX[:, 0]/tc, (EpX[:, 2]*tc), linestyle='solid',linewidth=0.75,  color='navy')

ax1.xaxis.set_label_coords(0.5,-0.24)
ax3.xaxis.set_label_coords(0.5,-0.24)
ax4.xaxis.set_label_coords(0.5,-0.24)
ax5.xaxis.set_label_coords(0.5,-0.24)
ax6.xaxis.set_label_coords(0.5,-0.24)
ax7.xaxis.set_label_coords(0.5,-0.24)
ax8.xaxis.set_label_coords(0.5,-0.24)
ax9.xaxis.set_label_coords(0.5,-0.24)

ax1.yaxis.set_label_coords(-0.25, 0.5)
ax3.yaxis.set_label_coords(-0.25, 0.5)
ax4.yaxis.set_label_coords(-0.25, 0.5)
ax5.yaxis.set_label_coords(-0.25, 0.5)
ax6.yaxis.set_label_coords(-0.25, 0.5)
ax7.yaxis.set_label_coords(-0.25, 0.5)
ax8.yaxis.set_label_coords(-0.25, 0.5)
ax9.yaxis.set_label_coords(-0.25, 0.5)

fig1.savefig('rayleighcollapse_emissions.pdf', bbox_inches='tight',pad_inches=0.035)
sys.exit(0)
