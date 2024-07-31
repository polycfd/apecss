import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset']='stix'
plt.rcParams['font.size']=25

color_names = list(mcolors.XKCD_COLORS)

cm = 1/2.54

# File designed to recover and copy/paste results to reproduce Ida's cavitation onset test case
# DOI : https://doi.org/10.1063/1.3265547

######### Step 1 : Recovering data ##################################################################################################################

inttype_list = ["NI", "IC", "QA"]

dic_2_bubbles = {}
dic_n_bubbles = {}

working_path = os.getcwd()

for inttype in inttype_list :
    if inttype not in list(dic_2_bubbles.keys()) :
        dic_2_bubbles[inttype] = {}
    
    if inttype not in list(dic_n_bubbles.keys()) :
        dic_n_bubbles[inttype] = {}

    inttype_path = os.path.join(working_path, inttype)
    inttype_path = os.path.join(inttype_path, "results")

    for file in os.listdir(inttype_path) :
        file_path = os.path.join(inttype_path, file)
        file_results = open(file_path, "r")
        lines = file_results.readlines()
        file_results.close()

        first_line = lines[0].split(" ")
        count = int(first_line[0])
        png = float(first_line[5])
        size = float(first_line[7])
        cluster = float(first_line[9])

        if cluster == 0 :
            if png not in list(dic_2_bubbles[inttype].keys()) :
                dic_2_bubbles[inttype][png] = {}
            if size not in list(dic_2_bubbles[inttype][png].keys()) :
                dic_2_bubbles[inttype][png][size] = []
            
            dic_data = dic_2_bubbles[inttype][png][size]
        
        else :
            if count not in list(dic_n_bubbles[inttype].keys()) :
                dic_n_bubbles[inttype][count] = []
            
            dic_data = dic_n_bubbles[inttype][count]
        
        second_line = lines[1].split(" ")
        for i in range(count) :
            init_radius = float(second_line[i+1])
            # list format : for each bubble, [R0, t_list, R_list, Pt_list]
            dic_data.append([init_radius, [], [], []])
        
        for line in lines[3:] :
            data = line.split(" ")
            t = float(data[0])
            for i in range(count) :
                r = float(data[1 + i])
                pt = float(data[1 + count + i])
                dic_data[i][1].append(t)
                dic_data[i][2].append(r)
                dic_data[i][3].append(pt)

######### Step 2 : Plotting results #################################################################################################################

######### Initial parameters ####################

T = 10.0e-06
P0 = 0.1013e06

######### Pressure time history & radius evolution without interaction ############################

nrow = 1
ncol = 2

fig, axs = plt.subplots(nrow, ncol, figsize=((ncol*20*cm, nrow*12.5*cm)))
plt.subplots_adjust(wspace=0.5*cm, hspace=0.5*cm)

axs[0].set_title("Pressure time history (" + r"$T$ = " + "{:.1f} ".format(T*1.0e06) + r"$\mu$s)")
axs[0].set_xlabel(r"t ($\mu$s)", fontsize=27.5)
axs[0].set_xlim(xmin=0.0, xmax=60.0)
axs[0].set_ylabel(r"$p_{\infty}$/$p_{0}$ (-)", fontsize=27.5)
axs[0].grid()

t_list = np.array(dic_2_bubbles["NI"][-25325][15.0][0][1]) * 1.0e6
p_list = np.array(dic_2_bubbles["NI"][-25325][15.0][0][3]) / P0
axs[0].plot(t_list, p_list, color="black", linewidth=2.5)

axs[1].set_title("Evolution of radius without interaction")
axs[1].set_xlabel(r"t ($\mu$s)", fontsize=27.5)
axs[1].set_xlim(xmin=0.0, xmax=60.0)
axs[1].set_ylabel(r"$R$ ($\mu$m)", fontsize=27.5)
axs[1].grid()

t_list = np.array(dic_2_bubbles["NI"][-25325][15.0][0][1]) * 1.0e6
r_list = np.array(dic_2_bubbles["NI"][-25325][15.0][0][2]) * 1.0e6
axs[1].plot(t_list, r_list, color="blue", label=r"$R_{1,0}$ = 2.0 $\mu$m", linewidth=2.5)

t_list = np.array(dic_2_bubbles["NI"][-25325][15.0][1][1]) * 1.0e6
r_list = np.array(dic_2_bubbles["NI"][-25325][15.0][1][2]) * 1.0e6
axs[1].plot(t_list, r_list, color="magenta", linestyle="dashed", label=r"$R_{2,0}$ = 20.0 $\mu$m", linewidth=2.5)

axs[1].legend(loc="upper left", frameon=False)

fig.savefig("cavitationonset_pressurehistory_radiusevolutionNI.pdf", bbox_inches='tight',pad_inches=0.35)

######### Cavitation inception for one single bubble ##############################################

nrow = 1
ncol = 1

fig, ax = plt.subplots(nrow, ncol, figsize=((ncol*20*cm, nrow*12.5*cm)))
plt.subplots_adjust(wspace=0.5*cm, hspace=0.5*cm)

png_list = [-17221, -17725.5, -18353.2, -18770.3]

ax.set_title("Cavitation inception of a single bubble \n depending on " + r"$p_{ng}$/$p_{0}$ ($R_{0}$ = 2 $\mu$m)")
ax.set_xlabel(r"t ($\mu$s)", fontsize=27.5)
ax.set_xlim(xmin=10.0, xmax=60.0)
ax.set_ylabel(r"$R$ ($\mu$m)", fontsize=27.5)
ax.set_ylim(ymin=0.0, ymax=14.0)
ax.grid()

for png in png_list :
    t_list = np.array(dic_2_bubbles["NI"][png][15.0][0][1]) * 1.0e6
    r_list = np.array(dic_2_bubbles["NI"][png][15.0][0][2]) * 1.0e6

    ax.plot(t_list, r_list, color="blue", linewidth=2.5)

ax.text(25.0, 3.00, r"$-0.17$")
ax.text(31.0, 6.25, r"$-0.175$")
ax.text(28.0, 9.00, r"$-0.176$")
ax.text(21.0, 13.0, r"$-0.18$")

fig.savefig("cavitationonset_singlebubble.pdf", bbox_inches='tight',pad_inches=0.35)

######### Cavitation inception with interactions with varying distance between 2 bubbles ##########

nrow = 1
ncol = 2

fig, axs = plt.subplots(nrow, ncol, figsize=((ncol*20*cm, nrow*12.5*cm)))
plt.subplots_adjust(wspace=0.35*cm, hspace=0.5*cm)

dist_list = [10, 12, 12.1, 12.5, 15, 20]

axs[0].set_title(r"Incompressible interactions" + "\n" +r"($p_{ng}$/$p_{0}$ = -0.25)")
axs[0].set_xlabel(r"t ($\mu$s)", fontsize=27.5)
axs[0].set_xlim(xmin=0.0, xmax=60.0)
axs[0].set_ylabel(r"$R$ ($\mu$m)", fontsize=27.5)
axs[0].set_ylim(ymin=-5.0, ymax=80.0)
axs[0].grid()

for dist in dist_list :
    t_list = np.array(dic_2_bubbles["IC"][-25325][dist][0][1]) * 1.0e6
    r_list = np.array(dic_2_bubbles["IC"][-25325][dist][0][2]) * 1.0e6

    axs[0].plot(t_list, r_list, color="blue", linewidth=2.5)

t_list = np.array(dic_2_bubbles["NI"][-25325][15.0][0][1]) * 1.0e6
r_list = np.array(dic_2_bubbles["NI"][-25325][15.0][0][2]) * 1.0e6
axs[0].plot(t_list, r_list, color="blue", linewidth=2.5)

t_list = np.array(dic_2_bubbles["NI"][-25325][15.0][1][1]) * 1.0e6
r_list = np.array(dic_2_bubbles["NI"][-25325][15.0][1][2]) * 1.0e6
axs[0].plot(t_list, r_list, color="magenta", linestyle="dashed", linewidth=2.5)

axs[0].text(0.5, 5.0, r"$R_{1,0}$", color="blue")
axs[0].text(0.5, 25.0, r"$R_{2,0}$", color="magenta")

axs[0].text(36.5, 73.5, r"$\infty$", color="blue")
axs[0].text(42.5, 73.5, r"20", color="blue")
axs[0].text(54.0, 73.5, r"15", color="blue")
axs[0].text(50.0, 54.0, r"12.5", color="blue")
axs[0].text(45.0, 15.0, r"12.1", color="blue")
axs[0].text(37.0, 6.5, r"12", color="blue")
axs[0].text(25.0, -2.5, r"10", color="blue")

dist_list = [10, 11.9, 12, 15, 20]

axs[1].set_title(r"Quasi acoustic interactions" + "\n" +r"($p_{ng}$/$p_{0}$ = -0.25)")
axs[1].set_xlabel(r"t ($\mu$s)", fontsize=27.5)
axs[1].set_xlim(xmin=0.0, xmax=60.0)
# axs[1].set_ylabel(r"$R$ ($\mu$m)")
axs[1].set_ylim(ymin=-5.0, ymax=80.0)
axs[1].grid()

for dist in dist_list :
    t_list = np.array(dic_2_bubbles["QA"][-25325][dist][0][1]) * 1.0e6
    r_list = np.array(dic_2_bubbles["QA"][-25325][dist][0][2]) * 1.0e6

    axs[1].plot(t_list, r_list, color="blue", linewidth=2.5)

t_list = np.array(dic_2_bubbles["NI"][-25325][15.0][0][1]) * 1.0e6
r_list = np.array(dic_2_bubbles["NI"][-25325][15.0][0][2]) * 1.0e6
axs[1].plot(t_list, r_list, color="blue", linewidth=2.5)

t_list = np.array(dic_2_bubbles["NI"][-25325][15.0][1][1]) * 1.0e6
r_list = np.array(dic_2_bubbles["NI"][-25325][15.0][1][2]) * 1.0e6
axs[1].plot(t_list, r_list, color="magenta", linestyle="dashed", linewidth=2.5)

axs[1].text(0.5, 5.0, r"$R_{1,0}$", color="blue")
axs[1].text(0.5, 25.0, r"$R_{2,0}$", color="magenta")

axs[1].text(36.5, 73.5, r"$\infty$", color="blue")
axs[1].text(42.5, 73.5, r"20", color="blue")
axs[1].text(54.0, 73.5, r"15", color="blue")
axs[1].text(51.0, 36.0, r"12", color="blue")
axs[1].text(40.0, 10.0, r"11.9", color="blue")
axs[1].text(25.0, -2.5, r"10", color="blue")

fig.savefig("cavitationonset_varyingdistance.pdf", bbox_inches='tight',pad_inches=0.35)

######### Cavitation inception with interactions with varying pressure between 2 bubbles ##########

nrow = 1
ncol = 2

fig, axs = plt.subplots(nrow, ncol, figsize=((ncol*20*cm, nrow*12.5*cm)))
plt.subplots_adjust(wspace=0.35*cm, hspace=0.5*cm)

png_list = [-25325, -27351, -27958.8, -29377]

axs[0].set_title(r"Incompressible interactions" + "\n" +r"($\Delta x_{12}$ = 10[$R_{1,0}$ + $R_{2,0}$])")
axs[0].set_xlabel(r"t ($\mu$s)", fontsize=27.5)
axs[0].set_xlim(xmin=0.0, xmax=60.0)
axs[0].set_ylabel(r"$R_{1}$ ($\mu$m)", fontsize=27.5)
axs[0].set_ylim(ymin=0.0, ymax=40.0)
axs[0].grid()

for png in png_list :
    t_list = np.array(dic_2_bubbles["IC"][png][10.0][0][1]) * 1.0e6
    r_list = np.array(dic_2_bubbles["IC"][png][10.0][0][2]) * 1.0e6

    axs[0].plot(t_list, r_list, color="blue", linewidth=2.5)

axs[0].text(0.5, 3.0, r"$R_{1,0}$", color="blue")

axs[0].text(32.0, 37.5, r"-0.29", color="blue")
axs[0].text(45.0, 25.5, r"-0.276", color="blue")
axs[0].text(30.0, 8.0, r"-0.27", color="blue")
axs[0].text(23.0, 0.5, r"-0.25", color="blue")

png_list = [-25325, -27351, -27958.8, -27654.9, -29377]

axs[1].set_title(r"Quasi acoustic interactions" + "\n" +r"($\Delta x_{12}$ = 10[$R_{1,0}$ + $R_{2,0}$])")
axs[1].set_xlabel(r"t ($\mu$s)", fontsize=27.5)
axs[1].set_xlim(xmin=0.0, xmax=60.0)
# axs[1].set_ylabel(r"$R$ ($\mu$m)")
axs[1].set_ylim(ymin=0.0, ymax=40.0)
axs[1].grid()

for png in png_list :
    t_list = np.array(dic_2_bubbles["QA"][png][10.0][0][1]) * 1.0e6
    r_list = np.array(dic_2_bubbles["QA"][png][10.0][0][2]) * 1.0e6

    axs[1].plot(t_list, r_list, color="blue", linewidth=2.5)

axs[1].text(0.5, 3.0, r"$R_{1,0}$", color="blue")

axs[1].text(30.5, 37.5, r"-0.29", color="blue")
axs[1].text(48.5, 37.5, r"-0.276", color="blue")
axs[1].text(50.0, 28.0, r"-0.273", color="blue")
axs[1].text(36.5, 8.0, r"-0.27", color="blue")
axs[1].text(25.0, 0.5, r"-0.25", color="blue")

fig.savefig("cavitationonset_varyingpressure.pdf", bbox_inches='tight',pad_inches=0.35)

######### Cavitation inception with monodispersed simple distributions ############################

nrow = 2
ncol = 2

fig, axs = plt.subplots(nrow, ncol, figsize=((ncol*20*cm, nrow*12.5*cm)), sharex=True)
plt.subplots_adjust(wspace=0.35*cm, hspace=0.25*cm)

dic_color = {1 : "black", 2 : "red", 3 : "magenta", 4 : "blue", 8 : "green"}
nbubble_list = [1, 2, 3, 4, 8]
dic_shape = {1 : "single bubble", 2 : "Line of 2 bubbles", 3 : "3 bubbles-regular triangle", 4 : "4 bubbles-regular tetragon", 8 : "8 bubbles-regular hexaedron"}
dic_linestyle = {1 : "solid", 2 : "dashed", 3 : "-.", 4 : "dotted", 8 : "dashdot"}
dic_marker = {1 : "*", 2 : "s", 3 : "X", 4 : "^", 8 : "D"}

for i in range(2) :
    for j in range(2) :
        axs[i, j].grid()
        axs[i, j].set_xlim(xmin=10.0, xmax=60.0)

axs[0, 0].set_title(r"Incompressible interactions" + "\n" +r"($p_{ng}$/$p_{0}$ = -0.25, $\Delta x_{12}$ = 20$R_{1,0}$)")
axs[0, 1].set_title(r"Quasi acoustic interactions" + "\n" +r"($p_{ng}$/$p_{0}$ = -0.25, $\Delta x_{12}$ = 20$R_{1,0}$)")

axs[1, 0].set_xlabel(r"t ($\mu$s)", fontsize=27.5)
axs[1, 1].set_xlabel(r"t ($\mu$s)", fontsize=27.5)

axs[0, 0].set_ylabel(r"$R_{1}$ ($\mu$m)", fontsize=27.5)
axs[1, 0].set_ylabel(r"$P_{1, \infty}$/$P_{0}$ (-)", fontsize=27.5)
axs[1, 0].set_ylim(ymin=-0.255, ymax=0.0)
axs[1, 1].set_ylim(ymin=-0.255, ymax=0.0)

for nbubble in nbubble_list :
    t_list = np.array(dic_n_bubbles["IC"][nbubble][0][1]) * 1.0e6
    r_list = np.array(dic_n_bubbles["IC"][nbubble][0][2]) * 1.0e6
    p_list = np.array(dic_n_bubbles["IC"][nbubble][0][3]) / P0

    axs[0, 0].plot(t_list, r_list, color=dic_color[nbubble], label=dic_shape[nbubble], marker=dic_marker[nbubble], markevery=600, markersize=10.0, linewidth=2.5)
    axs[1, 0].plot(t_list, p_list, color=dic_color[nbubble], marker=dic_marker[nbubble], markevery=600, markersize=10.0, linewidth=2.5)

    t_list = np.array(dic_n_bubbles["QA"][nbubble][0][1]) * 1.0e6
    r_list = np.array(dic_n_bubbles["QA"][nbubble][0][2]) * 1.0e6
    p_list = np.array(dic_n_bubbles["QA"][nbubble][0][3]) / P0

    axs[0, 1].plot(t_list, r_list, color=dic_color[nbubble], marker=dic_marker[nbubble], markevery=600, markersize=10.0, linewidth=2.5)
    axs[1, 1].plot(t_list, p_list, color=dic_color[nbubble], marker=dic_marker[nbubble], markevery=600, markersize=12.5, linewidth=2.5)

axs[0, 0].legend(bbox_to_anchor=(1.1, 1.175), loc="lower center", ncol=2, frameon=False)

fig.savefig("cavitationonset_monodispersedclusters.pdf", bbox_inches='tight',pad_inches=0.35)