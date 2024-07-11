import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset']='stix'
plt.rcParams['font.size']=12.5

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
axs[0].set_xlabel(r"t ($\mu$s)")
axs[0].set_xlim(xmin=0.0, xmax=60.0)
axs[0].set_ylabel(r"$p_{\infty}$/$p_{0}$ (-)")
axs[0].grid()

t_list = np.array(dic_2_bubbles["NI"][-25325][15.0][0][1]) * 1.0e6
p_list = np.array(dic_2_bubbles["NI"][-25325][15.0][0][3]) / P0
axs[0].plot(t_list, p_list, color="black")

axs[1].set_title("Evolution of radius without interaction")
axs[1].set_xlabel(r"t ($\mu$s)")
axs[1].set_xlim(xmin=0.0, xmax=60.0)
axs[1].set_ylabel(r"$R$ ($\mu$m)")
axs[1].grid()

t_list = np.array(dic_2_bubbles["NI"][-25325][15.0][0][1]) * 1.0e6
r_list = np.array(dic_2_bubbles["NI"][-25325][15.0][0][2]) * 1.0e6
axs[1].plot(t_list, r_list, color="blue", label=r"$R_{1,0}$ = 2.0 $\mu$m")

t_list = np.array(dic_2_bubbles["NI"][-25325][15.0][1][1]) * 1.0e6
r_list = np.array(dic_2_bubbles["NI"][-25325][15.0][1][2]) * 1.0e6
axs[1].plot(t_list, r_list, color="magenta", label=r"$R_{2,0}$ = 20.0 $\mu$m")

axs[1].legend(loc="upper left")

fig.savefig("cavitationonset_pressurehistory_radiusevolutionNI.pdf", bbox_inches='tight',pad_inches=0.35)

######### Cavitation inception for one single bubble ##############################################

nrow = 1
ncol = 1

fig, ax = plt.subplots(nrow, ncol, figsize=((ncol*20*cm, nrow*12.5*cm)))
plt.subplots_adjust(wspace=0.5*cm, hspace=0.5*cm)

png_list = [-17221, -17725.5, -18353.2, -18770.3]

ax.set_title("Cavitation inception of a single bubble depending on " + r"$p_{ng}$/$p_{0}$ ($R_{0}$ = 2 $\mu$m)")
ax.set_xlabel(r"t ($\mu$s)")
ax.set_xlim(xmin=10.0, xmax=60.0)
ax.set_ylabel(r"$R$ ($\mu$m)")
ax.set_ylim(ymin=0.0, ymax=14.0)
ax.grid()

for png in png_list :
    t_list = np.array(dic_2_bubbles["NI"][png][15.0][0][1]) * 1.0e6
    r_list = np.array(dic_2_bubbles["NI"][png][15.0][0][2]) * 1.0e6

    ax.plot(t_list, r_list, color="blue")

ax.text(25.0, 3.50, r"$-0.17$")
ax.text(31.0, 6.25, r"$-0.175$")
ax.text(28.0, 9.00, r"$-0.176$")
ax.text(22.5, 13.0, r"$-0.18$")

fig.savefig("cavitationonset_singlebubble.pdf", bbox_inches='tight',pad_inches=0.35)

######### Cavitation inception with interactions with varying distance between 2 bubbles ##########

nrow = 1
ncol = 2

######### Cavitation inception with interactions with varying pressure between 2 bubbles ##########

nrow = 1
ncol = 2

######### Cavitation inception with monodispersed simple distributions ############################

nrow = 2
ncol = 2