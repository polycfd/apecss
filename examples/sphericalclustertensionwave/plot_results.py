import os
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

fontsize = 15

plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset']='stix'
plt.rcParams['font.size']=fontsize

cm = 1/2.54

# File designed to plot the results for the cavitation onset test case with a spherical cluster

######### Step 1 : Retrieving data ##################################################################################################################
interaction_types = ["NI", "IC", "QA"]
cluster_types = ["mono"]

dic_bubbles = {}
dic_bubbles_loc = {}

for inttype in interaction_types :
    if inttype not in list(dic_bubbles.keys()) :
        dic_bubbles[inttype] = {}
    if inttype not in list(dic_bubbles_loc.keys()) :
        dic_bubbles_loc[inttype] = {}
    for cltype in cluster_types :
        if cltype not in list(dic_bubbles[inttype].keys()) :
            dic_bubbles[inttype][cltype] = {}
        if cltype not in list(dic_bubbles_loc[inttype].keys()) :
            dic_bubbles_loc[inttype][cltype] = {}

working_path = os.getcwd()
for inttype in interaction_types :
    interaction_path = os.path.join(working_path, inttype)
    interaction_path = os.path.join(interaction_path, "results")
    for file in os.listdir(interaction_path) :
        if "_loc.txt" in file :
            file_name = file.split("_loc.txt")[0]

            file_loc_path = os.path.join(interaction_path, file)
            file_loc = open(file_loc_path, "r")
            lines_loc = file_loc.readlines()
            file_loc.close()

            file_path = os.path.join(interaction_path, file_name + ".txt")
            file_res = open(file_path, "r")
            lines_res = file_res.readlines()
            file_res.close()

            firstline = lines_res[0].split(" ")
            count = int(firstline[0])
            p1 = float(firstline[5])

            dic_res = dic_bubbles[inttype]["mono"]
            dic_loc = dic_bubbles_loc[inttype]["mono"]
            
            if count not in list(dic_res.keys()) :
                dic_res[count] = {}
            if count not in list(dic_loc.keys()) :
                dic_loc[count] = {}
            
            if p1 not in list(dic_res[count].keys()) :
                dic_res[count][p1] = []
            if p1 not in list(dic_loc[count].keys()) :
                dic_loc[count][p1] = []
            
            dic_res[count][p1].append([])
            for i in range(count) :
                init_radius = float(lines_res[1].split(" ")[1 + i])
                x = float(lines_loc[1 + i].split(" ")[1])
                y = float(lines_loc[1 + i].split(" ")[2])
                z = float(lines_loc[1 + i].split(" ")[3])

                dic_loc[count][p1].append([x, y, z])
                dic_res[count][p1].append([init_radius, [], []])
            
            for line in lines_res[3:] :
                data = line.split(" ")
                t = float(data[0])
                dic_res[count][p1][0].append(t)

                for i in range(count) :
                    r = float(data[1 + i])
                    p = float(data[1 + count + i])
                    dic_res[count][p1][i+1][1].append(r)
                    dic_res[count][p1][i+1][2].append(p)

######### Step 2 : Define bubble distributions ######################################################################################################

# For monodispersed cluster, results are plotted based on bubbles locations
cluster_radius = 232.0e-6
interval_size = cluster_radius / 10

dic_loc_distrib_global = {}

for inttype in list(dic_bubbles_loc.keys()) :
    if inttype not in list(dic_loc_distrib_global.keys()) :
        dic_loc_distrib_global[inttype] = {}

    for count in list(dic_bubbles_loc[inttype]["mono"].keys()) :
        if count not in list(dic_loc_distrib_global[inttype].keys()) :
            dic_loc_distrib_global[inttype][count] = {}
        
        for p1 in list(dic_bubbles_loc[inttype]["mono"][count].keys()) :
            if p1 not in list(dic_loc_distrib_global[inttype][count].keys()) :
                dic_loc_distrib_global[inttype][count][p1] = {0.0 : [], 0.25 : [], 0.5 : [], 0.75 : [], 1.0 : []}

            for i in range(count) :
                radius_to_center = sqrt(dic_bubbles_loc[inttype]["mono"][count][p1][i][0]**2 + dic_bubbles_loc[inttype]["mono"][count][p1][i][1]**2 + dic_bubbles_loc[inttype]["mono"][count][p1][i][2]**2)
                if radius_to_center < 0.0 * cluster_radius + interval_size :
                    dic_loc_distrib_global[inttype][count][p1][0.0].append(i)
                elif 0.25 * cluster_radius - 0.5 * interval_size < radius_to_center < 0.25 * cluster_radius + 0.5 * interval_size :
                    dic_loc_distrib_global[inttype][count][p1][0.25].append(i)
                elif 0.5 * (cluster_radius - interval_size) < radius_to_center < 0.5 * (cluster_radius + interval_size) :
                    dic_loc_distrib_global[inttype][count][p1][0.5].append(i)
                elif 0.75 * cluster_radius - 0.5 * interval_size < radius_to_center < 0.75 * cluster_radius + 0.5 * interval_size :
                    dic_loc_distrib_global[inttype][count][p1][0.75].append(i)
                elif 1.0 * cluster_radius - interval_size < radius_to_center :
                    dic_loc_distrib_global[inttype][count][p1][1.0].append(i)

dic_loc_label = {1.0 : "{:.1f}".format(1.0 - interval_size/cluster_radius) + r" $\leq$ $r/R_{C}$ $\leq$ " + "{:.1f}".format(1.0),
                 0.75 : "{:.1f}".format(0.75 - interval_size/cluster_radius) + r" $\leq$ $r/R_{C}$ $\leq$ " + "{:.1f}".format(0.75 + interval_size/cluster_radius),
                 0.5 : "{:.1f}".format(0.5 - interval_size/cluster_radius) + r" $\leq$ $r/R_{C}$ $\leq$ " + "{:.1f}".format(0.5 + interval_size/cluster_radius),
                 0.25 : "{:.1f}".format(0.25 - interval_size/cluster_radius) + r" $\leq$ $r/R_{C}$ $\leq$ " + "{:.1f}".format(0.25 + interval_size/cluster_radius),
                 0.0 : "{:.1f}".format(0.0) + r" $\leq$ $r/R_{C}$ $\leq$ " + "{:.1f}".format(interval_size/cluster_radius)}

######### Step 3 : Plot results #####################################################################################################################

#### General parameters #####
rho_l = 1000.0
r_ref = 2.0e-6
p0 = 1.0e5

######### Radius and minimum distance between bubbles in both configurations ######################
p1 = -3.0e4

for cluster in cluster_types :
    dic_loc = dic_bubbles_loc["NI"][cluster]

    for count in list(dic_loc.keys()) :
        min_dist = cluster_radius
        loc_list = dic_loc[count][p1]
        for i in range(count) :
            x_i = loc_list[i][0]
            y_i = loc_list[i][1]
            z_i = loc_list[i][2]
            min_dist_bubble = cluster_radius

            for j in range(count) :
                if (i != j) :
                    x_j = loc_list[j][0]
                    y_j = loc_list[j][1]
                    z_j = loc_list[j][2]
                    dist = sqrt((x_i - x_j)**2 + (y_i - y_j)**2 + (z_i - z_j)**2)
                    if dist < min_dist_bubble :
                        min_dist_bubble = dist
            
            if min_dist_bubble < min_dist :
                min_dist = min_dist_bubble
            
        print("{} cluster with N={} bubbles, minimum distance between bubbles equals {:.2f} microns".format(cluster, count, min_dist*1.0e6))

######### Averaged radius evolution ###############
count = 250
p1 = -3.0e4
t_i = r_ref * sqrt(rho_l / (p0 - p1))

nrow = 1
ncol = 3
fig, axs = plt.subplots(nrow,ncol,figsize=(ncol*17.5*cm, nrow*12.5*cm))
for i in range(ncol) :
    axs[i].set_xlabel(r"$ t$ [$\mu$s]", fontsize=27.5)
    if i == 0 : axs[i].set_ylabel(r"$<R> / R_{0}$", fontsize=27.5)
    axs[i].set_xlim(xmin=0.0, xmax=15.0)
    axs[i].tick_params(axis="both", labelsize=25)
    axs[i].grid()

plt.subplots_adjust(wspace=0.35*cm)

axs[0].set_title(r"No interactions", fontsize=27.5, y=1.025)
axs[1].set_title(r"Incompressible interactions", fontsize=27.5, y=1.025)
axs[2].set_title(r"Quasi-acoustic interactions", fontsize=27.5, y=1.025)

dic_color_loc = {0.0 : "blue", 0.25 : "magenta", 0.5 : "red", 0.75 : "green", 1.0 : "black"}
dic_lines_loc = {0.0 : "dotted", 0.25 : "dashed", 0.5 : "dashed", 0.75 : "dashed", 1.0 : "solid"}

for k in [0.0, 0.5, 1.0] :
    # No interactions
    index_list = dic_loc_distrib_global["NI"][count][p1][k]
    t_list = 1.0e06 * np.array(dic_bubbles["NI"]["mono"][count][p1][0])

    r_list_0 = np.array(dic_bubbles["NI"]["mono"][count][p1][index_list[0] + 1][1])
    init_r_0 = dic_bubbles["NI"]["mono"][count][p1][index_list[0] + 1][0]
    avg_radius = r_list_0 / init_r_0

    for i in range(1, len(index_list)) :
        index = index_list[i]

        r_list = np.array(dic_bubbles["NI"]["mono"][count][p1][index + 1][1])
        init_r = dic_bubbles["NI"]["mono"][count][p1][index + 1][0]
        avg_radius = avg_radius + np.array(r_list) / init_r

    avg_radius = (1 / len(index_list)) * avg_radius
    
    axs[0].plot(t_list, avg_radius, linewidth=3.0, linestyle=dic_lines_loc[k], color=dic_color_loc[k])

    # Incompressible interactions
    index_list = dic_loc_distrib_global["IC"][count][p1][k]
    t_list = 1.0e06 * np.array(dic_bubbles["IC"]["mono"][count][p1][0])

    r_list_0 = np.array(dic_bubbles["IC"]["mono"][count][p1][index_list[0] + 1][1])
    init_r_0 = dic_bubbles["IC"]["mono"][count][p1][index_list[0] + 1][0]
    avg_radius = r_list_0 / init_r_0

    for i in range(1, len(index_list)) :
        index = index_list[i]

        r_list = np.array(dic_bubbles["IC"]["mono"][count][p1][index + 1][1])
        init_r = dic_bubbles["IC"]["mono"][count][p1][index + 1][0]
        avg_radius = avg_radius + np.array(r_list) / init_r

    avg_radius = (1 / len(index_list)) * avg_radius
    
    axs[1].plot(t_list, avg_radius, linewidth=3.0, linestyle=dic_lines_loc[k], color=dic_color_loc[k], label=dic_loc_label[k])

    # Quasi-acoustic interactions
    index_list = dic_loc_distrib_global["QA"][count][p1][k]
    t_list = 1.0e06 * np.array(dic_bubbles["QA"]["mono"][count][p1][0])

    r_list_0 = np.array(dic_bubbles["QA"]["mono"][count][p1][index_list[0] + 1][1])
    init_r_0 = dic_bubbles["QA"]["mono"][count][p1][index_list[0] + 1][0]
    avg_radius = r_list_0 / init_r_0

    for i in range(1, len(index_list)) :
        index = index_list[i]

        r_list = np.array(dic_bubbles["QA"]["mono"][count][p1][index + 1][1])
        init_r = dic_bubbles["QA"]["mono"][count][p1][index + 1][0]
        avg_radius = avg_radius + np.array(r_list) / init_r

    avg_radius = (1 / len(index_list)) * avg_radius
    
    axs[2].plot(t_list, avg_radius, linewidth=3.0, linestyle=dic_lines_loc[k], color=dic_color_loc[k])


axs[1].legend(bbox_to_anchor=(0.5, 1.175), loc="center", fontsize=26.5, ncol=3, frameon=False)
fig.savefig("sphericalclustertensionwave_radiusevolution.pdf", bbox_inches='tight',pad_inches=0.35)

######### Averaged pressure infinity evolution  ###

count = 250
p1 = -3.0e4
p0 = 1.0e5
t_i = r_ref * sqrt(rho_l / (p0 - p1))

nrow = 1
ncol = 3
fig, axs = plt.subplots(nrow,ncol,figsize=(ncol*17.5*cm, nrow*12.5*cm))
for i in range(ncol) :
    axs[i].set_xlabel(r"$ t$ [$\mu$s]", fontsize=27.5)
    if i == 0 : axs[i].set_ylabel(r"$<p_{\infty}> / p_{0}$", fontsize=27.5)
    axs[i].set_xlim(xmin=0.0, xmax=15.0)
    axs[i].tick_params(axis="both", labelsize=25)
    axs[i].grid()

plt.subplots_adjust(wspace=0.35*cm)

axs[0].set_title(r"No interactions", fontsize=27.5, y=1.025)
axs[1].set_title(r"Incompressible interactions", fontsize=27.5, y=1.025)
axs[2].set_title(r"Quasi-acoustic interactions", fontsize=27.5, y=1.025)

dic_color_loc = {0.0 : "blue", 0.25 : "magenta", 0.5 : "red", 0.75 : "green", 1.0 : "black"}
dic_lines_loc = {0.0 : "dotted", 0.25 : "dashed", 0.5 : "dashed", 0.75 : "dashed", 1.0 : "solid"}

zm_IC = axs[1].inset_axes([0.35, 0.35, 0.60, 0.60])
zm_IC.grid()
zm_IC.set_xlim(xmin=0.0,xmax=5.0)
zm_IC.set_ylim(ymin=-0.5,ymax=2.0)
zm_IC.set_xticks([0, 2, 4])
zm_IC.set_yticks([-0.267, 0.0, 1.0, 2.0])
zm_IC.set_yticklabels([r"$p_{\mathrm{C}}$", 0, 1, 2])
zm_IC.tick_params(axis="both", labelsize=25)

zm_QA = axs[2].inset_axes([0.30, 0.4, 0.65, 0.55])
zm_QA.grid()
zm_QA.set_xlim(xmin=0.0,xmax=5.0)
zm_QA.set_ylim(ymin=-0.5,ymax=2.0)
zm_QA.set_xticks([0, 2, 4])
zm_QA.set_yticks([-0.267, 0.0, 1.0, 2.0])
zm_QA.set_yticklabels([r"$p_{\mathrm{C}}$", 0, 1, 2])
zm_QA.tick_params(axis="both", labelsize=25)

for k in [0.0, 0.5, 1.0] :
    # No interactions
    index_list = dic_loc_distrib_global["NI"][count][p1][k]
    t_list = 1.0e06 * np.array(dic_bubbles["NI"]["mono"][count][p1][0])

    p_list_0 = np.array(dic_bubbles["NI"]["mono"][count][p1][index_list[0] + 1][2])
    avg_pressure = p_list_0 / p0

    for i in range(1, len(index_list)) :
        index = index_list[i]

        p_list = np.array(dic_bubbles["NI"]["mono"][count][p1][index + 1][2])
        avg_pressure = avg_pressure + np.array(p_list) / p0

    avg_pressure = (1 / len(index_list)) * avg_pressure
    
    axs[0].plot(t_list, avg_pressure, linewidth=3.0, linestyle=dic_lines_loc[k], color=dic_color_loc[k])

    # Incompressible interactions
    index_list = dic_loc_distrib_global["IC"][count][p1][k]
    t_list = 1.0e06 * np.array(dic_bubbles["IC"]["mono"][count][p1][0])

    p_list_0 = np.array(dic_bubbles["IC"]["mono"][count][p1][index_list[0] + 1][2])
    avg_pressure = p_list_0 / p0

    for i in range(1, len(index_list)) :
        index = index_list[i]

        p_list = np.array(dic_bubbles["IC"]["mono"][count][p1][index + 1][2])
        avg_pressure = avg_pressure + np.array(p_list) / p0

    avg_pressure = (1 / len(index_list)) * avg_pressure
    
    axs[1].plot(t_list, avg_pressure, linewidth=3.0, linestyle=dic_lines_loc[k], color=dic_color_loc[k], label=dic_loc_label[k])

    zm_IC.plot(t_list, avg_pressure, linewidth=3.0, linestyle=dic_lines_loc[k], color=dic_color_loc[k])

    # Quasi-acoustic interactions
    index_list = dic_loc_distrib_global["QA"][count][p1][k]
    t_list = 1.0e06 * np.array(dic_bubbles["QA"]["mono"][count][p1][0])

    p_list_0 = np.array(dic_bubbles["QA"]["mono"][count][p1][index_list[0] + 1][2])
    avg_pressure = p_list_0 / p0

    for i in range(1, len(index_list)) :
        index = index_list[i]

        p_list = np.array(dic_bubbles["QA"]["mono"][count][p1][index + 1][2])
        avg_pressure = avg_pressure + np.array(p_list) / p0

    avg_pressure = (1 / len(index_list)) * avg_pressure
    
    axs[2].plot(t_list, avg_pressure, linewidth=3.0, linestyle=dic_lines_loc[k], color=dic_color_loc[k])

    zm_QA.plot(t_list, avg_pressure, linewidth=3.0, linestyle=dic_lines_loc[k], color=dic_color_loc[k])

axs[1].legend(bbox_to_anchor=(0.5, 1.175), loc="center", fontsize=26.5, ncol=3, frameon=False)
fig.savefig("sphericalclustertensionwave_pressureevolution.pdf", bbox_inches='tight',pad_inches=0.35)