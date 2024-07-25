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
interaction_types = ["IC"]
cluster_types = ["mono", "poly"]

dic_bubbles = {}
dic_bubbles_loc = {}

for inttype in interaction_types :
    if inttype not in list(dic_bubbles.keys()) :
        dic_bubbles[inttype] = {}
    for cltype in cluster_types :
        if cltype not in list(dic_bubbles[inttype].keys()) :
            dic_bubbles[inttype][cltype] = {}
        if cltype not in list(dic_bubbles_loc.keys()) :
            dic_bubbles_loc[cltype] = {}

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
            png = float(firstline[5])

            if "mono" in file_name :
                dic_res = dic_bubbles[inttype]["mono"]
                dic_loc = dic_bubbles_loc["mono"]
            else :
                dic_res = dic_bubbles[inttype]["poly"]
                dic_loc = dic_bubbles_loc["poly"]
            
            if count not in list(dic_res.keys()) :
                dic_res[count] = {}
            if count not in list(dic_loc.keys()) :
                dic_loc[count] = {}
            
            if png not in list(dic_res[count].keys()) :
                dic_res[count][png] = []
            if png not in list(dic_loc[count].keys()) :
                dic_loc[count][png] = []
            
            dic_res[count][png].append([])
            for i in range(count) :
                init_radius = float(lines_res[1].split(" ")[1 + i])
                x = float(lines_loc[1 + i].split(" ")[1])
                y = float(lines_loc[1 + i].split(" ")[2])
                z = float(lines_loc[1 + i].split(" ")[3])

                dic_loc[count][png].append([x, y, z])
                dic_res[count][png].append([init_radius, [], []])
            
            for line in lines_res[3:] :
                data = line.split(" ")
                t = float(data[0])
                dic_res[count][png][0].append(t)

                for i in range(count) :
                    r = float(data[1 + i])
                    p = float(data[1 + count + i])
                    dic_res[count][png][i+1][1].append(r)
                    dic_res[count][png][i+1][2].append(p)

######### Step 2 : Define bubble distributions ######################################################################################################

# For monodispersed cluster, results are plotted based on bubbles locations
cluster_radius = 2.5e-03
interval_size = 0.20e-03

dic_loc_distrib_global = {}
dic_loc_distrib = {0.0 : [], 0.5 : [], 1.0 : []}

for count in list(dic_bubbles_loc["mono"].keys()) :
    if count not in list(dic_loc_distrib_global.keys()) :
        dic_loc_distrib_global[count] = {}
    
    for png in list(dic_bubbles_loc["mono"][count].keys()) :
        if png not in list(dic_loc_distrib_global[count].keys()) :
            dic_loc_distrib_global[count][png] = dic_loc_distrib

        for i in range(count) :
            radius_to_center = sqrt(dic_bubbles_loc["mono"][count][png][i][0]**2 + dic_bubbles_loc["mono"][count][png][i][1]**2 + dic_bubbles_loc["mono"][count][png][i][2]**2)
            if radius_to_center < 0.0 * cluster_radius + interval_size :
                dic_loc_distrib_global[count][png][0.0].append(i)
            elif 0.5 * (cluster_radius - interval_size) < radius_to_center < 0.5 * (cluster_radius + interval_size) :
                dic_loc_distrib_global[count][png][0.5].append(i)
            elif 1.0 * cluster_radius - interval_size < radius_to_center :
                dic_loc_distrib_global[count][png][1.0].append(i)

# For polydispersed cluster, results are plotted based on bubbles initial radii
n_quantiles = 5
quantile_name_dic = {4 : "quartile", 5 : "quintile"}
quantile_name = quantile_name_dic[n_quantiles]
quantiles_list = []
for i in range(1, n_quantiles + 1) :
    quantiles_list.append(i / n_quantiles)
quantiles_list = np.array(quantiles_list)

dic_polydisperse = {}

for inttype in interaction_types :
    dic_polydisperse[inttype] = {}
    for count in list(dic_bubbles[inttype]["poly"].keys()) :
        if count not in list(dic_polydisperse[inttype].keys()) :
            dic_polydisperse[inttype][count] = {}
        for png in list(dic_bubbles[inttype]["poly"][count].keys()) :
            if png not in list(dic_polydisperse[inttype][count].keys()) :
                dic_polydisperse[inttype][count][png] = {"quantiles" : []}
            for q in quantiles_list :
                dic_polydisperse[inttype][count][png][q] = []
            
            initial_radius_list = []
            for i in range(count) :
                initial_radius_list.append(dic_bubbles[inttype]["poly"][count][png][i][0])
            
            quantiles = np.quantile(np.array(initial_radius_list), quantiles_list)
            for q in quantiles : dic_polydisperse[inttype][count][png]["quantiles"].append(q)

            for i in range(count) :
                initial_radius = dic_bubbles[inttype]["poly"][count][png][i][0]

                index = 0
                while initial_radius > dic_polydisperse[inttype][count][png]["quantiles"][index] :
                    index += 1
                dic_polydisperse[inttype][count][png][quantiles_list[index]].append(i)

######### Step 3 : Plot results #####################################################################################################################

######### Averaged radius versus time evolution for monodispersed cluster #########################
count = 2500
png = -25325

nrow = 1
ncol = 2
fig, axs = plt.subplots(nrow,ncol,figsize=(ncol*17.5*cm, nrow*12.5*cm))
for i in range(ncol) :
    axs[i].set_xlabel(r"t ($\mu$s)", fontsize=15)
    axs[i].set_ylabel(r"$<R> / R_{0}$ (-)", fontsize=15)
    axs[i].set_xlim(xmin=10.0, xmax=60.0)
    axs[i].grid()

axs[0].set_title(r"Incompressible interactions" + "\n" +r"($N = 2500$, $R_{0}=10.0 \ \mu m$, $p_{ng} / p_{0} = -0.25$)", fontsize=15)
axs[1].set_title(r"Quasi acoustic interactions" + "\n" +r"($N = 2500$, $R_{0}=10.0 \ \mu m$, $p_{ng} / p_{0} = -0.25$)", fontsize=15)

dic_color_loc = {0.0 : "blue", 0.5 : "red", 1.0 : "black"}
dic_lines_loc = {0.0 : "dotted", 0.5 : "dashed", 1.0 : "solid"}

for k in dic_loc_distrib_global[count][png].keys() :
    index_list = dic_loc_distrib_global[count][png][k]

    # Incompressible interactions
    t_list = np.array(dic_bubbles["IC"]["mono"][count][png][0])

    r_list_0 = np.array(dic_bubbles["IC"]["mono"][count][png][index_list[0] + 1][1])
    init_r_0 = dic_bubbles["IC"]["mono"][count][png][index_list[0] + 1][0]
    avg_radius = r_list_0 / init_r_0

    for i in range(1, len(index_list)) :
        index = index_list[i]

        r_list = np.array(dic_bubbles["IC"]["mono"][count][png][index + 1][1])
        init_r = dic_bubbles["IC"]["mono"][count][png][index + 1][0]
        avg_radius = avg_radius + np.array(r_list) / init_r

    avg_radius = (1 / len(index_list)) * avg_radius
    
    axs[0].plot(t_list*1.0e6, avg_radius, linewidth=1.5, linestyle=dic_lines_loc[k], color=dic_color_loc[k], label=r"$r/R_{C} \approx$" + " {:.1f}".format(k))

    # Quasi acoustic interactions

axs[0].legend(bbox_to_anchor=(1.15, 1.175), loc="center", fontsize=15, ncol=3, frameon=False)
fig.savefig("sphericalclustercavitationonset_mono_radiusevolution.pdf", bbox_inches='tight',pad_inches=0.35)

######### Averaged radius versus time evolution for polydispersed cluster #########################

count = 2500
png = -25325

nrow = 1
ncol = 2
fig, axs = plt.subplots(nrow,ncol,figsize=(ncol*17.5*cm, nrow*12.5*cm))
for i in range(ncol) :
    axs[i].set_xlabel(r"t ($\mu$s)", fontsize=15)
    axs[i].set_ylabel(r"$<R/R_{0}>$ (-)", fontsize=15)
    axs[i].set_xlim(xmin=10.0, xmax=60.0)
    axs[i].grid()

axs[0].set_title(r"Incompressible interactions" + "\n" +r"($N = 2500$, $R_{0, ref}=10.0 \ \mu m$, $\overline{m} = 0$, $\varsigma = 0.7$, $p_{ng} / p_{0} = -0.25$)", fontsize=15)
axs[1].set_title(r"Quasi acoustic interactions" + "\n" +r"($N = 2500$, $R_{0, ref}=10.0 \ \mu m$, $\overline{m} = 0$, $\varsigma = 0.7$, $p_{ng} / p_{0} = -0.25$)", fontsize=15)

color_list = ["black", "magenta", "blue", "green", "red"]
dic_color_q = {}
dic_label_q = {}
for i in range(n_quantiles) :
    quantile = quantiles_list[i]
    dic_color_q[quantile] = color_list[i]

    dic_label_q[quantile] = "{}th-{}".format(i+1, quantile_name)
    if i == 0 :
        dic_label_q[quantile] = "1st-{}".format(quantile_name)
    elif i == 1 :
        dic_label_q[quantile] = "2nd-{}".format(quantile_name)
    elif i == 2 :
        dic_label_q[quantile] = "3rd-{}".format(quantile_name)

# Incompressible interactions
for q in list(dic_polydisperse["IC"][count][png].keys())[1:] :
    index_list = dic_polydisperse["IC"][count][png][q]

    t_list = np.array(dic_bubbles["IC"]["poly"][count][png][0])

    r_list_0 = np.array(dic_bubbles["IC"]["poly"][count][png][index_list[0] + 1][1])
    init_r_0 = dic_bubbles["IC"]["poly"][count][png][index_list[0] + 1][0]
    avg_radius = r_list_0 / init_r_0

    for i in range(1, len(index_list)) :
        index = index_list[i]

        r_list = np.array(dic_bubbles["IC"]["poly"][count][png][index + 1][1])
        init_r = dic_bubbles["IC"]["poly"][count][png][index + 1][0]
        avg_radius = avg_radius + np.array(r_list) / init_r

    avg_radius = (1 / len(index_list)) * avg_radius
    
    axs[0].plot(t_list*1.0e6, avg_radius, linewidth=1.5, linestyle="solid", color=dic_color_q[q], label=dic_label_q[q])

axs[0].legend(bbox_to_anchor=(1.1, 1.175), loc="center", fontsize=15, ncol=n_quantiles, frameon=False)
fig.savefig("sphericalclustercavitationonset_poly_radiusevolution.pdf", bbox_inches='tight',pad_inches=0.35)