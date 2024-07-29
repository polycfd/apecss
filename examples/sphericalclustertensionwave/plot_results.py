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
            p1 = float(firstline[5])

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
dic_loc_distrib = {0.0 : [], 0.5 : [], 1.0 : []}

for count in list(dic_bubbles_loc["mono"].keys()) :
    if count not in list(dic_loc_distrib_global.keys()) :
        dic_loc_distrib_global[count] = {}
    
    for p1 in list(dic_bubbles_loc["mono"][count].keys()) :
        if p1 not in list(dic_loc_distrib_global[count].keys()) :
            dic_loc_distrib_global[count][p1] = dic_loc_distrib

        for i in range(count) :
            radius_to_center = sqrt(dic_bubbles_loc["mono"][count][p1][i][0]**2 + dic_bubbles_loc["mono"][count][p1][i][1]**2 + dic_bubbles_loc["mono"][count][p1][i][2]**2)
            if radius_to_center < 0.0 * cluster_radius + interval_size :
                dic_loc_distrib_global[count][p1][0.0].append(i)
            elif 0.5 * (cluster_radius - interval_size) < radius_to_center < 0.5 * (cluster_radius + interval_size) :
                dic_loc_distrib_global[count][p1][0.5].append(i)
            elif 1.0 * cluster_radius - interval_size < radius_to_center :
                dic_loc_distrib_global[count][p1][1.0].append(i)

dic_loc_label = {1.0 : "{:.1f}".format(1.0 - interval_size/cluster_radius) + r" $\leq$ $r/R_{C}$ $\leq$ " + "{:.1f}".format(1.0),
                 0.5 : "{:.1f}".format(0.5 - interval_size/cluster_radius) + r" $\leq$ $r/R_{C}$ $\leq$ " + "{:.1f}".format(0.5 + interval_size/cluster_radius),
                 0.0 : "{:.1f}".format(0.0) + r" $\leq$ $r/R_{C}$ $\leq$ " + "{:.1f}".format(interval_size/cluster_radius)}

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
        for p1 in list(dic_bubbles[inttype]["poly"][count].keys()) :
            if p1 not in list(dic_polydisperse[inttype][count].keys()) :
                dic_polydisperse[inttype][count][p1] = {"quantiles" : []}
            for q in quantiles_list :
                dic_polydisperse[inttype][count][p1][q] = []
            
            initial_radius_list = []
            for i in range(count) :
                initial_radius_list.append(dic_bubbles[inttype]["poly"][count][p1][i][0])
            
            quantiles = np.quantile(np.array(initial_radius_list), quantiles_list)
            for q in quantiles : dic_polydisperse[inttype][count][p1]["quantiles"].append(q)

            for i in range(count) :
                initial_radius = dic_bubbles[inttype]["poly"][count][p1][i][0]

                index = 0
                while initial_radius > dic_polydisperse[inttype][count][p1]["quantiles"][index] :
                    index += 1
                dic_polydisperse[inttype][count][p1][quantiles_list[index]].append(i)

######### Step 3 : Plot results #####################################################################################################################

#### General parameters #####
rho_l = 1000.0
r_ref = 2.0e-6
p0 = 1.0e5
mu = 0.0
sigma = 0.7

######### Radius and minimum distance between bubbles in both configurations ######################
p1 = -3.0e4

for cluster in cluster_types :
    dic_loc = dic_bubbles_loc[cluster]

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
            
        print("{} cluster with N={} bubbles, minimum distance between bubbles equal {:.2f} microns".format(cluster, count, min_dist*1.0e6))

fig, ax = plt.subplots(1, 1, figsize=(17.5*cm, 12.5*cm))
ax.set_xlabel(r"$<R> / R_{0,ref}$ (-)", fontsize=15)
ax.grid()

for count in list(dic_bubbles["NI"]["poly"].keys()) :
    init_r_list = []
    for i in range(count) :
        init_r_list.append(dic_bubbles["NI"]["poly"][count][p1][i][0] / r_ref)
    count, bins, ignored = ax.hist(init_r_list, 100, density=True, align='mid', label="{} bubbles".format(count))

x = np.linspace(min(bins), max(bins), 10000)
pdf = (np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2)) / (x * sigma * np.sqrt(2 * np.pi)))
ax.plot(x, pdf, color="red", linewidth=1.5, linestyle="solid")

ax.legend(loc="upper right")
fig.savefig("sphericalclustertensionwave_poly_radiidistribution.pdf", bbox_inches='tight',pad_inches=0.35)

######### Averaged radius versus time evolution for monodispersed cluster #########################
count = 250
p1 = -3.0e4
t_i = r_ref * sqrt(rho_l / (p0 - p1))

nrow = 1
ncol = 3
fig, axs = plt.subplots(nrow,ncol,figsize=(ncol*17.5*cm, nrow*12.5*cm))
for i in range(ncol) :
    axs[i].set_xlabel(r"$ t / t_{i}$ (-)", fontsize=15)
    axs[i].set_ylabel(r"$<R> / R_{0}$ (-)", fontsize=15)
    # axs[i].set_xlim(xmin=10.0, xmax=60.0)
    axs[i].grid()

axs[0].set_title(r"No interaction" + "\n" +r"($N = 250$, $R_{0}=2.0 \ \mu $m, $p_{1} = -3.0 \times 10^{4} $ Pa)", fontsize=15)
axs[1].set_title(r"Incompressible interactions" + "\n" +r"($N = 250$, $R_{0}=2.0 \ \mu $m, $p_{1} = -3.0 \times 10^{4} $ Pa)", fontsize=15)
axs[2].set_title(r"Quasi acoustic interactions" + "\n" +r"($N = 250$, $R_{0}=2.0 \ \mu $m, $p_{1} = -3.0 \times 10^{4} $ Pa)", fontsize=15)

dic_color_loc = {0.0 : "blue", 0.5 : "red", 1.0 : "black"}
dic_lines_loc = {0.0 : "dotted", 0.5 : "dashed", 1.0 : "solid"}

for k in dic_loc_distrib_global[count][p1].keys() :
    index_list = dic_loc_distrib_global[count][p1][k]

    # No interaction
    t_list = (1/t_i) * np.array(dic_bubbles["NI"]["mono"][count][p1][0])

    r_list_0 = np.array(dic_bubbles["NI"]["mono"][count][p1][index_list[0] + 1][1])
    init_r_0 = dic_bubbles["NI"]["mono"][count][p1][index_list[0] + 1][0]
    avg_radius = r_list_0 / init_r_0

    for i in range(1, len(index_list)) :
        index = index_list[i]

        r_list = np.array(dic_bubbles["NI"]["mono"][count][p1][index + 1][1])
        init_r = dic_bubbles["NI"]["mono"][count][p1][index + 1][0]
        avg_radius = avg_radius + np.array(r_list) / init_r

    avg_radius = (1 / len(index_list)) * avg_radius
    
    axs[0].plot(t_list, avg_radius, linewidth=1.5, linestyle=dic_lines_loc[k], color=dic_color_loc[k])

    # Incompressible interactions
    t_list = (1/t_i) * np.array(dic_bubbles["IC"]["mono"][count][p1][0])

    r_list_0 = np.array(dic_bubbles["IC"]["mono"][count][p1][index_list[0] + 1][1])
    init_r_0 = dic_bubbles["IC"]["mono"][count][p1][index_list[0] + 1][0]
    avg_radius = r_list_0 / init_r_0

    for i in range(1, len(index_list)) :
        index = index_list[i]

        r_list = np.array(dic_bubbles["IC"]["mono"][count][p1][index + 1][1])
        init_r = dic_bubbles["IC"]["mono"][count][p1][index + 1][0]
        avg_radius = avg_radius + np.array(r_list) / init_r

    avg_radius = (1 / len(index_list)) * avg_radius
    
    axs[1].plot(t_list, avg_radius, linewidth=1.5, linestyle=dic_lines_loc[k], color=dic_color_loc[k], label=dic_loc_label[k])

    # # Quasi acoustic interactions
    # t_list = (1/t_i) * np.array(dic_bubbles["QA"]["mono"][count][p1][0])

    # r_list_0 = np.array(dic_bubbles["QA"]["mono"][count][p1][index_list[0] + 1][1])
    # init_r_0 = dic_bubbles["QA"]["mono"][count][p1][index_list[0] + 1][0]
    # avg_radius = r_list_0 / init_r_0

    # for i in range(1, len(index_list)) :
    #     index = index_list[i]

    #     r_list = np.array(dic_bubbles["QA"]["mono"][count][p1][index + 1][1])
    #     init_r = dic_bubbles["QA"]["mono"][count][p1][index + 1][0]
    #     avg_radius = avg_radius + np.array(r_list) / init_r

    # avg_radius = (1 / len(index_list)) * avg_radius
    
    # axs[2].plot(t_list, avg_radius, linewidth=1.5, linestyle=dic_lines_loc[k], color=dic_color_loc[k])


axs[1].legend(bbox_to_anchor=(0.5, 1.175), loc="center", fontsize=15, ncol=3, frameon=False)
fig.savefig("sphericalclustertensionwave_mono_radiusevolution.pdf", bbox_inches='tight',pad_inches=0.35)

######### Averaged radius versus time evolution for polydispersed cluster #########################

count = 250
p1 = -3.0e4
t_i = r_ref * sqrt(rho_l / (p0 - p1))

nrow = 1
ncol = 3
fig, axs = plt.subplots(nrow,ncol,figsize=(ncol*17.5*cm, nrow*12.5*cm))
for i in range(ncol) :
    axs[i].set_xlabel(r"$ t / t_{i}$ (-)", fontsize=15)
    axs[i].set_ylabel(r"$<R/R_{0}>$ (-)", fontsize=15)
    # axs[i].set_xlim(xmin=10.0, xmax=60.0)
    axs[i].grid()

axs[0].set_title(r"No interaction" + "\n" +r"($N = 250$, $R_{0, ref}=10.0 \ \mu $m, $\overline{m} = 0$, $\varsigma = 0.7$, $p_{1} = -3.0 \times 10^{4} $ Pa)", fontsize=15)
axs[1].set_title(r"Incompressible interactions" + "\n" +r"($N = 250$, $R_{0, ref}=10.0 \ \mu $m, $\overline{m} = 0$, $\varsigma = 0.7$, $p_{1} = -3.0 \times 10^{4} $ Pa)", fontsize=15)
axs[2].set_title(r"Quasi acoustic interactions" + "\n" +r"($N = 250$, $R_{0, ref}=10.0 \ \mu $m, $\overline{m} = 0$, $\varsigma = 0.7$, $p_{1} = -3.0 \times 10^{4} $ Pa)", fontsize=15)

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

# No interaction
for q in list(dic_polydisperse["NI"][count][p1].keys())[1:] :
    index_list = dic_polydisperse["NI"][count][p1][q]

    t_list = (1/t_i) * np.array(dic_bubbles["NI"]["poly"][count][p1][0])

    r_list_0 = np.array(dic_bubbles["NI"]["poly"][count][p1][index_list[0] + 1][1])
    init_r_0 = dic_bubbles["NI"]["poly"][count][p1][index_list[0] + 1][0]
    avg_radius = r_list_0 / init_r_0

    for i in range(1, len(index_list)) :
        index = index_list[i]

        r_list = np.array(dic_bubbles["NI"]["poly"][count][p1][index + 1][1])
        init_r = dic_bubbles["NI"]["poly"][count][p1][index + 1][0]
        avg_radius = avg_radius + np.array(r_list) / init_r

    avg_radius = (1 / len(index_list)) * avg_radius
    
    axs[0].plot(t_list, avg_radius, linewidth=1.5, linestyle="solid", color=dic_color_q[q])

# Incompressible interactions
for q in list(dic_polydisperse["IC"][count][p1].keys())[1:] :
    index_list = dic_polydisperse["IC"][count][p1][q]

    t_list = (1/t_i) * np.array(dic_bubbles["IC"]["poly"][count][p1][0])

    r_list_0 = np.array(dic_bubbles["IC"]["poly"][count][p1][index_list[0] + 1][1])
    init_r_0 = dic_bubbles["IC"]["poly"][count][p1][index_list[0] + 1][0]
    avg_radius = r_list_0 / init_r_0

    for i in range(1, len(index_list)) :
        index = index_list[i]

        r_list = np.array(dic_bubbles["IC"]["poly"][count][p1][index + 1][1])
        init_r = dic_bubbles["IC"]["poly"][count][p1][index + 1][0]
        avg_radius = avg_radius + np.array(r_list) / init_r

    avg_radius = (1 / len(index_list)) * avg_radius
    
    axs[1].plot(t_list, avg_radius, linewidth=1.5, linestyle="solid", color=dic_color_q[q], label=dic_label_q[q])

axs[1].legend(bbox_to_anchor=(0.5, 1.175), loc="center", fontsize=15, ncol=n_quantiles, frameon=False)
fig.savefig("sphericalclustertensionwave_poly_radiusevolution.pdf", bbox_inches='tight',pad_inches=0.35)

######### Averaged pressure infinity versus time evolution for monodispersed clusters #############

count = 250
p1 = -3.0e4
p0 = 1.0e5
t_i = r_ref * sqrt(rho_l / (p0 - p1))

nrow = 1
ncol = 3
fig, axs = plt.subplots(nrow,ncol,figsize=(ncol*17.5*cm, nrow*12.5*cm))
for i in range(ncol) :
    axs[i].set_xlabel(r"$ t / t_{i}$ (-)", fontsize=15)
    axs[i].set_ylabel(r"$<p_{\infty}> / p_{0}$ (-)", fontsize=15)
    # axs[i].set_xlim(xmin=10.0, xmax=60.0)
    axs[i].grid()

axs[0].set_title(r"No interaction" + "\n" +r"($N = 250$, $R_{0}=2.0 \ \mu $m, $p_{1} = -3.0 \times 10^{4} $ Pa)", fontsize=15)
axs[1].set_title(r"Incompressible interactions" + "\n" +r"($N = 250$, $R_{0}=2.0 \ \mu $m, $p_{1} = -3.0 \times 10^{4} $ Pa)", fontsize=15)
axs[2].set_title(r"Quasi acoustic interactions" + "\n" +r"($N = 250$, $R_{0}=2.0 \ \mu $m, $p_{1} = -3.0 \times 10^{4} $ Pa)", fontsize=15)

dic_color_loc = {0.0 : "blue", 0.5 : "red", 1.0 : "black"}
dic_lines_loc = {0.0 : "dotted", 0.5 : "dashed", 1.0 : "solid"}

for k in dic_loc_distrib_global[count][p1].keys() :
    index_list = dic_loc_distrib_global[count][p1][k]

    # No interaction
    t_list = (1/t_i) * np.array(dic_bubbles["NI"]["mono"][count][p1][0])

    p_list_0 = np.array(dic_bubbles["NI"]["mono"][count][p1][index_list[0] + 1][2])
    avg_pressure = p_list_0 / p0

    for i in range(1, len(index_list)) :
        index = index_list[i]

        p_list = np.array(dic_bubbles["NI"]["mono"][count][p1][index + 1][2])
        avg_pressure = avg_pressure + np.array(p_list) / p0

    avg_pressure = (1 / len(index_list)) * avg_pressure
    
    axs[0].plot(t_list, avg_pressure, linewidth=1.5, linestyle=dic_lines_loc[k], color=dic_color_loc[k])

    # Incompressible interactions
    t_list = (1/t_i) * np.array(dic_bubbles["IC"]["mono"][count][p1][0])

    p_list_0 = np.array(dic_bubbles["IC"]["mono"][count][p1][index_list[0] + 1][2])
    avg_pressure = p_list_0 / p0

    for i in range(1, len(index_list)) :
        index = index_list[i]

        p_list = np.array(dic_bubbles["IC"]["mono"][count][p1][index + 1][2])
        avg_pressure = avg_pressure + np.array(p_list) / p0

    avg_pressure = (1 / len(index_list)) * avg_pressure
    
    axs[1].plot(t_list, avg_pressure, linewidth=1.5, linestyle=dic_lines_loc[k], color=dic_color_loc[k], label=dic_loc_label[k])

    # # Quasi acoustic interactions
    # t_list = (1/t_i) * np.array(dic_bubbles["QA"]["mono"][count][p1][0])

    # p_list_0 = np.array(dic_bubbles["QA"]["mono"][count][p1][index_list[0] + 1][2])
    # avg_pressure = p_list_0 / p0

    # for i in range(1, len(index_list)) :
    #     index = index_list[i]

    #     p_list = np.array(dic_bubbles["QA"]["mono"][count][p1][index + 1][2])
    #     avg_pressure = avg_pressure + np.array(p_list) / p0

    # avg_pressure = (1 / len(index_list)) * avg_pressure
    
    # axs[2].plot(t_list, avg_pressure, linewidth=1.5, linestyle=dic_lines_loc[k], color=dic_color_loc[k])

axs[1].legend(bbox_to_anchor=(0.5, 1.175), loc="center", fontsize=15, ncol=3, frameon=False)
fig.savefig("sphericalclustertensionwave_mono_pressureevolution.pdf", bbox_inches='tight',pad_inches=0.35)

######### Averaged pressure infinity versus time evolution for monodispersed clusters #############

count = 250
p1 = -3.0e4
p0 = 1.0e5
t_i = r_ref * sqrt(rho_l / (p0 - p1))

nrow = 1
ncol = 3
fig, axs = plt.subplots(nrow,ncol,figsize=(ncol*17.5*cm, nrow*12.5*cm))
for i in range(ncol) :
    axs[i].set_xlabel(r"$ t / t_{i}$ (-)", fontsize=15)
    axs[i].set_ylabel(r"$<p_{\infty}>/p_{0}$ (-)", fontsize=15)
    # axs[i].set_xlim(xmin=10.0, xmax=60.0)
    axs[i].grid()

axs[0].set_title(r"No interaction" + "\n" +r"($N = 250$, $R_{0, ref}=10.0 \ \mu $m, $\overline{m} = 0$, $\varsigma = 0.7$, $p_{1} = -3.0 \times 10^{4} $ Pa)", fontsize=15)
axs[1].set_title(r"Incompressible interactions" + "\n" +r"($N = 250$, $R_{0, ref}=10.0 \ \mu $m, $\overline{m} = 0$, $\varsigma = 0.7$, $p_{1} = -3.0 \times 10^{4} $ Pa)", fontsize=15)
axs[2].set_title(r"Quasi acoustic interactions" + "\n" +r"($N = 250$, $R_{0, ref}=10.0 \ \mu $m, $\overline{m} = 0$, $\varsigma = 0.7$, $p_{1} = -3.0 \times 10^{4} $ Pa)", fontsize=15)

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

# No interaction
for q in list(dic_polydisperse["NI"][count][p1].keys())[1:] :
    index_list = dic_polydisperse["NI"][count][p1][q]

    t_list = (1/t_i) * np.array(dic_bubbles["NI"]["poly"][count][p1][0])

    p_list_0 = np.array(dic_bubbles["NI"]["poly"][count][p1][index_list[0] + 1][2])
    avg_pressure = p_list_0 / p0

    for i in range(1, len(index_list)) :
        index = index_list[i]

        p_list = np.array(dic_bubbles["NI"]["poly"][count][p1][index + 1][2])
        avg_pressure = avg_pressure + np.array(p_list) / p0

    avg_pressure = (1 / len(index_list)) * avg_pressure
    
    axs[0].plot(t_list, avg_pressure, linewidth=1.5, linestyle="solid", color=dic_color_q[q])

# Incompressible interactions
for q in list(dic_polydisperse["IC"][count][p1].keys())[1:] :
    index_list = dic_polydisperse["IC"][count][p1][q]

    t_list = (1/t_i) * np.array(dic_bubbles["IC"]["poly"][count][p1][0])

    p_list_0 = np.array(dic_bubbles["IC"]["poly"][count][p1][index_list[0] + 1][2])
    avg_pressure = p_list_0 / p0

    for i in range(1, len(index_list)) :
        index = index_list[i]

        p_list = np.array(dic_bubbles["IC"]["poly"][count][p1][index + 1][2])
        avg_pressure = avg_pressure + np.array(p_list) / p0

    avg_pressure = (1 / len(index_list)) * avg_pressure
    
    axs[1].plot(t_list, avg_pressure, linewidth=1.5, linestyle="solid", color=dic_color_q[q], label=dic_label_q[q])

axs[1].legend(bbox_to_anchor=(0.5, 1.175), loc="center", fontsize=15, ncol=n_quantiles, frameon=False)
fig.savefig("sphericalclustertensionwave_poly_pressureevolution.pdf", bbox_inches='tight',pad_inches=0.35)