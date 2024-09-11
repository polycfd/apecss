import os
import numpy as np
import matplotlib.pyplot as plt
from math import pi
from matplotlib.patches import Circle

fontsize = 20.5

plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset']='stix'
plt.rcParams['font.size']=fontsize

cm = 1/2.54

# This file is designed to gather data in order to reproduce the average radius evolution
# This exemple is a reproduction from the test from Fig.5 in 2021 Fan's article
# DOI : 10.1121/10.0008905

######### Step 1 : Retrieving data ##################################################################################################################
interaction_types = ["IC", "QA"]
excitation_w = [0.5, 0.9, 1.0, 1.5]
ratio_D_R0 = [400]

dic_bubbly_screen = {}
dic_radii = {}
dic_radii_index = {}
index = 1301
for intype in interaction_types :
    dic_bubbly_screen[intype] = {}
    dic_radii[intype] = {}
    dic_radii_index[intype] = {}
    for w in excitation_w :
        dic_bubbly_screen[intype][w] = {}
        dic_radii[intype][w] = {}
        dic_radii_index[intype][w] = {}
        for r in ratio_D_R0 :
            dic_bubbly_screen[intype][w][r] = np.zeros((51,51),dtype=float)
            dic_radii[intype][w][r] = [[] for i in range(51*51 + 1)]
            dic_radii_index[intype][w][r] = [[], []]

for intype in interaction_types :
    path = os.path.join(os.getcwd(), intype)
    path = os.path.join(path, "results")
    for file in os.listdir(path) :
        if "_extremum" in file :
            data = open(os.path.join(path, file), "r")
            lines = data.readlines()
            data.close()

            firstline = lines[0].split(" ")
            w0 = float(firstline[1])
            fa = float(firstline[3])
            pa = float(firstline[5])
            ratio = int(float(firstline[7]))

            R0 = float(lines[2].split(" ")[4])
            D = ratio * R0

            w = float("{:.1f}".format(2 * pi * fa / w0))

            if pa == 10**2 :
                for line in lines[2:] :
                    x = float(line.split(" ")[1])
                    y = float(line.split(" ")[2])
                    r_amp = (float(line.split(" ")[6]) - float(line.split(" ")[5])) / (2 * float(line.split(" ")[4]))

                    i = int((y/D)) + 25
                    j = int((x/D)) + 25

                    dic_bubbly_screen[intype][w][ratio][i][j] = r_amp
        
        if "_radii" in file :
            print(file)
            data = open(os.path.join(path, file), "r")
            lines = data.readlines()
            data.close()

            firstline = lines[0].split(" ")
            w0 = float(firstline[1])
            fa = float(firstline[3])
            pa = float(firstline[5])
            ratio = int(float(firstline[7]))

            R0 = float(lines[2].split(" ")[1])
            D = ratio * R0

            w = float("{:.1f}".format(2 * pi * fa / w0))

            if pa == 10**2 :
                # for i in range(2, len(lines), 25) :
                #     line = lines[i]
                #     t = float(line.split(" ")[0])
                #     dic_radii[intype][w][ratio][0].append(t)

                #     for i in range(51*51) :
                #         r = (float(line.split(" ")[i + 1]) - R0) / R0
                #         dic_radii[intype][w][ratio][i + 1].append(r)
                
                for i in range(2, len(lines)) :
                    line = lines[i]
                    t = float(line.split(" ")[0])
                    dic_radii_index[intype][w][ratio][0].append(t)
                    # r = (float(line.split(" ")[index]) - R0) / R0
                    r = float(line.split(" ")[index])
                    dic_radii_index[intype][w][ratio][1].append(r)

######### Functions #################################################################################################################################

def plot_oscillation_distribution(fig, axs, row, col, nbubbles_x, inttype,  ratio_w, ratio_d, order=0, clear=False) :
    space = int(0.5 * (nbubbles_x - 1))
    X = [-space + j for j in range(nbubbles_x)]
    Y = [-space + i for i in range(nbubbles_x)]
    X, Y = np.meshgrid(X, Y)

    axs[row, col].set_title(r"$\omega/\omega_{0}$ = " + "{:.1f}".format(ratio_w))
    cset = plt.pcolormesh(X, Y, dic_bubbly_screen[inttype][ratio_w][ratio_d]*10**(order))

    # Clearing options to not show pcolormesh (only used for the last plotted distribution on each figure)
    if clear : 
        axs[row, col].clear()
        axs[row, col].set_xlabel(r"x/D")
        axs[row, col].set_xlim(xmin=-(space + 0.5),xmax=(space + 0.5))
        axs[row, col].set_ylabel(r"y/D")
        axs[row, col].set_ylim(ymin=-(space + 0.5),ymax=(space + 0.5))
        axs[row, col].set_title(r"$\omega/\omega_{0}$ = " + "{:.1f}".format(ratio_w))
        axs[row, col].set_aspect("equal")

    clb = fig.colorbar(cset, ax=axs[row, col], pad=0.025, shrink=0.665)
    if order > 0 : 
        clb.ax.set_title(r"$|r'| \times 10^{-order}$".replace("order", str(order)), fontsize=fontsize+1.0)
    else :
        clb.ax.set_title(r"$|r'|$", fontsize=fontsize+1.0)
    colors = cset.cmap(cset.norm(cset.get_array()))
    for i in range(nbubbles_x) :
        for j in range(nbubbles_x) :
            # index = nbubbles_x * i + j
            index_color = colors[j][i]
            x = j - space
            y = i - space
            circle = Circle((x, y), 0.4, facecolor=(index_color[0], index_color[1], index_color[2]))
            axs[row, col].add_patch(circle)

def identify_oscillation_amplitude(t_list, r_list) :
    # t_list and r_list must be 1D-arrays
    # Determine the evolution of the oscillation amplitude by decomposing the entry signal
    r_mean = np.mean(r_list)
    r_mmean_list = r_list - r_mean

    index = 0
    index_positive_list = []
    while index < len(t_list) :
        if r_mmean_list[index] > 0 :
            index_positive_list.append(index)
        index += 1
    
    periods_index = []
    list_index = 1
    periods_index.append([index_positive_list[0]])
    while list_index < len(index_positive_list) :
        if index_positive_list[list_index] - index_positive_list[list_index - 1] > 1 :
            # A step larger than one in index means a new period has been reached
            periods_index[-1].append(index_positive_list[list_index] - 1)
            periods_index.append([index_positive_list[list_index]])
        list_index += 1
    
    if len(periods_index[-1]) == 1 :
        periods_index = periods_index[:len(periods_index)-1]
    
    t_amp_list = []
    r_amp_list = []
    for i in range(len(periods_index)) :
        start_index = periods_index[i][0]
        stop_index = periods_index[i][1]

        mean_index = int(0.5 * (start_index + stop_index))
        t_amp_list.append(t_list[mean_index])
        
        r_slice_list = r_list[start_index : stop_index + 1]
        r_max = np.max(r_slice_list)
        r_min = np.min(r_slice_list)
        r_amp_list.append(0.5 * (r_max - r_min))
    
    return t_amp_list, r_amp_list

######### Step 2 : Plot results #####################################################################################################################

######### D/R0 = 400 : Comparison between QA and IC ###############################################

nrow = 2
ncol = 4
fig, axs = plt.subplots(nrow,ncol,figsize=(55*cm, 27.5*cm))
for i in range(nrow) :
    for j in range(ncol) :
        axs[i,j].set_xlabel(r"y/D")
        axs[i,j].set_xlim(xmin=-25.5,xmax=25.5)
        axs[i,j].set_ylabel(r"z/D")
        axs[i,j].set_ylim(ymin=-25.5,ymax=25.5)
        axs[i,j].set_aspect("equal")
    
### First row : QA ###
plt.figtext(0.5, 0.875, r"(a) QA", fontsize=25.0, horizontalalignment="center")
plot_oscillation_distribution(fig, axs, 0, 0, 51, "QA", 0.5, 400, order=4)
plot_oscillation_distribution(fig, axs, 0, 1, 51, "QA", 0.9, 400, order=3)
plot_oscillation_distribution(fig, axs, 0, 2, 51, "QA", 1.0, 400, order=2)
plot_oscillation_distribution(fig, axs, 0, 3, 51, "QA", 1.5, 400, order=4)

### Second row : IC ###
plt.figtext(0.5, 0.475, r"(b) IC", fontsize=25.0, horizontalalignment="center")
plot_oscillation_distribution(fig, axs, 1, 0, 51, "IC", 0.5, 400, order=4)
plot_oscillation_distribution(fig, axs, 1, 1, 51, "IC", 0.9, 400, order=3)
plot_oscillation_distribution(fig, axs, 1, 2, 51, "IC", 1.0, 400, order=3)
plot_oscillation_distribution(fig, axs, 1, 3, 51, "IC", 1.5, 400, order=4, clear=True)

fig.subplots_adjust(hspace=0.01*cm, wspace=0.95*cm)
fig.savefig("bubblyscreen_ComparisonQAIC.pdf", bbox_inches="tight",pad_inches=0.035)

######### D/R0 = 400 : Radius evolution ###########################################################

fig, ax = plt.subplots(1, 1, figsize=(27.5*cm, 12.5*cm))
ax.set_xlabel(r"$t$ [$\mu$s]")
ax.set_xlim(xmin=0.0, xmax=60.0)
ax.set_ylabel(r"$|r'|$")
ax.set_yscale("log")
ax.set_ylim(ymin=10**(-4), ymax=5.0e-2)
ax.grid()

### QA ###
t_list, r_list = identify_oscillation_amplitude(np.array(dic_radii_index["QA"][1.0][400][0]), np.array(dic_radii_index["QA"][1.0][400][1]))
ax.plot(np.array(t_list)*1.0e6, np.array(r_list)/R0, linewidth=2.0, color="black", linestyle="dashed", label="QA")
# ax.plot(np.array(dic_radii_index["QA"][1.0][400][0])*1.0e6, np.array(dic_radii_index["QA"][1.0][400][1]), linewidth=2.0, color="black", linestyle="dashed", label="QA")
# ax.plot(np.array(dic_radii["QA"][1.0][400][0])*1.0e6, np.array(dic_radii["QA"][1.0][400][1301]), linewidth=2.0, color="black", linestyle="dashed", label="QA")
# ax.plot(np.array(dic_radii["QA"][1.0][400][0])*1.0e6, np.array(dic_radii["QA"][1.0][400][1]), linewidth=1.5, marker="s", markersize=5.0, markevery=1500, color="black", linestyle="dashed", label="QA, corner")

### IC ###
t_list, r_list = identify_oscillation_amplitude(np.array(dic_radii_index["IC"][1.0][400][0]), np.array(dic_radii_index["IC"][1.0][400][1]))
ax.plot(np.array(t_list)*1.0e6, np.array(r_list)/R0, linewidth=2.0, color="red", linestyle="solid", label="IC")
# ax.plot(np.array(dic_radii_index["IC"][1.0][400][0])*1.0e6, np.array(dic_radii_index["IC"][1.0][400][1]), linewidth=2.0, color="red", linestyle="solid", label="IC")
# ax.plot(np.array(dic_radii["IC"][1.0][400][0])*1.0e6, np.array(dic_radii["IC"][1.0][400][1301]), linewidth=2.0, color="red", linestyle="solid", label="IC")
# ax.plot(np.array(dic_radii["IC"][1.0][400][0])*1.0e6, np.array(dic_radii["IC"][1.0][400][1]), linewidth=1.5, marker="s", markersize=5.0, markevery=1500, color="red", linestyle="solid", label="IC, corner")

ax.legend(bbox_to_anchor=(0.5, 1.05), loc="center", frameon=False, ncol=2)
fig.savefig("bubblyscreen_radiievolution.pdf", bbox_inches="tight",pad_inches=0.035)