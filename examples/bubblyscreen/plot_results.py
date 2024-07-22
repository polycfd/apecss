import os
import numpy as np
import matplotlib.pyplot as plt
from math import pi
from matplotlib.patches import Circle

fontsize = 12

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
for intype in interaction_types :
    dic_bubbly_screen[intype] = {}
    for w in excitation_w :
        dic_bubbly_screen[intype][w] = {}
        for r in ratio_D_R0 :
            dic_bubbly_screen[intype][w][r] = np.zeros((51,51),dtype=float)

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

    clb = fig.colorbar(cset, ax=axs[row, col], shrink=0.7)
    if order > 0 : 
        clb.ax.set_title(r"$|r'| \times 10^{-order}$".replace("order", str(order)))
    else :
        clb.ax.set_title(r"$|r'|$")
    colors = cset.cmap(cset.norm(cset.get_array()))
    for i in range(nbubbles_x) :
        for j in range(nbubbles_x) :
            index = nbubbles_x * i + j
            index_color = colors[index]
            x = j - space
            y = i - space
            circle = Circle((x, y), 0.4, facecolor=(index_color[0], index_color[1], index_color[2]))
            axs[row, col].add_patch(circle)

######### Step 2 : Plot results #####################################################################################################################

######### D/R0 = 400 : Comparison between QA and IC ###############################################

nrow = 2
ncol = 4
fig, axs = plt.subplots(nrow,ncol,figsize=(55*cm, 27.5*cm))
for i in range(nrow) :
    for j in range(ncol) :
        axs[i,j].set_xlabel(r"x/D")
        axs[i,j].set_xlim(xmin=-25.5,xmax=25.5)
        axs[i,j].set_ylabel(r"y/D")
        axs[i,j].set_ylim(ymin=-25.5,ymax=25.5)
        axs[i,j].set_aspect("equal")
    
### First row : QA ###
plt.figtext(0.5, 0.85, r"QA, $D/R_{0}$ = 400", fontsize=16.5, horizontalalignment="center")
plot_oscillation_distribution(fig, axs, 0, 0, 51, "QA", 0.5, 400, order=4)
plot_oscillation_distribution(fig, axs, 0, 1, 51, "QA", 0.9, 400, order=3)
plot_oscillation_distribution(fig, axs, 0, 2, 51, "QA", 1.0, 400)
plot_oscillation_distribution(fig, axs, 0, 3, 51, "QA", 1.5, 400, order=4)

### Second row : IC ###
plt.figtext(0.5, 0.5, r"IC, $D/R_{0}$ = 400", fontsize=16.5, horizontalalignment="center")
plot_oscillation_distribution(fig, axs, 1, 0, 51, "IC", 0.5, 400, order=4)
plot_oscillation_distribution(fig, axs, 1, 1, 51, "IC", 0.9, 400, order=3)
plot_oscillation_distribution(fig, axs, 1, 2, 51, "IC", 1.0, 400, order=3)
plot_oscillation_distribution(fig, axs, 1, 3, 51, "IC", 1.5, 400, order=4, clear=True)

fig.subplots_adjust(hspace=0.05*cm, wspace=0.65*cm)
fig.savefig("bubblyscreen_ComparisonQAIC.pdf", bbox_inches="tight",pad_inches=0.035)