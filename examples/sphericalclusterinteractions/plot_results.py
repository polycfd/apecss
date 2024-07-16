import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from math import sqrt
from scipy.interpolate import interp1d

plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset']='stix'
plt.rcParams['font.size']=10

color_names = list(mcolors.XKCD_COLORS)

cm = 1/2.54

# The goal of this scrip is to plot a 3D visualization of the computed cluster
# Later, the interest is to vizualise the temporal evolution of the cluster

######### Step 0 : Retrieving initial data ##########################################################################################################

# Count the number of bubbles simulated
count = 0
listdir = os.listdir()
for f in listdir:
    if "Bubble_" in f and os.path.isdir(f):
        count += 1

# Retrieve data for each bubble
Bubbles = []
for i in range(count):
    path = os.path.join(os.getcwd(), "Bubble_{}".format(i))
    for file in os.listdir(path): 
        if "KellerMiksis" in file : 
            Bubbles.append(np.genfromtxt("Bubble_{}/".format(i) + file, delimiter=" "))

# Retrieve location data for each bubble and determine cluster radius
file_loc = open("bubble_loc.txt","r")
loc_lines = file_loc.readlines()

Bubbles_loc = []
R_c = 0
for i in range(count):
    line = loc_lines[i+1].split(" ")
    Bubbles_loc.append([float(line[1]), float(line[2]), float(line[3])])
    R_b = sqrt((float(line[1])**2)+ (float(line[2])**2)+ (float(line[3])**2))
    if R_c < R_b :
        R_c = R_b

file_loc.close()

######### Step 1 : Functions ########################################################################################################################

### signal sampling (for signal with non uniform timestep) ########################################
# return two lists : t_signal and y_signal (with uniform timestep)
def _sample(t_signal_nuni, y_signal_nuni, delta_t_uni_exp, subref) :
    # t_signal_nuni (list/1D array) represents computation time
    # y_signal_nuni (list/1D array) represents the function y(t) we need to sample
    # delta_t_uni_exp (float/int) is the power wanted for the time step of the uniform signal
    # subref (bool) indicates if you want to refine the uniform signal obtained

    # Need to check the size of both y_signal_nuni and t_signal_nuni
    if len(t_signal_nuni) != len(y_signal_nuni) : raise NameError("No matching size for sampling")

    delta_t_uni = 10.0**(-delta_t_uni_exp)
    N_points_signal_nuni = len(t_signal_nuni)

    t_start = t_signal_nuni[0]
    t_end = t_signal_nuni[-1]
    TS = t_end - t_start
    N_points_signal_uni = int(TS/delta_t_uni + 1)

    t_signal_uni = np.linspace(t_start, t_end, N_points_signal_uni)
    t_signal_nuni_round = [0]*N_points_signal_nuni

    for i in range(0, N_points_signal_nuni):
        t_signal_nuni_round[i] = round(t_signal_nuni[i], delta_t_uni_exp)

    y_signal_nuni_round_fw = [0]*N_points_signal_nuni
    amp = -1e100
    for i in range(1, N_points_signal_nuni):
        if y_signal_nuni[i] >= 0:
            if t_signal_nuni_round[i] == t_signal_nuni_round[i-1]:
                if y_signal_nuni[i] > amp:
                    amp = y_signal_nuni[i]
                    y_signal_nuni_round_fw[i] = amp
            else:
                amp = y_signal_nuni[i]
                y_signal_nuni_round_fw[i] = amp
                amp = -1e100
                
    y_signal_nuni_round_bw = [0]*N_points_signal_nuni
    amp = -1e100
    for i in range(N_points_signal_nuni-2,-1, -1):
        if y_signal_nuni[i] >= 0:
            if t_signal_nuni_round[i] == t_signal_nuni_round[i+1]:
                if y_signal_nuni[i] > amp:
                    amp = y_signal_nuni[i]
                    y_signal_nuni_round_bw[i] = amp
            else:
                amp = y_signal_nuni[i]
                y_signal_nuni_round_bw[i] = amp
                amp = -1e100
            
    y_signal_nuni_round = [0]*N_points_signal_nuni
    for i in range(0, N_points_signal_nuni):
        y_signal_nuni_round[i] = max(y_signal_nuni_round_fw[i], y_signal_nuni_round_bw[i])

    buff = [y_signal_nuni_round[0]]
    for i in range(1, N_points_signal_nuni):
        if t_signal_nuni_round[i] == t_signal_nuni_round[i-1]:
            buff.append(y_signal_nuni_round[i])
        else:
            maxBuff = max(buff)
            for j in range(i-len(buff),i):
                y_signal_nuni_round[j] = maxBuff
            buff = [y_signal_nuni_round[i]]
            
    i = 1
    while i < N_points_signal_nuni:
        if i == len(t_signal_nuni_round) - 1:
            break
        if t_signal_nuni_round[i] == t_signal_nuni_round[i-1]:
            t_signal_nuni_round.pop(i)
            y_signal_nuni_round.pop(i)
            i -= 1
        i += 1

    fint = interp1d(t_signal_nuni_round, y_signal_nuni_round, kind='nearest')
    y_signal_uni = fint(t_signal_uni)

    t_signal = t_signal_uni
    y_signal = y_signal_uni

    if (subref == 'true'):
        ref_factor = 10
        N_points_signal_ref = (N_points_signal_uni - 1)*ref_factor + 1

        t_signal_uni_ref = np.linspace(0.0, TS, N_points_signal_ref)
        y_signal_uni_ref = []
        for i in range(1, N_points_signal_uni):
            t_ref = 0.0
            delta_t = t_signal_uni[i] - t_signal_uni[i-1]
            y = y_signal_uni[i]
            y_old = y_signal_uni[i-1]
            for j in range(0, ref_factor):
                t_ref = delta_t/float(ref_factor)*float(j)
                y_ref = 0.5*((y_old-y)*np.cos(np.pi/delta_t*t_ref) + y + y_old)
                y_signal_uni_ref.append(y_ref)
        y_signal_uni_ref.append(y_signal_uni[-1])
        t_signal = t_signal_uni_ref
        y_signal = y_signal_uni_ref
    
    return t_signal, y_signal

def plot_cluster(ax, Bubbles, Bubbles_loc, t=0, figtitle ="Bubble cluster 3D visualisation",grid="True", show="False", legend="False", radius_vis="True") :
    # # fig (figure) is a figure object where to plot the cluster
    # ax is the ax where you want to plot the 3D cluster
    # Bubbles is a list containing all data gathered by APECSS on every bubble computed
    # Bubbles_loc is a list containing the center location of each bubble
    # t represents the time at which the cluster should be displayed
    # figtitle (str) is the title of the figure
    # grid (bool) indicates if grid shall be displayed or not
    # show (bool) indicates if the plot should be shown or not (high cost for a big cluster)
    # legend (bool) indicates if the legend should be displayed or not
    # radius_vis indicates if the displayed radius should be different than the true radius (to properly see even small bubbles)
    # ax = fig.add_subplot(projection='3d')

    u = np.linspace(0, 2*np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    dic_color = {}
    j = 0

    max_l = np.max([np.max(np.array(Bubbles_loc)[:, 0]), np.max(np.array(Bubbles_loc)[:, 1]), np.max(np.array(Bubbles_loc)[:, 2])])*1e6

    time_index_list = [0 for i in range(count)]
    if t > 0 :
        for i in range(count) :
            time_list = Bubbles[i][:, 1]
            index = 0
            while time_list[index] < t and index < len(time_list) :
                index += 1
            time_index_list[i] = index

    for i in range(count) :
        r = Bubbles[i][:, 3][0]*1e6
        if r not in list(dic_color.keys()) :
            dic_color[r] = [1, color_names[j], "", 0, 0]
            j += 1
        else :
            dic_color[r][0] += 1

    sorted_radius = sorted(list(dic_color.keys()))

    for j in range(len(sorted_radius)) :
        key = sorted_radius[j]
        dic_color[key][2] = r"$R_{0}=$" + "{:.2f}".format(key) + r" $\mu$m" + " ({})".format(dic_color[key][0])
        if 0.1 * max_l < key < max_l : dic_color[key][4] = key
        else : dic_color[key][4] = max_l / (7.5 * (len(sorted_radius) - j))
    
    radius_vis_list = []
    for i in range(count) :
        if t == 0 : 
            r = Bubbles[i][:, 3][0]*1e6
            radius_vis_list.append(r)
        else :
            r = Bubbles[i][:, 3][time_index_list[i]]*1e6
            radius_vis_list.append(r)
                
    x_b, y_b, z_b = [], [], []
    for i in range(count):
        r = Bubbles[i][:, 3][0]*1e6

        x_b.append(Bubbles_loc[i][0]*1e6)
        y_b.append(Bubbles_loc[i][1]*1e6)
        z_b.append(Bubbles_loc[i][2]*1e6)

        if (radius_vis == "True") : r_visu = dic_color[r][4]
        else : r_visu = r

        r_visu = radius_vis_list[i]

        x = r_visu * np.outer(np.cos(u), np.sin(v)) + x_b[-1]
        y = r_visu * np.outer(np.sin(u), np.sin(v)) + y_b[-1]
        z = r_visu * np.outer(np.ones(np.size(u)), np.cos(v)) + z_b[-1]

        # Display bubbles with colors depending on their location in the cluster
        interval_size = 0.20e-03
        cluster_radius = 2.5e-03
        color = "grey"

        radius_center = sqrt(Bubbles_loc[i][0]**2 + Bubbles_loc[i][1]**2 + Bubbles_loc[i][2]**2)
        if radius_center < interval_size :
            color = "blue"
        elif radius_center > cluster_radius - interval_size :
            color = "black"
        elif 0.5 * (cluster_radius - interval_size) < radius_center < 0.5 * (cluster_radius + interval_size) :
            color = "red"

        if dic_color[r][3] == 0 :
            c = ax.plot_surface(x, y, z, rstride=5, cstride=5, color=color, label=dic_color[r][2])
            dic_color[r][3] = 1
        else :
            c = ax.plot_surface(x, y, z, rstride=5, cstride=5, color=color)

        c._facecolors2d=c._facecolor3d
        c._edgecolors2d=c._edgecolor3d

    if figtitle != "" : ax.set_title(figtitle)
    ax.set_xlabel(r"$x$ [$\mu$m]")
    ax.set_ylabel(r"$y$ [$\mu$m]")
    ax.set_zlabel(r"$z$ [$\mu$m]")
    if legend != "False" : ax.legend()

    if grid : plt.grid()
    if show == "True" : plt.show()

def savefig(fig, filename="Cluster", format="pdf", path=None) :
    # fig is a figure object you want to save
    # filename (str) is the future filename
    # format (str) is the format in which you want to save the figure (pdf by default)
    # path (str) is the absolute path in which you want to store the figure
    if path == None : savepath = ""
    elif os.path.isdir(path) : savepath = path
    else : savepath = ""

    fig.savefig(os.path.join(savepath, filename + "." + format), bbox_inches='tight',pad_inches=0.35)

######### Step 2 : Plot / Savefig ###################################################################################################################

## Radii evolution depending on the location in the cluster ########################################################################################
dic_radii = {0.0 : [], 0.5 : [], 1.0 : []}
interval_size = 0.20e-03

cluster_radius = 2.5e-03
fa = 50e03

max_radius_to_center = 0.0

for i in range(count) :
    radius_to_center = sqrt(Bubbles_loc[i][0]**2 + Bubbles_loc[i][1]**2 + Bubbles_loc[i][2]**2)
    if radius_to_center < 0.0 * cluster_radius + interval_size :
        dic_radii[0.0].append(i)
    elif 0.5 * (cluster_radius - interval_size) < radius_to_center < 0.5 * (cluster_radius + interval_size) :
        dic_radii[0.5].append(i)
    elif 1.0 * cluster_radius - interval_size < radius_to_center :
        dic_radii[1.0].append(i)

    if radius_to_center > max_radius_to_center :
        max_radius_to_center = radius_to_center

dic_color_key = {0.0 : "blue", 0.5 : "red", 1.0 : "black"}
dic_style_key = {0.0 : "dotted", 0.5 : "dashed", 1.0 : "solid"}

fig = plt.figure(figsize=(20*cm, 15*cm))
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel(r"$t^{*}$", fontsize=15)
ax.set_ylabel(r"$<R>/R_{0}$", fontsize=15)
ax.set_title("Spherical cluster with " + r"$N=$" + "{} bubbles".format(count), fontsize=15)
ax.grid()

dic_radius_evolution = {1.0 : [], 0.5 : [], 0.0 : []}

for k in list(dic_radii.keys()) :
    index_list = dic_radii[k]
    for i in range(len(index_list)) :
        index = index_list[i]
        t_uni, r_uni = _sample(Bubbles[index][:, 1], Bubbles[index][:, 3], 8, False)
        if i == 0 :
            avg_radii = np.array(r_uni)
        else :
            avg_radii = avg_radii + np.array(r_uni)
        if i == len(index_list) - 1 :
            dic_radius_evolution[k].append(np.array(t_uni))
    dic_radius_evolution[k].append(avg_radii / len(index_list))

label_edge_front = 0
label_edge_back = 0

for index in dic_radii[1.0] :
    if -cluster_radius < Bubbles_loc[index][0] < -cluster_radius + interval_size :
        label_edge_front = index
    elif cluster_radius - interval_size < Bubbles_loc[index][0] < cluster_radius :
        label_edge_back = index


# for k in list(dic_radii.keys()) :
#     index_list = dic_radii[k]
#     if len(index_list) > 1 :
#         for i in index_list[1:] :
#             ax.plot(Bubbles[i][:, 1]*fa, Bubbles[i][:, 3]/Bubbles[i][:, 3][0], linestyle=dic_style_key[k], color=dic_color_key[k], linewidth=1.5)
#         ax.plot(Bubbles[index_list[0]][:, 1]*fa, Bubbles[index_list[0]][:, 3]/Bubbles[i][:, 3][0], linestyle=dic_style_key[k], color=dic_color_key[k], linewidth=1.5, label=r"$r/R_{c} \approx$"+"{:.1f}".format(k))
#     else :
#         for i in index_list :
#             ax.plot(Bubbles[i][:, 1]*fa, Bubbles[i][:, 3]/Bubbles[i][:, 3][0], linestyle=dic_style_key[k], color=dic_color_key[k], linewidth=1.5, label=r"$r/R_{c} \approx$"+"{:.1f}".format(k))

ax.plot(Bubbles[label_edge_front][:, 1]*fa, Bubbles[label_edge_front][:, 3]/Bubbles[label_edge_front][:, 3][0], linestyle="solid", color="black", linewidth=1.5, label=r"$r/R_{c} \approx 1.0$ (front)")
ax.plot(Bubbles[label_edge_back][:, 1]*fa, Bubbles[label_edge_back][:, 3]/Bubbles[label_edge_back][:, 3][0], linestyle="solid", color="grey", linewidth=1.5, label=r"$r/R_{c} \approx 1.0$ (back)") 

for k in list(dic_radius_evolution.keys())[1:] :
    ax.plot(dic_radius_evolution[k][0]*fa, dic_radius_evolution[k][1]/dic_radius_evolution[k][1][0], linestyle=dic_style_key[k], color=dic_color_key[k], linewidth=1.5, label=r"$r/R_{c} \approx$"+"{:.1f}".format(k))

ax.legend(loc="upper right")
savefig(fig, filename="sphericalclusterinteractions_radiievolution")

## Cluster evolution vizualisation ##################################################################################################################
plt.rcParams['font.size']=15

fig = plt.figure(figsize=(4*15*cm,15*cm))
plt.subplots_adjust(wspace=0, hspace=0)

t_star = 0
ax = fig.add_subplot(1, 4, 1, projection='3d')
plot_cluster(ax, Bubbles, Bubbles_loc, t=t_star/fa, figtitle="", grid="False", radius_vis="False")
ax.view_init(elev=90, azim=-90)
ax.set_aspect('equal')
ax.set_box_aspect(None, zoom=1.5)
ax.set_axis_off()
ax.set_title(r"$t^{*}$ = " + "{:.2f}".format(t_star))

t_star = 0.75
ax = fig.add_subplot(1, 4, 2, projection='3d')
plot_cluster(ax, Bubbles, Bubbles_loc, t=t_star/fa, figtitle="", grid="False", radius_vis="False")
ax.view_init(elev=90, azim=-90)
ax.set_aspect('equal')
ax.set_box_aspect(None, zoom=1.5)
ax.set_axis_off()
ax.set_title(r"$t^{*}$ = " + "{:.2f}".format(t_star))

t_star = 1.0
ax = fig.add_subplot(1, 4, 3, projection='3d')
plot_cluster(ax, Bubbles, Bubbles_loc, t=t_star/fa, figtitle="", grid="False", radius_vis="False")
ax.view_init(elev=90, azim=-90)
ax.set_aspect('equal')
ax.set_box_aspect(None, zoom=1.5)
ax.set_axis_off()
ax.set_title(r"$t^{*}$ = " + "{:.2f}".format(t_star))

t_star = 1.25
ax = fig.add_subplot(1, 4, 4, projection='3d')
plot_cluster(ax, Bubbles, Bubbles_loc, t=t_star/fa, figtitle="", grid="False", radius_vis="False")
ax.view_init(elev=90, azim=-90)
ax.set_aspect('equal')
ax.set_box_aspect(None, zoom=1.5)
ax.set_axis_off()
ax.set_title(r"$t^{*}$ = " + "{:.2f}".format(t_star))

# plot_cluster(fig, Bubbles, Bubbles_loc)
savefig(fig, filename="sphericalclusterinteractions_clusterevolution", format="png")