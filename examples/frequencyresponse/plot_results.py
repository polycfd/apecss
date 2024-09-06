import os
import numpy as np
import matplotlib.pyplot as plt

######### Graphical parameters ####################################################################
fontsize = 25

plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset']='stix'
plt.rcParams['font.size']=fontsize

cm = 1/2.54

######### Step 1 : Retrieving data ################################################################

#### Original paper results #####################
paper_results_NI = {}
paper_results_IC_3 = {}
paper_results_IC_4 = {}

working_path = os.getcwd()

# No interactions (NI)
paper_path = os.path.join(working_path, "articleresults")
paper_path = os.path.join(paper_path, "NI")

for file in os.listdir(paper_path) :
    file_path = os.path.join(paper_path, file)

    text = open(file_path, "r")
    lines = text.readlines()
    text.close()

    r_init = float(lines[1].split(" : ")[-1])
    f_list = []
    r_list = []
    for i in range(3, len(lines)) :
        line = lines[i].split(" ")
        f_list.append(float(line[0].replace(",", ".")))
        r_list.append(float(line[1].replace(",", ".")))
    
    paper_results_NI[r_init] = [f_list, r_list]

# Incompressible interactions (IC) with 3 bubbles
paper_path = os.path.join(working_path, "articleresults")
paper_path = os.path.join(paper_path, "IC")
paper_path = os.path.join(paper_path, "3MBs")

for file in os.listdir(paper_path) :
    file_path = os.path.join(paper_path, file)

    text = open(file_path, "r")
    lines = text.readlines()
    text.close()

    r_init = float(lines[1].split(" : ")[-1])
    f_list = []
    r_list = []
    for i in range(3, len(lines)) :
        line = lines[i].split(" ")
        f_list.append(float(line[0].replace(",", ".")))
        r_list.append(float(line[1].replace(",", ".")))
    
    paper_results_IC_3[r_init] = [f_list, r_list]

# Incompressible interactions (IC) with 4 bubbles
paper_path = os.path.join(working_path, "articleresults")
paper_path = os.path.join(paper_path, "IC")
paper_path = os.path.join(paper_path, "4MBs")

for file in os.listdir(paper_path) :
    file_path = os.path.join(paper_path, file)

    text = open(file_path, "r")
    lines = text.readlines()
    text.close()

    r_init = float(lines[1].split(" : ")[-1])
    f_list = []
    r_list = []
    for i in range(3, len(lines)) :
        line = lines[i].split(" ")
        f_list.append(float(line[0].replace(",", ".")))
        r_list.append(float(line[1].replace(",", ".")))
    
    paper_results_IC_4[r_init] = [f_list, r_list]

#### Computation results ########################
# Key distribution : Inttype->ODE->nBubbles->pressure->distance
dic_results = {}

inttype_list = ["NI", "IC", "QA"]
for inttype in inttype_list :
    dic_results[inttype] = {}
    inttype_path = os.path.join(working_path, inttype)
    inttype_path = os.path.join(inttype_path, "results")

    for file in os.listdir(inttype_path) :
        # Finding computation results
        file_path = os.path.join(inttype_path, file)
        file = open(file_path, "r")
        lines = file.readlines()
        file.close()

        # Updating the keys in the dictionnary
        firstline = lines[0].split(" ")
        nBubbles = int(firstline[1])
        pa = float(firstline[3])
        dist = float(firstline[5])
        ode = int(firstline[7])

        # Small correction for ode number (only 0 and 1 are accepted values)
        if (ode != 0 and ode != 1) : ode = 0

        if ode not in list(dic_results[inttype].keys()) :
            dic_results[inttype][ode] = {}
        if nBubbles not in list(dic_results[inttype][ode].keys()) :
            dic_results[inttype][ode][nBubbles] = {}
        if pa not in list(dic_results[inttype][ode][nBubbles].keys()) :
            dic_results[inttype][ode][nBubbles][pa] = {}
        if dist not in list(dic_results[inttype][ode][nBubbles][pa].keys()) :
            dic_results[inttype][ode][nBubbles][pa][dist] = []
        
        # Adding computation results to the dictionnary
        for i in range(nBubbles) :
            # Format of the results in the dictionnary : for each bubble, [r_init, [f_list], [r_list]]
            r_init = float(lines[2].split(" ")[1+i].split(";")[0])
            dic_results[inttype][ode][nBubbles][pa][dist].append([r_init, [], []])
        
        for i in range(2, len(lines)) :
            data = lines[i].split(" ")
            f = float(data[0])
            for j in range(nBubbles) :
                r = float(lines[i].split(" ")[1+j].split(";")[1])
                dic_results[inttype][ode][nBubbles][pa][dist][j][1].append(f)
                dic_results[inttype][ode][nBubbles][pa][dist][j][2].append(r)

######### Step 2 : Validation plots ###############################################################

#### No interactions ####
n_row = 1
n_col = 2

fig, axs = plt.subplots(n_row, n_col, figsize=(n_col*25.0*cm, n_row*12.5*cm), sharey=True)
for i in range(n_col) :
    axs[i].set_xlabel(r"$f$ [MHz]")
    axs[i].set_xlim(xmin=0.5, xmax=10.0)
    axs[i].grid()
axs[0].set_ylabel(r"$R_{\mathrm{max}}/R_{0}$")
fig.subplots_adjust(hspace=0.0, wspace=0.05*cm)

axs[0].set_title(r"full Keller-Miksis ODE")
axs[1].set_title(r"partial Keller-Miksis ODE")

data_list = dic_results["NI"][0][4][1.2e5][5e-6]
axs[0].plot(np.array(data_list[0][1])*1e-6, np.array(data_list[0][2])/data_list[0][0], color="blue", linewidth=2.5, linestyle="solid", label=r"$R_{0}=1.0$ $\mu$m")
axs[0].plot(np.array(data_list[1][1])*1e-6, np.array(data_list[1][2])/data_list[1][0], color="magenta", linewidth=2.5, linestyle="solid", label=r"$R_{0}=0.8$ $\mu$m")
axs[0].plot(np.array(data_list[2][1])*1e-6, np.array(data_list[2][2])/data_list[2][0], color="red", linewidth=2.5, linestyle="solid", label=r"$R_{0}=0.5$ $\mu$m")
axs[0].plot(np.array(data_list[3][1])*1e-6, np.array(data_list[3][2])/data_list[3][0], color="black", linewidth=2.5, linestyle="solid", label=r"$R_{0}=1.5$ $\mu$m")

data_list = dic_results["NI"][1][4][1.2e5][5e-6]
axs[1].plot(np.array(data_list[0][1])*1e-6, np.array(data_list[0][2])/data_list[0][0], color="blue", linewidth=2.5, linestyle="solid")
axs[1].plot(np.array(data_list[1][1])*1e-6, np.array(data_list[1][2])/data_list[1][0], color="magenta", linewidth=2.5, linestyle="solid")
axs[1].plot(np.array(data_list[2][1])*1e-6, np.array(data_list[2][2])/data_list[2][0], color="red", linewidth=2.5, linestyle="solid")
axs[1].plot(np.array(data_list[3][1])*1e-6, np.array(data_list[3][2])/data_list[3][0], color="black", linewidth=2.5, linestyle="solid")

for i in range(n_col) :
    axs[i].plot(paper_results_NI[1.0][0], paper_results_NI[1.0][1], color="blue", linewidth=2.5, linestyle="dashed")
    axs[i].plot(paper_results_NI[0.8][0], paper_results_NI[0.8][1], color="magenta", linewidth=2.5, linestyle="dashed")
    axs[i].plot(paper_results_NI[0.5][0], paper_results_NI[0.5][1], color="red", linewidth=2.5, linestyle="dashed")
    axs[i].plot(paper_results_NI[1.5][0], paper_results_NI[1.5][1], color="black", linewidth=2.5, linestyle="dashed")

axs[0].legend(loc="upper right", frameon=False)
fig.savefig("frequencyresponse_comparison_NI.pdf", bbox_inches='tight',pad_inches=0.35)

#### Incompressible interactions ####
n_row = 2
n_col = 2

fig, axs = plt.subplots(n_row, n_col, figsize=(n_col*25.0*cm, n_row*12.5*cm), sharey=True, sharex=True)
for i in range(n_row) :
    for j in range(n_col) :
        axs[i, j].grid()
        if i == 1 :
            axs[i, j].set_xlabel(r"$f$ [MHz]")
            axs[i, j].set_xlim(xmin=0.5, xmax=10.0)
        if j == 0 :
            axs[i, j].set_ylabel(r"$R_{\mathrm{max}}/R_{0}$")
fig.subplots_adjust(hspace=0.1*cm, wspace=0.05*cm)

axs[0, 0].set_title(r"full Keller-Miksis ODE")
axs[0, 1].set_title(r"partial Keller-Miksis ODE")

# 3 MBs
data_list = dic_results["IC"][0][3][1.2e5][5e-6]
axs[0, 0].plot(np.array(data_list[0][1])*1e-6, np.array(data_list[0][2])/data_list[0][0], color="blue", linewidth=2.5, linestyle="solid", label=r"$R_{0}=1.0$ $\mu$m")
axs[0, 0].plot(np.array(data_list[1][1])*1e-6, np.array(data_list[1][2])/data_list[1][0], color="magenta", linewidth=2.5, linestyle="solid", label=r"$R_{0}=0.8$ $\mu$m")
axs[0, 0].plot(np.array(data_list[2][1])*1e-6, np.array(data_list[2][2])/data_list[2][0], color="red", linewidth=2.5, linestyle="solid", label=r"$R_{0}=0.5$ $\mu$m")

data_list = dic_results["IC"][1][3][1.2e5][5e-6]
axs[0, 1].plot(np.array(data_list[0][1])*1e-6, np.array(data_list[0][2])/data_list[0][0], color="blue", linewidth=2.5, linestyle="solid")
axs[0, 1].plot(np.array(data_list[1][1])*1e-6, np.array(data_list[1][2])/data_list[1][0], color="magenta", linewidth=2.5, linestyle="solid")
axs[0, 1].plot(np.array(data_list[2][1])*1e-6, np.array(data_list[2][2])/data_list[2][0], color="red", linewidth=2.5, linestyle="solid")

# 4 MBs
data_list = dic_results["IC"][0][4][1.2e5][5e-6]
axs[1, 0].plot(np.array(data_list[0][1])*1e-6, np.array(data_list[0][2])/data_list[0][0], color="blue", linewidth=2.5, linestyle="solid", label=r"$R_{0}=1.0$ $\mu$m")
axs[1, 0].plot(np.array(data_list[1][1])*1e-6, np.array(data_list[1][2])/data_list[1][0], color="magenta", linewidth=2.5, linestyle="solid", label=r"$R_{0}=0.8$ $\mu$m")
axs[1, 0].plot(np.array(data_list[2][1])*1e-6, np.array(data_list[2][2])/data_list[2][0], color="red", linewidth=2.5, linestyle="solid", label=r"$R_{0}=0.5$ $\mu$m")
axs[1, 0].plot(np.array(data_list[3][1])*1e-6, np.array(data_list[3][2])/data_list[3][0], color="black", linewidth=2.5, linestyle="solid", label=r"$R_{0}=1.5$ $\mu$m")

data_list = dic_results["IC"][1][4][1.2e5][5e-6]
axs[1, 1].plot(np.array(data_list[0][1])*1e-6, np.array(data_list[0][2])/data_list[0][0], color="blue", linewidth=2.5, linestyle="solid", label=r"$R_{0}=1.0$ $\mu$m")
axs[1, 1].plot(np.array(data_list[1][1])*1e-6, np.array(data_list[1][2])/data_list[1][0], color="magenta", linewidth=2.5, linestyle="solid", label=r"$R_{0}=0.8$ $\mu$m")
axs[1, 1].plot(np.array(data_list[2][1])*1e-6, np.array(data_list[2][2])/data_list[2][0], color="red", linewidth=2.5, linestyle="solid", label=r"$R_{0}=0.5$ $\mu$m")
axs[1, 1].plot(np.array(data_list[3][1])*1e-6, np.array(data_list[3][2])/data_list[3][0], color="black", linewidth=2.5, linestyle="solid", label=r"$R_{0}=1.5$ $\mu$m")

for j in range(n_col) :
    axs[0, j].plot(paper_results_IC_3[1.0][0], paper_results_IC_3[1.0][1], color="blue", linewidth=2.5, linestyle="dashed")
    axs[0, j].plot(paper_results_IC_3[0.8][0], paper_results_IC_3[0.8][1], color="magenta", linewidth=2.5, linestyle="dashed")
    axs[0, j].plot(paper_results_IC_3[0.5][0], paper_results_IC_3[0.5][1], color="red", linewidth=2.5, linestyle="dashed")

    axs[1, j].plot(paper_results_IC_4[1.0][0], paper_results_IC_4[1.0][1], color="blue", linewidth=2.5, linestyle="dashed")
    axs[1, j].plot(paper_results_IC_4[0.8][0], paper_results_IC_4[0.8][1], color="magenta", linewidth=2.5, linestyle="dashed")
    axs[1, j].plot(paper_results_IC_4[0.5][0], paper_results_IC_4[0.5][1], color="red", linewidth=2.5, linestyle="dashed")
    axs[1, j].plot(paper_results_IC_4[1.5][0], paper_results_IC_4[1.5][1], color="black", linewidth=2.5, linestyle="dashed")

axs[0, 0].legend(loc="upper right", frameon=False)
axs[1, 0].legend(loc="upper right", frameon=False)
fig.savefig("frequencyresponse_comparison_IC.pdf", bbox_inches='tight',pad_inches=0.35)

######### Step 3 : Article plots ##################################################################

######## Reproduction of Haghi's results #################################
n_row = 3
n_col = 2

fig, axs = plt.subplots(n_row, n_col, figsize=(n_col*25.0*cm, n_row*12.5*cm), sharey=True, sharex=True)
for i in range(n_row) :
    for j in range(n_col) :
        axs[i, j].grid()
        if i == 2 :
            axs[i, j].set_xlabel(r"$f$ [MHz]")
            axs[i, j].set_xlim(xmin=0.5, xmax=10.0)
        if j == 0 :
            axs[i, j].set_ylabel(r"$R_{\mathrm{max}}/R_{0}$")
fig.subplots_adjust(hspace=0.1*cm, wspace=0.05*cm)

axs[0, 0].set_title(r"(a) 3MBs")
axs[0, 1].set_title(r"(b) 4MBs")

# 3 MBs
data_list = dic_results["NI"][0][4][1.2e5][5e-6]
axs[0, 0].plot(np.array(data_list[0][1])*1e-6, np.array(data_list[0][2])/data_list[0][0], color="blue", linewidth=2.5, linestyle="solid", label=r"$R_{0}=1.0$ $\mu$m")
axs[0, 0].plot(np.array(data_list[1][1])*1e-6, np.array(data_list[1][2])/data_list[1][0], color="magenta", linewidth=2.5, linestyle="solid", label=r"$R_{0}=0.8$ $\mu$m")
axs[0, 0].plot(np.array(data_list[2][1])*1e-6, np.array(data_list[2][2])/data_list[2][0], color="red", linewidth=2.5, linestyle="solid", label=r"$R_{0}=0.5$ $\mu$m")
axs[0, 0].text(3.0, 2.25,r"No interactions")

data_list = dic_results["IC"][0][3][1.2e5][5e-6]
axs[1, 0].plot(np.array(data_list[0][1])*1e-6, np.array(data_list[0][2])/data_list[0][0], color="blue", linewidth=2.5, linestyle="solid")
axs[1, 0].plot(np.array(data_list[1][1])*1e-6, np.array(data_list[1][2])/data_list[1][0], color="magenta", linewidth=2.5, linestyle="solid")
axs[1, 0].plot(np.array(data_list[2][1])*1e-6, np.array(data_list[2][2])/data_list[2][0], color="red", linewidth=2.5, linestyle="solid")
axs[1, 0].text(3.0, 2.25,r"Incompressible interactions")

data_list = dic_results["QA"][0][3][1.2e5][5e-6]
axs[2, 0].plot(np.array(data_list[0][1])*1e-6, np.array(data_list[0][2])/data_list[0][0], color="blue", linewidth=2.5, linestyle="solid")
axs[2, 0].plot(np.array(data_list[1][1])*1e-6, np.array(data_list[1][2])/data_list[1][0], color="magenta", linewidth=2.5, linestyle="solid")
axs[2, 0].plot(np.array(data_list[2][1])*1e-6, np.array(data_list[2][2])/data_list[2][0], color="red", linewidth=2.5, linestyle="solid")
axs[2, 0].text(3.0, 2.25,r"Quasi-acoustic interactions")

# 4 MBs
data_list = dic_results["NI"][0][4][1.2e5][5e-6]
axs[0, 1].plot(np.array(data_list[0][1])*1e-6, np.array(data_list[0][2])/data_list[0][0], color="blue", linewidth=2.5, linestyle="solid", label=r"$R_{0}=1.0$ $\mu$m")
axs[0, 1].plot(np.array(data_list[1][1])*1e-6, np.array(data_list[1][2])/data_list[1][0], color="magenta", linewidth=2.5, linestyle="solid", label=r"$R_{0}=0.8$ $\mu$m")
axs[0, 1].plot(np.array(data_list[2][1])*1e-6, np.array(data_list[2][2])/data_list[2][0], color="red", linewidth=2.5, linestyle="solid", label=r"$R_{0}=0.5$ $\mu$m")
axs[0, 1].plot(np.array(data_list[3][1])*1e-6, np.array(data_list[3][2])/data_list[3][0], color="black", linewidth=2.5, linestyle="solid", label=r"$R_{0}=1.5$ $\mu$m")
axs[0, 1].text(3.0, 2.25,r"No interactions")

data_list = dic_results["IC"][0][4][1.2e5][5e-6]
axs[1, 1].plot(np.array(data_list[0][1])*1e-6, np.array(data_list[0][2])/data_list[0][0], color="blue", linewidth=2.5, linestyle="solid")
axs[1, 1].plot(np.array(data_list[1][1])*1e-6, np.array(data_list[1][2])/data_list[1][0], color="magenta", linewidth=2.5, linestyle="solid")
axs[1, 1].plot(np.array(data_list[2][1])*1e-6, np.array(data_list[2][2])/data_list[2][0], color="red", linewidth=2.5, linestyle="solid")
axs[1, 1].plot(np.array(data_list[3][1])*1e-6, np.array(data_list[3][2])/data_list[3][0], color="black", linewidth=2.5, linestyle="solid")
axs[1, 1].text(3.0, 2.25,r"Incompressible interactions")

data_list = dic_results["QA"][0][4][1.2e5][5e-6]
axs[2, 1].plot(np.array(data_list[0][1])*1e-6, np.array(data_list[0][2])/data_list[0][0], color="blue", linewidth=2.5, linestyle="solid")
axs[2, 1].plot(np.array(data_list[1][1])*1e-6, np.array(data_list[1][2])/data_list[1][0], color="magenta", linewidth=2.5, linestyle="solid")
axs[2, 1].plot(np.array(data_list[2][1])*1e-6, np.array(data_list[2][2])/data_list[2][0], color="red", linewidth=2.5, linestyle="solid")
axs[2, 1].plot(np.array(data_list[3][1])*1e-6, np.array(data_list[3][2])/data_list[3][0], color="black", linewidth=2.5, linestyle="solid")
axs[2, 1].text(3.0, 2.25,r"Quasi-acoustic interactions")

# Reference values
axs[0, 0].plot(paper_results_NI[1.0][0], paper_results_NI[1.0][1], color="blue", linewidth=2.5, linestyle="dashed")
axs[0, 0].plot(paper_results_NI[0.8][0], paper_results_NI[0.8][1], color="magenta", linewidth=2.5, linestyle="dashed")
axs[0, 0].plot(paper_results_NI[0.5][0], paper_results_NI[0.5][1], color="red", linewidth=2.5, linestyle="dashed")

axs[0, 1].plot(paper_results_NI[1.0][0], paper_results_NI[1.0][1], color="blue", linewidth=2.5, linestyle="dashed")
axs[0, 1].plot(paper_results_NI[0.8][0], paper_results_NI[0.8][1], color="magenta", linewidth=2.5, linestyle="dashed")
axs[0, 1].plot(paper_results_NI[0.5][0], paper_results_NI[0.5][1], color="red", linewidth=2.5, linestyle="dashed")
axs[0, 1].plot(paper_results_NI[1.5][0], paper_results_NI[1.5][1], color="black", linewidth=2.5, linestyle="dashed")

axs[1, 0].plot(paper_results_IC_3[1.0][0], paper_results_IC_3[1.0][1], color="blue", linewidth=2.5, linestyle="dashed")
axs[1, 0].plot(paper_results_IC_3[0.8][0], paper_results_IC_3[0.8][1], color="magenta", linewidth=2.5, linestyle="dashed")
axs[1, 0].plot(paper_results_IC_3[0.5][0], paper_results_IC_3[0.5][1], color="red", linewidth=2.5, linestyle="dashed")

axs[1, 1].plot(paper_results_IC_4[1.0][0], paper_results_IC_4[1.0][1], color="blue", linewidth=2.5, linestyle="dashed")
axs[1, 1].plot(paper_results_IC_4[0.8][0], paper_results_IC_4[0.8][1], color="magenta", linewidth=2.5, linestyle="dashed")
axs[1, 1].plot(paper_results_IC_4[0.5][0], paper_results_IC_4[0.5][1], color="red", linewidth=2.5, linestyle="dashed")
axs[1, 1].plot(paper_results_IC_4[1.5][0], paper_results_IC_4[1.5][1], color="black", linewidth=2.5, linestyle="dashed")

axs[0, 0].legend(loc="upper right", frameon=False)
axs[0, 1].legend(loc="upper right", frameon=False)
fig.savefig("frequencyresponse.pdf", bbox_inches='tight',pad_inches=0.35)