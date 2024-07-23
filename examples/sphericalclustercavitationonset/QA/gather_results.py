import os

# File designed to recover data for plotting results in cavitation onset case

count_core = 0
for file in os.listdir() :
    if "onset_results_" in file :
        count_core += 1

lines_onset = []
for i in range(count_core) :
    file_onset = open("onset_results_{}.txt".format(i), "r")
    lines = file_onset.readlines()
    file_onset.close()
    lines_onset.append(lines)

firstline = lines_onset[0][0].split(" ")
count = int(firstline[0])
png = float(firstline[5])
cl_distrib = int(firstline[7])

if cl_distrib == 0 :
    # monodispersed spherical cluster
    file_name = "mono_{}_{:.4E}.txt".format(count, png)

else :
    # polydispersed spherical cluster
    file_name = "poly_{}_{:.4E}.txt".format(count, png)

working_path = os.getcwd()
results_path = os.path.join(working_path, "results")
file_results = open(os.path.join(results_path, file_name), "w")

results_onset = open(os.path.join(results_path, file_name), "w")
# Header
for line in lines_onset[0][:3] :
    results_onset.write(line)

# Combining with the good format the values located in the files associated with the cores used for computation
for t in range(3, len(lines_onset[0])) :
    time = lines_onset[0][t].split(" ")[0]
    r_list = []
    p_list = []

    for i in range(len(lines_onset)) :
        data = lines_onset[i][t].split(" ")
        nbubble_core = int((len(data) - 1) / 2)
        r_list_file = data[1:nbubble_core+1]
        p_list_file = data[nbubble_core+1:2*nbubble_core+1]

        for r in r_list_file :
            r_list.append(r)
        
        for p in p_list_file :
            if "\n" in p :
                p = p.split("\n")[0]
            p_list.append(p)
    
    results_onset.write(time)
    for rad in r_list :
        results_onset.write(" {}".format(rad))
    for pres in p_list :
        results_onset.write(" {}".format(pres))
    results_onset.write("\n")

results_onset.close()

file_loc = open("bubble_loc.txt", "r")
lines_loc = file_loc.readlines()
file_loc.close()

if cl_distrib == 0 :
    # monodispersed spherical cluster
    file_name_loc = "mono_{}_{:.4E}_loc.txt".format(count, png)

else :
    # polydispersed spherical cluster
    file_name_loc = "poly_{}_{:.4E}_loc.txt".format(count, png)

file_loc_results = open(os.path.join(results_path, file_name_loc), "w")

for line in lines_loc :
    file_loc_results.write(line)

file_loc_results.close()