import os

# File designed to recover data for plotting results in cavitation onset case

count_core = 0
for file in os.listdir() :
    if "tension_results_" in file :
        count_core += 1

lines_tension = []
for i in range(count_core) :
    file_tension = open("tension_results_{}.txt".format(i), "r")
    lines = file_tension.readlines()
    file_tension.close()
    lines_tension.append(lines)

firstline = lines_tension[0][0].split(" ")
count = int(firstline[0])
p1 = float(firstline[5])

file_name = "mono_{}_{:.4E}.txt".format(count, p1)

working_path = os.getcwd()
results_path = os.path.join(working_path, "results")
if not os.path.exists(results_path) :
    os.makedirs(results_path)
file_results = open(os.path.join(results_path, file_name), "w")

results_tension = open(os.path.join(results_path, file_name), "w")
# Header
for line in lines_tension[0][:3] :
    results_tension.write(line)

# Combining with the good format the values located in the files associated with the cores used for computation
for t in range(3, len(lines_tension[0])) :
    time = lines_tension[0][t].split(" ")[0]
    r_list = []
    p_list = []

    for i in range(len(lines_tension)) :
        data = lines_tension[i][t].split(" ")
        nbubble_core = int((len(data) - 1) / 2)
        r_list_file = data[1:nbubble_core+1]
        p_list_file = data[nbubble_core+1:2*nbubble_core+1]

        for r in r_list_file :
            r_list.append(r)
        
        for p in p_list_file :
            if "\n" in p :
                p = p.split("\n")[0]
            p_list.append(p)
    
    results_tension.write(time)
    for rad in r_list :
        results_tension.write(" {}".format(rad))
    for pres in p_list :
        results_tension.write(" {}".format(pres))
    results_tension.write("\n")

results_tension.close()

file_loc = open("bubble_loc.txt", "r")
lines_loc = file_loc.readlines()
file_loc.close()

file_name_loc = "mono_{}_{:.4E}_loc.txt".format(count, p1)

file_loc_results = open(os.path.join(results_path, file_name_loc), "w")

for line in lines_loc :
    file_loc_results.write(line)

file_loc_results.close()