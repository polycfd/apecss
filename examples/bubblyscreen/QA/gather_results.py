import os

# File designed to recover data for plotting results in bubbly screen case

file_radii = open("bubblyscreen_radii.txt", "r")
lines_radii = file_radii.readlines()
file_radii.close()

count_core = 0
for file in os.listdir() :
    if "bubblyscreen_extremum_" in file :
        count_core += 1

lines_extremum = []
for i in range(count_core) :
    file_extremum = open("bubblyscreen_extremum_{}.txt".format(i), "r")
    lines = file_extremum.readlines()
    file_extremum.close()
    lines_extremum.append(lines)

firstline = lines_radii[0].split(" ")
f = float(firstline[3])
p = float(firstline[5])
d = float(firstline[7])

file_name = "bubblyscreen_{:.3E}_{:.2E}_{:.1f}".format(f, p, d)

working_path = os.getcwd()
results_path = os.path.join(working_path, "results")

results_radii = open(os.path.join(results_path, file_name + "_radii.txt"), "w")
for line in lines_radii :
    results_radii.write(line)
results_radii.close()

results_extremum = open(os.path.join(results_path, file_name + "_extremum.txt"), "w")
for i in range(len(lines_extremum)) :
    if i == 0 :
        for line in lines_extremum[i] :
            results_extremum.write(line)
    else :
        for line in lines_extremum[i][2:] :
            results_extremum.write(line)
results_extremum.close()