import os

# File designed to recover data for plotting results in bubbly screen case

file_radii = open("bubblyscreen_radii.txt", "r")
lines_radii = file_radii.readlines()
file_radii.close()

file_extremum = open("bubblyscreen_extremum.txt", "r")
lines_extremum = file_extremum.readlines()
file_extremum.close()

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
for line in lines_extremum :
    results_extremum.write(line)
results_extremum.close()