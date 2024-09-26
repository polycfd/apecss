import os

# File designed to recover data for plotting results in cavitation onset case

file = open("tension_results.txt", "r")
lines = file.readlines()
file.close()

count = int(lines[0].split(" ")[0])
p1 = float(lines[0].split(" ")[5])

file_name = "mono_{}_{:.4E}.txt".format(count, p1)

working_path = os.getcwd()
results_path = os.path.join(working_path, "results")
if not os.path.exists(results_path) :
    os.makedirs(results_path)
file_results = open(os.path.join(results_path, file_name), "w")

index = 0
while index < len(lines) and "nan" not in lines[index] :
    file_results.write(lines[index])
    index += 1

file_results.close()

file_loc = open("bubble_loc.txt", "r")
lines_loc = file_loc.readlines()
file_loc.close()

file_name_loc = "mono_{}_{:.4E}_loc.txt".format(count, p1)

file_loc_results = open(os.path.join(results_path, file_name_loc), "w")

for line in lines_loc :
    file_loc_results.write(line)

file_loc_results.close()