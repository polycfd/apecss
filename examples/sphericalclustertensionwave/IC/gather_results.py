import os

# File designed to recover data for plotting results in cavitation onset case

file = open("tension_results.txt", "r")
lines = file.readlines()
file.close()

count = int(lines[0].split(" ")[0])
p1 = float(lines[0].split(" ")[5])
cl_distrib = int(lines[0].split(" ")[7])

if cl_distrib == 0 :
    # monodispersed spherical cluster
    file_name = "mono_{}_{:.4E}.txt".format(count, p1)

else :
    # polydispersed spherical cluster
    file_name = "poly_{}_{:.4E}.txt".format(count, p1)

working_path = os.getcwd()
results_path = os.path.join(working_path, "results")
file_results = open(os.path.join(results_path, file_name), "w")

# for line in lines :
#     if "nan" not in line :
#         file_results.write(line)

index = 0
while index < len(lines) and "nan" not in lines[index] :
    file_results.write(lines[index])
    index += 1

file_results.close()

file_loc = open("bubble_loc.txt", "r")
lines_loc = file_loc.readlines()
file_loc.close()

if cl_distrib == 0 :
    # monodispersed spherical cluster
    file_name_loc = "mono_{}_{:.4E}_loc.txt".format(count, p1)

else :
    # polydispersed spherical cluster
    file_name_loc = "poly_{}_{:.4E}_loc.txt".format(count, p1)

file_loc_results = open(os.path.join(results_path, file_name_loc), "w")

for line in lines_loc :
    file_loc_results.write(line)

file_loc_results.close()