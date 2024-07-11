import os

# File designed to recover data for plotting results in cavitation onset case

file = open("Ida2009_results.txt", "r")
lines = file.readlines()
file.close()

count = int(lines[0].split(" ")[0])
png = float(lines[0].split(" ")[5])
cl_size = float(lines[0].split(" ")[7])
cl_distrib = int(lines[0].split(" ")[9])

if cl_distrib == 0 :
    # Two bubbles of different size
    file_name = "{}_{:.4E}_{:.2f}.txt".format(count, png, cl_size)

else :
    # Monodispersed system
    file_name = "{}_{:.4E}.txt".format(count, png)

working_path = os.getcwd()
results_path = os.path.join(working_path, "results")
file_results = open(os.path.join(results_path, file_name), "w")

for line in lines :
    file_results.write(line)

file_results.close()