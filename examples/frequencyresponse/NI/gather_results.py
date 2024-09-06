import os

######### File designed to retrieve and gather in the "results" folder the data obtained with frequency response computations
 
working_path = os.getcwd()

results_file = open(os.path.join(working_path, "max_radius.txt"), "r")
lines = results_file.readlines()
results_file.close()

firstline = lines[0].split(" ")
nBubbles = int(firstline[1])
pa = float(firstline[5])
dist = float(firstline[7])
ode = int(firstline[9])

gather_path = os.path.join(working_path, "results")

file_name = "{}_p{:.5E}_d{:.2E}_{}.txt".format(nBubbles, pa, dist, ode)

data_file = open(os.path.join(gather_path, file_name), "w")
data_file.write("nb {} p(Pa) {:.5E} dist {:.2E} ODE {}\n".format(nBubbles, pa, dist, ode))
data_file.write("#f(Hz) R0(m);R(m)\n")

for line in lines :
    if line != "" :
        data = line.split(" ")
        f = float(data[3])
        data_file.write("{:.5E}".format(f))

        r_list = data[11:]
        for r in r_list :
            data_file.write(" {}".format(r))

data_file.close()