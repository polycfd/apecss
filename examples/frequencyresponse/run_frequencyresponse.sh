set -e

######################## No interactions (NI) #####################################################

cd NI/build
./compile.sh
cd ..

for ((f=500000;f<=10000000;f+=3000))
do
    ./build/frequencyresponse_apecss -options run.apecss -nb 4 -freq $f -amp 120e3
done
python3 gather_results.py
rm -rf max_radius.txt

for ((f=500000;f<=10000000;f+=3000))
do
    ./build/frequencyresponse_apecss -options run.apecss -nb 4 -freq $f -amp 120e3 -ode 1
done
python3 gather_results.py
rm -rf max_radius.txt

cd ..

######################## Incompressible interactions (IC) #########################################

cd IC/build
./compile.sh
cd ..

for ((f=500000;f<=10000000;f+=3000))
do
    ./build/frequencyresponse_apecss -options run.apecss -nb 3 -freq $f -amp 120e3 -dt_inter 1e-9
done
python3 gather_results.py
rm -rf max_radius.txt

for ((f=500000;f<=10000000;f+=3000))
do
    ./build/frequencyresponse_apecss -options run.apecss -nb 3 -freq $f -amp 120e3 -ode 1 -dt_inter 1e-9
done
python3 gather_results.py
rm -rf max_radius.txt

for ((f=500000;f<=10000000;f+=3000))
do
    ./build/frequencyresponse_apecss -options run.apecss -nb 4 -freq $f -amp 120e3 -dt_inter 1e-9
done
python3 gather_results.py
rm -rf max_radius.txt

for ((f=500000;f<=10000000;f+=3000))
do
    ./build/frequencyresponse_apecss -options run.apecss -nb 4 -freq $f -amp 120e3 -ode 1 -dt_inter 1e-9
done
python3 gather_results.py
rm -rf max_radius.txt

cd ..

######################## Cleaning #################################################################

cd NI
for ((c=0;c<=4;c++))
do
    rm -rf Bubble_$c
done
cd ..

cd IC
for ((c=0;c<=4;c++))
do
    rm -rf Bubble_$c
done
cd ..