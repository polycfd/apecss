######### No Interaction computations ###############################################################################################################

cd NI/build
./compile.sh
cd ..

######### Pressure time history & Radius evolution without interaction ############################
./build/cavitationonset_apecss -options run.apecss -amp -25325 -tend 60.0e-6 -clsize 15
python3 gather_results.py

cd ..

echo ""
echo "No interaction test cases passed"
echo ""

######### Incompressible computations ###############################################################################################################

cd IC/build
./compile.sh
cd ..

######### Cavitation inception with interactions with varying distance between 2 bubbles ##########
dist_list=(10 11 12 12.05 12.1 12.5 15 20)
for d in "${dist_list[@]}"
do
    ./build/cavitationonset_apecss -options run.apecss -amp -25325 -tend 60.0e-6 -clsize $d
    python3 gather_results.py
done

######### Cavitation inception with interactions with varying pressure between 2 bubbles ##########
png_list=(-25325 -27351 -27654.9 -27958.8 -29377)
for png in "${png_list[@]}"
do
    ./build/cavitationonset_apecss -options run.apecss -amp $png -tend 60.0e-6 -clsize 10
    python3 gather_results.py
done

cd ..

echo ""
echo "Incompressible test cases passed"
echo ""

######### Quasi-acoustic computations ###############################################################################################################

cd QA/build
./compile.sh
cd ..

######### Cavitation inception with interactions with varying distance between 2 bubbles ##########
dist_list=(10 11.9 12 15 20)
for d in "${dist_list[@]}"
do
    ./build/cavitationonset_apecss -options run.apecss -amp -25325 -tend 60.0e-6 -clsize $d -dt_inter 1.0e-09
    python3 gather_results.py
done

######### Cavitation inception with interactions with varying pressure between 2 bubbles ##########
png_list=(-25325 -27351 -27654.9 -27958.8 -29377)
for png in "${png_list[@]}"
do
    ./build/cavitationonset_apecss -options run.apecss -amp $png -tend 60.0e-6 -clsize 10 -dt_inter 1.0e-09
    python3 gather_results.py
done

cd ..

echo ""
echo "Quasi-acoustic test cases passed"
echo ""

######### Plotting results ##########################################################################################################################

python3 plot_results.py

######### Cleaning ##################################################################################################################################

cd NI
rm -rf "Ida2009_results.txt"
for ((c=0; c<2; c++))
do
    rm -rf Bubble_$c
done
cd ..

cd IC
rm -rf "Ida2009_results.txt"
for ((c=0; c<8; c++))
do
    rm -rf Bubble_$c
done
cd ..

cd QA
rm -rf "Ida2009_results.txt"
for ((c=0; c<8; c++))
do
    rm -rf Bubble_$c
done
cd ..

echo ""
echo "Cleaning completed"
echo ""