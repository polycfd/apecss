######### Parameters ################################################################################################################################

ncores=17
nbubbles=2500

######### No interaction computations ###############################################################################################################

cd NI/build
./compile.sh
cd ..

######### Monodispersed system ####################################################################
./build/sphericalclustercavitationonset_apecss -options run.apecss -amp -25325 -tend 60.0e-6 -cldistrib 0
python3 gather_results.py

######### Polydispersed system ####################################################################
./build/sphericalclustercavitationonset_apecss -options run.apecss -amp -25325 -tend 60.0e-6 -cldistrib 1
python3 gather_results.py

cd ..

echo ""
echo "No interaction test cases passed"
echo ""

######### Incompressible computations ###############################################################################################################

cd IC/build
./compile.sh
cd ..

######### Monodispersed system ####################################################################
./build/sphericalclustercavitationonset_apecss -options run.apecss -amp -25325 -tend 60.0e-6 -cldistrib 0
python3 gather_results.py

######### Polydispersed system ####################################################################
./build/sphericalclustercavitationonset_apecss -options run.apecss -amp -25325 -tend 60.0e-6 -cldistrib 1
python3 gather_results.py

cd ..

echo ""
echo "Incompressible test cases passed"
echo ""

######### Incompressible computations ###############################################################################################################

######### Plot results ##############################################################################################################################

python3 plot_results.py

######### Cleaning ##################################################################################################################################

cd NI
rm -rf "onset_results.txt" "bubble_loc.txt"
for ((c=0; c<$nbubbles; c++))
do
    rm -rf Bubble_$c
done
cd ..

cd IC
rm -rf "onset_results.txt" "bubble_loc.txt"
for ((c=0; c<$nbubbles; c++))
do
    rm -rf Bubble_$c
done
cd ..

cd QA
rm -rf "bubble_loc.txt"
for ((c=0; c<$ncores; c++))
do
    rm -rf onset_results_$c.txt
done
for ((c=0; c<$nbubbles; c++))
do
    rm -rf Bubble_$c
done
cd ..

echo ""
echo "Cleaning completed"
echo ""