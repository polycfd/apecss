######### Parameters ################################################################################################################################

nbubbles=2500

######### Incompressible computations ###############################################################################################################

cd IC/build
./compile.sh
cd ..

######### Monodispersed system ####################################################################
./build/sphericalclustercavitationonset_apecss -options run.apecss -amp -25325 -tend 25.0e-6 -cldistrib 0
python3 gather_results.py

######### Polydispersed system ####################################################################
./build/sphericalclustercavitationonset_apecss -options run.apecss -amp -25325 -tend 25.0e-6 -cldistrib 1
python3 gather_results.py

cd ..

echo ""
echo "Incompressible test cases passed"
echo ""

######### Incompressible computations ###############################################################################################################

######### Plot results ##############################################################################################################################

######### Cleaning ##################################################################################################################################

cd IC
# rm -rf "onset_results.txt" "bubble_loc.txt"
for ((c=0; c<$nbubbles; c++))
do
    rm -rf Bubble_$c
done
cd ..

echo ""
echo "Cleaning completed"
echo ""