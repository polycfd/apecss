######### Parameters ################################################################################################################################

ncore=15
nbubbles=2601

######### Incompressible computations ###############################################################################################################

cd IC/build
./compile.sh
cd ..

f_list=(1.631e06 2.934e06 3.262e06 4.893e06)
for f in "${f_list[@]}"
do
    ./build/bubblyscreen_apecss -options run.apecss -freq $f -amp 100.0 -tend 10.0e-6
    python3 gather_results.py
done

cd ..

echo ""
echo "Incompressible test cases passed"
echo ""

######### Quasi acoustic computations ###############################################################################################################

# cd QA/build
# ./compile.sh
# cd ..

# cd ..

# echo ""
# echo "Quasi acoustic test cases passed"
# echo ""

######### Plotting results ##########################################################################################################################

python3 plot_results.py

######### Cleaning ##################################################################################################################################

cd IC
for ((c=0; c<$nbubbles; c++))
do
    rm -rf Bubble_$c
done
rm -rf bubblyscreen_radii.txt bubblyscreen_extremum.txt
cd ..

# cd QA
# for ((c=0; c<$nbubbles; c++))
# do
#     rm -rf Bubble_$c
# done
# rm -rf bubblyscreen_radii.txt
# for ((c=0; c<$ncore; c++))
# do
#     rm -rf bubblyscreen_extremum_$c.txt
# done
# cd ..

echo ""
echo "Cleaning completed"
echo ""