######### Test case #################################################################################################################################

# cd IC/build
# ./compile.sh
# cd ..
# nbubble=8
# png=-25325
# ./build/cavitationonset_apecss -options run.apecss -amp $png -tend 60.0e-6 -nbb $nbubble -cldistrib 1 -clsize 10 -inttype 1
# python3 plot_temp.py
# for ((c=0; c<$nbubble; c++))
# do
#     rm -rf Bubble_$c
# done
 
# cd ..

######### No Interaction computations ###############################################################################################################

cd NI/build
./compile.sh
cd ..

######### Pressure time history & Radius evolution without interaction ############################
./build/cavitationonset_apecss -options run.apecss -amp -25325 -tend 60.0e-6 -nbb 2 -cldistrib 0 -clsize 15 -inttype 0
python3 gather_results.py

######### Cavitation inception pressure treshold for one single bubble ############################
png_list=(-17221 -17725.5 -18353.2 -18770.3)
for png in "${png_list[@]}"
do
    ./build/cavitationonset_apecss -options run.apecss -amp $png -tend 60.0e-6 -nbb 2 -cldistrib 0 -clsize 15 -inttype 0
    python3 gather_results.py
done

cd ..

echo ""
echo "No interaction test cases passed"
echo ""

######### Incompressible computations ###############################################################################################################

cd IC/build
./compile.sh
cd ..

######### Cavitation inception with interactions with varying distance between 2 bubbles ##########

######### Cavitation inception with interactions with varying pressure between 2 bubbles ##########

######### Cavitation inception with monodispersed simple distributions ############################ 

cd ..

echo ""
echo "Incompressible test cases passed"
echo ""

######### Quasi acoustic computations ###############################################################################################################

cd QA/build
./compile.sh
cd ..

######### Cavitation inception with interactions with varying distance between 2 bubbles ##########

######### Cavitation inception with interactions with varying pressure between 2 bubbles ##########

######### Cavitation inception with monodispersed simple distributions ############################

cd ..

echo ""
echo "Quasi acoustic test cases passed"
echo ""

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