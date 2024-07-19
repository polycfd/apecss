cd build
./compile.sh
cd ..
mpiexec -n 15 ./build/sphericalclustercavitationonset_apecss -options run.apecss -amp -25325 -tend 60.0e-06

# rm -rf bubble_loc.txt
# for ((c=0; c<2500; c++))
# do
#     rm -rf Bubble_$c
# done