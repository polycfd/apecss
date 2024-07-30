nbubbles=2500
ncores=20

cd build
./compile.sh
cd ..
mpiexec -n $ncores ./build/sphericalclustercavitationonset_apecss -options run.apecss -amp -25325 -tend 60.0e-06 -cldistrib 1
python3 gather_results.py

# rm -rf bubble_loc.txt
# for ((c=0; c<$ncores; c++))
# do
#     rm -rf onset_results_$c.txt
# done
# for ((c=0; c<$nbubbles; c++))
# do
#     rm -rf Bubble_$c
# done