nbubbles=250
ncores=10

cd build
./compile.sh
cd ..
mpiexec -n $ncores ./build/sphericalclustertensionwave_apecss -options run.apecss -amp -3.0e4 -tend 25.0e-06 -cldistrib 0
python3 gather_results.py

rm -rf bubble_loc.txt
for ((c=0; c<$ncores; c++))
do
    rm -rf tension_results_$c.txt
done
for ((c=0; c<$nbubbles; c++))
do
    rm -rf Bubble_$c
done