cd build
./compile.sh
cd ..
mpiexec -n 15 ./build/sphericalclusterinteractions_apecss -options run.apecss -freq 50.0e03 -amp -2.5e5 -tend 5.0e-07

rm -rf bubble_loc.txt
for ((c=0; c<2500; c++))
do
    rm -rf Bubble_$c
done