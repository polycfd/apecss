cd build
./compile.sh
cd ..
mpiexec -n 15 ./build/bubblyscreen_apecss -options run.apecss -freq 3.262e6 -amp 100.0 -tend 1.0e-6

for ((c=0; c<2601; c++))
do
    rm -rf Bubble_$c
done