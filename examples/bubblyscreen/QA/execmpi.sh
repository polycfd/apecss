ncore=15
nbubble=2601

cd build
./compile.sh
cd ..
mpiexec -n $ncore ./build/bubblyscreen_apecss -options run.apecss -freq 1.631e6 -amp 100.0 -tend 35.0e-6 -coeff 1.0
python3 gather_results.py

for ((c=0; c<$nbubble; c++))
do
    rm -rf Bubble_$c
done

rm -rf bubblyscreen_radii.txt
for ((c=0; c<$ncore; c++))
do
    rm -rf bubblyscreen_extremum_$c.txt
done