ncore=15
nbubble=2601

cd build
./compile.sh
cd ..
mpiexec -n $ncore ./build/bubblyscreen_apecss -options run.apecss -freq 3.262e6 -amp 100.0 -tend 1.0e-6
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