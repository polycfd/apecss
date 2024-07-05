cd build
./compile.sh
cd ..
mpiexec -n 5 ./build/parallelinteraction_apecss -options run.apecss -freq 15.7e3 -amp -120e3 -tend 7.5e-4

#mpiexec -n 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valout.txt ./build/parallelinteraction_apecss -options run.apecss -freq 15.7e3 -amp -120e3 -tend 7.5e-8
#valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valout.txt ./build/parallelinteraction_apecss -options run.apecss -freq 15.7e3 -amp -120e3 -tend 7.5e-8
#valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valout.txt mpiexec -n 2 ./build/parallelinteraction_apecss -options run.apecss -freq 15.7e3 -amp -120e3 -tend 7.5e-8

for ((c=0; c<25; c++))
do
    rm -rf Bubble_$c
done