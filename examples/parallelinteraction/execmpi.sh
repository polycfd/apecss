cd build
./compile.sh
cd ..
mpiexec -n 7 ./build/parallelinteraction_apecss -options run.apecss -freq 15.7e3 -amp -120e3 -tend 7.5e-4
