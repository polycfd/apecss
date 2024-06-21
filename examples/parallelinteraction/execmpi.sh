cd build
./compile.sh
cd ..
/usr/lib64/mpich/bin/mpiexec -np 4 ./build/parallelinteraction_apecss -options run.apecss -freq 15.7e3 -amp -120e3 -tend 7.5e-4
