#!/bin/bash
set -e
cd acousticemitter/build
./compile.sh
cd ..
./build/acousticemitter_apecss -options planar.apecss -tend 0.5 -fa 1e3 -dpa 1e6
python3 plot_result_planar.py
rm EmissionsTime_*
./build/acousticemitter_apecss -options spherical.apecss -tend 4.0e-3 -fa 1500 -dpa 1
python3 plot_result_spherical.py
rm EmissionsTime_*
cd ../binaryinteraction/build
./compile.sh
cd ..
./build/binaryinteraction_apecss -options run.apecss -freq 15.7e3 -amp -120e3 -tend 7.5e-4
python3 plot_result.py
rm -r Bubble_0 Bubble_1
./cavitationonset/run_cavitationonset.sh
cd ../gastemperature/build
./compile.sh
cd ..
./build/gastemperature_apecss -options run.apecss -freq 20e3 -amp 70e3 -tend 0.001
python3 plot_result.py
rm KellerMiksis_R1.300e-04_fa2.000e+04_pa7.000e+04.txt
cd ../laserinducedcavitation/build
./compile.sh
cd ..
./build/lic_apecss -options run.apecss -tend 12e-6
python3 plot_result.py
rm EmissionsTime_6.49* Gilmore_R1.330e-06.txt
cd ../rayleighcollapse/build/
./compile.sh
cd ..
./build/rayleigh_apecss -options simple.apecss -tend 1
python3 plot_result_simple.py
rm RP_R1.000e+00.txt
./build/rayleigh_apecss -options cylindrical.apecss -tend 5e-4
python3 plot_result_cylindrical.py
rm Gilmore_R2.110e-03.txt
./build/rayleigh_apecss -options emissions.apecss -tend 0.2
python3 plot_result_emissions.py
rm EmissionsSpace_* Gilmore_R1.000e+00.txt
./build/rayleigh_apecss -options shock.apecss -tend 0.11
python3 plot_result_shock.py
rm EmissionsSpace_* Gilmore_R1.000e+00.txt
cd ../ultrasound/build
./compile.sh
cd ..
./build/ultrasound_apecss -options lipidcoated_simple.apecss -freq 2.9e6 -amp 130e3 -tend 2e-6
python3 plot_result_lipidcoated_simple.py
rm RPAR_R9.750e-07_fa2.900e+06_pa1.300e+05.txt
./build/ultrasound_apecss -options lipidcoated_emissions.apecss -tend 11e-6 -freq 1e6 -amp 600e3
python3 plot_result_lipidcoated_emissions.py
rm EmissionsNode_* Gilmore_R1.000e-06_fa1.000e+06_pa6.000e+05.txt
./build/ultrasound_apecss -options sonolum_emissions.apecss -tend 40e-6 -freq 23.5e3 -amp 145e3
python3 plot_result_sonolum.py
rm EmissionsNode_* Gilmore_R5.000e-06_fa2.350e+04_pa1.450e+05.txt
./build/ultrasound_apecss -options kelvinvoigt.apecss -tend 6e-6 -freq 1e6 -amp 3e6
python3 plot_result_kelvinvoigt.py
rm KellerMiksis_R1.000e-06_fa1.000e+06_pa3.000e+06.txt
./build/ultrasound_apecss -options zener.apecss -tend 5e-6 -freq 1e6 -amp 1e6
python3 plot_result_zener.py
rm Gilmore_R1.000e-06_fa1.000e+06_pa1.000e+06.txt
./build/ultrasound_apecss -options oldroydb.apecss -tend 3e-6 -freq 3e6 -amp 400e3
python3 plot_result_oldroydb.py
rm RP_R1.000e-06_fa3.000e+06_pa4.000e+05.txt
cd ../
./sphericalclustertensionwave/run_sphericalclustertensionwave.sh
