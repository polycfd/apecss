cd IC/build
./compile.sh
cd ..
nbubble=2
png=-25325
./build/cavitationonset_apecss -options run.apecss -amp $png -tend 60.0e-6 -nbb $nbubble -cldistrib 0 -clsize 20 -inttype 1
# for ((c=0; c<$nbubble; c++))
# do
#     rm -r Bubble_$c
# done
 
cd ..