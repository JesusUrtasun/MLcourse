#! /bin/bash

echo -n "\nRun ggHiggs.cc"

# For loop
for ((s = 1000; s <= 15000; s+=50))
do
./ggHiggs_partonic_Nspace -m 125 -E $s
done
printf "\n"