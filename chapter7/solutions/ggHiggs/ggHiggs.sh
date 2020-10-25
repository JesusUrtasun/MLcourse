#! /bin/bash

echo -n "\nRun ggHiggs.cc"

# For loop
for ((s = 1000; s <= 15000; s+=500))
do
./ggHiggs -m 125 -E $s
done
printf "\n"