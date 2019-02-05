#!/bin/bash --login

templist=$(seq 250.0 10 400.0)

for temp in ${templist}
do
    echo ${temp} $(./do_mc.o -t ${temp} | grep Average | awk '{print $7}')
done
