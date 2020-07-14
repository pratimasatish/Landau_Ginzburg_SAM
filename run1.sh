#!/bin/bash --login

templist=$(seq 320.0 5 420.0)

for temp in ${templist}
do
    echo ${temp} $(./do_mc.o -t ${temp} -bc_on 0 | grep Average | awk '{print $7}')
done
