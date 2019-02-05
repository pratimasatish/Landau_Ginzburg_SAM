#!/bin/bash --login

templist=$(seq 340.0 5 390.0)

for temp in ${templist}
do
    echo ${temp} $(./test.o -t ${temp} | grep Average | awk '{print $7}')
done
