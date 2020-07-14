#!/bin/bash --login

awk 'NR != 1' $1 > dump1.xyz
awk 'NR % 443 != 0' dump1.xyz > dump2.xyz
awk 'NR != 1' dump2.xyz > dump3.xyz
awk 'NR % 442 != 0' dump3.xyz > dump4.xyz
awk '{print $2,$3,$4}' dump4.xyz > dump5.xyz
rm dump1.xyz dump2.xyz dump3.xyz dump4.xyz
mv dump5.xyz $2 

