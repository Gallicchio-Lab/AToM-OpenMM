#!/bin/bash

for pair in <LEG1PAIRS>  ; do
    mols_t=$(echo $pair | tr '-' ' ' )
    read -ra mols <<< $mols_t
    leg1=${mols[0]}-${mols[1]}
    leg2=${mols[1]}-${mols[0]}
    cd <RECEPTOR>-$leg1
    #DGb = 18.33984 +- 0.3736325 DE = -0.4071684 +- 0.3738681  range: 150 310
    res1=`./analyze.sh 20`
    cd ..
    cd <RECEPTOR>-$leg2
    res2=`./analyze.sh 20`
    cd ..
    res=`echo $res1 $res2 | awk '{printf("DDGb= %6.2f +- %6.2f DDE =  %6.2f +- %6.2f \n", $3-$16, sqrt($5*$5+$18*$18), $8+$21, sqrt($10*$10+$23*$23) )}'`
    echo "$leg1 $res1"
    echo "$leg2 $res2"
    echo "$leg1 $res"
done
