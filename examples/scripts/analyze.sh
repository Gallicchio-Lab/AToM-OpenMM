#!/bin/bash

#usage:
#
# analyze.sh [ discard_low ] [ discard_high ]
#


many=1000000

if [ $# -eq 0 ] 
then
    discard_samples_low=0
    discard_samples_high=$many
else
    if [ $# -eq 1 ]
    then
	discard_samples_low=$1
	discard_samples_high=$many
    else
	discard_samples_low=$1
        discard_samples_high=$2
    fi
fi


datafilename=repl.cycle.state.temp.lambda1.lambda2.alpha.u0.w0.epot.epert.dat
rm -f ${datafilename} || exit 1


for repl_dir in r? r?? ; do
    repl=${repl_dir#r}
    awk -v low=$discard_samples_low -v high=$discard_samples_high -v repl=$repl 'FNR >= low && FNR <= high {print repl, FNR, $0 }' $repl_dir/*.out || exit 1
done >> $datafilename

maxsamples=$discard_samples_high
for repl_dir in r? r?? ; do
    repl=${repl_dir#r}
    nsamples=`wc $repl_dir/*.out | awk 'FNR ==1 {print $1}'` || exit 1
    if [ $maxsamples -gt $nsamples ] ; then
	maxsamples=$nsamples
    fi
done


rm -f result.log Rplots.pdf p-l*.dat || exit 1
  
R CMD BATCH uwham_analysis.R || exit 1
res=`grep 'DGb =' result.log` || exit 1

echo ${res} " range:" $discard_samples_low $maxsamples || exit 1
