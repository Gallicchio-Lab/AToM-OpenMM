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


rm -f result.log Rplots.pdf p-*.dat || exit 1

#get the jobname from the name of the current folder
jobname=$(basename $PWD)

if [ -f r0/${jobname}.out ] ; then
    
    R CMD BATCH -${jobname} -${discard_samples_low} -${discard_samples_high}  uwham_analysis.R || exit 1
    res=`grep -e '^DDGb =' uwham_analysis.Rout` || exit 1

    echo ${res}
else
    echo "One or more replica directories are missing" 
fi
