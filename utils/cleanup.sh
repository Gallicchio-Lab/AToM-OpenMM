#!/bin/bash
# Clean metafiles for all replicas in asynchronous Replica Exchang jobs
# Contributors: 
#    Junchao Xia <junchao.xia@temple.edu>

SCHRODINGER=$1
scripts=$2
jobname=$3
imp_version=$4
rb=$5 
re=$6
export MALLOC_CHECK_=0
export SCHRODINGER
for (( ir=$rb; ir<=$re; ir++ ))
do
   cd r$ir
   $SCHRODINGER/run $scripts/cleanup.py $jobname $imp_version && rm -f ${jobname}_*.{inp,err,log,trj,idx}
   echo "Finished the cleanup in r$ir"
   cd ../
done
