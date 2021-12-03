#!/bin/bash

for lig in <LIGS>  ; do
    ( cd <RECEPTOR>-${lig} && res=`./analyze.sh 20` && echo "$lig $res")
done
