#!/bin/bash

for pair in <LIGPAIRS>  ; do
    ( cd <RECEPTOR>-${pair} && res=`./analyze.sh 20` && echo "$pair $res")
done
