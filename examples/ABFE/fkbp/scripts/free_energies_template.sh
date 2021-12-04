#!/bin/bash

rcpt=<RECEPTOR>
for lig in <LIGS>  ; do
    ( cd ${rcpt}-${lig} && res=`./analyze.sh 20` && echo "${rcpt}-${lig} ${res}")
done
