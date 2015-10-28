# Clean metafiles for asynchronous Replica Exchang jobs
"""A module to clean metfiles for asynchronous RE jobs
See documentation in doc/ directory.

Contributors: 

Junchao Xia <junchao.xia@temple.edu>
Emilio Gallicchio <emilio.gallicchio@gmail.com>

"""


import os
import sys
import math
import getopt
import re,glob
from schrodinger import structure, structureutil 

# Parse arguments
items = sys.argv[1:]
if(len(items) < 1):
    print "Please specify basename"
    sys.exit(-1)

basename = items[0]
imp_version = items[1]

# get all cycle number 
cycles = []
if imp_version == "academic" :
    dms = "%s_*.dms" % basename
    outdms_lig = basename + "_lig.tar"
    outdms_rcpt =  basename + "_rcpt.tar"
    dms_files = glob.glob("%s_lig_*.dms" % basename)
    print 1,dms_files
    to_cycle = re.compile(basename + r"_lig_(\d+).dms")
    for f in dms_files:
        c = re.match(to_cycle, f).group(1)
        cycles.append(int(c))
        cycles.sort()
else:
    maeg = "%s_*.maegz" % basename
    outmae = basename + ".maegz"
    mae_files = glob.glob("%s_*.maegz" % basename)
    print 1,mae_files
    to_cycle = re.compile(basename + r"_(\d+).maegz")
    for f in mae_files:
        c = re.match(to_cycle, f).group(1)
        cycles.append(int(c))
        cycles.sort()

#leave the last two alone
cycles.pop()
cycles.pop()


# concatenate mae or tar files
if imp_version == "academic" :
    for c in cycles:
        #construct mae file name
        file_lig = "%s_lig_%d.dms" % (basename,c)
        file_rcpt = "%s_rcpt_%d.dms" % (basename,c)
        print "dms structure files %s %s" % (file_lig,file_rcpt)
        try:
            tarcom_lig = "tar -r --file=" + outdms_lig + " " + file_lig 
            tarcom_rcpt = "tar -r --file=" + outdms_rcpt + " " + file_rcpt
            os.system(tarcom_lig)
            os.system(tarcom_rcpt)
            os.remove(file_lig)
            os.remove(file_rcpt)
        except:
            print "Warning: Cannot open dms structure file %s %s" % (file_lig, file_rcpt)
else:
    for c in cycles:
        #construct mae file name
        file = "%s_%d.maegz" % (basename,c)
        print file
        try:
            ct1 = structure.StructureReader(file).next()
            ct1.append(outmae)
            os.remove(file)
        except:
            print "Warning: Cannot open Maestro structure file %s" % file

#concatenate out files
outfile = "%s.out" % basename
fout = open(outfile,"a")
for r in cycles:
        #construct file name
        file = "%s_%d.out" % (basename,r)
        print file
        try:
            f = open(file,"r")
            fout.write( f.read() )
            f.close()
            os.remove(file)
        except:
            print "Warning: Cannot open output file %s" % file
fout.close()

#delete .rst files
for r in cycles:
        #construct file name
        file = "%s_%d.rst" % (basename,r)
        print file
        try:
            os.remove(file)
        except:
            print "Warning: Cannot open output file %s" % file


