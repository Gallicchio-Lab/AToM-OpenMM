import os
import sys
import math
import getopt
import re,glob

"""
example:
catdcdpath=`/bin/ls -1 -d /usr/local/lib/vmd/plugins/LINUXAMD64/bin/catdcd*`
export PATH=$PATH:${catdcdpath}

for i in `seq 0 4` ; do 
 cd r$i
 python cleanup.py basename
 cd ..
done
"""


# Parse arguments
items = sys.argv[1:]
if(len(items) < 1):
    print "Please specify basename"
    sys.exit(-1)

basename = items[0]

out_files = glob.glob("%s_*.out" % basename)
to_cycle = re.compile(basename + r"_(\d+).out")
cycles = []
for f in out_files:
    try:
        c = re.match(to_cycle, f).group(1)
        cycles.append(int(c))
    except:
        pass

if len(cycles) < 3:
    print("Nothing to clean")
    sys.exit(0)

cycles.sort()
#leaves the last two cycles alone
cycles.pop()
cycles.pop()

# concatenate dcd files
dcdfiles = ""
outdcd = "%s.dcd" % basename
for c in cycles:
    dcdfiles += " %s_%d.dcd" % (basename, c)
if os.path.exists(outdcd):
    command = "catdcd -o .tmp.dcd %s %s" % (outdcd, dcdfiles)
else:
    command = "catdcd -o .tmp.dcd %s" % (dcdfiles)

try:
    print(command)
    os.system(command)
    command = "mv .tmp.dcd %s && rm %s" % (outdcd, dcdfiles)
    print(command)
    os.system(command)
except:
    sys.exit("Error archiving dcd files.")

#concatenate out files
outfiles = ""
outout =  "%s.out" % basename
for c in cycles:
    outfiles += " %s_%d.out" % (basename, c)
if os.path.exists(outout):
    command = "cat %s %s > .tmp.out " % (outout,outfiles)
else:
    command = "cat %s > .tmp.out " % (outfiles)
    
try:
    print(command)
    os.system(command)
    command = "mv .tmp.out %s && rm %s "  % (outout, outfiles)
    print(command)
    os.system(command)
except:
    sys.exit("Error archiving out files.")

#delete .err, .py. .log .dms .pdb files
files = ""
for c in cycles:
    files += " %s_%d.err" % (basename, c)
command = "rm %s" % files
print(command)
os.system(command)

files = ""
for c in cycles:
    files += " %s_%d.py" % (basename, c)
command = "rm %s" % files
print(command)
os.system(command)

files = ""
for c in cycles:
    files += " %s_%d.log" % (basename, c)
command = "rm %s" % files
print(command)
os.system(command)

files = ""
for c in cycles:
    files += " %s_%d.dms" % (basename, c)
command = "rm %s" % files
print(command)
os.system(command)

files = ""
for c in cycles:
    files += " %s_%d.pdb" % (basename, c)
command = "rm %s" % files
print(command)
os.system(command)
