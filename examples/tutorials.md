ASyncRE-OpenMM Tutorials 
========================

VERSION 0.2.0

Temperature Replica Exchange with `SSH` Transport
-------------------------------------------------

```
cd temperature-RE-ssh-transport
```

Set the `nodefile` and the `runopenmm` launch script as described below. Then, assuming ASyncRE is installed under `$HOME/src/asyncre-openmm`, do:

```
./runopenmm $HOME/src/async_re-openmm/tempt_async_re.py t4l_asyncre.cntl
```

Do:
```
catdcdpath=`/bin/ls -1 -d /usr/local/lib/vmd/plugins/LINUXAMD64/bin/catdcd*`
export PATH=$PATH:${catdcdpath}

for i in `seq 0 4` ; do 
 cd r$i
 python cleanup.py t4l
 cd ..
done
```
to assemble the output and trajectory files. Temperatures and potential energies are written to the `t4l.out` in each replica directory `r0`, `r1`, etc. Similarly, `t4l.dcd` holds the replica trajectory.  


Temperature Replica Exchange with the `LOCAL_OPENMM` Transport
--------------------------------------------------------------

```
cd temperature-RE-local-transport
```

Set the `nodefile` and the `runopenmm` launch script as described below. Then, assuming ASyncRE is installed under `$HOME/src/asyncre-openmm`, do:

```
./runopenmm $HOME/src/async_re-openmm/tempt_async_re.py t4l_asyncre.cntl
```

Temperatures and potential energies are written to the `t4l.out` in each replica directory `r0`, `r1`, etc. Similarly, `t4l.dcd` holds the replica trajectory.  


The `nodefile`
-------------

The SSH ASyncRE system dispatches jobs to GPU computing devices on the machines listed in the `nodefile`. The format of each line is as follows:

```
<machine name> , <OpenCL platform id>:<OpenCL GPU device id> , <num GPUs> , <platform name> , <remote username>, <remote temp directory>
```

For example, when using 2 GPUs on 2 servers named `compute-server-1` and `compute-server-2`, we will set the nodefile as:

```
compute-server-1,0:0,1,OpenCL,,/tmp
compute-server-2,1:0,1,OpenCL,,/tmp
```

This assumes that on the first machine the OpenCL platform for the GPUs is the first platform (0), and for the second server it is the second OpenCL platform (1). On both platforms we are using the first GPU (0). If you are running this tutorial on your desktop with only one GPU, a `nodefile` such as the following will probably work: 

```
locahost,0:0,1,OpenCL,,/tmp
```

The `<num GPUs>` setting is for MD threads running on multiple GPUs, that are not currently supported. The only allowed value is 1. In this tutorial, the RE simulation of each complex employs 5 replicas. ASyncRE assumes that there are more replicas than computing devices. In this case, it is not recommended to use more than 2 GPUs.

With the 'local' transport, see below, jobs are dispatched only to the local machine. In this case, the machine identifier is ignored.  

The `runopenmm` launch script
-----------------------------

Edit the `runopenmm` launch script to reflect your environment. For example, if OpenMM is installed in `$HOME/local/openmm-7.3.1` and the corresponding python bindings are installed in the Conda environment under `$HOME/miniconda2`, the `runopenmm` should look like:

```
#!/bin/bash
openmm_dir=$HOME/local/openmm-7.3.1
pythondir=$HOME/miniconda2
export LD_LIBRARY_PATH=${openmm_dir}/lib:${openmm_dir}/lib/plugins:$LD_LIBRARY_PATH
${pythondir}/bin/python "$@"
```




