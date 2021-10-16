ASyncRE-OpenMM
==============

ASynchronous Replica Exchange for OpenMM (ASyncRE-OpenMM) is an extensible Python package enabling asynchronous parallel replica exchange molecular simulations with OpenMM. 

This version of ASyncRE is tailored to Alchemical Transfer Method alchemical calculations using the [ATMetaForce plugin](https://github.com/Gallicchio-Lab/openmm-atmmetaforce-plugin)

Credits
-------

This software is developed and maintained by the [Emilio Gallicchio's lab](http://www.compmolbiophysbc.org) with support from current and past grants from the National Science Foundation (ACI 1440665 and CHE 1750511).

Authors:

Emilio Gallicchio <egallicchio@brooklyn.cuny.edu>

Baofeng Zhang BZhang@brooklyn.cuny.edu

Rajat Pal <rajatfor2014@gmail.com>

The asynchronous replica exchange method was first implemented in the [AsyncRE](https://github.com/ComputationalBiophysicsCollaborative/AsyncRE) package for the IMPACT program.

Citation
--------

[Asynchronous Replica Exchange Software for Grid and Heterogeneous Computing](http://www.compmolbiophysbc.org/publications#asyncre_software_2015)

Installation & Usage
--------------------

It is recommended that the installation is performed in a personal python environment (`conda`, `miniconda`, or similar). ASyncRE requires the `configobj` and `numpy` python modules. 

```
git clone https://github.com/Gallicchio-Lab/async_re-openmm.git
cd async_re-openmm
python setup.py install
```

See [examples](examples/) for examples and tutorials.


