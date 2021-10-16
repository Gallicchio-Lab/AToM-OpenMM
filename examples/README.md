ASyncRE-OpenMM Examples
-----------------------

### Contents


 - [ABFE](ABFE/): absolute binding free energy between the TEMOA host and the G1 guest from the [SAMPL8 GDCC](https://github.com/samplchallenges/SAMPL8/tree/master/host_guest/GDCC) challenge.
 - RBFE: (upcoming) relative binding free energy between the G1 and G4 guests to the TEMOA host from the [SAMPL8 GDCC](https://github.com/samplchallenges/SAMPL8/tree/master/host_guest/GDCC) challenge.

### Software requirements

Tu run the examples we assume that the [OpenMM](http://openmm.org) libraries (version 7.5.0 or newer) have been installed in a location where python can find them. Presumably you did so under a `conda` environment. If so activate the environment to run the examples.

We also assume that the [ATMMetaForce OpenMM plugin](https://github.com/Gallicchio-Lab/openmm-atmmetaforce-plugin) is installed and available within the same OpenMM environment and that python can find the corresponding python bindings.

A sample `runopenmm` script is provided in the  [scripts/]( https://github.com/Gallicchio-Lab/async_re-openmm/tree/master/examples/scripts) to get you started if needed. It requires you to define the environment variable `OPENMM_DIR` pointing to the OpenMM installation on your system. It also assumes that the python OpenMM bindings have been stored under the same folder and that you are using python 3.7.

Free energy analysis requires R with th UWHAM R module. To install the UWHAM module run `install.packages("UWHAM")` under R.

