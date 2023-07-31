# randoutbreaksim
A C implementation of a disease propagation branching process. It consists in a
fast and memory efficient implementation of the simulation process, with extensions
of the model beyond a branching process. The performance improvement is due to the
native implementation, avoidance of any indexing/hashing/searching (only direct access
is used throughout the algorithm), minimised memory operations and on-the-fly computation
of parameters. Small memory usage is achieved through the simulation algorithm that
decouples the simulation process from data processing. Memory usage for the simulation
process scales with the number of layers such that it is kept to a minimum. Data processing
functions can be defined by the user.

Also now included is a non branching model more suitable for the simulation of significant
level of propagation within a finite population. It shares most of its code with the
branching model.

Depends on [librngstream](https://github.com/pldrouin/librngstream) and on [libgsl](https://www.gnu.org/software/gsl/).
The branching process is detailed in the included
[report](doc/DRDC-RDDC-2023-R025_DOCUMENT-PDFA.pdf).
