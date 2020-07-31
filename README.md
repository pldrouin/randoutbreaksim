# branchsimc
A C implementation of the branchsim simulation. It currently focusses only on the simulation functionalities, so it is not a replacement for the original package. It is a faster and more memory efficient implementation of the simulation process. Speed gains are due to the native implementation, minimised memory operations and on-the-fly computation of parameters. Reduction of memory usage is due to the different simulation algorithm that decouples the simulation process from data processing. Memory usage for the simulation process scales with the number of layers such that it is kept to a minimum. Data processing functions can be defined by the user. Speed gains of approximately three orders of magnitude can be expected (depends on data processing and simulation parameters).

Depends only on libgsl. See src/bin/main.c for usage example.

No API documentation yet.
