# Evolving structures in complex systems

This repository contains the code to reproduce the Figures and numerical results
from the paper: _Evolving structures in complex systems_ by Hugo Cisneros, Josef
Sivic and Tomas Mikolov.


## Build the C library

The C library provided implements a general cellular automaton simulator and all
the steps for computing the metrics we discuss in the paper. 

The library can be built with
```make all```
which will create a binary in `bin/automaton`.

### Generate automata

Automata rules are encoded in mapping files with the following format: 

Possible transitions are enumerated for a 3x3 neighborhood and `N` states in the
way described below, starting from the top-left cell and incrementing in a
row-first manner a __base-N counter__ with __8 cells__.

![Incrementing](./figures/increment.png)

Mapping files just enumerate the resulting state of the middle cell for the
corresponding neighborhood state.

Mapping files with the 3-states rules reported in the paper can be obtained at
the following
[link](https://drive.google.com/uc?id=1fymRRN-Yeig560CkXrLTfpl879YLP_UF&export=download).
The unzipped `maps` directory contains the subdirectories `train` and `test`
that correspond to the training and testing sets used in the paper.

Once the `maps` directory is placed at the root of this directory, running 

```./scripts/generate_results_from_maps.sh```

will simulate all automata from the maps and compute the various metrics we
discuss in the paper (this might take a while to complete). The script calls the
executable `bin/automaton` with some options for each `.map` file.


## Python scripts

