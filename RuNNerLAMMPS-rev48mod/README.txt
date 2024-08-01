####################################################################################################
#
# README.txt for RuNNerLAMMPS
#
# author:   Andreas Singraber
# date:     2013-02-28
# email:    andreas.singraber@univie.ac.at
#
# Copyright (c) 2013 Andreas Singraber
#
####################################################################################################

RuNNerLAMMPS is a C/C++ implementation of mode 3 (prediction) of the RuNNer neural network (NN) code
in the MD package LAMMPS. It combines the ability of a previously fitted NN to predict energies and
forces with the usability and features of a large MD simulation package. It also allows to run
massively parallelized simulations with the NN potential using MPI.

####################################################################################################
# BUILDING LAMMPS WITH RuNNer PAIR STYLE
####################################################################################################

First, download the LAMMPS tarball from

http://lammps.sandia.gov/index.html

and extract its contents. RuNNerLAMMPS was tested with LAMMPS version 22Feb13 which can be obtained
following this direct link:

    wget http://lammps.sandia.gov/tars/lammps-22Feb13.tar.gz
    tar -xzvf lammps-22Feb13.tar.gz

Then copy "pair_runner.cpp" and "pair_runner.h" from the "source" directory to the LAMMPS "src"
directory. 

    cp source/pair_runner.* lammps-22Feb13/src/

Next, change to the LAMMPS Makefile directory:

    cd lammps-22Feb13/src/MAKE/

There you find a lot of preconfigured Makefiles which can be adapted to fit your compiler and
library settings. For more information on how to build LAMMPS visit the website: 

http://lammps.sandia.gov/doc/Section_start.html#start_2

(If you don't have a fftw library you can turn it off by setting "FFT_INC = -DFFT_NONE" and
leaving "FFT_LIB = " blank.)

After editing your Makefile change back to the "lammps-22Feb13/src/" directory and type 

    make <target>

where <target> is the suffix of the Makefile you want to use, e.g. make openmpi.  If no errors
occurred during compilation you should find the executable with the name lmp_<target> (e.g.
lmp_openmpi).

####################################################################################################
# TESTING LAMMPS
####################################################################################################

To test if the executable works correctly copy it to one of the directories in the "examples"
directory and execute:

    mpirun -np 4 ./lmp_<target> < sim.lmp

You should see LAMMPS running an MD simulations, data is produced in the "traj" directory. 

####################################################################################################
# INFORMATION ABOUT RuNNer PAIR STYLE 
####################################################################################################

The file "sim.lmp" in one of the "examples" directories contains all commands executed by LAMMPS
during runtime. See the LAMMPS website for a detailed description of all commands: 

http://lammps.sandia.gov/doc/Section_commands.html#cmd_5

For some command templates see the file "templates.lmp".

The source files "pair_runner.*" introduce the neural network potential in the context of a pair
style named "runner". LAMMPS comes with a large collection of these pair styles, e.g. for a LJ or
Tersoff potential, for more information see:

http://lammps.sandia.gov/doc/pair_style.html

The setup of a pair style is done by issuing two commands: pair_style and pair_coeff. For the
RuNNer pair style the correct syntax is:

    pair_style runner keyword value ...
    pair_coeff * * <maxrcut>

    .) One or more keyword/value pairs may be listed. 
    .) keyword = dir OR showew OR resetew OR maxew

    dir = directory which contains the following NN files (default: RuNNer/):
        input.nn
        scaling.data
        weights.???.data
    showew = print extrapolation warnings (yes OR no, default: yes). 
    showewsum  = print extrapolation warning summary every this many timesteps (0=off, default: 0). 
    resetew = reset extrapolation warning counter every timestep (each MPI thread counts
        extrapolation warnings individually for its local atoms) (yes OR no, default: no).
    maxew = maximum number of extrapolation warnings allowed before the simulation is stopped
        (default: 0).

    <maxrcut> is the largest cutoff used by the symmetry functions in units of Angstrom (units
    are converted internally from Angstrom <-> Bohr and eV <-> Ha).

ATTENTION: The numbering of elements in your LAMMPS data file (containing the atomic starting
positions) has to be consistent with the order used in RuNNer. RuNNer sorts elements in order of
increasing nuclear charge, i.e. 
    
    RuNNer  nuclear charge    LAMMPS element index
    S       16                1
    Cu      29                2 

See the "tools" directory for converters from RuNNer and VASP format to LAMMPS data format.

Detailed comments on the example simulations can be found in the file "sim.lmp".

####################################################################################################
# LIMITATIONS AND IMPLEMENTED FEATURES 
####################################################################################################

Currently the LAMMPS implementation of the NN potential supports only a basic set of keywords in the
"input.nn" file.

Limitations:

    .) only short-range NN
    .) only global number of hidden layers
    .) only global node number for hidden layers
    .) only certain global activation functions implemented:
        - l (linear)
        - t (tanh)
    .) only certain symmetry functions implemented:
        - type 2
        - type 3

Implemented keywords:

    .) number_of_elements
    .) elements
    .) symfunction_short
    .) element_symfunction_short
    .) global_symfunction_short
    .) global_hidden_layers_short
    .) global_nodes_short
    .) global_activation_short
    .) normalize_nodes
    .) scale_symmetry_functions, scale_min_short, scale_max_short
    .) center_symmetry_functions

####################################################################################################
# PARALLELIZATION 
####################################################################################################

LAMMPS uses a spatial decomposition of the simulation domain to distribute work to MPI processes.
Each process holds a set of positions of local atoms (and ghost atoms) which are used by the
"runner" pair style to calculate energy contributions and forces.

If you intend to run a large simulation with the NN potential you may be interested in how many
cores C to use for a given system size of N atoms. As a general rule of thumb 

    N / C ~ 5-10 atoms/core

should give a reasonable parallel efficiency of about 70-80% on a typical cluster setup. This
behaviour has been tested on the VSC-2 (Vienna Scientific Cluster) which consists of nodes with 16
cores each (2x AMD Opteron 6132HE) connected via Infiniband network.

Example: Two-component system copper sulfide (Cu2S) with ~70 symmetry functions, 2 hidden layers
         with 30 nodes each:

         144  atoms: 16  cores N/C = 9.0 - parallel efficiency 79%
                     32  cores N/C = 4.5 - parallel efficiency 68%

         1152 atoms: 120 cores N/C = 9.6 - parallel efficiency 79%
                     232 cores N/C = 5.0 - parallel efficiency 66%

         3888 atoms: 400 cores N/C = 9.7 - parallel efficiency 84%
                     768 cores N/C = 5.1 - parallel efficiency 66%

####################################################################################################
# TOOLS 
####################################################################################################

The "tools" directory contains useful scripts for file conversion:

    .) runner2lammps.sh ... converts a RuNNer configuration (input.data) into a LAMMPS data file.
    .) poscar2lammps.sh ... converts a VASP configuration (POSCAR) into a LAMMPS data file.

USAGE: ./<tool>.sh < <input_file>
