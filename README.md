# Nonequilibrium Free-Energy Calculations of Fluids using LAMMPS

What is this?
------------
This repository contains a set of source codes and input scripts that allows to perform nonequilibrium free-energy calculations of fluid-phase systems using [LAMMPS](http://lammps.sandia.gov/) code. Details about the code implementation, capabilities and obtained results with these methods will be found in:

["Nonequilibrium Free-Energy Calculations of Fluids using LAMMPS"  
Rodolfo Paula Leite, Maurice de Koning  
Computational Material Science 159, 316-326 (2019)  
DOI: 10.1016/j.commatsci.2018.12.029](https://www.sciencedirect.com/science/article/pii/S0927025618308085)

What are the repository contents?
--------------
[`doc`](doc): This directory contains an updated user manual.

[`examples`](examples): This directory contains input scripts to run MD simulations.

source codes: The main directory contains our LAMMPS source codes.

[`README`](README.md): A brief overview of the distribution.

What is new in these source codes?
--------------
| Code Name                       | Already exists? |  Modification |
| :---                            |     :---:      |     :---      |
|angle.cpp / .h                   | yes            | Added extract method and reinit function.|
|angle_harmonic.cpp / .h          | yes            | Added a fscale variable (which multiplies only the forces) and the extract() method to adapt parameter over the time.|
|bond_harmonic.cpp / .h           | yes            | Added a fscale variable (which multiplies only the forces) and the extract() method to adapt parameter over the time.|
|fix_adapt.cpp / .h               | yes            | Added a fscale keyword for kspace pppm and pppm/tip4p styles to adapt only the forces during MD simulations. Also added the possibility to adapt angle potential (such as angle_harmonic that is used in the examples).|
|kspace.cpp / .h                  | yes            | Added a fscale variable in extract() method.                          |
|meam_force.cpp                   | yes            | Added a fscale variable into meam_force function.                     |
|meam.h                           | yes            | Added a fscale variable into meam_force function.                     |
|pair_lj_cut_coul_long.cpp / .h   | yes            | Added a fscale variable (which multiplies only the forces) in the extract() method.                           |
|pair_lj_cut_tip4p_long.cpp / .h  | yes            | Added a fscale variable (which multiplies only the forces) in the extract() method.                         |
|pair_meam.cpp / .h               | yes            | Added a fscale variable (which multiplies only the forces) and the extract() method to adapt parameter over the time.                         |
|pair_sw.cpp / .h                 | yes            | Added a fscale variable (which multiplies only the forces) and the extract() method to adapt parameter over the time.                         |
|pair_ufm.cpp / .h                | yes            | Added citation information inside the code.                          |
|pair_ufm_rw.cpp / .h             | no             | A new pair style which is based on LAMMPS implementation of TIP4P water model.                          |
|pppm_tip4p.cpp / .h              | yes            | Added an option to set fscale variable (which multiplies only the forces) when pppm_tip4p styles is invoked.                          |
|pppm.cpp / .h                    | yes            | Added an option to set fscale variable (which multiplies only the forces) when pppm styles is invoked.                         |

How to install?
--------------
To install these source codes, please follow the steps below:

1) Go to LAMMPS webpage (https://lammps.sandia.gov/download.html) and download the source code at your local machine.

2) Clone this repository to a subdirectory named `USER-FFE` inside the `src/` directory of your LAMMPS installation:
```
cd <your-local-lammps>/src/
git clone https://github.com/plrodolfo/FluidFreeEnergyforLAMMPS.git USER-FFE
```
WARNING: The next step will overwrite some native source codes in your src folder. However, our source codes have compatibility with old versions. Make a backup of your src folder if you want.

3) Choose some machine file (e.g. Makefile.mpi) and build LAMMPS using the following commands:

i) make yes-manybody

ii) make yes-kspace

iii) make yes-molecule

iv) make yes-rigid

v) make yes-user-meamc

vi) make yes-user-ffe

vii) make mpi

NOTE: Steps (i-vi) are necessary to install the required packages to reproduce the results presented in our [paper](https://www.sciencedirect.com/science/article/pii/S0927025618308085).

4) If LAMMPS was successully built, an executable called "lmp_mpi" will be created in the src directory. Otherwise, an error message is reported. For futher details, please visit the ["Build LAMMPS"](https://lammps.sandia.gov/doc/Build.html) section on user documentation.

How to use these codes?
--------------
Instructions of how to use these codes can be found inside each LAMMPS input script in the [`examples`](examples) directory and in our [paper](https://www.sciencedirect.com/science/article/pii/S0927025618308085).

Current compatibility:
--------------
lammps-22Aug18

Contact:
--------------
Rodolfo Paula Leite - pl.rodolfo@gmail.com
