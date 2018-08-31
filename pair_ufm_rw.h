/* -*- c++ -*- ----------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 http://lammps.sandia.gov, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov
 
 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.
 
 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
 Based on pair_tip4p_cut.ccp/h by Pavel Elkind (Gothenburg University)
 Contributing author:
            Rodolfo Paula Leite (Unicamp/Brazil) - pl.rodolfo@gmail.com
            Maurice de Koning (Unicamp/Brazil) - dekoning@ifi.unicamp.br
 ------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(ufm/rw,PairUFMRW)

#else

#ifndef LMP_PAIR_UFMRW_H
#define LMP_PAIR_UFMRW_H

#include "pair.h"

namespace LAMMPS_NS {

class PairUFMRW : public Pair {
 public:
  PairUFMRW(class LAMMPS *);
  virtual ~PairUFMRW();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);
  double memory_usage();
  void *extract(const char *, int &);

 protected:
  double cut_com_global;
  double cut_com,cut_comsq;
    
  int typeH,typeO;             // atom types of UF/Rw water H and O atoms
  double epsilon,sigma;        // energy and distance scales for UF/Rw
  double uf1,uf2,uf3,uf4;      // UFM parameters
  double **fscale;             // scale parameter used to perform TI

  int nmax;                    // info on molecule's center-of-mass
  int **hneigh;                // 0,1 = indices of 2 H associated with O
                               // 2 = 0 if site loc not yet computed, 1 if yes
  double **com;                // locations of the molecule center-of-mass
  
  void allocate();
  void compute_com(double *, double *, double *, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:
 
 E: UF/Rw hydrogen is missing
 
 The UF/Rw pairwise computation failed to find the correct H atom
 within a water molecule.
 
 E: UF/Rw hydrogen has incorrect atom type
 
 The UF/Rw pairwise computation found an H atom whose type does not
 agree with the specified H type.
 
 E: Illegal ... command
 
 Self-explanatory.  Check the input script syntax and compare to the
 documentation for the command.  You can use -echo screen as a
 command-line option when running LAMMPS to see the offending line.
 
 E: Incorrect args for pair coefficients
 
 Self-explanatory.  Check the input script or data file.
 
 E: Pair style ufm/rw requires atom IDs
 
 This is a requirement to use this potential.
 
 E: Pair style ufm/rw requires newton pair on
 
 See the newton command.  This is a restriction to use this potential.
 
 E: Pair style ufm/rw requires atom attribute q
 
 The atom style defined does not have this attribute.
 
 E: Must use a bond style with UF/Rw potential
 
UF/Rw potentials assume bond lengths in water are constrained
 by a fix shake command.
 
 E: Must use an angle style with UF/Rw potential
 
UF/Rw potentials assume angles in water are constrained by a fix shake
 command.
 
 */
