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

/* -----------------------------------------------------------------------
Based on pair_tip4p_cut.ccp by Pavel Elkind (Gothenburg University)
Contributing author:
            Rodolfo Paula Leite (Unicamp/Brazil) - pl.rodolfo@gmail.com
            Maurice de Koning (Unicamp/Brazil) - dekoning@ifi.unicamp.br
 ------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_ufm_rw.h"
#include "atom.h"
#include "angle.h"
#include "bond.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "domain.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "citeme.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/*
static const char cite_pair_ufm_rw[] =
"UF/Rw pair style:\n\n"
"@article{PaulaLeite2018,\n"
"  author={Paula Leite, Rodolfo and de Koning, Maurice},\n"
"  title={Nonequilibrium free-energy calculations of fluids using LAMMPS},\n"
"  journal={Computational Materials Science},\n"
"  volume={XXX},\n"
"  pages={XX--XX},\n"
"  year={2018},\n"
"  publisher={Elsevier}\n"
"}\n\n";
 */

/* ---------------------------------------------------------------------- */

PairUFMRW::PairUFMRW(LAMMPS *lmp) : Pair(lmp)
{

//  if (lmp->citeme) lmp->citeme->add(cite_pair_ufm_rw);

  single_enable = 0;

  nmax = 0;
  hneigh = NULL;
  com = NULL;
}

/* ---------------------------------------------------------------------- */

PairUFMRW::~PairUFMRW()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(hneigh);
    memory->destroy(com);
    memory->destroy(fscale);
  }
}

/* ---------------------------------------------------------------------- */

void PairUFMRW::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fcom;
  double rsq, expuf, factor;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double *x1,*x2,*xH1,*xH2;
  int iH1,iH2,jH1,jH2;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  tagint *tag = atom->tag;
  double *m = atom->mass;
  double M = m[1] + 2 * m[2];
  int *type = atom->type;
    
  if (atom->nmax > nmax) {
      nmax = atom->nmax;
      memory->destroy(hneigh);
      memory->create(hneigh,nmax,3,"pair:hneigh");
      memory->destroy(com);
      memory->create(com,nmax,3,"pair:com");
  }
    
  if (neighbor->ago == 0)
      for (i = 0; i < nall; i++) hneigh[i][0] = -1;
  for (i = 0; i < nall; i++) hneigh[i][2] = 0;
    
  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;

  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
      
    if (itype == typeO) {
        if (hneigh[i][0] < 0) {
            iH1 = atom->map(tag[i] + 1);
            iH2 = atom->map(tag[i] + 2);
            if (iH1 == -1 || iH2 == -1)
                error->one(FLERR,"UF/Rw hydrogen is missing");
            if (atom->type[iH1] != typeH || atom->type[iH2] != typeH)
                error->one(FLERR,"UF/Rw hydrogen has incorrect atom type");
            // set iH1,iH2 to closest image to O
            iH1 = domain->closest_image(i,iH1);
            iH2 = domain->closest_image(i,iH2);
            compute_com(x[i],x[iH1],x[iH2],com[i]);
            hneigh[i][0] = iH1;
            hneigh[i][1] = iH2;
            hneigh[i][2] = 1;
        } else {
            iH1 = hneigh[i][0];
            iH2 = hneigh[i][1];
            if (hneigh[i][2] == 0) {
                hneigh[i][2] = 1;
                compute_com(x[i],x[iH1],x[iH2],com[i]);
            }
        }
        x1 = com[i];
    }
      
    jlist = firstneigh[i];
    jnum = numneigh[i];
      
    for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        factor = special_lj[sbmask(j)];
        j &= NEIGHMASK;
                
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        jtype = type[j];
        
        if (rsq < cut_comsq) {
            if (itype == typeO && jtype == typeO) {

                if (jtype == typeO) {
                    if (hneigh[j][0] < 0) {
                        jH1 = atom->map(tag[j] + 1);
                        jH2 = atom->map(tag[j] + 2);
                        if (jH1 == -1 || jH2 == -1)
                            error->one(FLERR,"UF/Rw hydrogen is missing");
                        if (atom->type[jH1] != typeH || atom->type[jH2] != typeH)
                            error->one(FLERR,"UF/Rw hydrogen has incorrect atom type");
                        // set jH1,jH2 to closest image to O
                        jH1 = domain->closest_image(j,jH1);
                        jH2 = domain->closest_image(j,jH2);
                        compute_com(x[j],x[jH1],x[jH2],com[j]);
                        hneigh[j][0] = jH1;
                        hneigh[j][1] = jH2;
                        hneigh[j][2] = 1;
                    } else {
                        jH1 = hneigh[j][0];
                        jH2 = hneigh[j][1];
                        if (hneigh[j][2] == 0) {
                            hneigh[j][2] = 1;
                            compute_com(x[j],x[jH1],x[jH2],com[j]);
                        }
                    }
                    x2 = com[j];
                }
                delx = x1[0] - x2[0];
                dely = x1[1] - x2[1];
                delz = x1[2] - x2[2];
                rsq = delx*delx + dely*dely + delz*delz;
            }
      // Uhlenbeck-Ford interaction based on molecule center-of-mass.
        
             if (rsq < cut_comsq && itype == typeO && jtype == typeO) {
                 expuf = exp(- rsq * uf2);
                 fcom = factor * fscale[itype][jtype] * uf1 * expuf /(1.0 - expuf);
        
                 f[i][0] += delx*fcom*m[typeO]/M;
                 f[i][1] += dely*fcom*m[typeO]/M;
                 f[i][2] += delz*fcom*m[typeO]/M;
        
                 f[iH1][0] += delx*fcom*m[typeH]/M;
                 f[iH1][1] += dely*fcom*m[typeH]/M;
                 f[iH1][2] += delz*fcom*m[typeH]/M;
        
                 f[iH2][0] += delx*fcom*m[typeH]/M;
                 f[iH2][1] += dely*fcom*m[typeH]/M;
                 f[iH2][2] += delz*fcom*m[typeH]/M;
        
                 if (newton_pair || j < nlocal) {
                     f[j][0] -= delx*fcom*m[typeO]/M;
                     f[j][1] -= dely*fcom*m[typeO]/M;
                     f[j][2] -= delz*fcom*m[typeO]/M;
            
                     f[jH1][0] -= delx*fcom*m[typeH]/M;
                     f[jH1][1] -= dely*fcom*m[typeH]/M;
                     f[jH1][2] -= delz*fcom*m[typeH]/M;
            
                     f[jH2][0] -= delx*fcom*m[typeH]/M;
                     f[jH2][1] -= dely*fcom*m[typeH]/M;
                     f[jH2][2] -= delz*fcom*m[typeH]/M;
                 }

                 if (eflag) {
                     evdwl = -uf3 * log(1.0 - expuf);
                     evdwl *= factor;
                 }

                 if (evflag) ev_tally(i,j,nlocal,newton_pair,
                                      evdwl,0.0,fcom,delx,dely,delz);
             }
        }
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairUFMRW::allocate()
{
  allocated = 1;

  int n = atom->ntypes;
    
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(fscale,n+1,n+1,"pair:fscale");
    
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairUFMRW::settings(int narg, char **arg)
{
  if (narg != 5) error->all(FLERR,"Illegal pair_style command");

  typeO = force->inumeric(FLERR,arg[0]);
  typeH = force->inumeric(FLERR,arg[1]);
  epsilon = force->numeric(FLERR,arg[2]);
  sigma = force->numeric(FLERR,arg[3]);
  cut_com = force->numeric(FLERR,arg[4]);
  cut_comsq = cut_com * cut_com;
    
  if (sigma <= 0.0) error->all(FLERR,"UF sigma parameter must be > 0.0");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairUFMRW::coeff(int narg, char **arg)
{
  if (narg != 2)
      error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      setflag[i][j] = 1;
      fscale[i][j] = 1.0;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init specific to this pair style
 ------------------------------------------------------------------------- */

void PairUFMRW::init_style()
{
    if (atom->tag_enable == 0)
        error->all(FLERR,"Pair style ufm/rw requires atom IDs");
    if (!force->newton_pair)
        error->all(FLERR,"Pair style ufm/rw requires newton pair on");
    if (force->bond == NULL)
        error->all(FLERR,"Must use a bond style with UF/Rw potential");
    if (force->angle == NULL)
        error->all(FLERR,"Must use an angle style with UF/Rw potential");
    
    neighbor->request(this,instance_me);
    
    // set the uf parameters
    
    uf1 = 2.0 * epsilon / pow(sigma,2.0);
    uf2 = 1.0 / pow(sigma,2.0);
    uf3 = epsilon;
    uf4 = sigma;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairUFMRW::init_one(int i, int j)
{
  fscale[j][i] = fscale[i][j];
  return cut_com;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairUFMRW::write_restart(FILE *fp)
{
  write_restart_settings(fp);
  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++)
      fwrite(&setflag[i][j],sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairUFMRW::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairUFMRW::write_restart_settings(FILE *fp)
{
    fwrite(&typeO,sizeof(int),1,fp);
    fwrite(&typeH,sizeof(int),1,fp);
    fwrite(&epsilon,sizeof(double),1,fp);
    fwrite(&sigma,sizeof(double),1,fp);
    fwrite(&cut_com,sizeof(double),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairUFMRW::read_restart_settings(FILE *fp)
{
    if (comm->me == 0) {
        fread(&typeO,sizeof(int),1,fp);
        fread(&typeH,sizeof(int),1,fp);
        fread(&epsilon,sizeof(double),1,fp);
        fread(&sigma,sizeof(double),1,fp);
        fread(&cut_com,sizeof(double),1,fp);
    }
    
    MPI_Bcast(&typeO,1,MPI_INT,0,world);
    MPI_Bcast(&typeH,1,MPI_INT,0,world);
    MPI_Bcast(&epsilon,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&sigma,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&cut_com,1,MPI_DOUBLE,0,world);
    
    cut_comsq = cut_com * cut_com;
}

/* ----------------------------------------------------------------------
 compute the center-of-mass position xCOM between O atom and 2 H atoms
 return it as xCOM
 ------------------------------------------------------------------------- */

void PairUFMRW::compute_com(double *xO, double *xH1, double *xH2, double *xCOM)
{
    double delx1 = xH1[0] - xO[0];
    double dely1 = xH1[1] - xO[1];
    double delz1 = xH1[2] - xO[2];
    
    double delx2 = xH2[0] - xO[0];
    double dely2 = xH2[1] - xO[1];
    double delz2 = xH2[2] - xO[2];
    
    double *m = atom->mass;
    double M = m[typeO] + 2 * m[typeH];
    
    xCOM[0] = xO[0] + m[typeH] * (delx1 + delx2) / M;
    xCOM[1] = xO[1] + m[typeH] * (dely1 + dely2) / M;
    xCOM[2] = xO[2] + m[typeH] * (delz1 + delz2) / M;
}

/* ----------------------------------------------------------------------
 memory usage of hneigh
 ------------------------------------------------------------------------- */

double PairUFMRW::memory_usage()
{
    double bytes = maxeatom * sizeof(double);
    bytes += maxvatom*6 * sizeof(double);
    bytes += 2 * nmax * sizeof(double);
    return bytes;
}

/* ---------------------------------------------------------------------- */

void *PairUFMRW::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"epsilon") == 0) return (void *) &epsilon;
  if (strcmp(str,"sigma") == 0) return (void *) &sigma;
  dim = 2;
  if (strcmp(str,"fscale") == 0) return (void *) fscale;
  return NULL;
}
