/* -------------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------
   Contributing authors:
             Rodrigo Freitas (UC Berkeley) - rodrigof@berkeley.edu
             Mark Asta (UC Berkeley) - mdasta@berkeley.edu
             Maurice de Koning (Unicamp/Brazil) - dekoning@ifi.unicamp.br
------------------------------------------------------------------------- */

#include <stdlib.h>
#include <string.h>
#include "fix_ti_spring_mod.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "respa.h"
#include "memory.h"
#include "error.h"
#include "citeme.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

static const char cite_fix_ti_spring[] =
  "ti/spring command:\n\n"
  "@article{freitas2016,\n"
  "  author={Freitas, Rodrigo and Asta, Mark and de Koning, Maurice},\n"
  "  title={Nonequilibrium free-energy calculation of solids using LAMMPS},\n"
  "  journal={Computational Materials Science},\n"
  "  volume={112},\n"
  "  pages={333--341},\n"
  "  year={2016},\n"
  "  publisher={Elsevier}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

FixTISpringMod::FixTISpringMod(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_ti_spring);

  if (narg != 5 )
    error->all(FLERR,"Illegal fix ti/spring command");

  // Flags.
  restart_peratom = 1;
  scalar_flag = 1;
  global_freq = 1;
  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  // disallow resetting the time step, while this fix is defined
  time_depend = 1;

  // Spring constant.
  k = force->numeric(FLERR,arg[3]);
  if (k <= 0.0) error->all(FLERR,"Illegal fix ti/spring/mod command");

  // Perform initial allocation of atom-based array
  // Register with Atom class
  xoriginal = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  // xoriginal = initial unwrapped positions of atoms

  double **x = atom->x;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) domain->unmap(x[i],image[i],xoriginal[i]);
    else xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;
  }

  // Time variables.
  t0 = update->ntimestep;  // timestep of original/starting coordinates
  lambda = force->numeric(FLERR,arg[4]); // Coupling parameter
  if ((lambda < 0) || (lambda > 1))
    error->all(FLERR,"Illegal fix ti/spring/mod command");

  espring = 0.0;
}

/* ---------------------------------------------------------------------- */

FixTISpringMod::~FixTISpringMod()
{
  // unregister callbacks to this fix from Atom class
  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  // delete locally stored array
  memory->destroy(xoriginal);
}

/* ---------------------------------------------------------------------- */

int FixTISpringMod::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTISpringMod::init()
{
  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixTISpringMod::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixTISpringMod::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTISpringMod::post_force(int vflag)
{

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  double dx, dy, dz;
  double unwrap[3];

  espring = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - xoriginal[i][0];
      dy = unwrap[1] - xoriginal[i][1];
      dz = unwrap[2] - xoriginal[i][2];
      f[i][0] = (1-lambda) * f[i][0] + lambda * (-k*dx);
      f[i][1] = (1-lambda) * f[i][1] + lambda * (-k*dy);
      f[i][2] = (1-lambda) * f[i][2] + lambda * (-k*dz);
      espring += k * (dx*dx + dy*dy + dz*dz);
    }

  espring *= 0.5;
}

/* ---------------------------------------------------------------------- */

void FixTISpringMod::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTISpringMod::min_post_force(int vflag)
{
  post_force(vflag);
}


/* ----------------------------------------------------------------------
   energy of stretched springs
------------------------------------------------------------------------- */

double FixTISpringMod::compute_scalar()
{
  double all;
  MPI_Allreduce(&espring,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}


/* ----------------------------------------------------------------------
     memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixTISpringMod::memory_usage()
{
  double bytes = atom->nmax*3 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
     allocate atom-based array
------------------------------------------------------------------------- */

void FixTISpringMod::grow_arrays(int nmax)
{
  memory->grow(xoriginal,nmax,3,"fix_ti/spring/mod:xoriginal");
}

/* ----------------------------------------------------------------------
     copy values within local atom-based array
------------------------------------------------------------------------- */

void FixTISpringMod::copy_arrays(int i, int j, int delflag)
{
  xoriginal[j][0] = xoriginal[i][0];
  xoriginal[j][1] = xoriginal[i][1];
  xoriginal[j][2] = xoriginal[i][2];
}

/* ----------------------------------------------------------------------
    pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixTISpringMod::pack_exchange(int i, double *buf)
{
  buf[0] = xoriginal[i][0];
  buf[1] = xoriginal[i][1];
  buf[2] = xoriginal[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
    unpack values in local atom-based array from exchange with another proc
 ------------------------------------------------------------------------- */

int FixTISpringMod::unpack_exchange(int nlocal, double *buf)
{
  xoriginal[nlocal][0] = buf[0];
  xoriginal[nlocal][1] = buf[1];
  xoriginal[nlocal][2] = buf[2];
  return 3;
}

/* ----------------------------------------------------------------------
    pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixTISpringMod::pack_restart(int i, double *buf)
{
  buf[0] = 4;
  buf[1] = xoriginal[i][0];
  buf[2] = xoriginal[i][1];
  buf[3] = xoriginal[i][2];
  return 4;
}

/* ----------------------------------------------------------------------
    unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixTISpringMod::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  xoriginal[nlocal][0] = extra[nlocal][m++];
  xoriginal[nlocal][1] = extra[nlocal][m++];
  xoriginal[nlocal][2] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
     maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixTISpringMod::maxsize_restart()
{
  return 4;
}

/* ----------------------------------------------------------------------
     size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixTISpringMod::size_restart(int nlocal)
{
  return 4;
}
