/* ----------------------------------------------------------------------
  Flavio Giuliani 1aug2024 flavio.giuliani@uniroma1.it
------------------------------------------------------------------------- */

#include "pair_hybrid_scaled.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "memory.h"

#include "neighbor.h"
#include "neigh_request.h"

#include "respa.h"
#include "suffix.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairHybridScaled::PairHybridScaled(LAMMPS *lmp) :
    PairHybrid(lmp), fsum(NULL), tsum(NULL), scaleval(NULL), scaleidx(NULL)
{
  nmaxfsum = -1;
}

/* ---------------------------------------------------------------------- */

PairHybridScaled::~PairHybridScaled()
{
  memory->destroy(fsum);
  memory->destroy(tsum);
  delete[] scaleval;
  delete[] scaleidx;
}

/* ----------------------------------------------------------------------
  call each sub-style's compute() or compute_outer() function
  accumulate sub-style global/peratom energy/virial in hybrid
  for global vflag = VIRIAL_PAIR:
    each sub-style computes own virial[6]
    sum sub-style virial[6] to hybrid's virial[6]
  for global vflag = VIRIAL_FDOTR:
    call sub-style with adjusted vflag to prevent it calling
      virial_fdotr_compute()
    hybrid calls virial_fdotr_compute() on final accumulated f
------------------------------------------------------------------------- */

void PairHybridScaled::compute(int eflag, int vflag)
{
  int i, j, m, n;

  // update scale values from variables where needed

  const int nvars = scalevars.size();
  if (nvars > 0) {
    auto vals = new double[nvars];
    for (int k = 0; k < nvars; ++k) {
      int m = input->variable->find(const_cast<char*>(scalevars[k].c_str()));
      if (m < 0)
        error->all(FLERR, "Variable not found when updating scale factors");
      vals[k] = input->variable->compute_equal(m);
    }
    for (int k = 0; k < nstyles; ++k) {
      if (scaleidx[k] >= 0) scaleval[k] = vals[scaleidx[k]];
    }
    delete[] vals;
  }

  // if no_virial_fdotr_compute is set and global component of
  //   incoming vflag = 2, then
  // reset vflag as if global component were 1
  // necessary since one or more sub-styles cannot compute virial as F dot r

  if (no_virial_fdotr_compute && vflag % 4 == 2) vflag = 1 + vflag/4 * 4;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = eflag_global = vflag_global =
         eflag_atom = vflag_atom = 0;

  // grow fsum array if needed, and copy existing forces (usually 0.0) to it.

  if (atom->nmax > nmaxfsum) {
    memory->destroy(fsum);
    if (atom->torque_flag) memory->destroy(tsum);
    nmaxfsum = atom->nmax;
    memory->create(fsum, nmaxfsum, 3, "pair:fsum");
    if (atom->torque_flag) memory->create(tsum, nmaxfsum, 3, "pair:tsum");
  }
  const int nall = atom->nlocal + atom->nghost;
  auto f = atom->f;
  auto t = atom->torque;
  for (i = 0; i < nall; ++i) {
    fsum[i][0] = f[i][0];
    fsum[i][1] = f[i][1];
    fsum[i][2] = f[i][2];
    if (atom->torque_flag) {
      tsum[i][0] = t[i][0];
      tsum[i][1] = t[i][1];
      tsum[i][2] = t[i][2];
    }
  }

  // check if global component of incoming vflag = 2
  // if so, reset vflag passed to substyle as if it were 0
  // necessary so substyle will not invoke virial_fdotr_compute()

  int vflag_substyle;
  if (vflag % 4 == 2) vflag_substyle = vflag/4 * 4;
  else vflag_substyle = vflag;

  double *saved_special = save_special();

  // check if we are running with r-RESPA using the hybrid keyword

  Respa *respa = NULL;
  respaflag = 0;
  if (strstr(update->integrate_style,"respa")) {
    respa = (Respa *) update->integrate;
    if (respa->nhybrid_styles > 0) respaflag = 1;
  }

  for (m = 0; m < nstyles; m++) {

    // clear forces and torques

    memset(&f[0][0], 0, nall * 3 * sizeof(double));
    if (atom->torque_flag) memset(&t[0][0], 0, nall * 3 * sizeof(double));

    set_special(m);

    if (!respaflag || (respaflag && respa->hybrid_compute[m])) {

      // invoke compute() unless compute flag is turned off or
      // outerflag is set and sub-style has a compute_outer() method

      if (styles[m]->compute_flag == 0) continue;
      if (outerflag && styles[m]->respa_enable)
        styles[m]->compute_outer(eflag, vflag_substyle);
      else
        styles[m]->compute(eflag, vflag_substyle);
    }

    // add scaled forces to global sum
    const double scale = scaleval[m];
    for (i = 0; i < nall; ++i) {
      fsum[i][0] += scale * f[i][0];
      fsum[i][1] += scale * f[i][1];
      fsum[i][2] += scale * f[i][2];
      if (atom->torque_flag) {
        tsum[i][0] += scale * t[i][0];
        tsum[i][1] += scale * t[i][1];
        tsum[i][2] += scale * t[i][2];
      }
    }

    restore_special(saved_special);

    // jump to next sub-style if r-RESPA does not want global accumulated data

    if (respaflag && !respa->tally_global) continue;

    if (eflag_global) {
      eng_vdwl += scale * styles[m]->eng_vdwl;
      eng_coul += scale * styles[m]->eng_coul;
    }
    if (vflag_global) {
      for (n = 0; n < 6; n++) virial[n] += scale * styles[m]->virial[n];
    }
    if (eflag_atom) {
      n = atom->nlocal;
      if (force->newton_pair) n += atom->nghost;
      double *eatom_substyle = styles[m]->eatom;
      for (i = 0; i < n; i++) eatom[i] += scale * eatom_substyle[i];
    }
    if (vflag_atom) {
      n = atom->nlocal;
      if (force->newton_pair) n += atom->nghost;
      double **vatom_substyle = styles[m]->vatom;
      for (i = 0; i < n; i++)
        for (j = 0; j < 6; j++) vatom[i][j] += scale * vatom_substyle[i][j];
    }
  }

  // copy accumulated scaled forces to original force array

  for (i = 0; i < nall; ++i) {
    f[i][0] = fsum[i][0];
    f[i][1] = fsum[i][1];
    f[i][2] = fsum[i][2];
    if (atom->torque_flag) {
      t[i][0] = tsum[i][0];
      t[i][1] = tsum[i][1];
      t[i][2] = tsum[i][2];
    }
  }
  delete[] saved_special;

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   create one pair style for each arg in list
------------------------------------------------------------------------- */

void PairHybridScaled::settings(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR, "Illegal pair_style command");

  if (atom->avec->forceclearflag)
    error->all(FLERR, "Atom style is not compatible with pair_style hybrid/scaled");

  // delete old lists, since cannot just change settings

  if (nstyles > 0) {
    for (int m = 0; m < nstyles; m++) { delete styles[m]; }
    delete[] styles;
    for (int m = 0; m < nstyles; m++) delete [] keywords[m];
    delete [] keywords;
    delete[] scaleval;
    delete[] scaleidx;
    scalevars.clear();
  }

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cutghost);
    memory->destroy(nmap);
    memory->destroy(map);
  }
  allocated = 0;

  // allocate list of sub-styles as big as possibly needed if no extra args

  styles = new Pair *[narg];
  cutmax_style = new double[narg];
  memset(cutmax_style, 0.0, narg * sizeof(double));
  keywords = new char *[narg];
  multiple = new int[narg];

  special_lj = new double *[narg];
  special_coul = new double *[narg];
  compute_tally = new int[narg];

  scaleval = new double[narg];
  scaleidx = new int[narg];
  scalevars.reserve(narg);

  // allocate each sub-style
  // allocate uses suffix, but don't store suffix version in keywords,
  //   else syntax in coeff() will not match
  // call settings() with set of args that are not pair style names
  // use force->pair_map to determine which args these are

  int iarg, jarg, dummy;

  iarg = 0;
  nstyles = 0;
  while (iarg < narg - 1) {

    // first process scale factor or variable
    // idx < 0 indicates constant value otherwise index in variable name list

    double val = 0.0;
    int idx = -1;
    if (strstr(arg[iarg], "v_")) {
      for (std::size_t i = 0; i < scalevars.size(); ++i) {
        if (scalevars[i] == arg[iarg] + 2) {
          idx = i;
          break;
        }
      }
      if (idx < 0) {
        idx = scalevars.size();
        scalevars.emplace_back(arg[iarg] + 2);
      }
    } else {
      val = force->numeric(FLERR,arg[iarg]);
    }
    scaleval[nstyles] = val;
    scaleidx[nstyles] = idx;
    ++iarg;

    if (strcmp(arg[iarg],"hybrid") == 0)
      error->all(FLERR,"Pair style hybrid cannot have hybrid as an argument");
    if (strcmp(arg[iarg],"none") == 0)
      error->all(FLERR,"Pair style hybrid cannot have none as an argument");

    styles[nstyles] = force->new_pair(arg[iarg], 1, dummy);
    force->store_style(keywords[nstyles],arg[iarg],0);
    special_lj[nstyles] = special_coul[nstyles] = NULL;
    compute_tally[nstyles] = 1;

    // determine list of arguments for pair style settings
    // by looking for the next known pair style name.

    jarg = iarg + 1;
    while ((jarg < narg) && !force->pair_map->count(arg[jarg]))
      jarg++;

    // decrement to account for scale factor except when last argument

    if (jarg < narg) --jarg;

    styles[nstyles]->settings(jarg - iarg - 1, arg + iarg + 1);
    iarg = jarg;
    nstyles++;
  }

  // multiple[i] = 1 to M if sub-style used multiple times, else 0

  for (int i = 0; i < nstyles; i++) {
    int count = 0;
    for (int j = 0; j < nstyles; j++) {
      if (strcmp(keywords[j], keywords[i]) == 0) count++;
      if (j == i) multiple[i] = count;
    }
    if (count == 1) multiple[i] = 0;
  }

  // set pair flags from sub-style flags

  flags();
}

/* ----------------------------------------------------------------------
   call sub-style to compute single interaction
   error if sub-style does not support single() call
   since overlay could have multiple sub-styles, sum results explicitly
------------------------------------------------------------------------- */

double PairHybridScaled::single(int i, int j, int itype, int jtype, double rsq,
                                double factor_coul, double factor_lj,
                                double &fforce)
{
  if (nmap[itype][jtype] == 0) error->one(FLERR, "Invoked pair single on pair style none");

  // update scale values from variables where needed

  const int nvars = scalevars.size();
  if (nvars > 0) {
    auto vals = new double[nvars];
    for (int k = 0; k < nvars; ++k) {
      int m = input->variable->find(const_cast<char*>(scalevars[k].c_str()));
      if (m < 0)
        error->all(FLERR, "Variable not found when updating scale factors");
      vals[k] = input->variable->compute_equal(m);
    }
    for (int k = 0; k < nstyles; ++k) {
      if (scaleidx[k] >= 0) scaleval[k] = vals[scaleidx[k]];
    }
    delete[] vals;
  }

  double fone;
  fforce = 0.0;
  double esum = 0.0;

  for (int m = 0; m < nmap[itype][jtype]; m++) {
    auto pstyle = styles[map[itype][jtype][m]];
    if (rsq < pstyle->cutsq[itype][jtype]) {
      if (pstyle->single_enable == 0)
        error->one(FLERR, "Pair hybrid sub-style does not support single call");

      if ((special_lj[map[itype][jtype][m]] != NULL) ||
          (special_coul[map[itype][jtype][m]] != NULL))
        error->one(FLERR, "Pair hybrid single() does not support"
                  " per sub-style special_bond");

      double scale = scaleval[map[itype][jtype][m]];
      esum += scale * pstyle->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fone);
      fforce += scale * fone;
    }
  }

  if (single_extra) copy_svector(itype, jtype);
  return esum;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairHybridScaled::coeff(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds(FLERR, arg[0], atom->ntypes, ilo, ihi);
  force->bounds(FLERR, arg[1], atom->ntypes, jlo, jhi);

  // 3rd arg = pair sub-style name
  // 4th arg = pair sub-style index if name used multiple times
  // allow for "none" as valid sub-style name

  int multflag = 0;
  int m;

  for (m = 0; m < nstyles; m++) {
    multflag = 0;
    if (strcmp(arg[2], keywords[m]) == 0) {
      if (multiple[m]) {
        multflag = 1;
        if (narg < 4) error->all(FLERR, "Incorrect args for pair coefficients");
        if (!isdigit(arg[3][0]))
          error->all(FLERR,"Incorrect args for pair coefficients");
        int index = force->inumeric(FLERR,arg[3]);
        if (index == multiple[m]) break;
        else continue;
      } else break;
    }
  }

  int none = 0;
  if (m == nstyles) {
    if (strcmp(arg[2], "none") == 0) none = 1;
    else error->all(FLERR, "Pair coeff for hybrid has invalid style");
  }

  // move 1st/2nd args to 2nd/3rd args
  // if multflag: move 1st/2nd args to 3rd/4th args
  // just copy ptrs, since arg[] points into original input line

  arg[2 + multflag] = arg[1];
  arg[1 + multflag] = arg[0];

  // ensure that one_coeff flag is honored

  if (!none && styles[m]->one_coeff)
    if ((strcmp(arg[0], "*") != 0) || (strcmp(arg[1], "*") != 0))
      error->all(FLERR, "Incorrect args for pair coefficients");

  // invoke sub-style coeff() starting with 1st remaining arg

  if (!none) styles[m]->coeff(narg - 1 - multflag, &arg[1 + multflag]);

  // set setflag and which type pairs map to which sub-style
  // if sub-style is none: set hybrid subflag, wipe out map
  // else: set hybrid setflag & map only if substyle setflag is set
  //       if sub-style is new for type pair, add as multiple mapping
  //       if sub-style exists for type pair, don't add, just update coeffs

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      if (none) {
        setflag[i][j] = 1;
        nmap[i][j] = 0;
        count++;
      } else if (styles[m]->setflag[i][j]) {
        int k;
        for (k = 0; k < nmap[i][j]; k++)
          if (map[i][j][k] == m) break;
        if (k == nmap[i][j]) map[i][j][nmap[i][j]++] = m;
        setflag[i][j] = 1;
        count++;
      }
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairHybridScaled::write_restart(FILE *fp)
{
  PairHybrid::write_restart(fp);

  fwrite(scaleval, sizeof(double), nstyles, fp);
  fwrite(scaleidx, sizeof(int), nstyles, fp);

  int n = scalevars.size();
  fwrite(&n, sizeof(int), 1, fp);
  for (auto &var : scalevars) {
    n = var.size() + 1;
    fwrite(&n, sizeof(int), 1, fp);
    fwrite(var.c_str(), sizeof(char), n, fp);
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairHybridScaled::read_restart(FILE *fp)
{
  PairHybrid::read_restart(fp);

  delete[] scaleval;
  delete[] scaleidx;
  scalevars.clear();
  scaleval = new double[nstyles];
  scaleidx = new int[nstyles];

  int n, me = comm->me;
  if (me == 0) {
    fread(scaleval,sizeof(double),nstyles,fp);
    fread(scaleidx,sizeof(int),nstyles,fp);
  }
  MPI_Bcast(scaleval, nstyles, MPI_DOUBLE, 0, world);
  MPI_Bcast(scaleidx, nstyles, MPI_INT, 0, world);

  char *tmp;
  if (me == 0) fread(&n,sizeof(int),1,fp);
  MPI_Bcast(&n, 1, MPI_INT, 0, world);
  scalevars.resize(n);
  for (auto &scale : scalevars) {
    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    tmp = new char[n];
    if (me == 0) fread(tmp,sizeof(char),n,fp);
    MPI_Bcast(tmp, n, MPI_CHAR, 0, world);
    scale = tmp;
    delete[] tmp;
  }
}

/* ----------------------------------------------------------------------
   we need to handle Pair::svector special for hybrid/scaled
------------------------------------------------------------------------- */

void PairHybridScaled::init_svector()
{
  // single_extra = list all sub-style single_extra
  // allocate svector

  single_extra = 0;
  for (int m = 0; m < nstyles; m++) single_extra += styles[m]->single_extra;

  if (single_extra) {
    delete[] svector;
    svector = new double[single_extra];
  }
}

/* ----------------------------------------------------------------------
   we need to handle Pair::svector special for hybrid/scaled
------------------------------------------------------------------------- */

void PairHybridScaled::copy_svector(int itype, int jtype)
{
  int n = 0;
  Pair *this_style = NULL;

  // fill svector array.
  // copy data from active styles and use 0.0 for inactive ones
  for (int m = 0; m < nstyles; m++) {
    for (int k = 0; k < nmap[itype][jtype]; ++k) {
      if (m == map[itype][jtype][k]) {
        this_style = styles[m];
      } else {
        this_style = NULL;
      }
    }
    for (int l = 0; l < styles[m]->single_extra; ++l) {
      if (this_style) {
        svector[n++] = this_style->svector[l];
      } else {
        svector[n++] = 0.0;
      }
    }
  }
}
