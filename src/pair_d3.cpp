////////////////////////////////////////////////////////////////////////////////
//
// D3 pair style for LAMMPS
// https://www.chemie.uni-bonn.de/grimme/de/software/dft-d3
// ----------------------------
//
// authors: Flavio Giuliani, Riccardo Piombo
// date   : 2024-03
// email  : flavio.giuliani@uniroma1.it
// email  : riccardo.piombo@uniroma1.it
//
////////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "pair_d3.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "citeme.h"

using namespace LAMMPS_NS;

static const char cite_d3[] =
   "Pair style D3:\n\n"
   "@Article{pair_d3_2010,\n"
   "author = {Grimme, Stefan and Antony, Jens and Ehrlich, Stephan and Krieg, Helge},\n"
   "title = {A consistent and accurate ab initio parametrization of density functional dispersion correction (DFT-D) for the 94 elements H-Pu},\n"
   "journal = {The Journal of Chemical Physics},\n"
   "volume = {132},\n"
   "number = {15},\n"
   "pages = {154104},\n"
   "year = {2010},\n"
    "month = {04},\n"
    "issn = {0021-9606},\n"
    "doi = {10.1063/1.3382344},\n"
    "url = {https://doi.org/10.1063/1.3382344},\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */

PairD3::PairD3(LAMMPS *lmp) : Pair(lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_d3);
  manybody_flag = 1;
  //one_coeff = 1;
  restartinfo = 0;
  single_enable = 0;
  // set comm size needed by this Pair
  comm_forward = 1;
  //comm_reverse = 1;

  CN = NULL; // null pointer

  maxshort = 100;
  neighshort = NULL;
}

/* ---------------------------------------------------------------------- */

PairD3::~PairD3()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);

    memory->destroy(CN);

    memory->destroy(Z);
    memory->destroy(rcov);
    memory->destroy(Q);
    memory->destroy(r0);

    memory->destroy(neighshort);

    nmax=0;
  }
}

/* ---------------------------------------------------------------------- */

void PairD3::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,kk,inum,jnum,jnum1,knum,itype,jtype,ktype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double r,rr,CNtmp,rsq,factor_lj, c6,c8;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double delr1[3],delr2[3],delr3[3],rsq1,rsq2,rsq3,fj[3],fk[3];
  double fxtmp,fytmp,fztmp, c6_threebody[3],c9;

  // Inizializza l'energia di van der Waals
  evdwl = 0.0;

  // Configura i calcoli di energia e forza se richiesti
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  // grow energy and fp arrays if necessary
  // need to be atom->nmax in length
  if (atom->nmax > nmax) {
    memory->destroy(CN);
    nmax = atom->nmax;
    memory->create(CN,nmax,"pair:CN");

    if (debugging_mode >= 1) {
      if(comm->me == 0) {
          if(screen) fprintf(screen, "PairD3: CN reallocated to %d\n",nmax);
          if(logfile) fprintf(logfile, "PairD3: CN reallocated to %d\n",nmax);
      }
      if(screen) fprintf(screen, "PairD3: comm->me=%d. I have nlocal = %d nghost = %d\n", comm->me,atom->nlocal,atom->nghost);
      if(logfile) fprintf(logfile, "PairD3: comm->me=%d. I have nlocal = %d nghost = %d\n", comm->me,atom->nlocal,atom->nghost);
    }
  }


  // Estrai i dati necessari dalla struttura atomica e dalla lista delle interazioni
  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero out and compute Coordination Number for each atom

  if (newton_pair) {
    for (i = 0; i < nall; i++) CN[i] = 0.0;
  } else
    for (i = 0; i < nlocal; i++) CN[i] = 0.0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];


    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cut_CN_sq) {

        jtype = type[j];

        r = sqrt(rsq);

        rr = rcov[itype][jtype]/r; // note: k2 is included in rcov
        CNtmp = 1./( 1.+exp(-k1*(rr-1.)) );
        CN[i] += CNtmp;

        if (newton_pair || j < nlocal) {
          CN[j] += CNtmp;
        }

      }
    }
  }

  // communicate CNs between processors !!

  if (newton_pair) comm->reverse_comm_pair(this);
  //if (newton_pair) comm->forward_comm_pair(this);

  // print Coordination Numbers, C6(i,i) and C8(i,i)
  if (debugging_mode >= 1) {

    if(screen) fprintf(screen, "PairD3: comm->me = %d. SUMMARY: i , x,y,z(Ang.), Z,rcov(Ang.),r0(Ang.), CN , C6(eV*Ang^6),C8(eV*Ang^8), C6(Ry*Bohr^6),C8(Ry*Bohr^8)\n", comm->me);
    if(logfile) fprintf(logfile, "PairD3: comm->me = %d. SUMMARY: i , x,y,z(Ang.), Z,rcov(Ang.),r0(Ang.), CN , C6(eV*Ang^6),C8(eV*Ang^8), C6(Ry*Bohr^6),C8(Ry*Bohr^8)\n", comm->me);
    double RyBohr6_to_eVang6 = RY_TO_EV*pow(AU_TO_ANG,6.);
    double eVang6_to_RyBohr6 = 1./RyBohr6_to_eVang6;
    double RyBohr8_to_eVang8 = RY_TO_EV*pow(AU_TO_ANG,8.);
    double eVang8_to_RyBohr8 = 1./RyBohr8_to_eVang8;

    for (i = 0; i < nlocal; i++) {

      c6 = interpolate_c6(i, i, itype, itype);
      c8 = c6*Q[itype][itype];

      if(screen) fprintf(screen, "PairD3: %d  %f %f %f  %d %f %f  %f  %f %f  %f %f\n",
        i, x[i][0],x[i][1],x[i][2], Z[itype],0.5*rcov[itype][itype],0.5*r0[itype][itype], CN[i], c6,c8, c6*eVang6_to_RyBohr6,c8*eVang8_to_RyBohr8);
      if(logfile) fprintf(logfile, "PairD3: %d  %f %f %f  %d %f %f  %f  %f %f  %f %f\n",
        i, x[i][0],x[i][1],x[i][2], Z[itype],0.5*rcov[itype][itype],0.5*r0[itype][itype], CN[i], c6,c8, c6*eVang6_to_RyBohr6,c8*eVang8_to_RyBohr8);

    }
  }

  // compute 2-body & 3-body forces and energies on each atom

  // Ciclo sugli atomi i
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    // accumulate forces on atom i from: each j in 2body, each j,k in 3body
    fxtmp = fytmp = fztmp = 0.0;

    int numshort=0;

    // Ciclo sui primi vicini dell'atomo i
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      //factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      // Calcola la distanza ij tra gli atomi
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      // Se la distanza è inferiore alla soglia di cutoff

      if (rsq < cutsq[itype][jtype]) {

        // build on the fly a neighbour list for the 3-body interaction
        if (include_threebody == 1 && rsq < cut_threebody_sq) {
          neighshort[numshort++] = j;
          if (numshort >= maxshort) {
            maxshort += maxshort/2;
            memory->grow(neighshort,maxshort,"pair:neighshort");
          }
        }

        twobody(i,j,itype,jtype,rsq,fpair,eflag,evdwl);

        fxtmp += delx*fpair; // accumulate forces on atom i
        fytmp += dely*fpair;
        fztmp += delz*fpair;

        // Apply Newton's 3rd law to compute Fji = -Fij.
        // (only if half-list or if j is not a ghost atom)
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair; // sum force on atom j from i
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }

    // three-body interaction

    jnum1 = numshort - 1;
    if (debugging_mode >= 1 && comm->me == 0) {
      if(screen) fprintf(screen, "PairD3: numshort = %d\n",numshort);
      if(logfile) fprintf(logfile, "PairD3: numshort = %d\n",numshort);
    }

    for (jj = 0; jj < jnum1; jj++) {
      if(include_threebody == 0) continue; // skip 3-body interaction

      j = neighshort[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;
      jtype = type[j];

      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];

      c6_threebody[0] = interpolate_c6(i, j, itype, jtype);

      // at fixed i, accumulate forces on atom j from each k in 3body terms
      double fjxtmp,fjytmp,fjztmp;
      fjxtmp = fjytmp = fjztmp = 0.0;

      for (kk = jj+1; kk < numshort; kk++) {
        k = neighshort[kk];
        k &= NEIGHMASK;
        ktype = type[k];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

        delr3[0] = x[k][0] - x[j][0];
        delr3[1] = x[k][1] - x[j][1];
        delr3[2] = x[k][2] - x[j][2];
        rsq3 = delr3[0]*delr3[0] + delr3[1]*delr3[1] + delr3[2]*delr3[2];

        // all 3 distances rij.rik.rjk must be within the cutoff
        if (rsq3 >= cut_threebody_sq) continue;

        c6_threebody[1] = interpolate_c6(i, k, itype, ktype);
        c6_threebody[2] = interpolate_c6(j, k, jtype, ktype);
        c9 = -1.0*pow(c6_threebody[0]*c6_threebody[1]*c6_threebody[2], 0.5);

        threebody(itype,jtype,ktype,c9,rsq1,rsq2,rsq3,delr1,delr2,delr3,fj,fk,eflag,evdwl);

        fxtmp -= fj[0] + fk[0]; // accumulate forces on atom i for each j,k
        fytmp -= fj[1] + fk[1];
        fztmp -= fj[2] + fk[2];
        fjxtmp += fj[0]; // accumulate forces on atom j for each k
        fjytmp += fj[1];
        fjztmp += fj[2];
        f[k][0] += fk[0]; // sum forces on atom k from i,j
        f[k][1] += fk[1];
        f[k][2] += fk[2];

        if (evflag) ev_tally3(i,j,k,evdwl,0.0,fj,fk,delr1,delr2);
      }

      f[j][0] += fjxtmp; // sum accumulated forces on atom j
      f[j][1] += fjytmp;
      f[j][2] += fjztmp;

    }

    f[i][0] += fxtmp; // sum accumulated forces on atom i
    f[i][1] += fytmp;
    f[i][2] += fztmp;

  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairD3::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");

  memory->create(Z,n+1,"pair:Z");

  memory->create(rcov,n+1,n+1,"pair:rcov");
  memory->create(Q,n+1,n+1,"pair:Q");
  memory->create(r0,n+1,n+1,"pair:r0");

  memory->create(neighshort,maxshort,"pair:neighshort");

}

/* ----------------------------------------------------------------------
   global settings - questi sono i numeretti dopo pair style
------------------------------------------------------------------------- */

void PairD3::settings(int narg, char **arg)
{
  if (narg < 5 || narg > 7) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);
  sR6 = force->numeric(FLERR,arg[1]);
  s8 =  force->numeric(FLERR,arg[2]);
  double cut_3 = force->numeric(FLERR,arg[3]);
  cut_threebody_sq = cut_3*cut_3;
  double cut_CN = force->numeric(FLERR,arg[4]);
  cut_CN_sq=cut_CN*cut_CN;

  include_threebody = 1; // default
  if (narg >= 6) include_threebody = force->numeric(FLERR,arg[5]);
  if (include_threebody < 0 || include_threebody > 1)
    error->all(FLERR,"Illegal include_threebody in pair_style. Options are 0,1");

  debugging_mode = 0; // default
  if (narg >= 7) debugging_mode = force->numeric(FLERR,arg[6]);
  if (debugging_mode < 0 || debugging_mode > 2)
    error->all(FLERR,"Illegal debugging_mode in pair_style. Options are 0,1,2");

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs - questi sono i numeretti dopo pair coeff
------------------------------------------------------------------------- */
//MODFICATO
void PairD3::coeff(int narg, char **arg)
{
  if (narg != 2)
    error->all(FLERR,"Incorrect args for pair coefficients");

  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);

  int Z_one = force->numeric(FLERR,arg[1]);
  if(Z_one<1 || Z_one>94)
    error->all(FLERR,"Invalid value for pair coefficient 'Z'");
  if(Z_one!=51 || atom->ntypes>1)
    error->all(FLERR,"Currently implemented only monospecies Sb (Z=51)");

  int k, count = 0;
  for (int i = ilo; i <= ihi; i++) {
    Z[i] = Z_one;
    rcov[i][i] = (rcov_Z[Z[i]]+rcov_Z[Z[i]])*AU_TO_ANG;  // for CN computation
    Q[i][i] = 3.0*SQUARE(sqrtQ_Z[Z[i]]*AU_TO_ANG);         // for c8 computation
    k = get_pair_index(Z[i],Z[i],1,94);
    r0[i][i] = r0_ZZ[k]*AU_TO_ANG;
    if(debugging_mode>=1 && comm->me == 0) {
        if(screen) fprintf(screen,
          "PairD3: diagonal coefficients (Angstrom): type=%d \t rcov=%f \t Q=%f \t r0=%f\n",
          i,rcov[i][i], Q[i][i], r0[i][i]
        );
        if(logfile) fprintf(logfile,
          "PairD3: diagonal coefficients (Angstrom): type=%d \t rcov=%f \t Q=%f \t r0=%f\n",
          i,rcov[i][i], Q[i][i], r0[i][i]
        );
    }
    setflag[i][i] = 1;
    cut[i][i] = cut_global;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");

}

/* ----------------------------------------------------------------------
   initialuzation of pair style
--------------------------------------------------------------------------*/
void PairD3::init_style() {

    int irequest = neighbor->request(this,instance_me);

    if(setup_completed == 0) {
        if(comm->me == 0) {
            if(screen) fprintf(screen, "PairD3: Starting pair style setup...\n");
            if(logfile) fprintf(logfile, "PairD3: Starting pair style setup...\n");
        }

        init_vars();

        setup_completed = 1;

        if(comm->me == 0) {
            if(screen) fprintf(screen, "PairD3: pair style setup completed. include_threebody = %d , debugging_mode = %d\n",include_threebody,debugging_mode);
            if(logfile) fprintf(logfile, "PairD3: pair style setup completed. include_threebody = %d , debugging_mode = %d\n",include_threebody,debugging_mode);

        }
    }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairD3::init_one(int i, int j)
{
  if (setflag[i][i] == 0 || setflag[j][j]==0)
    error->all(FLERR,"Not all diagonal pair coeffs are set");

  cut[i][j] = cut[j][i] = cut_global;
  // Combination rules: Qij=sqrt(Qi*Qj), rcov_ij = rcov_i+rcov_j
  rcov[i][j] = (rcov_Z[Z[i]]+rcov_Z[Z[j]])*AU_TO_ANG;  // for CN computation
  Q[i][j] = 3.0*sqrtQ_Z[Z[i]]*sqrtQ_Z[Z[j]]*SQUARE(AU_TO_ANG);         // for c8 computation
  int k = get_pair_index(Z[i],Z[j],1,94);
  r0[i][j] = r0_ZZ[k];

  rcov[j][i] = rcov[i][j];
  Q[j][i] = Q[i][j];
  r0[j][i] = r0[i][j];

  setflag[i][j] = setflag[j][i] = 1;

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairD3::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++) {
    fwrite(&setflag[i][j],sizeof(int),1,fp);
    fwrite(&Z[i],sizeof(int),1,fp);
  }
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairD3::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  size_t foo;

  for (i = 1; i <= atom->ntypes; i++) {
    if (me == 0) foo=fread(&setflag[i][i],sizeof(int),1,fp);
    MPI_Bcast(&setflag[i][i],1,MPI_INT,0,world);
    if (setflag[i][i]) {
      if (me == 0) {
        foo=fread(&Z[i],sizeof(double),1,fp);
      }
      MPI_Bcast(&Z[i],1,MPI_DOUBLE,0,world);
    }
  }

  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) foo=fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          foo=fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairD3::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&sR6,sizeof(double),1,fp);
  fwrite(&s8,sizeof(double),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairD3::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    size_t foo;
    foo=fread(&cut_global,sizeof(double),1,fp);
    foo=fread(&sR6,sizeof(double),1,fp);
    foo=fread(&s8,sizeof(double),1,fp);
    foo=fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&sR6,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&s8,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */
/*
double PairD3::single(int i, int j, int itype, int jtype, double rsq,
			 double factor_coul, double factor_lj,
			 double &fforce)
{
  double fdamp6, fdamp6_deriv, K6, fdamp8, fdamp8_deriv, K8;
  double r,r2inv,r6inv,r8inv,r10inv, phi;

  r2inv = 1.0/rsq;
  r6inv = r2inv*r2inv*r2inv;
  r8inv = r6inv*r2inv;
  r10inv = r8inv*r2inv;
  r = sqrt(rsq);

  K6 = 6.0*pow(r/rr6[itype][jtype],-alpha6);
  K8 = 6.0*pow(r/rr8[itype][jtype],-alpha8);
  fdamp6 = 1.0 / (1.0 + K6);
  fdamp8 = 1.0 / (1.0 + K8);

  fdamp6_deriv = fdamp6*fdamp6*K6*alpha6; // derivative of f_damping_6, times r
  fdamp8_deriv = fdamp8*fdamp8*K8*alpha8; // derivative of f_damping_8, times r

  fforce = s6*c6[itype][jtype]*(fdamp6_deriv-fdamp6)*r8inv;
  fforce += s8*c8[itype][jtype]*(fdamp8_deriv-fdamp8)*r10inv;
  fforce *= factor_lj;

  phi = -s6*c6[itype][jtype]*fdamp6*r6inv - s8*c8[itype][jtype]*fdamp8*r8inv;
  return factor_lj*phi;
}
*/
/* ----------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */


int PairD3::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = CN[i];
  return m;
}

void PairD3::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    CN[j] += buf[m++];
  }
}

/*
int PairD3::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = CN[j];
  }
  return m;
}

void PairD3::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) CN[i] = buf[m++];
}
*/
////////////////////////////////////////////////////////////////////////////
// New member functions of PairD3
////////////////////////////////////////////////////////////////////////////
void PairD3::init_vars() {
    setup_completed = 0;

    alpha6 = 14.0;
    alpha8 = 16.0;
    alpha9 = 16.0;
    s6 = 1.0;
    sR8 = 1.0;
    sR9 = 4./3.;

    k1=16.0;
    //k2=4./3.; // not used: it's contained in rcov
    k3=4.0;

}



// 2D interpolation of Coordination Number between references of each type in the pair
double PairD3::interpolate_c6(int i, int j, int itype, int jtype)
{
  double CNi=CN[i], CNj=CN[j], Di,Dj, dist2, weight, z=0., w=0.;
  double mindist2=999999.9;
  int imindist2=-1;
  if(debugging_mode >= 2){
    if(screen) fprintf(screen,
      "PairD3: interpolate_c6: i=%d \t j=%d \t itype=%d \t jtype=%d \t CNi=%f \t CNj=%f\n",
      i,j,itype,jtype,CNi,CNj
    );
    if(logfile) fprintf(logfile,
      "PairD3: interpolate_c6: i=%d \t j=%d \t itype=%d \t jtype=%d \t CNi=%f \t CNj=%f\n",
      i,j,itype,jtype,CNi,CNj
    );
  }
  // NOTE: notation itype=0,..,n-1
  for (int iref = 0; iref < num_ref[itype-1][jtype-1]; iref++)
  {
    Di = CNi-ref_c6[itype-1][jtype-1][iref].CNa;
    Dj = CNj-ref_c6[itype-1][jtype-1][iref].CNb;
    dist2 = Di*Di + Dj*Dj;
    weight = exp( -k3*dist2 );
    w += weight;
    z += weight*ref_c6[itype-1][jtype-1][iref].c6;
    if(dist2<mindist2) {
      mindist2 = dist2;
      imindist2 = iref;
    }

    if(debugging_mode >= 2){
      if(screen) fprintf(screen,
        "PairD3: interpolate_c6:   iref=%d \t Di=%f \t Dj=%f \t weight=%g\n",
        iref,Di,Dj,weight
      );
      if(logfile) fprintf(logfile,
        "PairD3: interpolate_c6:   iref=%d \t Di=%f \t Dj=%f \t weight=%g\n",
        iref,Di,Dj,weight
      );
    }
  }

  if(imindist2==-1) error->all(FLERR,"Coordination Number is too distant from reference!");

  /*
    if(w==0.) {
      if(screen) fprintf(screen,
        "PairD3: interpolate_c6: i=%d \t j=%d \t CNi=%f \t CNj=%f \t z=%f \t w=%f\n",
        i,j,CNi,CNj,z,w
      );
      if(logfile) fprintf(logfile,
        "PairD3: interpolate_c6: i=%d \t j=%d \t CNi=%f \t CNj=%f \t z=%f \t w=%f\n",
        i,j,CNi,CNj,z,w
      );
      error->all(FLERR,"Coordination Number is too far from reference data");
    }
  */

  if(w>1e-6) return z/w;
  else       {
    if(debugging_mode >= 2){
      if(screen) fprintf(screen,
        "PairD3: interpolate_c6:   imindist2=%d \t c6=%f\n",
        imindist2,ref_c6[itype-1][jtype-1][imindist2].c6
      );
      if(logfile) fprintf(logfile,
        "PairD3: interpolate_c6:   imindist2=%d \t c6=%f\n",
        imindist2,ref_c6[itype-1][jtype-1][imindist2].c6
      );
    }
    return ref_c6[itype-1][jtype-1][imindist2].c6; // take the closest one
  }
}

int PairD3::get_pair_index(int i, int j, int min, int max)
{
  int ii,jj, k=0;
  if(i<min || i>max) error->all(FLERR,"Invalid value for i in pair");
  if(j<min || j>max) error->all(FLERR,"Invalid value for j in pair");
  for(ii=min;ii<=max;ii++)
    for(jj=min;jj<=ii;jj++)
    {
      k++;
      if(ii==i && jj==j) return k;
    }
  error->all(FLERR,"This loop should have been broken");
  return -1;
}

void PairD3::twobody(int i, int j, int itype, int jtype,
                    double rsq, double &fpair, int eflag, double &evdwl)
{
  double r,r2inv,r6inv,r8inv,r10inv;
  double c6,c8, fdamp6,fdamp6_deriv,k6, fdamp8,fdamp8_deriv,k8;
  r2inv = 1.0/rsq;
  r6inv = r2inv*r2inv*r2inv;
  r8inv = r6inv*r2inv;
  r10inv = r8inv*r2inv;

  // Calcola la correzione di van der Waals D3
  r = sqrt(rsq);

  k6 = 6.0*pow((sR6*r0[itype][jtype])/r, alpha6);
  k8 = 6.0*pow((sR8*r0[itype][jtype])/r, alpha8);
  fdamp6 = 1.0 / (1.0 + k6);
  fdamp8 = 1.0 / (1.0 + k8);

  fdamp6_deriv = fdamp6*fdamp6*k6*alpha6; // derivative of f_damping_6, times r
  fdamp8_deriv = fdamp8*fdamp8*k8*alpha8; // derivative of f_damping_8, times r

  c6 = interpolate_c6(i, j, itype, jtype);
  c8 = c6*Q[itype][jtype];

  fpair  = s6*c6*(fdamp6_deriv-fdamp6)*r8inv;
  fpair += s8*c8*(fdamp8_deriv-fdamp8)*r10inv;
  //fpair *= factor_lj;

  if (eflag) {
    evdwl = -s6*c6*fdamp6*r6inv - s8*c8*fdamp8*r8inv;
    //evdwl *= factor_lj;

    //if(screen) fprintf(screen, "PairD3: c6 = %f c8 = %f \t evdwl = %f\n", c6, c8, evdwl);
    //if(logfile) fprintf(logfile, "PairD3: c6 = %f c8 = %f \t evdwl = %f\n", c6, c8, evdwl);
  }
  return;
}

void PairD3::threebody(int itype, int jtype, int ktype,
                       double c9,
                       double rsq1, double rsq2, double rsq3,
                       double *delr1, double *delr2, double *delr3,
                       double *fj, double *fk, int eflag, double &evdwl)
{
  double r123, rinv123, tmp;
  double r1, r2, r3;
  double cos12,cos13,cos23;
  double rgeom, r0geom, k9, fdamp9;

  // Convention:
  // delr1 = r_j - r_i
  // delr2 = r_k - r_i
  // delr3 = r_k - r_j

  r1 = sqrt(rsq1);
  r2 = sqrt(rsq2);
  r3 = sqrt(rsq3);

  r123=r1*r2*r3;
  rinv123 = 1.0/r123;

  // r_ji * r_ki
  cos12 =  (delr1[0]*delr2[0] + delr1[1]*delr2[1] + delr1[2]*delr2[2]) / (r1*r2);
  // r_ji * r_jk
  cos13 = -(delr1[0]*delr3[0] + delr1[1]*delr3[1] + delr1[2]*delr3[2]) / (r1*r3);
  // r_ki * r_kj
  cos23 =  (delr2[0]*delr3[0] + delr2[1]*delr3[1] + delr2[2]*delr3[2]) / (r2*r3);

  // damping function
  rgeom = pow(r123, 1./3.);
  r0geom = pow(r0[itype][jtype]*r0[itype][ktype]*r0[ktype][jtype], 1./3.);
  k9 = 6.0*pow((sR9*r0geom)/rgeom, alpha9);
  fdamp9 = 1.0 / (1.0 + k9);

  // DA COMPLETARE!
  fj[0] = 0.;
  fj[1] = 0.;
  fj[2] = 0.;

  fk[0] = 0.;
  fk[1] = 0.;
  fk[2] = 0.;

/* // SW reference:
  fj[0] = delr1[0]*(frad1+csfac1)-delr2[0]*facang12;
  fj[1] = delr1[1]*(frad1+csfac1)-delr2[1]*facang12;
  fj[2] = delr1[2]*(frad1+csfac1)-delr2[2]*facang12;

  csfac2 = rinvsq2*csfacang;

  fk[0] = delr2[0]*(frad2+csfac2)-delr1[0]*facang12;
  fk[1] = delr2[1]*(frad2+csfac2)-delr1[1]*facang12;
  fk[2] = delr2[2]*(frad2+csfac2)-delr1[2]*facang12;
*/
  if (eflag) evdwl = - fdamp9 * c9 * (3.0*cos12*cos13*cos23 + 1.0) * CUBE(rinv123);

  if(debugging_mode >= 1 && comm->me==0){
    if(screen) fprintf(screen,"PairD3: evdwl3 = %g\n",- fdamp9 * c9 * (3.0*cos12*cos13*cos23 + 1.0) * CUBE(rinv123));
    if(logfile) fprintf(logfile,"PairD3: evdwl3 = %g\n",- fdamp9 * c9 * (3.0*cos12*cos13*cos23 + 1.0) * CUBE(rinv123));
  }

  return;
}
