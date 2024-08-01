/* -*- c++ -*- ----------------------------------------------------------
   Flavio Giuliani 1aug2024 flavio.giuliani@uniroma1.it
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(hybrid/scaled,PairHybridScaled);

#else

#ifndef LMP_PAIR_HYBRID_SCALED_H
#define LMP_PAIR_HYBRID_SCALED_H

#include "pair_hybrid.h"

#include <string>
#include <vector>

namespace LAMMPS_NS {

class PairHybridScaled : public PairHybrid {
 public:
  PairHybridScaled(class LAMMPS *);
  ~PairHybridScaled() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;

  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;

  void init_svector(); // override;
  void copy_svector(int, int); // override;

 protected:
  double **fsum, **tsum;
  double *scaleval;
  int *scaleidx;
  std::vector<std::string> scalevars;
  int nmaxfsum;
};

}    // namespace LAMMPS_NS

#endif
#endif
