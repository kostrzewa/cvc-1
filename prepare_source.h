#ifndef _PREPARE_SOURCE_H
#define _PREPARE_SOURCE_H

#include "types.h"

namespace cvc {

int prepare_volume_source ( double * const s, unsigned int const V);

int init_eo_spincolor_pointsource_propagator(double *s_even, double *s_odd, int global_source_coords[4], int isc, double*gauge_field, int sign, int have_source, double *work0);

int fini_eo_propagator(double *p_even, double *p_odd, double *r_even,   double *r_odd , double*gauge_field, int sign, double *work0);

int init_eo_sequential_source(double *s_even, double *s_odd, double *p_even, double *p_odd, int tseq, double*gauge_field, int sign, int pseq[3], int gseq, double *work0);


int init_clover_eo_spincolor_pointsource_propagator(double *s_even, double *s_odd, int global_source_coords[4], int isc, double *gauge_field, double*mzzinv, int have_source, double *work0);

int fini_clover_eo_propagator(double *p_even, double *p_odd, double *r_even, double *r_odd , double*gauge_field, double*mzzinv, double *work0);

int init_clover_eo_sequential_source ( double * const s_even, double * const s_odd, double * const p_even, double * const p_odd, int const tseq, double * const gauge_field, double * const mzzinv, int const pseq[3], int const gseq, double * const work0);


int check_vvdagger_locality (double** V, int numV, int gcoords[4], char*tag, double **sp);

int init_sequential_source ( double * const s, double * const p, int const tseq, int const pseq[3], int const gseq);

int init_coherent_sequential_source(double *s, double **p, int tseq, int ncoh, int pseq[3], int gseq);

int init_timeslice_source_oet ( double ** const s, int const tsrc, int * const momentum, int const spin_dilution, int const color_dilution, int const init);
  
  /**
   * @brief Generate oet timeslice source
   * Takes a random spinor, restricts it to a single time slice, multiplies
   * with a gamma structure of choice and the momentum phase, if any. 
   * @param s spinor, assumed to have 24*VOLUME elements
   * @param ran spinor, assumed to contain random elements of some type
   *            in all 24*VOLUME elements
   * @param gamma_id gamma structure to be used at the source
   * @param t_src source time slice
   * @param momentum
   *
   * @return 0 on success, -1 on memory failure 
   */
int prepare_gamma_timeslice_oet(double * const s,
      double const * const ran,
      int const gamma_id,
      int const t_src,
      const mom_t & momentum);
  
  int prepare_gamma_timeslice_oet_ext_momentum(
      double * const s,
      double const * const ran,
      int const gamma_id,
      int const t_src,
      ::cvc::complex const * const phases,
      const mom_t & momentum);

int init_sequential_fht_source ( double * const s, double * const p, int const pseq[3], int const gseq);

}  /* end of namespace cvc */

#endif
