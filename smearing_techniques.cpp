/******************************************
 * smearing_techniques.cpp
 *
 * originally from
 *   smearing_techniques.cc
 *   Author: Marc Wagner
 *   Date: September 2007
 *
 * February 2010
 * ported to cvc_parallel, parallelized
 *
 * Fri Dec  9 09:28:18 CET 2016
 * ported to cvc
 ******************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_MPI
#include <mpi.h>  
#endif
#include "cvc_complex.h"
#include "cvc_linalg.h"
#include "global.h"
#include "mpi_init.h"
#include "cvc_geometry.h"
#include "cvc_utils.h"
#include "smearing_techniques.h"
#include "debug_printf.hpp"
#include "types.h"
#include "loop_tools.h"

namespace cvc {

/************************************************************
 * Performs a number of APE smearing steps
 *
 ************************************************************/
int APE_Smearing(double * const smeared_gauge_field, double const APE_smearing_alpha, int const APE_smearing_niter) {
  const unsigned int gf_items = 72*VOLUME;
  const size_t gf_bytes = gf_items * sizeof(double);

  int iter;
  double *smeared_gauge_field_old = NULL;
  alloc_gauge_field(&smeared_gauge_field_old, VOLUMEPLUSRAND);

#pragma omp parallel
  {
    double M1[18], M2[18];
    int index;
    int index_mx_1, index_mx_2, index_mx_3;
    int index_px_1, index_px_2, index_px_3;
    int index_my_1, index_my_2, index_my_3;
    int index_py_1, index_py_2, index_py_3;
    int index_mz_1, index_mz_2, index_mz_3;
    int index_pz_1, index_pz_2, index_pz_3;
    double *U = NULL;
  
    for(int iter=0; iter<APE_smearing_niter; iter++) {
    
#pragma omp single
      {
        memcpy((void*)smeared_gauge_field_old, (void*)smeared_gauge_field, gf_bytes);
#ifdef HAVE_MPI
        xchange_gauge_field(smeared_gauge_field_old);
#endif
      } // there is an implicit barrier here and there's an implicit barrier at the
        // end of the volume loop below, so we should not have any synchro problems

      FOR_IN_PARALLEL(idx, 0, VOLUME){
  
        /************************
         * Links in x-direction.
         ************************/
        index = _GGI(idx, 1);
  
        index_my_1 = _GGI(g_idn[idx][2], 2);
        index_my_2 = _GGI(g_idn[idx][2], 1);
        index_my_3 = _GGI(g_idn[g_iup[idx][1]][2], 2);
  
        index_py_1 = _GGI(idx, 2);
        index_py_2 = _GGI(g_iup[idx][2], 1);
        index_py_3 = _GGI(g_iup[idx][1], 2);
  
        index_mz_1 = _GGI(g_idn[idx][3], 3);
        index_mz_2 = _GGI(g_idn[idx][3], 1);
        index_mz_3 = _GGI(g_idn[g_iup[idx][1]][3], 3);
  
        index_pz_1 = _GGI(idx, 3);
        index_pz_2 = _GGI(g_iup[idx][3], 1);
        index_pz_3 = _GGI(g_iup[idx][1], 3);
  
  
        U = smeared_gauge_field + index;
        _cm_eq_zero(U);
  
        /* negative y-direction */
        _cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_my_2, smeared_gauge_field_old + index_my_3);
  
        _cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_my_1, M1);
        _cm_pl_eq_cm(U, M2);
  
        /* positive y-direction */
        _cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_py_2, smeared_gauge_field_old + index_py_3);
  
        _cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_py_1, M1);
        _cm_pl_eq_cm(U, M2);
  
        /* negative z-direction */
        _cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_mz_2, smeared_gauge_field_old + index_mz_3);
  
        _cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_mz_1, M1);
        _cm_pl_eq_cm(U, M2);
  
        /* positive z-direction */
        _cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_pz_2, smeared_gauge_field_old + index_pz_3);
  
        _cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_pz_1, M1);
        _cm_pl_eq_cm(U, M2);
  
        _cm_ti_eq_re(U, APE_smearing_alpha);
  
        /* center */
        _cm_pl_eq_cm(U, smeared_gauge_field_old + index);
  
        /* Projection to SU(3). */
        cm_proj(U);
  
  
        /***********************
         * Links in y-direction.
         ***********************/
  
        index = _GGI(idx, 2);
  
        index_mx_1 = _GGI(g_idn[idx][1], 1);
        index_mx_2 = _GGI(g_idn[idx][1], 2);
        index_mx_3 = _GGI(g_idn[g_iup[idx][2]][1], 1);
  
        index_px_1 = _GGI(idx, 1);
        index_px_2 = _GGI(g_iup[idx][1], 2);
        index_px_3 = _GGI(g_iup[idx][2], 1);
  
        index_mz_1 = _GGI(g_idn[idx][3], 3);
        index_mz_2 = _GGI(g_idn[idx][3], 2);
        index_mz_3 = _GGI(g_idn[g_iup[idx][2]][3], 3);
  
        index_pz_1 = _GGI(idx, 3);
        index_pz_2 = _GGI(g_iup[idx][3], 2);
        index_pz_3 = _GGI(g_iup[idx][2], 3);
  
        U = smeared_gauge_field + index;
        _cm_eq_zero(U);
  
        /* negative x-direction */
        _cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_mx_2, smeared_gauge_field_old + index_mx_3);
        _cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_mx_1, M1);
        _cm_pl_eq_cm(U, M2);
  
        /* positive x-direction */
        _cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_px_2, smeared_gauge_field_old + index_px_3);
        _cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_px_1, M1);
        _cm_pl_eq_cm(U, M2);
  
        /* negative z-direction */
        _cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_mz_2, smeared_gauge_field_old + index_mz_3);
        _cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_mz_1, M1);
        _cm_pl_eq_cm(U, M2);
  
        /* positive z-direction */
        _cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_pz_2, smeared_gauge_field_old + index_pz_3);
        _cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_pz_1, M1);
        _cm_pl_eq_cm(U, M2);
  
        _cm_ti_eq_re(U, APE_smearing_alpha);
  
        /* center */
        _cm_pl_eq_cm(U, smeared_gauge_field_old + index);
  
        /* Projection to SU(3). */
        cm_proj(U);
  
        /**************************
         * Links in z-direction.
         **************************/
  
        index = _GGI(idx, 3);
  
        index_mx_1 = _GGI(g_idn[idx][1], 1);
        index_mx_2 = _GGI(g_idn[idx][1], 3);
        index_mx_3 = _GGI(g_idn[g_iup[idx][3]][1], 1);
  
        index_px_1 = _GGI(idx, 1);
        index_px_2 = _GGI(g_iup[idx][1], 3);
        index_px_3 = _GGI(g_iup[idx][3], 1);
  
        index_my_1 = _GGI(g_idn[idx][2], 2);
        index_my_2 = _GGI(g_idn[idx][2], 3);
        index_my_3 = _GGI(g_idn[g_iup[idx][3]][2], 2);
  
        index_py_1 = _GGI(idx, 2);
        index_py_2 = _GGI(g_iup[idx][2], 3);
        index_py_3 = _GGI(g_iup[idx][3], 2);
  
        U = smeared_gauge_field + index;
        _cm_eq_zero(U);
  
        /* negative x-direction */
        _cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_mx_2, smeared_gauge_field_old + index_mx_3);
        _cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_mx_1, M1);
        _cm_pl_eq_cm(U, M2);
  
        /* positive x-direction */
        _cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_px_2, smeared_gauge_field_old + index_px_3);
        _cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_px_1, M1);
        _cm_pl_eq_cm(U, M2);
  
        /* negative y-direction */
        _cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_my_2, smeared_gauge_field_old + index_my_3);
        _cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_my_1, M1);
        _cm_pl_eq_cm(U, M2);
  
        /* positive y-direction */
        _cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_py_2, smeared_gauge_field_old + index_py_3);
        _cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_py_1, M1);
        _cm_pl_eq_cm(U, M2);
  
        _cm_ti_eq_re(U, APE_smearing_alpha);
  
        /* center */
        _cm_pl_eq_cm(U, smeared_gauge_field_old + index);
  
        /* Projection to SU(3). */
        cm_proj(U);
      }  /* end of loop on ix */
    }  /* end of loop on number of smearing steps */
#ifdef HAVE_OPENMP
  }  /* end of parallel region */
#endif

  free(smeared_gauge_field_old);
#ifdef HAVE_MPI
  xchange_gauge_field(smeared_gauge_field);
#endif
  return(0);
}  /* end of APE_Smearing */


/*****************************************************
 *
 * Performs a number of Jacobi smearing steps
 *
 * smeared_gauge_field  = gauge field used for smearing (in)
 * psi = quark spinor (in/out)
 * N, kappa = Jacobi smearing parameters (in)
 *
 *****************************************************/
int Jacobi_Smearing(double * smeared_gauge_field, double * const psi, int const N, double const kappa) {
  const size_t sf_items = _GSI(VOLUME);
  const size_t sf_bytes = sf_items * sizeof(double);
  const double norm = 1.0 / (1.0 + 6.0*kappa);

  if ( N == 0 || kappa == 0.0 ) {
    if ( g_cart_id == 0 ) fprintf(stdout, "# [Jacobi_Smearing] Warning, nothing to smear\n");
    return ( 0 );
  }

  if ( smeared_gauge_field == NULL || psi == NULL ) {
    fprintf(stderr, "[Jacobi_Smearing] Error, input / output fields are NULL \n");
    return(4);
  }


  int i1;
  double *psi_old = (double*)malloc( _GSI(VOLUME+RAND)*sizeof(double));
  if(psi_old == NULL) {
    fprintf(stderr, "[Jacobi_Smearing] Error from malloc\n");
    return(3);
  }

#ifdef HAVE_MPI
  xchange_gauge_field (smeared_gauge_field);
#endif

  /* loop on smearing steps */
  for(i1 = 0; i1 < N; i1++) {

    memcpy((void*)psi_old, (void*)psi, sf_bytes);
#ifdef HAVE_MPI
    xchange_field(psi_old);
#endif

#ifdef HAVE_OPENMP
#pragma omp parallel
{
#endif

  unsigned int idx, idy;
  unsigned int index_s, index_s_mx, index_s_px, index_s_my, index_s_py, index_s_mz, index_s_pz, index_g_mx;
  unsigned int index_g_px, index_g_my, index_g_py, index_g_mz, index_g_pz; 
  double *s=NULL, spinor[24];

#ifdef HAVE_OPENMP
#pragma omp for
#endif
    for(idx = 0; idx < VOLUME; idx++) {
      index_s = _GSI(idx);

      index_s_mx = _GSI(g_idn[idx][1]);
      index_s_px = _GSI(g_iup[idx][1]);
      index_s_my = _GSI(g_idn[idx][2]);
      index_s_py = _GSI(g_iup[idx][2]);
      index_s_mz = _GSI(g_idn[idx][3]);
      index_s_pz = _GSI(g_iup[idx][3]);

      idy = idx;
      index_g_mx = _GGI(g_idn[idy][1], 1);
      index_g_px = _GGI(idy, 1);
      index_g_my = _GGI(g_idn[idy][2], 2);
      index_g_py = _GGI(idy, 2);
      index_g_mz = _GGI(g_idn[idy][3], 3);
      index_g_pz = _GGI(idy, 3);

      s = psi + _GSI(idy);
      _fv_eq_zero(s);

      /* negative x-direction */
      _fv_eq_cm_dag_ti_fv(spinor, smeared_gauge_field + index_g_mx, psi_old + index_s_mx);
      _fv_pl_eq_fv(s, spinor);

      /* positive x-direction */
      _fv_eq_cm_ti_fv(spinor, smeared_gauge_field + index_g_px, psi_old + index_s_px);
      _fv_pl_eq_fv(s, spinor);

      /* negative y-direction */
      _fv_eq_cm_dag_ti_fv(spinor, smeared_gauge_field + index_g_my, psi_old + index_s_my);
      _fv_pl_eq_fv(s, spinor);

      /* positive y-direction */
      _fv_eq_cm_ti_fv(spinor, smeared_gauge_field + index_g_py, psi_old + index_s_py);
      _fv_pl_eq_fv(s, spinor);

      /* negative z-direction */
      _fv_eq_cm_dag_ti_fv(spinor, smeared_gauge_field + index_g_mz, psi_old + index_s_mz);
      _fv_pl_eq_fv(s, spinor);

      /* positive z-direction */
      _fv_eq_cm_ti_fv(spinor, smeared_gauge_field + index_g_pz, psi_old + index_s_pz);
      _fv_pl_eq_fv(s, spinor);


      /* Put everything together; normalization. */
      _fv_ti_eq_re(s, kappa);
      _fv_pl_eq_fv(s, psi_old + index_s);
      _fv_ti_eq_re(s, norm);
    }  /* end of loop on ix */
#ifdef HAVE_OPENMP
}  /* end of parallel region */
#endif
  }  /* end of loop on smearing steps */

  free(psi_old);
  return(0);
}  /* end of Jacobi_Smearing */


int Momentum_Smearing(double const * const smeared_gauge_field,
    mom_t const & momentum, double const mom_scale_factor,
    double * const psi, int const N, double const kappa)
{
  size_t const gf_bytes = 72*VOLUMEPLUSRAND*sizeof(double);
  double * smeared_gauge_field_with_phase = (double*)malloc(gf_bytes);
  if( (void*)smeared_gauge_field_with_phase == NULL ){
    error_printf("malloc failure in [MomentumSmearing]");
    EXIT(CVC_EXIT_MALLOC_FAILURE);
  }

  const double TWO_MPI = 2.0 * M_PI;
  ::cvc::complex smear_phase[3];
  smear_phase[0].re = cos( mom_scale_factor * TWO_MPI * momentum.x / (double)LX_global );
  smear_phase[0].im = sin( mom_scale_factor * TWO_MPI * momentum.x / (double)LX_global );
  smear_phase[1].re = cos( mom_scale_factor * TWO_MPI * momentum.y / (double)LY_global );
  smear_phase[1].im = sin( mom_scale_factor * TWO_MPI * momentum.y / (double)LY_global );
  smear_phase[2].re = cos( mom_scale_factor * TWO_MPI * momentum.z / (double)LZ_global );
  smear_phase[2].im = sin( mom_scale_factor * TWO_MPI * momentum.z / (double)LZ_global );

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {
    // multiply gauge field by momentum smearing phase 
    FOR_IN_PARALLEL(ix, 0, VOLUME)
    {
      _cm_eq_cm_ti_co(smeared_gauge_field_with_phase + _GGI(ix,1),
                      smeared_gauge_field + _GGI(ix,1),
                      &smear_phase[0]);
      _cm_eq_cm_ti_co(smeared_gauge_field_with_phase + _GGI(ix,2),
                      smeared_gauge_field + _GGI(ix,2),
                      &smear_phase[1]);
      _cm_eq_cm_ti_co(smeared_gauge_field_with_phase + _GGI(ix,3),
                      smeared_gauge_field + _GGI(ix,3),
                      &smear_phase[2]);
    }
  } // omp parallel closing brace
  
  Jacobi_Smearing(smeared_gauge_field_with_phase, psi, N, kappa);
  free(smeared_gauge_field_with_phase);
  return 0;
}

}  /* end of namespace cvc */
