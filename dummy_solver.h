#ifndef _DUMMY_SOLVER_H
#define _DUMMY_SOLVER_H

#include "tmLQCD.h"

#if ( defined DUMMY_SOLVER ) && ( defined GPU_DIRECT_SOLVER )
#error exclusive XXX_SOLVER definitions
#endif

#if ( defined DUMMY_SOLVER ) || ( ! defined HAVE_TMLQCD_LIBWRAPPER ) 
#  define _TMLQCD_INVERT dummy_solver
namespace cvc {
// int dummy_solver ( double * const propagator, double * const source, const int op_id, int write_prop);
int dummy_solver ( double * const propagator, double * const source, int const op_id );
}  /* end of namespace cvc */
#elif ( defined GPU_DIRECT_SOLVER )
#  define _TMLQCD_INVERT(PROP, SOURCE, OP_ID) \
  invert_quda_direct((PROP), (SOURCE), (OP_ID));
#else
#  define _TMLQCD_INVERT(PROP, SOURCE, OP_ID) \
  tmLQCD_invert((PROP), (SOURCE), (OP_ID), 0);
#endif

#endif

