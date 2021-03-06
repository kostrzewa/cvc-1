#ifndef _PROPAGATOR_IO_H
#define _PROPAGATOR_IO_H


#include "dml.h"

#ifdef __cplusplus
extern "C"
{
#endif

#include "lime.h"

#ifdef __cplusplus
}
#endif


namespace cvc {
int write_propagator(double * const s, char const * const filename, const int append, const int prec);

int write_propagator_format(char const * const filename, const int prec, const int no_flavours);

int write_propagator_type(const int type, char const * const filename);

int get_propagator_type(char const * const filename);

#ifdef HAVE_LIBLEMON
int write_lemon_spinor(double * const s, char const * const filename, const int append, const int prec);
#else
int write_lime_spinor(double * const s, char const * const filename, const int append, const int prec);
#endif

int write_checksum(char const * const filename, DML_Checksum *checksum); 

/* int write_binary_spinor_data(double * const s, LimeWriter * limewriter,
  const int prec, DML_Checksum * ans); */

#ifdef HAVE_LIBLEMON
int read_lime_spinor(double * const s, char * filename, const int position);
#else
int read_lime_spinor(double * const s, char const * const filename, const int position);
#endif
int read_cmi(double *v, char const * const filename);
int write_binary_spinor_data_timeslice(double * const s, LimeWriter * limewriter,
  const int prec, int timeslice, DML_Checksum * ans);


int write_lime_spinor_timeslice(double * const s, char const * const filename,
   const int prec, int timeslice, DML_Checksum *checksum);

int prepare_propagator(int timeslice, int iread, int is_mms, int no_mass, double sign, double mass, int isave, double *work, int pos, double*gauge_field);

int prepare_propagator2(int *source_coords, int iread, int sign, void *work, int pos, int format, size_t prec_out);

int rotate_propagator_ETMC_UKQCD (double *spinor, long unsigned int V);

int read_lime_spinor_single(float * const s, char const * const filename, const int position);
int read_binary_spinor_data_timeslice(double * const s, int timeslice, LimeReader * limereader, const int prec, DML_Checksum *ans);
int read_lime_spinor_timeslice(double * const s, int timeslice, char const * const filename, const int position, DML_Checksum*checksum);

int write_source_type(const int type, char const * const filename);

int deflator_read_evecs ( double *const evecs, unsigned int const nev, char * const fileformat, char const * const filename , int const prec);

}  // end of namespace cvc
#endif
