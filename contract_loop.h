#ifndef _CONTRACT_LOOP_H
#define _CONTRACT_LOOP_H
/****************************************************
 * contract_loop.h
 *
 * PURPOSE:
 * DONE:
 * TODO:
 ****************************************************/
namespace cvc {

int contract_local_loop_stochastic ( double *** const loop, double * const source, double * const prop, int const momentum_number, int ( * const momentum_list)[3] );

int contract_loop_write_to_h5_file (double *** const loop, void * file, char*tag, int const momentum_number, int const nc, int const io_proc );

int loop_read_from_h5_file (double *** const loop, void * file, char*tag, int const momentum_number, int const nc, int const io_proc );

int loop_get_momentum_list_from_h5_file ( int (*momentum_list)[3], void * file, int const momentum_number, int const io_proc );

int loop_write_momentum_list_to_h5_file ( int (*momentum_list)[3], void * file, int const momentum_number, int const io_proc );


}  /* end of namespace cvc */
#endif
