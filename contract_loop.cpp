/****************************************************
 * contract_loop.cpp
 *
 * PURPOSE:
 * DONE:
 * TODO:
 ****************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef HAVE_MPI
#  include <mpi.h>
#endif
#ifdef HAVE_OPENMP
#  include <omp.h>
#endif
#include <getopt.h>

#ifdef HAVE_LHPC_AFF
#include "lhpc-aff.h"
#endif

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

#include "cvc_complex.h"
#include "iblas.h"
#include "ilinalg.h"
#include "cvc_linalg.h"
#include "global.h"
#include "cvc_geometry.h"
#include "cvc_utils.h"
#include "mpi_init.h"
#include "table_init_d.h"
#include "table_init_z.h"
#include "project.h"
#include "Q_phi.h"
#include "Q_clover_phi.h"
#include "contract_cvc_tensor.h"
#include "scalar_products.h"
#include "contract_loop_inline.h"

#define MAX_SUBGROUP_NUMBER 20

namespace cvc {

/***********************************************************
 * local loop contractions
 *
 * contraction with spin dilution
 *
 * OUT: loop must hold 
 *   (a) T x momentum_number x [4 x 4] complex values if momentum projection is done
 *       ( momentum_number > 0 && momentum_list != NULL )
 *   (b) T x VOL3 x [4x4] complex values, if momentum projetion is not done
 *
 * IN: source = stochastic source, full VOLUME spinor field
 * IN: prop   = stochastic propagator, full VOLUME spinor field
 * IN: momentum_number = number of entries in momentum_list
 * IN: momentum_list = list of int[3] with momenta
 *
 ***********************************************************/
int contract_local_loop_stochastic ( double *** const loop, double * const source, double * const prop, int const momentum_number, int ( * const momentum_list)[3] ) {

  double * loop_x = init_1level_dtable ( 32 * VOLUME );
  if ( loop_x == NULL ) {
    fprintf ( stderr, "[contract_local_loop_stochastic] Error from init_1level_dtable %s %d\n", __FILE__, __LINE__ );
    return ( 1 );
  }

  if ( 
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
  for ( int ix = 0; ix < VOLUME; ix++ ) {
      unsigned int const offset = _GSI ( ix );
      double * const source_ = source + offset;
      double * const prop_   = prop   + offset;
      _contract_loop_x_spin_diluted ( loop_x+16*ix, source_ , prop_ );
    }  /* end of loop on volume */
  }  /* end of loop on timeslices */

  if ( momentum_number > 0 && momentum_list != NULL ) {
    /***********************************************************
     * in that case we do a momentum projection
     ***********************************************************/
  
    for ( int it = 0; it < T; it++ ) {

      unsigned int const offset = 32 * it * VOL3;

      /***********************************************************
       * momentum projection for a timeslice
       * VOL3 x 4x4 -> momentum_number x 4x4
       *
       * NOTE: momentum_proction2 contains the MPI-reduction in
       * spatial dimensions
       ***********************************************************/
      int exitstatus = momentum_projection2 ( loop_x + offset, loop[it][0], 16, momentum_number, momentum_list, NULL );
      if ( exitstatus == 0 ) {
        fprintf ( stderr, "[contract_local_loop_stochastic] Error from momentum_projection2, status was %d %s %d\n", exitstatus, __FILE__, __LINE__ );
        return ( 2 );
      }

    }  /* end of loop on timeslices */
  } else { 
    /***********************************************************
     * in that case we just copy the results from loop_x
     ***********************************************************/
    memcpy ( loop[0][0], loop_x, 32*VOLUME*sizeof(double) );
  }  /* end of if momentum projection else */

  fini_1level_dtable ( &loop_x );

 return ( 0 );
}  /* end of contract_local_loop_stochastic */

/***************************************************************************/
/***************************************************************************/

#ifdef HAVE_HDF5

/***************************************************************************
 * write time-momentum-dependent loop to HDF5 file
 *
 * p runs slower than t, i.e. c_tp[momentum][time]
 *
 * NOTE: data it c_tp MUST have been MPI_Reduced already in 
 *       timeslice g_ts_comm;
 *       here we only MPI_Gather in time ray g_tr_comm
 ***************************************************************************/

int contract_loop_write_to_h5_file (double *** const loop, void * file, char*tag, int const momentum_number, int const nc, int const io_proc ) {

  if ( io_proc > 0 ) {

    double ratime, retime;

    double * zbuffer = NULL;

    char * filename = (char *)file;

    ratime = _GET_TIME;

    if ( io_proc == 2 ) {
      /***************************************************************************
       * recvbuf for MPI_Gather; io_proc 2 is receive process,
       * zbuffer is significant
       ***************************************************************************/
      zbuffer = init_1level_dtable ( 32 * T_global * momentum_number );
      if( zbuffer == NULL ) {
        fprintf(stderr, "[contract_loop_write_to_h5_file] Error from init_1level_dtable %s %d\n", __FILE__, __LINE__);
        return(1);
      }
    }  /* end of if io_proc == 2 else */

    /***************************************************************************
     * loops from contract_local_loop_stochastic are already in form
     * T x momentum_number x 16
     *
     * ERGO: we do not need to reorder
     ***************************************************************************/
    int mitems = momentum_number * 32 * T;
#ifdef HAVE_MPI
    /***************************************************************************
     * io_proc's 1 and 2 gather the data to g_tr_id = 0 into zbuffer
     ***************************************************************************/
#  if (defined PARALLELTX) || (defined PARALLELTXY) || (defined PARALLELTXYZ) 
    int exitstatus = MPI_Gather ( loop[0][0], mitems, MPI_DOUBLE, zbuffer, mitems, MPI_DOUBLE, 0, g_tr_comm);
    if(exitstatus != MPI_SUCCESS) {
      fprintf(stderr, "[contract_loop_write_to_h5_file] Error from MPI_Gather, status was %d %s %d\n", exitstatus, __FILE__, __LINE__);
      return(4);
    }
#  else
    int exitstatus = MPI_Gather ( loop[0][0], mitems, MPI_DOUBLE, zbuffer, mitems, MPI_DOUBLE, 0, g_cart_grid);
    if(exitstatus != MPI_SUCCESS) {
      fprintf(stderr, "[contract_loop_write_to_h5_file] Error from MPI_Gather, status was %d %s %d\n", exitstatus, __FILE__, __LINE__);
      return(5);
    }
#  endif

#else
    /***************************************************************************
     * just copy data into zbuffer
     * T = T_global
     ***************************************************************************/
    memcpy ( zbuffer, loop[0][0], mitems * sizeof(double) );
#endif

    /***************************************************************************
     * io_proc 2 is origin of Cartesian grid and does the write to disk
     ***************************************************************************/
    if(io_proc == 2) {
  
      /***************************************************************************
       * create or open file
       ***************************************************************************/

      hid_t   file_id;
      herr_t  status;

      struct stat fileStat;
      if ( stat( filename, &fileStat) < 0 ) {
        /* creat a new file */
  
        if ( g_verbose > 1 ) fprintf ( stdout, "# [contract_loop_write_to_h5_file] create new file\n" );
  
        unsigned flags = H5F_ACC_TRUNC; /* IN: File access flags. Allowable values are:
                                           H5F_ACC_TRUNC --- Truncate file, if it already exists, erasing all data previously stored in the file.
                                           H5F_ACC_EXCL  --- Fail if file already exists.
  
                                           H5F_ACC_TRUNC and H5F_ACC_EXCL are mutually exclusive; use exactly one.
                                           An additional flag, H5F_ACC_DEBUG, prints debug information.
                                           This flag can be combined with one of the above values using the bit-wise OR operator (`|'),
                                           but it is used only by HDF5 Library developers; it is neither tested nor supported for use in applications.  */
        hid_t fcpl_id = H5P_DEFAULT; /* IN: File creation property list identifier, used when modifying default file meta-data.
                                        Use H5P_DEFAULT to specify default file creation properties. */
        hid_t fapl_id = H5P_DEFAULT; /* IN: File access property list identifier. If parallel file access is desired,
                                        this is a collective call according to the communicator stored in the fapl_id.
                                        Use H5P_DEFAULT for default file access properties. */
  
        /*  hid_t H5Fcreate ( const char *name, unsigned flags, hid_t fcpl_id, hid_t fapl_id ) */
        file_id = H5Fcreate (         filename,          flags,       fcpl_id,       fapl_id );
  
      } else {
        /* open an existing file. */
        if ( g_verbose > 1 ) fprintf ( stdout, "# [contract_loop_write_to_h5_file] open existing file\n" );
  
        unsigned flags = H5F_ACC_RDWR;  /* IN: File access flags. Allowable values are:
                                           H5F_ACC_RDWR   --- Allow read and write access to file.
                                           H5F_ACC_RDONLY --- Allow read-only access to file.
  
                                           H5F_ACC_RDWR and H5F_ACC_RDONLY are mutually exclusive; use exactly one.
                                           An additional flag, H5F_ACC_DEBUG, prints debug information.
                                           This flag can be combined with one of the above values using the bit-wise OR operator (`|'),
                                           but it is used only by HDF5 Library developers; it is neither tested nor supported for use in applications. */
        hid_t fapl_id = H5P_DEFAULT;
        /*  hid_t H5Fopen ( const char *name, unsigned flags, hid_t fapl_id ) */
        file_id = H5Fopen (         filename,         flags,        fapl_id );
      }
  
      if ( g_verbose > 1 ) fprintf ( stdout, "# [contract_loop_write_to_h5_file] file_id = %ld\n", file_id );
  
      /***************************************************************************
       * H5 data space and data type
       ***************************************************************************/
      hid_t dtype_id = H5Tcopy( H5T_NATIVE_DOUBLE );
      status = H5Tset_order ( dtype_id, H5T_ORDER_LE );
      /* big_endian() ?  H5T_IEEE_F64BE : H5T_IEEE_F64LE; */
 
      /* shape of the output arary */ 
      hsize_t dims[3] = { T_global, momentum_number, 32 };
  
      /*
                 int rank                             IN: Number of dimensions of dataspace.
                 const hsize_t * current_dims         IN: Array specifying the size of each dimension.
                 const hsize_t * maximum_dims         IN: Array specifying the maximum size of each dimension.
                 hid_t H5Screate_simple( int rank, const hsize_t * current_dims, const hsize_t * maximum_dims )
       */
      hid_t space_id = H5Screate_simple(        3,                         dims,                          NULL);
  
      /***************************************************************************
       * some default settings for H5Dwrite
       ***************************************************************************/
      hid_t mem_type_id   = H5T_NATIVE_DOUBLE;
      hid_t mem_space_id  = H5S_ALL;
      hid_t file_space_id = H5S_ALL;
      hid_t xfer_plit_id  = H5P_DEFAULT;
      hid_t lcpl_id       = H5P_DEFAULT;
      hid_t dcpl_id       = H5P_DEFAULT;
      hid_t dapl_id       = H5P_DEFAULT;
      hid_t gcpl_id       = H5P_DEFAULT;
      hid_t gapl_id       = H5P_DEFAULT;
      /* size_t size_hint    = 0; */
  
      /***************************************************************************
       * create the target (sub-)group and all
       * groups in hierarchy above if they don't exist
       ***************************************************************************/
      hid_t grp_list[MAX_SUBGROUP_NUMBER];
      int grp_list_nmem = 0;
      char grp_name[400], grp_name_tmp[400];
      char * grp_ptr = NULL;
      char grp_sep[] = "/";
      strcpy ( grp_name, tag );
      strcpy ( grp_name_tmp, grp_name );
      if ( g_verbose > 1 ) fprintf ( stdout, "# [contract_loop_write_to_h5_file] full grp_name = %s\n", grp_name );
      grp_ptr = strtok ( grp_name_tmp, grp_sep );
  
      while ( grp_ptr != NULL ) {
        hid_t grp;
        hid_t loc_id = ( grp_list_nmem == 0 ) ? file_id : grp_list[grp_list_nmem-1];
        if ( g_verbose > 1 ) fprintf ( stdout, "# [contract_loop_write_to_h5_file] grp_ptr = %s\n", grp_ptr );
  
        grp = H5Gopen2( loc_id, grp_ptr, gapl_id );
        if ( grp < 0 ) {
          fprintf ( stderr, "[contract_loop_write_to_h5_file] Error from H5Gopen2 for group %s, status was %ld %s %d\n", grp_ptr, grp, __FILE__, __LINE__ );
          grp = H5Gcreate2 (       loc_id,         grp_ptr,       lcpl_id,       gcpl_id,       gapl_id );
          if ( grp < 0 ) {
            fprintf ( stderr, "[contract_loop_write_to_h5_file] Error from H5Gcreate2 for group %s, status was %ld %s %d\n", grp_ptr, grp, __FILE__, __LINE__ );
            return ( 6 );
          } else {
            if ( g_verbose > 1 ) fprintf ( stdout, "# [contract_loop_write_to_h5_file] created group %s %ld %s %d\n", grp_ptr, grp, __FILE__, __LINE__ );
          }
        } else {
          if ( g_verbose > 1 ) fprintf ( stdout, "# [contract_loop_write_to_h5_file] opened group %s %ld %s %d\n", grp_ptr, grp, __FILE__, __LINE__ );
        }
        grp_ptr = strtok(NULL, grp_sep );
  
        grp_list[grp_list_nmem] = grp;
        grp_list_nmem++;
      }  /* end of loop on sub-groups */
  
      /***************************************************************************
       * hdf5 id to write to
       * either file itself or current group identifier
       ***************************************************************************/
      hid_t loc_id = ( grp_list_nmem == 0 ) ? file_id : grp_list[grp_list_nmem - 1 ];

      /***************************************************************************
       * write data set
       ***************************************************************************/
      char name[] = "loop";
  
      /***************************************************************************
       * create a data set
       ***************************************************************************/
      /*
                   hid_t loc_id         IN: Location identifier
                   const char *name     IN: Dataset name
                   hid_t dtype_id       IN: Datatype identifier
                   hid_t space_id       IN: Dataspace identifier
                   hid_t lcpl_id        IN: Link creation property list
                   hid_t dcpl_id        IN: Dataset creation property list
                   hid_t dapl_id        IN: Dataset access property list
                   hid_t H5Dcreate2 ( hid_t loc_id, const char *name, hid_t dtype_id, hid_t space_id, hid_t lcpl_id, hid_t dcpl_id, hid_t dapl_id )
  
                   hid_t H5Dcreate ( hid_t loc_id, const char *name, hid_t dtype_id, hid_t space_id, hid_t lcpl_id, hid_t dcpl_id, hid_t dapl_id ) 
       */
      hid_t dataset_id = H5Dcreate (       loc_id,             name,       dtype_id,       space_id,       lcpl_id,       dcpl_id,       dapl_id );
  
      /***************************************************************************
       * write the current data set
       ***************************************************************************/
      /*
               hid_t dataset_id           IN: Identifier of the dataset to write to.
               hid_t mem_type_id          IN: Identifier of the memory datatype.
               hid_t mem_space_id         IN: Identifier of the memory dataspace.
               hid_t file_space_id        IN: Identifier of the dataset's dataspace in the file.
               hid_t xfer_plist_id        IN: Identifier of a transfer property list for this I/O operation.
               const void * buf           IN: Buffer with data to be written to the file.
        herr_t H5Dwrite ( hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id, const void * buf )
       */
      status = H5Dwrite (       dataset_id,       mem_type_id,       mem_space_id,       file_space_id,        xfer_plit_id,    zbuffer[0][0] );
  
      if( status < 0 ) {
        fprintf(stderr, "[contract_loop_write_to_h5_file] Error from H5Dwrite, status was %d %s %d\n", status, __FILE__, __LINE__);
        return(8);
      }
  
      /***************************************************************************
       * close the current data set
       ***************************************************************************/
      status = H5Dclose ( dataset_id );
      if( status < 0 ) {
        fprintf(stderr, "[contract_loop_write_to_h5_file] Error from H5Dclose, status was %d %s %d\n", status, __FILE__, __LINE__);
        return(9);
      }
  
      /***************************************************************************
       * close the data space
       ***************************************************************************/
      status = H5Sclose ( space_id );
      if( status < 0 ) {
        fprintf(stderr, "[contract_loop_write_to_h5_file] Error from H5Sclose, status was %d %s %d\n", status, __FILE__, __LINE__);
        return(10);
      }
  
      /***************************************************************************
       * close all (sub-)groups in reverse order
       ***************************************************************************/
      for ( int i = grp_list_nmem - 1; i>= 0; i-- ) {
        status = H5Gclose ( grp_list[i] );
        if( status < 0 ) {
          fprintf(stderr, "[contract_loop_write_to_h5_file] Error from H5Gclose, status was %d %s %d\n", status, __FILE__, __LINE__);
          return(11);
        } else {
          if ( g_verbose > 1 ) fprintf(stdout, "# [contract_loop_write_to_h5_file] closed group %ld %s %d\n", grp_list[i], __FILE__, __LINE__);
        }
      }
  
      /***************************************************************************
       * close the data type
       ***************************************************************************/
      status = H5Tclose ( dtype_id );
      if( status < 0 ) {
        fprintf(stderr, "[contract_loop_write_to_h5_file] Error from H5Tclose, status was %d %s %d\n", status, __FILE__, __LINE__);
        return(12);
      }
  
      /***************************************************************************
       * close the file
       ***************************************************************************/
      status = H5Fclose ( file_id );
      if( status < 0 ) {
        fprintf(stderr, "[contract_loop_write_to_h5_file] Error from H5Fclose, status was %d %s %d\n", status, __FILE__, __LINE__);
        return(13);
     } 
  
    }  /* if io_proc == 2 */

    fini_2level_dtable ( &zbuffer );

    retime = _GET_TIME;
    if ( io_proc == 2 ) fprintf(stdout, "# [contract_loop_write_to_h5_file] time for saving momentum space results = %e seconds\n", retime-ratime);

  }  /* end of of if io_proc > 0 */

  return(0);

}  /* end of contract_loop_write_to_h5_file */

#endif  /* of if defined HAVE_HDF5 */

/***********************************************************/
/***********************************************************/

}  /* end of namespace cvc */
