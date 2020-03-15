/****************************************************
 * loop_extract
 *
 * PURPOSE:
 * DONE:
 * TODO:
 ****************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
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

#define MAIN_PROGRAM

#include "cvc_complex.h"
#include "cvc_linalg.h"
#include "enums.hpp"
#include "cvc_global.h"
#include "cvc_geometry.h"
#include "cvc_utils.h"
#include "mpi_init.h"
#include "set_default.h"
#include "io.h"
#include "propagator_io.h"
#include "read_input_parser.h"
#include "contractions_io.h"
#include "Q_clover_phi.h"
#include "contract_cvc_tensor.h"
#include "prepare_source.h"
#include "prepare_propagator.h"
#include "project.h"
#include "table_init_z.h"
#include "table_init_d.h"
#include "dummy_solver.h"
#include "Q_phi.h"
#include "clover.h"
#include "contract_loop.h"
#include "ranlxd.h"

#define _OP_ID_UP 0
#define _OP_ID_DN 1
#define _OP_ID_ST 2

using namespace cvc;

void usage() {
  fprintf(stdout, "Code to analyse loop data\n");
  fprintf(stdout, "Usage:    [options]\n");
  fprintf(stdout, "Options:  -f <input filename> : input filename for cvc      [default cpff.input]\n");
  EXIT(0);
}

int main(int argc, char **argv) {
  
  const char outfile_prefix[] = "loop";

  /* const char fbwd_str[2][4] =  { "fwd", "bwd" }; */

  int conf_traj ;  /* For LegacyTraj */

  int c;
  int filename_set = 0;
  int exitstatus;
  int io_proc = -1;
  char filename[1200];

  struct timeval ta, tb;
  long unsigned int seconds, useconds;

  char output_filename[400];

  char data_tag[400];

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  while ((c = getopt(argc, argv, "h?f:")) != -1) {
    switch (c) {
    case 'f':
      strcpy(filename, optarg);
      filename_set=1;
      break;
    case 'h':
    case '?':
    default:
      usage();
      break;
    }
  }

  g_the_time = time(NULL);

  /* set the default values */
  if(filename_set==0) snprintf ( filename, 400, "%s.input", outfile_prefix );
  /* fprintf(stdout, "# [loop_extract] Reading input from file %s\n", filename); */
  read_input_parser(filename);

  /***************************************************************************
   * initialize MPI parameters for cvc
   ***************************************************************************/
  mpi_init(argc, argv);
#ifdef HAVE_MPI
  if (g_nproc > 1){
    fprintf(stderr,"# [loop_extract] The MPI part of the code is broken, please"
                   "use it with a single MPI task\n");
    EXIT(1);
  }
#endif

  /***************************************************************************
   * report git version
   ***************************************************************************/
  if ( g_cart_id == 0 ) {
    fprintf(stdout, "# [loop_extract] git version = %s\n", g_gitversion);
  }


  /***************************************************************************
   * set number of openmp threads
   ***************************************************************************/
  set_omp_number_threads ();

  /***************************************************************************
   * initialize geometry fields
   ***************************************************************************/
  if ( init_geometry() != 0 ) {
    fprintf(stderr, "[loop_extract] Error from init_geometry %s %d\n", __FILE__, __LINE__);
    EXIT(4);
  }

  geometry();

  /***************************************************************************
   * some additional xchange operations
   ***************************************************************************/
  mpi_init_xchange_contraction(2);

  /***************************************************************************
   * set io process
   ***************************************************************************/
  io_proc = get_io_proc ();
  if( io_proc < 0 ) {
    fprintf(stderr, "[loop_extract] Error, io proc must be ge 0 %s %d\n", __FILE__, __LINE__);
    EXIT(14);
  }
  fprintf(stdout, "# [loop_extract] proc%.4d has io proc id %d\n", g_cart_id, io_proc );
  if (g_LoopExtract_OutQSq > g_LoopExtract_InQSq){
    if (g_proc_id == 0){
     fprintf(stderr, "# [loop_extract] We do not have momenta  <=%d\n", g_LoopExtract_OutQSq);
     EXIT(1);
    }
  }
  if (g_proc_id == 0){
    fprintf(stdout, "# [loop_extract] Filtering value of Qsq %d\n", g_LoopExtract_OutQSq);
  }
  if (g_LoopExtract_FilterLoopTypesNumber == 0){
   if (g_proc_id == 0){
    fprintf(stderr, "# [loop_extract] Specify the type of loops you are interested in e.g.\n");
    fprintf(stderr, "# [loop_extract] filterlooptype = Scalar, Loops\n");
    EXIT(1);
   }
  }

  if (g_proc_id == 0){
   fprintf(stdout, "# [loop_extract] Following loops will be filtered\n");
    for (int i=0; i<g_LoopExtract_FilterLoopTypesNumber; ++i){
      switch( g_LoopExtract_FilterLoopTypes[i] ){
        case LOOP_EXTRACT_LOOP_TYPE_NAIVE:
          if (g_proc_id == 0) fprintf(stdout, "# [loop_extract] Naive\n");
          break;
        case LOOP_EXTRACT_LOOP_TYPE_SCALAR:
          if (g_proc_id == 0) fprintf(stdout, "# [loop_extract] Scalar\n");
          break;
        case LOOP_EXTRACT_LOOP_TYPE_DOP:
          if (g_proc_id == 0) fprintf(stdout, "# [loop_extract] dOp\n");
          break;
        case LOOP_EXTRACT_LOOP_TYPE_LOOPS:
          if (g_proc_id == 0) fprintf(stdout, "# [loop_extract] Loops\n");
          break;
        case LOOP_EXTRACT_LOOP_TYPE_LPSDW:
          if (g_proc_id == 0) fprintf(stdout, "# [loop_extract] LpsDw\n");
          break;
        case LOOP_EXTRACT_LOOP_TYPE_LOOPSCV:
          if (g_proc_id == 0) fprintf(stdout, "# [loop_extract] LoopsCV\n");
          break;
        case LOOP_EXTRACT_LOOP_TYPE_LPSDWCV:
          if (g_proc_id == 0) fprintf(stdout, "# [loop_extract] LpsDwCv\n");
          break;
        default:
          EXIT_WITH_MSG(122, "Invalid loop type passed!\n");
          break;
      }
    }
  }
  /***************************************************************************
  * On some ensembles (for example cB211.072.64 light) the HDFfile traj TAG
  * has not been set correctly, all the conf have 4 for this tag, we account
  * for this be allowing setting g_LoopExtract_LegacyTraj
  * **************************************************************************/
  if (g_LoopExtract_LegacyTraj == 1)
    conf_traj=4;
  else
    conf_traj=Nconf;

  /***************************************************************************
   * loop data filename
   ***************************************************************************/
  char accumulate;
  if (g_LoopExtract_NstochAccumulated == 1){
    accumulate='N';
  } else if (g_LoopExtract_NstochAccumulated == 0){
    accumulate='n';
  } else {
    if (g_proc_id == 0){
      fprintf(stderr, "# [loop_extract] invalid option for LoopExtract NstochAccumulated\n");
      EXIT(1);
    }
  }

  if (g_LoopExtract_FilenameBuild == 1){
    snprintf(filename, 1200, "%s/%s.%.4d_%s_Ns%.4d_step%.4d_Qsq%d.h5",
          g_LoopExtract_InPath, g_LoopExtract_FilenamePrefix,
          Nconf, g_LoopExtract_FilenameSuffix, g_LoopExtract_Nstoch,
          g_LoopExtract_NsaveStoch, g_LoopExtract_InQSq );
  }
  else{
    snprintf(filename, 1200, "%s/%s.%04d.h5", g_LoopExtract_InPath, g_LoopExtract_FilenamePrefix, Nconf);
  }

  if ( io_proc == 2 && g_verbose > 2 ) fprintf ( stdout, "# [loop_extract] loop filename = %s\n", filename );

  /***************************************************************************
   * count momenta and build momentum list
   ***************************************************************************/
  g_sink_momentum_number = 0;

  // we don't want to have to set the spatial size to its actual value because then
  // initialising the geometry takes needlessly long
  // Instead, we take the sqrt of InQSq, add one and then convert back to int
  const int max_mom_unit = (int)(sqrt( (double)g_LoopExtract_InQSq ) + 1); 
  for( int x1 = -max_mom_unit; x1 <= max_mom_unit; x1++ ) {
  for( int x2 = -max_mom_unit; x2 <= max_mom_unit; x2++ ) {
  for( int x3 = -max_mom_unit; x3 <= max_mom_unit; x3++ ) {
    int const qq = x1*x1 + x2*x2 + x3*x3;
    if ( qq <= g_LoopExtract_InQSq ) {
      g_sink_momentum_list[g_sink_momentum_number][0] = x1;
      g_sink_momentum_list[g_sink_momentum_number][1] = x2;
      g_sink_momentum_list[g_sink_momentum_number][2] = x3;
      g_sink_momentum_number++;
    }
  }}}
  if ( g_sink_momentum_number <= 0 ) {
    fprintf ( stderr, "[loop_extract] Error, momentum list is empty %s %d\n", __FILE__, __LINE__ );
    EXIT(1);
  } else {
    if (io_proc == 2 && g_verbose > 1 ) fprintf ( stdout, "# [loop_extract] number of momenta <= %3d is %3d\n", g_LoopExtract_InQSq, g_sink_momentum_number );
  }

  exitstatus = loop_get_momentum_list_from_h5_file ( g_sink_momentum_list, filename, g_sink_momentum_number, io_proc );
  if ( exitstatus != 0 ) {
    fprintf ( stderr, "[loop_extract] Error from loop_get_momentum_list_from_h5_file"
                      ", status was %d, %s %d\n", exitstatus, __FILE__, __LINE__ );
    EXIT(1);
  }


  int filtered_sink_momentum_list[MAX_MOMENTUM_NUMBER][3];
  int filtered_sink_momentum_index[MAX_MOMENTUM_NUMBER];
  int index=0;
  for (int i=0;i<g_sink_momentum_number; ++i)
  {
    int momentum_sq = g_sink_momentum_list[i][0]*g_sink_momentum_list[i][0] + 
                      g_sink_momentum_list[i][1]*g_sink_momentum_list[i][1] + 
                      g_sink_momentum_list[i][2]*g_sink_momentum_list[i][2];

    if (momentum_sq <= g_LoopExtract_OutQSq ){
      filtered_sink_momentum_list[index][0]=g_sink_momentum_list[i][0];
      filtered_sink_momentum_list[index][1]=g_sink_momentum_list[i][1];
      filtered_sink_momentum_list[index][2]=g_sink_momentum_list[i][2];
      filtered_sink_momentum_index[index]=i;
      index++;
    }
  }
  int filtered_sink_momentum_number = index;

  char *attribute_correlator;
//  get_attribute_from_h5_file (&attribute_correlator, filename, "Correlator-info", io_proc ) ;

  char *attribute_ensemble_info;
  get_attribute_from_h5_file (&attribute_ensemble_info, filename, "Ensemble-info", io_proc ) ;
  char attribute_nmoms[100];
  snprintf(attribute_nmoms,100, "%d",filtered_sink_momentum_number);
  char attribute_qsq[100];
  snprintf(attribute_qsq,100,"%d",(int)g_LoopExtract_OutQSq);


  /***************************************************************************
   * allocate memory for contractions
   ***************************************************************************/
  double **** loop = init_4level_dtable ( g_LoopExtract_Nstoch/g_LoopExtract_NsaveStoch, T, g_sink_momentum_number, 32 );
  if ( loop == NULL ) {
    fprintf(stderr, "[loop_extract] Error from init_4level_dtable %s %d\n", __FILE__, __LINE__ );;
    EXIT(48);
  }

  /***************************************************************************
   * allocate memory for filtered_contractions
   ***************************************************************************/
  double **** loop_filtered = init_4level_dtable ( g_LoopExtract_Nstoch/g_LoopExtract_NsaveStoch, T, filtered_sink_momentum_number, 32 );
  if ( loop_filtered == NULL ) {
    fprintf(stderr, "[loop_extract] Error from init_4level_dtable %s %d\n", __FILE__, __LINE__ );;
    EXIT(48);
  }

  if (( g_LoopExtract_SpinProject == 1) && (g_LoopExtract_SpinProjectGammaStructure_Number == 0)){
    fprintf(stderr, "[loop_extract] Spin projection is on, you must specify a list of gamma structures %s %d\n", __FILE__, __LINE__ );;
    EXIT(1);
  }

  int total_gamma;
  int index_gamma;
  if (g_LoopExtract_SpinProject == 0)
    total_gamma=1;
  else
    total_gamma=g_LoopExtract_SpinProjectGammaStructure_Number;
  /****************************************************************************************
 * Loop over the different gamma structures to be projected 
 * ****************************************************************************************/
  for (index_gamma=0; index_gamma<total_gamma;++index_gamma){

    if (g_LoopExtract_SpinProjectGammaStructure_List[index_gamma] != 4){
      snprintf(filename, 400, "%s.%.4d_%s_Ns%.4d_step%.4d_Qsq%d_gamma%d.h5",
              g_LoopExtract_OutputFilenamePrefix, Nconf, g_LoopExtract_OutputFilenameSuffix,
              g_LoopExtract_Nstoch, g_LoopExtract_NsaveStoch , (int)g_LoopExtract_OutQSq,
              g_LoopExtract_SpinProjectGammaStructure_List[index_gamma]  );
    } else {
      snprintf(filename, 400, "%s.%.4d_%s_Ns%.4d_step%.4d_Qsq%d.h5",
               g_LoopExtract_OutputFilenamePrefix, Nconf, g_LoopExtract_OutputFilenameSuffix,
               g_LoopExtract_Nstoch, g_LoopExtract_NsaveStoch, (int)g_LoopExtract_OutQSq);
    }

    if ( io_proc == 2 && g_verbose > 2 ) fprintf ( stdout, "# [loop_extract] loop filename = %s\n", filename );

    exitstatus = loop_write_momentum_list_to_h5_file ( filtered_sink_momentum_list, filename, filtered_sink_momentum_number, io_proc );
    if ( exitstatus != 0 ) {
      fprintf ( stderr, "[loop_extact] Error from loop_write_momentum_list_to_h5_file, status was %d %s %d\n", exitstatus, __FILE__, __LINE__ );
      EXIT(1);
    }
//    set_attribute_in_h5_file (attribute_correlator, filename, "Correlator-info", io_proc ) ;
    set_attribute_in_h5_file (attribute_ensemble_info, filename, "Ensemble-info", io_proc ) ;
    set_attribute_in_h5_file (attribute_nmoms, filename, "Nmoms", io_proc ) ;
    set_attribute_in_h5_file (attribute_qsq, filename, "Qsq", io_proc ) ;

//    free(attribute_correlator);
    free(attribute_ensemble_info);

    if ((io_proc == 2) && (g_LoopExtract_ASCII_Output == 1)) {
      for ( int imom = 0; imom < filtered_sink_momentum_number; imom++ ) {
        fprintf(stdout, " %3d  ( %3d, %3d, %3d)\n", 
                imom, filtered_sink_momentum_list[imom][0],
                filtered_sink_momentum_list[imom][1], filtered_sink_momentum_list[imom][2] );
      }
    } 

    /***************************************************************************
     * loop on different loop types
     ***************************************************************************/

    if (io_proc == 2) fprintf(stdout,"# [loop_extract] Loop number %d\n", g_LoopExtract_FilterLoopTypesNumber);
    for ( int iloop_type=0; iloop_type < g_LoopExtract_FilterLoopTypesNumber ; ++iloop_type){

    /***************************************************************************
     * loop on stochastic oet samples
     ***************************************************************************/
      for ( int isample = 0; isample < g_LoopExtract_Nstoch/g_LoopExtract_NsaveStoch; isample += 1 )
      {

        int const Nstoch = isample * g_LoopExtract_NsaveStoch + g_LoopExtract_FirstStochIndex;
        char loop_type[100];
        char loop_name[100];
        

        int inner_loop_length;

        switch( g_LoopExtract_FilterLoopTypes[iloop_type] ){
          case LOOP_EXTRACT_LOOP_TYPE_NAIVE:
            snprintf ( loop_type, 100, "%s", "Naive" );
            inner_loop_length=1;
            if (g_proc_id == 0) fprintf(stdout, "# [loop_extract] Naive loops\n");
            break;
          case LOOP_EXTRACT_LOOP_TYPE_SCALAR:
            snprintf ( loop_type, 100, "%s", "Scalar" );
            inner_loop_length=1;
            if (g_proc_id == 0) fprintf(stdout, "# [loop_extract] Scalar loops\n");
            break;
          case LOOP_EXTRACT_LOOP_TYPE_DOP:
            snprintf ( loop_type, 100, "%s", "dOp" );
            inner_loop_length=1;
            if (g_proc_id == 0) fprintf(stdout, "# [loop_extract] dOp loops\n");
            break;
          case LOOP_EXTRACT_LOOP_TYPE_LOOPS: 
            snprintf ( loop_type, 100, "%s", "Loops" );
            inner_loop_length=4;
            if (g_proc_id == 0) fprintf(stdout, "# [loop_extract] Loops loops\n");
            break;
          case LOOP_EXTRACT_LOOP_TYPE_LPSDW:
            snprintf ( loop_type, 100, "%s", "LpsDw" );
            inner_loop_length=4;
            if (g_proc_id == 0) fprintf(stdout, "# [loop_extract] LpsDw loops\n");
            break;
          case LOOP_EXTRACT_LOOP_TYPE_LOOPSCV:
            snprintf ( loop_type, 100, "%s", "LoopsCv" );
            inner_loop_length=4;
            if (g_proc_id == 0) fprintf(stdout, "# [loop_extract] LoopsCV loops\n");
            break;
          case LOOP_EXTRACT_LOOP_TYPE_LPSDWCV:
            snprintf ( loop_type, 100, "%s", "LpsDwCv" );
            inner_loop_length=4;
            if (g_proc_id == 0) fprintf(stdout, "# [loop_extract] LpsDwCv loops\n");
            break;
          default:
            EXIT_WITH_MSG(123, "[loop_extract]: Invalid loop type encountered!\n");
            break;
        }

        snprintf ( loop_name, 100, "%s", "loop" );

        for (int loop_direction=0; loop_direction < inner_loop_length ; ++loop_direction){

          char direction_part[100];
          snprintf(direction_part, 100, "dir_%02d", loop_direction);
     
          if (inner_loop_length > 1)
            snprintf(data_tag, 400, "/conf_%.4d/%cstoch_%.4d/%s/%s/%s",
                     conf_traj, accumulate, Nstoch, loop_type, direction_part, loop_name );
          else
            snprintf(data_tag, 400, "/conf_%.4d/%cstoch_%.4d/%s/%s",
                     conf_traj, accumulate, Nstoch, loop_type, loop_name );

          if ( io_proc == 2 && g_verbose > 2 ) fprintf( stdout, "# [loop_extract] data_tag = %s\n", data_tag);

          if (g_LoopExtract_FilenameBuild ==1){
            snprintf(filename, 1200, "%s/%s.%.4d_%s_Ns%.4d_step%.4d_Qsq%d.h5",
                   g_LoopExtract_InPath, g_LoopExtract_FilenamePrefix, Nconf,
                   g_LoopExtract_FilenameSuffix, g_LoopExtract_Nstoch, g_LoopExtract_NsaveStoch, g_LoopExtract_InQSq );
          }
          else{
            snprintf(filename, 1200, "%s/%s.%04d.h5", g_LoopExtract_InPath, g_LoopExtract_FilenamePrefix, Nconf);
          }

          exitstatus = loop_read_from_h5_file ( loop[isample], filename, data_tag, g_sink_momentum_number, 16, io_proc );
          if ( exitstatus != 0 ) {
            fprintf ( stderr, "[loop_extract] Error from loop_read_from_h5_file, status was %d %s %d\n", exitstatus, __FILE__, __LINE__ );
            EXIT(1);
          }
          if ( io_proc > 0 ) {
          /*****************************************************************
           * Collecting data for writing
           *****************************************************************/
            for ( int iproc = 0; iproc < g_nproc_t; iproc++ ) {

              if ( g_tr_id == iproc ) {
                char output_filename[400];
                FILE *ofs;
                if (io_proc== 2 && g_LoopExtract_ASCII_Output == 1){
                  snprintf ( output_filename, 400, "Nconf_%.4d.%s.%s_gamma%d", Nconf, loop_type, loop_name,g_LoopExtract_SpinProjectGammaStructure_List[index_gamma] );
                  ofs = fopen ( output_filename, "a" );
                  if ( ofs == NULL ) {
                    fprintf ( stderr, "[loop_extract] Error from fopen %s %d\n", __FILE__, __LINE__ );
                    EXIT(1);
                  }
                  fprintf ( ofs, "# [loop_extract] %s\n", data_tag );
                } /* end of g_LoopExtract_ASCII_Output */
                for ( int x0 = 0; x0 < T; x0++ ) {
                  int const y0 = x0 + g_proc_coords[0] * T;

                  for ( int imom = 0; imom < filtered_sink_momentum_number; imom++ ) {
                    double sp[32];

                    if (g_LoopExtract_SpinProject == 1){
                      _fm_eq_gamma_ti_fm(sp, g_LoopExtract_SpinProjectGammaBasis, g_LoopExtract_SpinProjectGammaStructure_List[index_gamma], loop[isample][x0][filtered_sink_momentum_index[imom]]);
                    }
                    else{
                      int index;
                      for (index=0; index<32;++index){
                        sp[index]=loop[isample][x0][filtered_sink_momentum_index[imom]][index];
                      } 
                    }

                    if (g_LoopExtract_SpinTrace == 1){

                    /* Taking the trace */

                      loop_filtered[isample][x0][imom][0] = sp[0]+sp[10]+sp[20]+sp[30];
                      loop_filtered[isample][x0][imom][1] = sp[1]+sp[11]+sp[21]+sp[31];
                      loop_filtered[isample][x0][imom][0] = loop_filtered[isample][x0][imom][0]/(double)g_LoopExtract_NsaveStoch ;
                      loop_filtered[isample][x0][imom][1] = loop_filtered[isample][x0][imom][1]/(double)g_LoopExtract_NsaveStoch ;

                    }
                    else{
                      for( int ic = 0; ic < 16; ic++ ) {
                        loop_filtered[isample][x0][imom][2*ic+0]=sp[2*ic+0];
                        loop_filtered[isample][x0][imom][2*ic+1]=sp[2*ic+1];
                        loop_filtered[isample][x0][imom][2*ic+0]=loop_filtered[isample][x0][imom][2*ic+0]/(double)g_LoopExtract_NsaveStoch;
                        loop_filtered[isample][x0][imom][2*ic+1]=loop_filtered[isample][x0][imom][2*ic+1]/(double)g_LoopExtract_NsaveStoch;
                      }   
                    }
                    if (io_proc== 2 && g_LoopExtract_ASCII_Output == 1){
                    /*****************************************************************
                     * write in ASCII format
                     *****************************************************************/
                      for( int ic = 0; ic < 16; ic++ ) {
                        fprintf ( ofs, " %3d %4d   %3d% 3d% 3d   %d %d  %25.16e %25.16e\n", Nstoch, y0, 
                         g_sink_momentum_list[filtered_sink_momentum_index[imom]][0], g_sink_momentum_list[filtered_sink_momentum_index[imom]][1], g_sink_momentum_list[filtered_sink_momentum_index[imom]][2],
                         ic/4, ic%4, sp[2*ic]/(double)g_LoopExtract_NsaveStoch, sp[2*ic+1]/(double)g_LoopExtract_NsaveStoch );
                      }  /* end of loop on components */
                    }  /* end of g_LoopExtract_ASCII_Output */
                  }  /* end of loop on momenta */
                }  /* end of loop on time slices */

                if (io_proc== 2 && g_LoopExtract_ASCII_Output == 1){
                  fclose ( ofs );
                }
              }  /* end of if g_tr_id == iproc */

#ifdef HAVE_MPI
              MPI_Barrier ( g_tr_comm );
#endif
            }  /* end of loop on procs in time direction */
          }  /* end of if io_proc > 0 */

        /*****************************************************************
         * Write data in h5 format with the correct tag
         *****************************************************************/

          if (inner_loop_length > 1)
            snprintf (data_tag, 400, "/conf_%.4d/%cstoch_%.4d/%s/%s",
                      Nconf, accumulate, Nstoch, loop_type, direction_part);
          else
            snprintf (data_tag, 400, "/conf_%.4d/%cstoch_%.4d/%s",
                      Nconf, accumulate, Nstoch, loop_type);

          if ( io_proc == 2 && g_verbose > 2 ) fprintf( stdout, "# [loop_extract] data_tag = %s\n", data_tag);

          if (g_LoopExtract_SpinProjectGammaStructure_List[index_gamma] != 4){
            snprintf ( filename, 400, "%s.%.4d_%s_Ns%.4d_step%.4d_Qsq%d_gamma%d.h5", g_LoopExtract_OutputFilenamePrefix, Nconf, g_LoopExtract_FilenameSuffix, g_LoopExtract_Nstoch, g_LoopExtract_NsaveStoch, (int)g_LoopExtract_OutQSq, g_LoopExtract_SpinProjectGammaStructure_List[index_gamma]);
          }
          else{
            snprintf ( filename, 400, "%s.%.4d_%s_Ns%.4d_step%.4d_Qsq%d.h5", g_LoopExtract_OutputFilenamePrefix, Nconf, g_LoopExtract_OutputFilenameSuffix, g_LoopExtract_Nstoch, g_LoopExtract_NsaveStoch, (int)g_LoopExtract_OutQSq);
          }
          if ( io_proc == 2 && g_verbose > 2 ) fprintf ( stdout, "# [loop_extract] loop filename = %s\n", filename );

          if (g_LoopExtract_SpinTrace == 1){
            exitstatus = contract_loop_write_to_h5_file ( loop_filtered[isample], filename, data_tag, filtered_sink_momentum_number, 1, io_proc );
          }
          else{
            exitstatus = contract_loop_write_to_h5_file ( loop_filtered[isample], filename, data_tag, filtered_sink_momentum_number, 16, io_proc );
          }

          if ( exitstatus != 0 ) {
            fprintf ( stderr, "[loop_extract] Error from loop_write_to_h5_file, status was %d %s %d\n", exitstatus, __FILE__, __LINE__ );
            EXIT(1);
          }

        }  /* end of loop for the different loop directions */

      } /* end of loop on oet samples */

    }  /* loop on different types */

  }/* end of loop on different gamma structures */
  /***************************************************************************
   * decallocate fields
   ***************************************************************************/
  fini_4level_dtable ( &loop );
  fini_4level_dtable ( &loop_filtered ); 

  /***************************************************************************
   * free the allocated memory, finalize
   ***************************************************************************/

  free_geometry();

#ifdef HAVE_MPI
  mpi_fini_xchange_contraction();
  mpi_fini_datatypes();
  MPI_Finalize();
#endif

  if(g_cart_id==0) {
    g_the_time = time(NULL);
    fprintf(stdout, "# [loop_extract] %s# [loop_extract] end of run\n", ctime(&g_the_time));
    fprintf(stderr, "# [loop_extract] %s# [loop_extract] end of run\n", ctime(&g_the_time));
  }

  return(0);

}
