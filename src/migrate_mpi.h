#ifndef MIGRATE_MPI_
#define MIGRATE_MPI_
/*
 (c) Peter Beerli 2013 Tallahassee FL
 
 Permission is hereby granted, free of charge, to any person obtaining
 a copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject
 to the following conditions:
 
 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
 ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/


#ifdef MPI
#ifdef POOCH
#include "poochmpi.h"
#else
#include <mpi.h>
#endif
#include <stdio.h>
#include <stdarg.h>
#include "migration.h"

#include "bayes.h"

#define MYINT  int
#define REPTAG 90
#define PRINTTAG 800
#define TEMPTAG 7000
#define BURNTAG 6000
#define ONDEMANDTAG 9000

typedef struct _filedb_fmt {
	FILE *file;
	long handle;
} filedb_fmt;


extern filedb_fmt filedb[30];
extern long filenum;

extern MPI_Comm comm_world;
extern MPI_Comm comm_workers;

#endif
extern int myID;
extern int myRepID;
extern int allID;
extern int color;
extern int numcpu;
extern int locidone;
extern int profiledone;

#define SLAVE 1

#define MIGMPI_END 0
#define MIGMPI_LIKE 1
#define MIGMPI_GRADIENT 2
#define MIGMPI_RESULT 3
#define MIGMPI_SUMFILE 4
#define MIGMPI_MIGHIST 5
#define MIGMPI_BAYESHIST 6
#define MIGMPI_SKYLINE 7
#define MIGMPI_TREESPACE 8
#define MIGMPI_PLOTPLANE 9
#define MIGMPI_HAPLOTYPING 10
#define MIGMPI_PARALIO 11
#define MIGMPI_ASSIGN 12
#define MIGMPI_SEQERROR 13
#ifdef MPI

#ifndef PTHREADS
#define tpool_t char
#else
#include "heating.h"
#endif

extern void broadcast_options (option_fmt * options);

extern void mpi_runloci_master (long loci, int *who, world_fmt *world, option_fmt *options, data_fmt *data, boolean options_readsum, boolean menu);

extern void mpi_runloci_worker (world_fmt ** universe, int usize,
                                    option_fmt * options, data_fmt * data,
                                    tpool_t * heating_pool, long maxreplicate,
                                    long *treefilepos, long *Gmax);

extern MYREAL mpi_likelihood_master (MYREAL *param, MYREAL *lparam,
                                         world_fmt * world, nr_fmt * nr,
                                         helper_fmt * helper, int *who);
extern void broadcast_options_master (option_fmt * options, data_fmt *data);
extern void broadcast_options_worker (option_fmt * options);

extern void broadcast_data_master (data_fmt * data, option_fmt * options);
extern void broadcast_data_worker (data_fmt * data, option_fmt * options, world_fmt *world);
extern void mpi_results_master (long sendtype, world_fmt * world,
                                    long maxrep,
                                    void (*unpack) (MYREAL *buffer,
                                                    world_fmt * world, long locus,
                                                    long maxrep, long numpop));
extern void mpi_results_worker (long bufsize, world_fmt * world, long maxrep,
                                    long (*pack) (MYREAL **buffer,
                                                  world_fmt * world, long locus,
                                                  long maxrep, long numpop));

extern long pack_result_buffer (MYREAL **buffer, world_fmt * world,
                                    long locus, long maxrep, long numpop);
extern void unpack_result_buffer (MYREAL *buffer, world_fmt * world,
                                      long locus, long maxrep, long numpop);
extern long pack_sumfile_buffer (MYREAL **buffer, world_fmt * world,
                                     long locus, long maxrep, long numpop);
extern void unpack_sumfile_buffer (MYREAL *buffer, world_fmt * world,
                                       long locus, long maxrep, long numpop);
extern long pack_treespace_buffer (MYREAL **buffer, world_fmt * world,
				   long locus, long maxrep, long numpop);
extern void unpack_treespace_buffer (MYREAL *buffer, world_fmt * world,
				     long locus, long maxrep, long numpop);

extern long pack_mighist_buffer (MYREAL **buffer, world_fmt * world,
                                     long locus, long maxrep, long numpop);
extern void unpack_mighist_buffer (MYREAL *buffer, world_fmt * world,
                                       long locus, long maxrep, long numpop);
extern long pack_skyline_buffer (MYREAL **buffer, world_fmt * world,
				 long locus, long maxrep, long numpop);
extern void unpack_skyline_buffer (MYREAL *buffer, world_fmt * world,
				   long locus, long maxrep, long numpop);

extern long unpack_single_bayes_buffer(MYREAL *buffer,bayes_fmt * bayes, world_fmt * world,long locus);

extern long pack_bayes_buffer (MYREAL **buffer, world_fmt * world,
                               long locus, long maxrep, long numpop);

extern void unpack_bayes_buffer (MYREAL *buffer, world_fmt * world,
                                 long locus, long maxrep, long numpop);

extern long pack_haplotypes_buffer (MYREAL **buffer, world_fmt * world,
                                 long locus, long maxrep, long numpop);
extern void unpack_haplotypes_buffer (MYREAL *buffer, world_fmt * world,
                                 long locus, long maxrep, long numpop);

extern long pack_assign_buffer(MYREAL **buffer, world_fmt * world,
			       long locus, long maxrep, long numpop);
extern void unpack_assign_buffer(MYREAL *buffer, world_fmt * world,
				 long locus, long maxrep, long numpop);

extern void unpack_seqerror_buffer(MYREAL *buffer, world_fmt * world,
			    long locus, long maxrep, long numpop);


extern void mpi_gradient_master (nr_fmt * nr, world_fmt * world, int *who);

extern void mpi_maximize_worker (world_fmt * world, option_fmt *options,long kind, long rep);

extern void mpi_startparam_master(world_fmt * world);
extern void mpi_startparam_worker(world_fmt * world);

extern void mpi_gmax_master(world_fmt * world, long *Gmax);
extern void mpi_gmax_worker(world_fmt * world);

extern void mpi_send_stop (world_fmt * world);
extern void mpi_results_stop (void);
extern void assignloci_worker (world_fmt * world, option_fmt *options, long * Gmax);
void setup_parameter0_mpi (world_fmt * world, nr_fmt * nr, long repkind,
                           long repstart, long repstop, long loci, long kind,
                           boolean multilocus);

extern void setup_filehandle_db(FILE *file, world_fmt *world, option_fmt *options, data_fmt *data);
extern long get_filehandle(void *file, world_fmt *world, option_fmt *options, data_fmt *data);
extern void mpi_fprintf(FILE *file, const char *fmt, ...);
extern void mpi_fprintf2(FILE *file, long filesize, const char *fmt, ...);
#ifdef PARALIO
extern void mpi_mdim_send(MPI_File *file, char *values, long size);
#else
extern void mpi_mdim_send(float *values, long size);
#endif
//extern void mpi_mdim_send(float *values, long size);
extern long get_filehandle2(void *file, world_fmt *world);
extern void send_divtime(const char *fmt, ...);

#ifdef SLOWNET
extern void
mpi_broadcast_results (world_fmt * world, long loci,
                       long (*pack) (MYREAL **buffer, world_fmt * world,
                                     long locus, long maxrep, long numpop),
                       void (*unpack) (MYREAL *buffer, world_fmt * world,
                                       long locus, long maxrep, long numpop));
extern void mpi_profiles_master (world_fmt * world, long nparam, int *profilewho);
extern void mpi_profiles_worker (world_fmt * world, long *gmaxptr);
extern void setup_parameter0_slowmpi (world_fmt * world, nr_fmt * nr, long repkind,
                               long repstart, long repstop, long loci,
                               long kind, boolean multilocus);
extern long  mpi_send_stop_mcmc_lociworker(long numcpu, long loci);
extern long  mpi_send_stop_mcmc_worker_orig(long numcpu, long loci, MPI_Comm *comm, MPI_Request *irequests, MPI_Status *istatus, long id);
#endif /*slownet */
extern long  mpi_send_stop_mcmc_replicateworker(long numcpu, long loci);
extern int get_replicant_color(int numcpu, long maxreplicate, long loci);
extern void request_data(long pop,long ind, long locus, long sublocus,long allelenum, world_fmt *world, data_fmt *data, option_fmt * options, site_fmt ***datapart);
//#include <sys/resource.h>
//extern void set_memory_limit(rlim_t softsize,rlim_t maxsize);
//extern void check_memory_limit();

#endif /*mpi */


#endif
