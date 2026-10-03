#ifndef JOINT_COMBINE_H
#define JOINT_COMBINE_H
/* Joint multi-locus combination of Bayesian posteriors, see joint_combine.c */
#include "migration.h"

/// one row of the Monte Carlo error table of the joint estimates
typedef struct
{
  char name[LINESIZE];
  double median, mcerr, blo, bhi, sd, ratio;
  double booterr;     /* -1: no bootstrap */
  boolean flagged;    /* too noisy: All is used instead of Joint */
} jc_mcerr_row;

extern void jc_set_run_options (const char *outfilename, boolean recover);
extern boolean jc_supported (world_fmt *world);
extern boolean jc_in_memory (world_fmt *world);
extern void jc_record_sample (world_fmt *world);
extern void jc_clear_locus (world_fmt *world, long locus);
extern long jc_pack_buffer (MYREAL **buffer, world_fmt *world, long locus, long maxrep, long numpop);
extern void jc_unpack_buffer (MYREAL *buffer, world_fmt *world, long locus, long maxrep, long numpop);
extern boolean jc_combine (world_fmt *world);
extern boolean jc_fill_histogram (world_fmt *world, bayeshistogram_fmt *hist);
extern void jc_print_note (world_fmt *world, FILE *out);
extern long jc_mcerr_rows (world_fmt *world, jc_mcerr_row **rows);
extern void jc_mcerr_info (world_fmt *world, long *nblock, boolean *byrep, double *bootess);
extern boolean jc_param_flagged (world_fmt *world, long p);
extern boolean jc_joint_scaling (world_fmt *world, double *logc, double *err);
#ifdef MPI
extern void jc_worker_service (world_fmt *world, long part);
#endif

#endif
