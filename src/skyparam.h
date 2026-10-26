#ifndef __SKYPARAM__
#define __SKYPARAM__
// skyline plot parametrized
#include "migration.h"
#include "tree.h"

void insert_time_boundaries(timelist_fmt *timevector, world_fmt *world);
#ifdef MPI
void skyline_param_record(float *temp, long *z, world_fmt * world);
#else
void skyline_param_record(char *temp, long *c, world_fmt * world);
#endif
MYREAL bayes_update_timeparam(world_fmt * world, boolean *success);
/* one row of the skyline=PARAM table (locus -1: all loci; t1 -1: open) */
typedef struct
{
  char label[64];
  long seg;
  double t0;
  double t1;
  long locus;
  double q025;
  double q50;
  double q975;
  double v025;   /* segment value (parameter times multiplier), -1 if none */
  double v50;
  double v975;
  double ess;
} skyline_row_fmt;
void skyline_param_alloc(world_fmt *world, long loci);
void skyline_param_sample(world_fmt *world);
void skyline_param_count(world_fmt *world, boolean accepted);
void skyline_locus_finish(world_fmt *world, long locus);
long skyline_param_bufsize(world_fmt *world);
long skyline_param_pack(world_fmt *world, long locus, MYREAL *buf);
long skyline_param_unpack(world_fmt *world, long locus, MYREAL *buf);
long skyline_param_summary(world_fmt *world, skyline_row_fmt **rows, double *trials, double *accept);
void skyline_param_print(FILE *file, world_fmt *world);
#endif /*SKYPARAM*/
