#ifndef  MARGINALLIKE_MIGRATE
#define MARGINALLIKE_MIGRATE
// marginal likelihood summaries in migrate
// MIT opensource license
// consolidation from other files to simplify future changes
// (c) Peter Beerli 2021
//

#include "migration.h"

extern void calculate_BF(world_fmt **universe, option_fmt *options);
extern MYREAL combine_scaling_factor(world_fmt *world);
#if defined(MPI) && !defined(PARALIO)
extern void      print_marginal_like(float *temp, long *z, world_fmt * world);
#else /*not MPI*/
extern void      print_marginal_like(char *temp, long *c, world_fmt * world);
#endif
extern double ss_locus_logml (world_fmt *world, long locus);
extern boolean ckpt_locus_skipped (option_fmt *options, long locus);
extern void ckpt_write_locus (world_fmt *world, option_fmt *options, long locus);
extern void ckpt_remove_all (option_fmt *options, long loci);
extern void ckpt_restore (world_fmt *world, option_fmt *options);
extern boolean ti_binned (world_fmt *world, long locus, double *ti, double *bti);
extern MYREAL sumbezier(long intervals, MYREAL x0, MYREAL y0, MYREAL x1, MYREAL y1, MYREAL x2, MYREAL y2, MYREAL *ratio);

#endif
