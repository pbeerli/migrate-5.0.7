#ifndef TABLES_MIGRATE// Tables printing in migrate
#define TABLES_MIGRATE
// MIT opensource license
// consolidation from other files to simplify future changes
// (c) Peter Beerli 2021
//

#include "migration.h"
#include "pretty.h"


extern void print_bayesfactor(world_fmt **universe, option_fmt * options);
extern void print_burnin_autostop(world_fmt * world);
extern void print_heatingreport(world_fmt **universe, option_fmt* options);
extern void print_mcmc_run_character(world_fmt * world);
extern void print_finaltree(world_fmt *world, option_fmt *options);

#endif
