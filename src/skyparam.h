#ifndef __SKYPARAM__
// skyline plot parametrized
#include "migration.h"
#include "tree.h"

void insert_time_boundaries(timelist_fmt *timevector, world_fmt *world);
void skyline_param_reader(world_fmt *world, long step, long locus, char **input);
MYREAL bayes_update_timeparam(world_fmt * world, boolean *success);
#endif /*SKYPARAM*/
