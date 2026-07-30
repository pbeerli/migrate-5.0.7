#ifndef WINDOWUPDATE_H
#define WINDOWUPDATE_H
/*------------------------------------------------------------------------
Joint local-M + windowed migration-history MCMC move.

See windowupdate.c for the full design rationale. Kept in its own
translation unit rather than folded into bayes.c because the eventual
statistical core (a CTMC bridge sampler, a bystander mean-field forward
pass, windowed forward-filtering/backward-sampling) is substantial new
numerical code with no existing counterpart elsewhere in migrate, not an
extension of anything already in bayes.c.
------------------------------------------------------------------------*/
#include "migration.h"

extern long windowed_joint_update(world_fmt *world);

/* Exposed for task-level testing/reuse; not part of the public move API. */
extern long window_pick(world_fmt *world, node **window, long window_size);
extern long window_touched_branches(node **window, long window_n,
                                     node **touched, long touched_alloc);

#endif /*WINDOWUPDATE_H*/
