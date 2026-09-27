#ifndef MLF_HASTINGS_H
#define MLF_HASTINGS_H
/*------------------------------------------------------------------------
Mittag-Leffler Hastings correction for genealogy moves, see mlf_hastings.c
------------------------------------------------------------------------*/
#include "migration.h"

extern boolean mlh_supported(world_fmt *world);
extern double mlh_log_correction(world_fmt *world, proposal_fmt *proposal,
                                 timelist_fmt *R);
#ifdef MIGRATE_MLH_VERIFY
extern void mlh_verify_record(world_fmt *world, double dlogp, double dlogq);
extern void mlh_verify_before(world_fmt *world, double s);
extern void mlh_verify_after(world_fmt *world, proposal_fmt *R2prop,
                             timelist_fmt *R2, double s);
#endif

#endif /*MLF_HASTINGS_H*/
