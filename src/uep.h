#ifndef UEP_H
#define UEP_H
/* unique event polymorphism UEP
 (c) Peter Beerli 2013
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
#ifdef UEP

#include "migration.h"
#include "tree.h"

extern void alloc_ueplike (MYREAL ***ueplike, long size,
                               MYREAL *****ueplikestore, long storesize,
                               long ****ueprootstore, long numpop, long loci);
extern void alloc_ueptime (ueptime_fmt ** ueptime, long size,
                               ueptime_fmt **** ueptimestore, long storesize,
                               long numpop, long loci);

extern void allocate_uep (node * p, world_fmt * world, char datatype,
                              boolean withtips);

extern void destroy_uep(world_fmt * world, long storesize);


extern void constrain_distance_uep (int **uep, long uepsites, MYREAL **distm,
                                        long tips);
extern void make_uep_values (world_fmt * world, data_fmt * data, long locus);
extern void update_uep (node * p, world_fmt * world);
extern void check_uep_root (node * p, world_fmt * world);
extern boolean is_success_pseudo_uep (proposal_fmt * proposal);
extern boolean uep_compare (int *target, int *origin, long uepsites);
extern void store_uep (world_fmt * world);
extern MYREAL pseudo_ueplikelihood (world_fmt * world,
                                        proposal_fmt * proposal);
extern MYREAL ueplikelihood (world_fmt * world);
extern void calc_ueplike (node * p, world_fmt * world, MYREAL **ueplike);
extern MYREAL collect_pseudo_ueplike (proposal_fmt * proposal,
                                          MYREAL **ueplike);
extern void collect_ueplike (node * p, node * d1, node * d2, long uepsite,
                                 world_fmt * world, MYREAL *ueplike);
extern void show_uep_store (world_fmt * world);
extern node *first_uep (node * p, node * nop, long uepsites);
extern node *first_uep2 (node * p, node * root, long uepsites);
extern void analyze_uep (world_fmt * world);
extern void alloc_uepanc (long ***uepanc, long size, long numpop, long loci);
extern void update_uepanc (world_fmt * world);
extern void setup_uep (world_fmt * world, option_fmt * options);
extern void pseudonu_twostate (proposal_fmt * proposal, ueparray_fmt * xx1,
                                   MYREAL *lx1, MYREAL v1, ueparray_fmt * xx2,
                                   MYREAL *lx2, MYREAL v2);
extern MYREAL pseudo_tl_uep (ueparray_fmt * xx1, ueparray_fmt * xx2,
                                 MYREAL v1, MYREAL v2, proposal_fmt * proposal,
                                 world_fmt * world);
extern void twostate_nuview (node * mother, world_fmt * world, const long locus);

extern void copy_uepx (proposal_fmt * proposal, ueparray_fmt xx1, ueparray_fmt xx2);

#endif


#endif /* __UEP__ */
