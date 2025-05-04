/*! \file=tree.h */
#ifndef TREE_INCLUDE
#define TREE_INCLUDE
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 T R E E B U I L D I N G   R O U T I N E S 
 
 Peter Beerli 1996, Seattle
 beerli@fsu.edu
 updated 2009,2016

(c) 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
(c) 2003-2016 Peter Beerli, Tallahassee FL
 
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

 
$Id: tree.h 2157 2013-04-16 22:27:40Z beerli $
-------------------------------------------------------*/
#include "migration.h"
void buildtree (world_fmt * world, option_fmt * options, data_fmt * data, long locus);
void allocate_lineages (timelist_fmt **timevector, const long offset, const long numpop);
void allocate_tip (world_fmt * world, option_fmt * options, node ** p, long pop, long locus, long a, long ind, char **tipnames);
void fix_times (world_fmt * world, option_fmt * options);
void first_smooth (world_fmt * world, long locus);
void smooth (const node * root, node * p, world_fmt * world, const long locus);
void set_all_dirty (const node * root, node * p, world_fmt * world, const long locus);
void set_dirty (node * p);

void set_pop (node * theNode, long pop, long actualpop);
void set_v (node * p);
void ltov (node * p);
void find_minmax_msat_allele (world_fmt * world, data_fmt * data, long locus, long *smallest, long *biggest);


void create_treetimelist (world_fmt * world, timelist_fmt ** ltl);
void increase_timelist2 (timelist_fmt ** timevector,long extender);
void construct_tymelist (world_fmt * world, timelist_fmt * timevector);
void add_partlineages (long numpop, timelist_fmt ** timevector, world_fmt * world);
node *add_migration (world_fmt *world, node * p, char event, long from, long to, MYREAL utime);

MYREAL treelikelihood (world_fmt * world);
MYREAL pseudotreelikelihood (world_fmt * world,proposal_fmt * proposal);

void pseudonuview (proposal_fmt * proposal, xarray_fmt *axx1, MYREAL **alx1, MYREAL v1, 
		   xarray_fmt *axx2, MYREAL **alx2, MYREAL v2);
MYREAL find_tipdate(char * id, long pop, world_fmt *world);
void find_tips (node * p, node ** nodelist, long *z);
void allocatetips (world_fmt * world, option_fmt * options, data_fmt * data, long locus);
void allocate_xseq(xarray_fmt *x, long sites, long categs);
void swap_tree (world_fmt * tthis, world_fmt * tthat);

void free_tree (node * p, world_fmt * world);
void calc_treelength (node * p, MYREAL *treelen);
void free_mignodelet (node * p, world_fmt * world);

void print_tree (world_fmt * world, long g, long *filepos);

#ifdef NEXUSTREE
void nexus_treesheader(world_fmt *world, option_fmt *options, data_fmt *data);
#endif

extern nuview_function * nuview;
extern prob_micro_function * prob_micro;
extern pseudonuview_function * pseudonuv;

#endif /*TREE_INCLUDE */
