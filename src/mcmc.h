#ifndef MCMC_INCLUDE
#define MCMC_INCLUDE
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 M C M C   R O U T I N E S 
 
 Markov Monte Carlo stuff: treechange, acceptance
 
 
 Peter Beerli 1996, Seattle
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2010 Peter Beerli, Tallahassee FL
 
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

$Id: mcmc.h 2157 2013-04-16 22:27:40Z beerli $
-------------------------------------------------------*/

extern long tree_update (world_fmt * world, long g, boolean assign);
extern void free_timevector (timelist_fmt * timevector);
extern void free_masterproposal (proposal_fmt *proposal);
extern void count_migrations (node * p, long *count);
extern void traverse_tagnew(node *theNode, node *origin);
extern void set_tree_down_dirty (node * p);
extern void set_tree_dirty (node * p);
extern MYREAL eventtime (proposal_fmt * proposal, long pop, vtlist * tentry, char *event);
extern int migrateb (proposal_fmt * proposal, node * up, char event, long from, long to);
extern int migrate (proposal_fmt * proposal, node * up, char event, long from, long to);
extern boolean acceptlike (world_fmt * world, proposal_fmt * proposal, long g, timelist_fmt * tyme, boolean assign);
extern long migration_from (long to, proposal_fmt * proposal);
extern void new_proposal (proposal_fmt ** proposal, timelist_fmt * tl, world_fmt * world);
extern void construct_localtimelist (timelist_fmt * timevector, proposal_fmt * proposal);
extern void chooseTarget (proposal_fmt * proposal, timelist_fmt * timevector, node ** bordernodes, long *bordernum);
extern void chooseOrigin (proposal_fmt * proposal);
extern void new_localtimelist (timelist_fmt ** ntl, timelist_fmt * otl, long numpop);
extern void allocate_nodelist(node ***nodelist, long *oldnode_elem, long elements);
extern void traverse_check(node *theNode);
#endif
