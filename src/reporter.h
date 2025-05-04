#ifndef REPORTER_INCLUDE
#define REPORTER_INCLUDE
/** \file reporter.h */
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 
 R E P O R T E R   R O U T I N E S 
                                                                                                               
 Peter Beerli 1999, Seattle
 beerli@fsu.edu
 
(c) 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
(c) 2003-2004 Peter Beerli, Tallahassee FL
 
 
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

$Id: reporter.h 2157 2013-04-16 22:27:40Z beerli $
-------------------------------------------------------*/

#include "migration.h"

extern void convergence_check (world_fmt * world, boolean progress);
extern void calc_chain_s(MYREAL *cs, MYREAL *cm, world_fmt *world, long replicate);
extern void convergence_check_bayes (world_fmt *world, long maxreplicate);
extern void chain_means (MYREAL *thischainmeans, world_fmt * world);
extern void convergence_progress(FILE *file, world_fmt *world);
extern MYREAL single_chain_var(world_fmt *world, unsigned long T, MYREAL *variance, MYREAL *autoc, MYREAL *effsample);
extern boolean max_ess(const MYREAL * ess, const  long n, const MYREAL minimum, MYREAL *miness);
extern void print_bayes_ess(FILE *file, world_fmt *world, MYREAL *autocorr, MYREAL *effsample);
extern void calculate_ess_frombayes(world_fmt *world, long T, MYREAL *params, long locus, MYREAL *autoc, MYREAL *effsample);
extern void collect_ess_values(world_fmt *world);
extern void collect_acceptance(world_fmt *world);
#endif
