#ifndef GROWTHINCLUDE
#define GROWTHINCLUDE
/*------------------------------------------------------
 inference of population parameters
 using a Metropolis-Hastings Monte Carlo algorithm
 -------------------------------------------------------
 growth routines
 
 Peter Beerli 2013, Tallahassee
 beerli@fsu.edu
 
 Copyright 2017 Peter Beerli
 
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
#include "migration.h"
extern void reset_growth(world_fmt * world);
#if defined(MPI) && !defined(PARALIO) /* */
extern void print_growth_record(float *temp, long *z, world_fmt *world);
#else
extern void  print_growth_record(char *temp, long *c, world_fmt * world);
#endif
extern void construct_locusgrowth_histogram(world_fmt *world, long locus, MYREAL *mini, MYREAL *maxi, double **results);
extern boolean init_growpop(worldoption_fmt * wopt, option_fmt *options, long numpop);
extern void init_growth(world_fmt * world, long numpop);
#endif
