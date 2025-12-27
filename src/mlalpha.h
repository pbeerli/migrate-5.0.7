#ifndef MLALPHAINCLUDE
#define MLALPHAINCLUDE
/*------------------------------------------------------
 inference of population parameters
 using a Metropolis-Hastings Monte Carlo algorithm
 -------------------------------------------------------
 mlalpha routines
 
 Peter Beerli 2025, Tallahassee
 beerli@fsu.edu
 
 Copyright 2025 Peter Beerli
 
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
extern void reset_mlalpha(world_fmt * world);
#if defined(MPI) && !defined(PARALIO) /* */
extern void print_mlalpha_record(float *temp, long *z, world_fmt *world);
#else
extern void  print_mlalpha_record(char *temp, long *c, world_fmt * world);
#endif
extern void construct_locusmlalpha_histogram(world_fmt *world, long locus, MYREAL *mini, MYREAL *maxi, double **results);
extern boolean init_mlalphapop(worldoption_fmt * wopt, option_fmt *options, long numpop);
extern void init_mlalpha(world_fmt * world, long numpop);
extern void print_parm_mlalphapops(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options, data_fmt *data);
extern void set_mlalpha(char **value, char **tmp, option_fmt *options);
#endif
