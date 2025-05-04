// Bayes update scheme
#ifndef BAYESUPDATE_
#define BAYESUPDATE_
/*
send questions concerning this software to:
Peter Beerli
beerli@fsu.edu

Copyright 1997-2002 Peter Beerli and Joseph Felsenstein 
Copyright 2003-2008 Peter Beerli
Copyright 2009-2012 Peter Beerli and Michal Palczewski

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies
or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

$Id: bayes.h 2093 2012-09-06 23:41:30Z beerli $
started in 2000
*/
#include "migration.h"
extern boolean acceptBayes (MYREAL newval, MYREAL oldval);
extern long bayes_update (world_fmt *world);
extern void bayes_free(world_fmt *world);
extern void bayes_fill(world_fmt *world, option_fmt *options);
extern void bayes_init(bayes_fmt *bayes, world_fmt *world, option_fmt *options);
extern void bayes_save(world_fmt *world, long step);
extern void bayes_stat(world_fmt *world, data_fmt *data);
extern long setup_bayes_map(longpair *map, world_fmt *world, long size);
extern void bayes_init_histogram(world_fmt * world, option_fmt * options);
extern void adjust_bayes_bins(world_fmt * world, long locus);
extern void calculate_credibility_interval(world_fmt * world, long locus);
extern void bayes_reset(world_fmt * world);
extern void bayes_check_and_fix_param(world_fmt *world);
extern long number_replicates(world_fmt * world);
extern long number_replicates2(option_fmt *options);
#ifdef ZNZ
extern void print_bayes_mdimfileheader(znzFile file, long interval, world_fmt *world, data_fmt *data);
#else
extern void print_bayes_mdimfileheader(FILE *file, long interval, world_fmt *world, data_fmt *data);
#endif
extern boolean bayes_accept (MYREAL newval, MYREAL oldval, MYREAL heat, MYREAL hastingsratio);
extern void bayes_print_accept(FILE * outfile, world_fmt * world);
extern void bayes_print_hyperprior(FILE * file,  world_fmt *world);
extern MYREAL probg_treetimes(world_fmt *world);
extern void recalc_timelist (world_fmt * world, MYREAL new_ratio, MYREAL old_ratio);
extern void bayes_smooth(double *x, long xelem, long el, boolean lastfirst, boolean boundary);
extern void bayes_set_param(MYREAL *param, MYREAL newparam, long which, char *custm2, long numpop);
extern MYREAL log_prior_ratio_uni(MYREAL newparam, 
			   MYREAL oldparam, 
			   bayes_fmt * bayes, 
				  long which);
extern MYREAL log_prior_ratio_all(world_fmt *world, MYREAL *newvals);
extern MYREAL calculate_prior(world_fmt *world);
extern MYREAL scaling_prior(world_fmt *world, long numparam, MYREAL val);
extern void calc_hpd_credibility(world_fmt *world,long locus, long numpop2, long numparam);
extern void check_min_max_param(MYREAL *param0, MYREAL minparam, MYREAL maxparam);
extern void destroy_global_function_arrays(void);
extern void construct_param_hist(world_fmt *world, long locus, long npa, long pa, long numbin, MYREAL *mini,
                          MYREAL *maxi, double **results,
                          long *total, MYREAL *themean, MYREAL *thestd);

extern logpriorratioptr * log_prior_ratio;
extern logpriorptr * log_prior;
extern logprior1ptr * log_prior_1;
extern profuncptr * propose_new;
extern hastratioptr * hastings_ratio;

extern MYREAL probg_treetimes_local(world_fmt* world, timelist_fmt * treetimes);
#endif
