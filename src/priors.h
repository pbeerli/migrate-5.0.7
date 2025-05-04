#ifndef PRIORS__
#define PRIORS__
/*------------------------------------------------------
Bayesian inference 
of migration rate  and effectice population size
using a Metropolis-Hastings Monte Carlo algorithm
-------------------------------------------------------
    Prior distributions   R O U T I N E S

send questions concerning this software to:
Peter Beerli
beerli@fsu.edu

Copyright 2013 Peter Beerli 

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

 $Id$
    */
/*! \file priors.h 

priors.h supplies definitions of prior distributions specified in priors.c
*/
#include "migration.h"
extern void is_priorkind(prior_fmt *p, char *priorkind);
extern MYREAL propose_uniform(MYREAL param, MYREAL minparam, MYREAL maxparam, MYREAL *r, MYREAL delta);
extern MYREAL propose_uni_newparam (MYREAL param, long which, world_fmt * world, MYREAL *r);
extern MYREAL propose_exp_newparam (MYREAL param, long which, world_fmt * world, MYREAL *r);
extern MYREAL propose_expb_newparam (MYREAL param,long which, world_fmt * world, MYREAL *r);
extern MYREAL propose_mult_newparam (MYREAL param,long which, world_fmt * world, MYREAL *r);
extern MYREAL propose_normal_newparam (MYREAL param,long which, world_fmt * world, MYREAL *r);
extern MYREAL propose_gamma_newparam (MYREAL param, long which, world_fmt * world, MYREAL *r);
extern MYREAL      log_prior_ratio_uni  (MYREAL newparam, MYREAL oldparam, bayes_fmt *bayes, long which);
extern MYREAL      log_prior_ratio_exp  (MYREAL newparam, MYREAL oldparam, bayes_fmt *bayes, long which);
extern MYREAL      log_prior_ratio_wexp (MYREAL newparam, MYREAL oldparam, bayes_fmt *bayes, long which);
extern MYREAL      log_prior_ratio_mult (MYREAL newparam, MYREAL oldparam, bayes_fmt *bayes, long which);
extern MYREAL      log_prior_ratio_normal (MYREAL newparam, MYREAL oldparam, bayes_fmt *bayes, long which);
extern MYREAL      log_prior_uni  (world_fmt * world, long numparam);//for heating
extern MYREAL      log_prior_exp  (world_fmt * world, long numparam);//for heating
extern MYREAL      log_prior_wexp (world_fmt * world, long numparam);//for heating
extern MYREAL      log_prior_mult (world_fmt * world, long numparam);//for heating
extern MYREAL      log_prior_normal (world_fmt * world, long numparam);//for heating
extern MYREAL      log_prior_uni1  (world_fmt * world, long numparam, MYREAL);
extern MYREAL      log_prior_exp1  (world_fmt * world, long numparam, MYREAL);
extern MYREAL      log_prior_wexp1 (world_fmt * world, long numparam, MYREAL);
extern MYREAL      log_prior_mult1 (world_fmt * world, long numparam, MYREAL);
extern MYREAL      log_prior_normal1 (world_fmt * world, long numparam, MYREAL);
extern MYREAL hastings_ratio_uni(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r,
				 bayes_fmt * bayes, long whichparam);
extern MYREAL hastings_ratio_exp(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r,
				 bayes_fmt * bayes, long whichparam);
extern MYREAL hastings_ratio_expb(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r,
				  bayes_fmt * bayes, long whichparam);
extern MYREAL hastings_ratio_mult(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r,
				  bayes_fmt * bayes, long whichparam);
extern MYREAL hastings_ratio_normal(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r,
				    bayes_fmt * bayes, long whichparam);





extern MYREAL find_beta_truncgamma(MYREAL mean, MYREAL alpha, MYREAL lower, MYREAL upper);
extern MYREAL hastings_ratio_gamma(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r,
				   bayes_fmt * bayes, long whichparam);
extern MYREAL log_prior_ratio_gamma(MYREAL newparam, MYREAL oldparam,  bayes_fmt * bayes,  long which);
extern MYREAL log_prior_gamma(world_fmt *world, long numparam);
extern MYREAL log_prior_gamma1(world_fmt *world, long numparam, MYREAL val);
extern MYREAL hastings_ratio_gamma(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r,
				   bayes_fmt * bayes, long whichparam);
extern void set_option_prior(prior_fmt **p, int type, MYREAL mini, MYREAL maxi, MYREAL mean, long bins);

extern prior_fmt * copy_option_prior(prior_fmt *pmodel, option_fmt *options);
extern prior_fmt * get_prior_list(prior_fmt **list, int type);
extern void check_bayes_priors(option_fmt *options, data_fmt * data, world_fmt *world);
extern long get_species_record(world_fmt * world,long which);
extern void hyper_normal(MYREAL *mean, MYREAL *std, MYREAL origmean, MYREAL origstd, MYREAL lower, MYREAL upper);
extern MYREAL normal_rand(MYREAL mean, MYREAL std);
extern double trunc_random_normal2(double mi, double ma, double mean, double std);
extern void init_hyperpriorrecord(hyper_fmt ** hyperp, long numparam);

extern prior_fmt * find_prior_menu(long from, long to, long priortype, option_fmt * options);
extern long set_all_default_priors(option_fmt *options);

#endif
