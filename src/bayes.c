/*------------------------------------------------------
 Bayesian inference
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm
 -------------------------------------------------------
 Bayesian   R O U T I N E S
 
 send questions concerning this software to:
 Peter Beerli
 beerli@fsu.edu
 
 Copyright 1997-2002 Peter Beerli and Joseph Felsenstein
 Copyright 2003-2008 Peter Beerli
 Copyright 2009-2018 Peter Beerli and Michal Palczewski
 
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
 
 $Id: bayes.c 2170 2013-09-19 12:08:27Z beerli $
 */
/*! \file bayes.c
 
 bayes.c contains functions to handle the Bayesian implementation of migrate
 it is triggered through the options->bayes_infer=TRUE.
 
 */
#include <assert.h>
#include <stdlib.h>
#include "definitions.h"
#include "migration.h"
#include "options.h"
#include "bayes.h"
#include "random.h"
#include "tools.h"
#include "sighandler.h"
#include "world.h"
#include "tree.h"
#include "slice.h"
#include "reporter.h"
#include "haplotype.h"
#include "mcmc.h"
#include "skyparam.h"
#include "speciate.h"
#include "correlation.h"
#include "autotune.h"
#include "priors.h"
#include "mittag_leffler.h"
#include "growth.h"
#include "savitzky_golay.h"
#include "kernel_smooth.h"
#ifdef PRETTY
#include "pretty.h"
#endif
#ifdef MPI
#include "migrate_mpi.h"
#else
extern int myID;
#endif

// counter for bayesallfile on a per locus basis
//extern long * mdimfilecount;
extern double page_height;

extern void reprecalc_world(world_fmt *world, long that);
#if defined(MPI) && !defined(PARALIO)
extern void print_marginal_like(float *temp, long *z, world_fmt * world);
#else
extern void print_marginal_like(char *temp, long *z, world_fmt * world);
#endif
extern void  print_marginal_order(char *buf, long *bufsize, world_fmt *world);

void  set_map_groups(long numpop, long t, long size, longpair *map, long lastold);
//boolean bayes_accept (MYREAL newval, MYREAL oldval, MYREAL heat, MYREAL hastingsratio);
void bayes_print_accept(FILE * file,  world_fmt *world);
//long bayes_update (world_fmt * world);

// function pointers for the bayes machinery

//MYREAL (**log_prior_ratio) (MYREAL, MYREAL, bayes_fmt *, long);
 logpriorratioptr * log_prior_ratio;
//MYREAL (*log_prior_ratiotheta) (MYREAL, MYREAL, bayes_fmt *, long);
//MYREAL (*log_prior_ratiomig) (MYREAL, MYREAL, bayes_fmt *, long);
//MYREAL (*log_prior_ratiospec) (MYREAL, MYREAL, bayes_fmt *, long);
//MYREAL (*log_prior_ratiorate) (MYREAL, MYREAL, bayes_fmt *, long);
//MYREAL (*log_prior_ratioall) (MYREAL, MYREAL, bayes_fmt *, long);

//MYREAL (**log_prior) (world_fmt *, long);//for heating
 logpriorptr * log_prior;
//MYREAL (*log_prior_theta) (world_fmt *, long);//for heating
//MYREAL (*log_prior_mig) (world_fmt *, long);//for heating
//MYREAL (*log_prior_spec) (world_fmt *, long);//for heating
//MYREAL (*log_prior_rate) (world_fmt *, long);//for heating
// for scaling multiple loci correctly
//
//MYREAL (**log_prior_1) (world_fmt *, long, MYREAL);
 logprior1ptr * log_prior_1;
// MYREAL (*log_prior_theta1) (world_fmt *, long, MYREAL);
// MYREAL (*log_prior_mig1) (world_fmt *, long, MYREAL);
// MYREAL (*log_prior_spec1) (world_fmt *, long, MYREAL);
// MYREAL (*log_prior_rate1) (world_fmt *, long, MYREAL);
//
//MYREAL (**propose_new) (MYREAL, long, bayes_fmt *, MYREAL *);
 profuncptr * propose_new;
// MYREAL (*propose_newtheta) (MYREAL, long, bayes_fmt *, MYREAL *);
// MYREAL (*propose_newmig) (MYREAL,  long, bayes_fmt *, MYREAL *);
// MYREAL (*propose_newspec) (MYREAL,  long, bayes_fmt *, MYREAL *);
// MYREAL (*propose_newrate) (MYREAL,  long, bayes_fmt *, MYREAL *);

//MYREAL (**hastings_ratio) (MYREAL , MYREAL, MYREAL, MYREAL, bayes_fmt *, long );
 hastratioptr * hastings_ratio;
// MYREAL (*hastings_ratiotheta) (MYREAL , MYREAL, MYREAL, MYREAL, bayes_fmt *, long );
// MYREAL (*hastings_ratiomig) (MYREAL , MYREAL, MYREAL, MYREAL, bayes_fmt *, long );
// MYREAL (*hastings_ratiospec) (MYREAL , MYREAL, MYREAL, MYREAL, bayes_fmt *, long );
// MYREAL (*hastings_ratiorate) (MYREAL , MYREAL, MYREAL, MYREAL, bayes_fmt *, long );


MYINLINE MYREAL probWait(long *lines,MYREAL *locallparam, long numpop, long numpop2,long skyz);
MYINLINE MYREAL skyprobWait(world_fmt *world, long *lines, MYREAL *locallparam, long numpop, long numpop2, long skyz);

void calculate_credibility_interval(world_fmt * world, long locus);
void calc_hpd_credibility(world_fmt *world,long locus, long numpop2, long numparam);
void bayes_combine_loci(world_fmt * world);
void print_locus_histogram_header(FILE *bayesfile, MYREAL *deltas, char *custm2, long numpop, long numparam, boolean usem, boolean mu, world_fmt *world);
void print_locus_histogram(FILE *bayesfile, world_fmt * world, long locus, long numparam);
void print_loci_histogram(FILE *bayesfile, world_fmt * world, long locus, long numparam);
void bayes_set_param(MYREAL *param, MYREAL newparam, long which, char *custm2, long numpop);
void adjust_bayes_min_max(world_fmt* world, MYREAL **mini, MYREAL **maxi, MYREAL **adjmini, MYREAL **adjmaxi);
void bayes_progress(world_fmt *world, long ten);
void init_global_function_arrays( long numparam);
void destroy_global_function_arrays(void);
void which_prior (prior_fmt *bayes_priors,  long numparam);
void adjust_llocalparam(MYREAL *locallparam, world_fmt * world,  long locus,  long skyz);
MYREAL probg_treetimesCOMPLEX(world_fmt *world);
MYREAL scaling_prior(world_fmt *world, long numparam, MYREAL val);
void calculate_plotter_priors(world_fmt *world);
void print_bayes_verbose(long which, world_fmt *world, MYREAL newparam, boolean success);
MYREAL uniform_proposal(long which, world_fmt * world, MYREAL *oldparam, boolean *success);
void traverse_adjust(node *theNode, MYREAL new_old_ratio);
void recalc_timelist_1 (world_fmt * world, MYREAL new_old_ratio);
void bayes_save_parameter(world_fmt *world, long pnum, long step);
MYINLINE  void select_prior_param(int selector, long i, bayes_fmt *bayes, prior_fmt *prior);
void construct_param_hist(world_fmt *world, long locus, long npa, long pa, long numbin, MYREAL *mini,
                          MYREAL *maxi, double **results,
                          long *total, MYREAL *themean, MYREAL *thestd);
void construct_locus_histogram(world_fmt *world, long locus, MYREAL *mini, MYREAL *maxi, double **results);
void find_bayes_min_max(world_fmt* world, MYREAL **mini, MYREAL **maxi, MYREAL **adjmaxi);
void print_param_order(char **buf, long *bufsize, long *allocbufsize, world_fmt *world, long numparam);
void print_bayes_credibility(FILE *file, MYREAL *cred, long numpop);
void     find_bayes_mode(world_fmt *world, long locus, long numparam);
//##
#ifdef ZNZ
void	  print_bayes_tofile(znzFile mdimfile, MYREAL *params, bayes_fmt *bayes, world_fmt *world, long numpop2, long mdiminterval, long locus, long replicate, long step);
#else
void	  print_bayes_tofile(FILE *mdimfile, MYREAL *params, bayes_fmt *bayes, world_fmt *world, long numpop2, long mdiminterval, long locus, long replicate, long step);
#endif
#ifdef DEBUG
extern void print_bf_values(world_fmt * world);
#endif

void recalc_timelist (world_fmt * world, MYREAL new_ratio , MYREAL old_ratio);
void handle_averages(char representant, MYREAL newparam, long which, long numpop, char *custm2,  MYREAL * param);

MYREAL probg_treetimes_local(world_fmt* world, timelist_fmt * treetimes);

/// \brief set param within min max
///
void check_min_max_param(MYREAL *param0, MYREAL minparam, MYREAL maxparam)
{
    if(*param0 > maxparam)
        *param0 = maxparam;
    else
    {
        if(*param0 < minparam)
            *param0 = minparam;
    }
}

void init_global_function_arrays( long numparam)
{
    // windows pointer mockery needed here ?
    log_prior_ratio = (logpriorratioptr *) mycalloc(numparam , sizeof( MYREAL (*)  (MYREAL, MYREAL, bayes_fmt *, long)));
    log_prior = (logpriorptr *) mycalloc(numparam , sizeof( MYREAL (*)  (world_fmt *, long)));
    log_prior_1 = (logprior1ptr * ) mycalloc(numparam , sizeof( MYREAL (*)  (world_fmt *, long, MYREAL)));
    //  propose_new = (MYREAL (*) (MYREAL, long, bayes_fmt *, MYREAL *)) calloc(numparam , sizeof( MYREAL (*)  (MYREAL, long, bayes_fmt *, MYREAL *)));
    propose_new = (profuncptr *) mycalloc(numparam,sizeof(profuncptr));
    hastings_ratio = (hastratioptr *) mycalloc(numparam , sizeof(hastratioptr));
}

void destroy_global_function_arrays()
{
  myfree(log_prior_ratio);
  myfree(log_prior);
  myfree(log_prior_1);
  myfree(propose_new);
  myfree(hastings_ratio);
}

/// \brief Decide which prior distribution for the THETA parameter to use
/// Decide which prior distribution to use for THETA: the functionpointers propose_newparam_x will hold
/// either the Exponential prior distribution or a Uniform prior distribution or .. other priors ..
/// each prior distribution has its own specific hastings ratio that will be calculated in the
/// function ptr hastings_ratio
void which_prior (prior_fmt *bayes_priors,  long numparam)
{
    long kind;
     long i;
    init_global_function_arrays(numparam);
    for (i=0; i<numparam; i++)
    {
        kind = bayes_priors[i].kind;
        switch (kind)
        {
            case UNIFORMPRIOR:
                log_prior_ratio[i] = (MYREAL (*) (MYREAL, MYREAL, bayes_fmt *, long)) log_prior_ratio_uni;
                log_prior[i] = (MYREAL (*) (world_fmt *, long)) log_prior_uni;
                log_prior_1[i] = (MYREAL (*) (world_fmt *, long, MYREAL)) log_prior_uni1;
                propose_new[i] = (MYREAL (*) (MYREAL, long, world_fmt *, MYREAL *)) propose_uni_newparam;
                hastings_ratio[i] = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, bayes_fmt *, long)) hastings_ratio_uni;
                break;
            case EXPPRIOR:
                log_prior_ratio[i] = (MYREAL (*) (MYREAL,  MYREAL, bayes_fmt *, long)) log_prior_ratio_exp;
                log_prior[i] = (MYREAL (*) (world_fmt *, long)) log_prior_exp;
                log_prior_1[i] = (MYREAL (*) (world_fmt *,  long, MYREAL)) log_prior_exp1;
                propose_new[i] = (MYREAL (*) (MYREAL,  long, world_fmt *, MYREAL * )) propose_exp_newparam;
                hastings_ratio[i] = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, bayes_fmt *, long )) hastings_ratio_exp;
                break;
            case WEXPPRIOR:
                log_prior_ratio[i] = (MYREAL (*) (MYREAL,  MYREAL, bayes_fmt *, long)) log_prior_ratio_wexp;
                log_prior[i] = (MYREAL (*) (world_fmt *, long)) log_prior_wexp;
                log_prior_1[i] = (MYREAL (*) (world_fmt *,  long, MYREAL)) log_prior_wexp1;
                propose_new[i] = (MYREAL (*) (MYREAL,  long, world_fmt *, MYREAL * )) propose_expb_newparam;
                hastings_ratio[i] = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, bayes_fmt *, long)) hastings_ratio_expb;
                break;
            case MULTPRIOR:
                log_prior_ratio[i] = (MYREAL (*) (MYREAL,  MYREAL,bayes_fmt *, long)) log_prior_ratio_mult;
                log_prior[i] = (MYREAL (*) (world_fmt * , long)) log_prior_mult;
                log_prior_1[i] = (MYREAL (*) (world_fmt * ,  long, MYREAL)) log_prior_mult1;
                propose_new[i] = (MYREAL (*) (MYREAL,  long, world_fmt *, MYREAL * )) propose_mult_newparam;
                hastings_ratio[i] = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, bayes_fmt *, long )) hastings_ratio_mult;
                break;
            case GAMMAPRIOR:
                log_prior_ratio[i] = (MYREAL (*) (MYREAL,  MYREAL, bayes_fmt *, long)) log_prior_ratio_gamma;
                log_prior[i] = (MYREAL (*) (world_fmt *, long)) log_prior_gamma;
                log_prior_1[i] = (MYREAL (*) (world_fmt *,  long, MYREAL)) log_prior_gamma1;
                propose_new[i] = (MYREAL (*) (MYREAL,  long, world_fmt *, MYREAL * )) propose_gamma_newparam;
                hastings_ratio[i] = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, bayes_fmt *, long)) hastings_ratio_gamma;
                break;
            case NORMALPRIOR:
                log_prior_ratio[i] = (MYREAL (*) (MYREAL,  MYREAL, bayes_fmt *, long)) log_prior_ratio_normal;
                log_prior[i] = (MYREAL (*) (world_fmt *, long)) log_prior_normal;
                log_prior_1[i] = (MYREAL (*) (world_fmt *,  long, MYREAL)) log_prior_normal1;
                propose_new[i] = (MYREAL (*) (MYREAL,  long, world_fmt *, MYREAL * )) propose_normal_newparam;
                hastings_ratio[i] = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, bayes_fmt *, long)) hastings_ratio_normal;
                break;
            default:
                error("Prior distribution hookup failed");
                //break;
        }
    }
}

/// \brief fixes start parameter that are outside of the bayesian mini/maxi values
///
/// Fixes parameter start values that are outside of the Bayesian minimum and maximum
/// values. Only called in a Bayesian analysis
void bayes_check_and_fix_param(world_fmt *world)
{
     long i;
    //    long frompop, topop;
    //MYREAL theta;
    const  long npp = world->numpop2 + ( long) world->bayes->mu + 2 * world->species_model_size + world->grownum;
    MYREAL *maxparam = world->bayes->maxparam;
    MYREAL *minparam = world->bayes->minparam;
    
    for (i=0; i< npp; i++)
    {
        check_min_max_param(&world->param0[i],minparam[i],maxparam[i]);
    }
}


void adjust_llocalparam(MYREAL *locallparam, world_fmt * world,  long locus,  long skyz)
{
     long msta;
     long msto;
     long r;
     long j;
     long numpop = world->numpop;
     long numpop2 = world->numpop2;
    boolean usem = world->options->usem;
    
    const MYREAL *geo = world->data->geo;
    const MYREAL *lgeo = world->data->lgeo;
    const MYREAL mu_rate = world->options->mu_rates[locus];
    const MYREAL lmu_rate = world->options->lmu_rates[locus];
    const MYREAL  inheritance = world->options->inheritance_scalars[locus];
    const MYREAL linheritance = log(inheritance);
    MYREAL *param0 = world->param0;
    MYREAL *localparam;
    MYREAL *sm;
    MYREAL *skyparam;
    MYREAL pk;
     long s = skyz * (numpop2);
    localparam = locallparam + numpop2;
    sm = localparam + numpop;
    skyparam = sm + numpop;
    for (r = 0; r < numpop; r++)
    {
        locallparam[r] = LOG2 - (log(skyparam[s+r]*param0[r]) + linheritance + lmu_rate);
        localparam[r] = -1. / (skyparam[s+r]*param0[r] * mu_rate * inheritance);
        msta = world->mstart[r];
        msto = world->mend[r];
        sm[r] = 0.0;
        for (j = msta; j < msto; j++)
        {
            if (param0[j] > 0.0)
	      {
		pk = usem ? param0[j] : param0[j] / param0[r];
		sm[r] -= skyparam[s+j] * geo[j] * pk / mu_rate;
		locallparam[j] = LOG(skyparam[s+j] * pk) + lgeo[j] - lmu_rate;
	      }
        }
    }
}

extern long m2mmm(long frompop, long topop, long numpop);
/*MYREAL probg_treetimes(world_fmt *world)
{
  for(i=1; i<T-1;i++)
    {
      tli1 = &tl[i-1];
      t0 = tli1->age;
      tli = &tl[i];
      t1 = tli->age; 
      k = tli->lineages;
      type = tli->eventnode->type;
      
      deltatime = t0 - t1; // should be negative
      sumprob += deltatime * probWait(k,locallparam, world->numpop, world->numpop2, skyz);//theta and M
      sumprob += deltatime * wait_event_species(world,tli,tli1->age,tli->age,TRUE,&eventprob);
    }
}
*/
//MYREAL probg_treetimesSIMPLE(world_fmt* world)
MYREAL probg_treetimes(world_fmt* world)
{
    const MYREAL *geo = world->data->geo;
    //const MYREAL *lgeo = world->data->lgeo;
    const  long numpop = world->numpop;
    const  long numpop2 = world->numpop2;
    //const long locus = world->locus;
    const  long npp = numpop2 + ( long) world->bayes->mu;
    const  long nppall = npp + 2 * world->species_model_size; //OK, do not add growth here!
    species_fmt *s;
    long i;
     long pop;
    MYREAL t0;
    MYREAL t1;
    double tx;
    MYREAL *param0 = world->param0;
    const long T = world->treetimes->T;
    vtlist *tl = world->treetimes->tl;
    vtlist *tli;
    vtlist *tli1;
    MYREAL mu;
    MYREAL sigma;
    MYREAL sumprob = 0.0;
    MYREAL deltatime = 0.0;
    MYREAL deltatime2 = 0.0;
    MYREAL waitprobcoal=0.0;
    MYREAL waitprobmig=0.0;
    MYREAL waitprob_spec = 0.0;
    MYREAL eventprob=0.0;
    long *k;
    char type;
    MYREAL pk;
    boolean usem = world->options->usem;
    long pop2;
    long kpop;
    long msta;
    long msto;
    double kpopmurate;
    double sum;
    const MYREAL mu_rate = world->options->mu_rates[world->locus];
    double mlalpha = world->mlalpha;
    //double mlinheritance = world->mlinheritance;
    //assert(mlinheritance==2.0);
    boolean has_mlalpha = mlalpha<1.0;
    //boolean hasnot_mlalpha = !has_mlalpha;
    double x;
    double g;
    long * growpops = world->options->growpops;
    double *growth=NULL;
    if (world->has_growth)
      {
	growth = world->growth;
      }

    //double alphapart = exp(LGAMMA(1.0+mlalpha));
    //const MYREAL lmu_rate = world->options->lmu_rates[world->locus];
    for(i=1; i<T-1;i++)
    {
      tli1 = &tl[i-1];
      t0 = tli1->age;
      tli = &tl[i];
      t1 = tli->age; 
      k = tli->lineages;// the lineages are filled in the actual timeslice: e.g. tli1=time=0=lasttip lineages are not
      // build up yet and may look like [4,5], if tli is a 'm' event then the lineages at the event are [5,5] and
      // after that, say, [6,4] [CHECK ON OTHER TIME INTERVALS -- not done yet]
      type = tli->eventnode->type;

      deltatime = t0 - t1; // should be negative
      deltatime2 = t1 - t0; // should be positive
      if(has_mlalpha)
	{
	  deltatime = -pow(deltatime2,mlalpha);
	  //fprintf(stderr,"%i> -(t1-t0)^a=-(%f)^%f=%f\n",myID,t1-t0,mlalpha,deltatime);
	}
      if(type == 't')
        {
	  //eventprob = 0.0;
	  //  waitprob = 0.0;
	  continue;
        }
      else
        {
	  waitprobcoal = 0.0;
	  waitprobmig = 0.0;
	  waitprob_spec = 0.0;
	  // coalescence with fixed or exp growing population size
	  for(pop=0;pop<numpop;pop++)
	    {
	      kpop = k[pop];
	      if (kpop>1)
		{
		  if (world->has_growth && growpops[pop]!=0 && (fabs(growth[growpops[pop]-1])>EPSILON))
		    {
		      x = mu_rate * param0[pop];
		      g = growth[growpops[pop]-1];
		      //		  waitprob += kpop * (kpop - 1) / (x * exp(-g * t1));
		      waitprobcoal +=  -kpop * (kpop - 1) * (exp(g * (t1)) - exp(g * (t0)))/(x * g);
		      if (isnan(waitprobcoal))
			waitprobcoal = (double) -HUGE;
		      assert(!isnan(waitprobcoal));
		    }
		  else
		    {
		      waitprobcoal += deltatime * kpop * (kpop - 1) / (mu_rate * param0[pop]);
		      assert(!isnan(waitprobcoal));
		    }
		}
	      msta = world->mstart[pop];
	      msto = world->mend[pop];
	      kpopmurate = kpop/mu_rate;

	      sum=0.0;
	      if (usem)
		{
		  for (pop2 = msta; pop2 < msto; pop2++)
		    {
		      pk = param0[pop2] * geo[pop2];
		      sum += kpopmurate * pk; 
		    }
		}
	      else
		{
		  for (pop2 = msta; pop2 < msto; pop2++)
		    {
		      pk = param0[pop2] / param0[pop]; 
		      sum += kpopmurate * geo[pop2] * pk;
		    }
		}
	      waitprobmig += sum;
	    }
	  waitprobmig *= deltatime;
	  assert(!isnan(waitprobmig));
	  if (world->has_speciation)
	    {
	      long sindex;
	      for(sindex=npp; sindex<nppall; sindex++)
		{
		  s = get_which_species_model(sindex, world->species_model, world->species_model_size);
		  if ((s == NULL) || (sindex != s->paramindex_mu))
		    {
		      continue;
		    }
		  else
		    {
		      mu = param0[s->paramindex_mu];
		      sigma = param0[s->paramindex_sigma];
		      kpop = k[s->to];
		      if (kpop>0)
			waitprob_spec += kpop * (*log_prob_wait_speciate)(t0,t1,mu,sigma,s);
		      assert(!isnan(waitprob_spec));
		    }
		}
	    }
	  long xpop = tli->eventnode->actualpop;
	  switch(type)
	    {
	    case 'i':
	      if (world->has_growth && growpops[xpop]!= 0 && (fabs(growth[growpops[xpop]-1])>EPSILON))
		{
		  x = mu_rate * param0[xpop];
		  g = growth[growpops[xpop]-1];
		  eventprob = LOG2 + g * (t1) - LOG(x);
		  assert(!isnan(eventprob));
		  //eventprob =  LOG2 + log((exp(g * t1)/(x * g));
		}
	      else
		{
		  x = mu_rate * param0[xpop];
		  eventprob = LOG2 - log(x);
		  assert(!isnan(eventprob));
		}
	      break;
	    case 'm':
	      if (usem)
		{
		  eventprob = log(param0[m2mmm(tli->eventnode->pop,
					     tli->eventnode->actualpop, (long) numpop)]);
		}
	      else
		{
		  eventprob = log(param0[m2mmm(tli->eventnode->pop,
					     tli->eventnode->actualpop, (long) numpop)]/param0[tli->eventnode->actualpop]);
		}
	      assert(!isnan(eventprob));
	      break;
	    case 'd':
	      s = get_fixed_species_model(tli->eventnode->pop,tli->eventnode->actualpop, world->species_model, world->species_model_size);
	      mu = param0[s->paramindex_mu];
	      sigma = param0[s->paramindex_sigma];
	      if (world->species_model_dist == NORMALSHORTCUT_DIST)
	      	tx = (t0+t1)/2.0;
	      else
		tx = t1;
	      eventprob = (*log_point_prob_speciate)(tx,mu,sigma,s);//-log(k[s->to]) no 1/ necessary;
	      assert(!isnan(eventprob));
	      break;
	    default:
	      error("point prob failed");
	    }
	}
      //sumprob += deltatime * (waitprob + waitprob_spec) + eventprob;
      if(has_mlalpha)
	{
	   //mittag-leffler
	  double pw = waitprobcoal + waitprobmig + waitprob_spec;
#ifdef WINDOWS
	  _Dcomplex pwc = {pw,0.0};
	  _Dcomplex mlfc = mittag_leffler(mlalpha, mlalpha, pwc);
#else
	  MYCOMPLEX mlfc = mittag_leffler(mlalpha, mlalpha, pw);
#endif
	  double mlfval = creal(mlfc);    
	  double mittag_result = LOG(deltatime2)*(mlalpha-1.0) + mlfval + eventprob;
	  sumprob += mittag_result;
	  assert(!isnan(sumprob));
	  //fprintf(stderr,"%i> locus %li sumprob=%f\n",myID,locus,sumprob);
	}
	  else
	{
	  //#ifdef DEBUGWEIRD
	  if (isinf(waitprobcoal) && isinf(eventprob))
	    {
	      waitprobcoal = 0.0;
	      eventprob = 0.0;
	    }
	  //#endif
	  sumprob += waitprobcoal + waitprobmig + waitprob_spec + eventprob;
	  assert(!isnan(sumprob));
	}
      //printf("%i> progtreetimes(): sumprob=%f waitcoal=%f + waitmig=%f + waitspec=%f + eventprob=%f\n",myID, sumprob,waitprobcoal,waitprobmig,waitprob_spec, eventprob);
    }
    assert(!isnan(sumprob));
    return sumprob;
}
//newlocal
// used for assignment calculation [see world.c: updating(), speciate.c newtree_update()]
MYREAL probg_treetimes_local(world_fmt* world, timelist_fmt * treetimes)
{
    const MYREAL *geo = world->data->geo;
    //const MYREAL *lgeo = world->data->lgeo;
    const  long numpop = world->numpop;
    const  long numpop2 = world->numpop2;
    //const long locus = world->locus;
    const  long npp = numpop2 + ( long) world->bayes->mu;
    const  long nppall = npp + 2 * world->species_model_size; //OK, do not add growth here!
    species_fmt *s;
    long i;
     long pop;
    MYREAL t0;
    MYREAL t1;
    double tx;
    MYREAL *param0 = world->param0;
    const long T = treetimes->T;
    vtlist *tl = treetimes->tl;
    vtlist *tli;
    vtlist *tli1;
    MYREAL mu;
    MYREAL sigma;
    MYREAL sumprob = 0.0;
    MYREAL deltatime = 0.0;
    MYREAL deltatime2 = 0.0;
    MYREAL waitprobcoal=0.0;
    MYREAL waitprobmig=0.0;
    MYREAL waitprob_spec = 0.0;
    MYREAL eventprob=0.0;
    long *k;
    char type;
    MYREAL pk;
    boolean usem = world->options->usem;
    long pop2;
    long kpop;
    long msta;
    long msto;
    double kpopmurate;
    double sum;
    const MYREAL mu_rate = world->options->mu_rates[world->locus];
    double mlalpha = world->mlalpha;
    //double mlinheritance = world->mlinheritance;
    //assert(mlinheritance==2.0);
    boolean has_mlalpha = mlalpha<1.0;
    //boolean hasnot_mlalpha = !has_mlalpha;
    double x;
    double g;
    long * growpops = world->options->growpops;
    double *growth=NULL;
    if (world->has_growth)
      {
	growth = world->growth;
      }

    //double alphapart = exp(LGAMMA(1.0+mlalpha));
    //const MYREAL lmu_rate = world->options->lmu_rates[world->locus];
    for(i=1; i<T-1;i++)
    {
      tli1 = &tl[i-1];
      t0 = tli1->age;
      tli = &tl[i];
      t1 = tli->age; 
      k = tli->lineages;// the lineages are filled in the actual timeslice: e.g. tli1=time=0=lasttip lineages are not
      // build up yet and may look like [4,5], if tli is a 'm' event then the lineages at the event are [5,5] and
      // after that, say, [6,4] [CHECK ON OTHER TIME INTERVALS -- not done yet]
      if (tli->eventnode == NULL)
	type = tli->type;
      else
	type = tli->eventnode->type;

      deltatime = t0 - t1; // should be negative
      deltatime2 = t1 - t0; // should be positive
      if(has_mlalpha)
	{
	  deltatime = -pow(deltatime2,mlalpha);
	  //fprintf(stderr,"%i> -(t1-t0)^a=-(%f)^%f=%f\n",myID,t1-t0,mlalpha,deltatime);
	}
      if(type == 't')
        {
	  //eventprob = 0.0;
	  //  waitprob = 0.0;
	  continue;
        }
      else
        {
	  waitprobcoal = 0.0;
	  waitprobmig = 0.0;
	  waitprob_spec = 0.0;
	  // coalescence with fixed or exp growing population size
	  for(pop=0;pop<numpop;pop++)
	    {
	      kpop = k[pop];
	      if (kpop>1)
		{
		  if (world->has_growth && growpops[pop]!=0 && (fabs(growth[growpops[pop]-1])>EPSILON))
		    {
		      x = mu_rate * param0[pop];
		      g = growth[growpops[pop]-1];
		      //		  waitprob += kpop * (kpop - 1) / (x * exp(-g * t1));
		      waitprobcoal +=  -kpop * (kpop - 1) * (exp(g * (t1)) - exp(g * (t0)))/(x * g);
		      if (isnan(waitprobcoal))
			waitprobcoal = (double) -HUGE;
		      assert(!isnan(waitprobcoal));
		    }
		  else
		    {
		      waitprobcoal += deltatime * kpop * (kpop - 1) / (mu_rate * param0[pop]);
		      assert(!isnan(waitprobcoal));
		    }
		}
	      msta = world->mstart[pop];
	      msto = world->mend[pop];
	      kpopmurate = kpop/mu_rate;

	      sum=0.0;
	      if (usem)
		{
		  for (pop2 = msta; pop2 < msto; pop2++)
		    {
		      pk = param0[pop2] * geo[pop2];
		      sum += kpopmurate * pk; 
		    }
		}
	      else
		{
		  for (pop2 = msta; pop2 < msto; pop2++)
		    {
		      pk = param0[pop2] / param0[pop]; 
		      sum += kpopmurate * geo[pop2] * pk;
		    }
		}
	      waitprobmig += sum;
	    }
	  waitprobmig *= deltatime;
	  assert(!isnan(waitprobmig));
	  if (world->has_speciation)
	    {
	      long sindex;
	      for(sindex=npp; sindex<nppall; sindex++)
		{
		  s = get_which_species_model(sindex, world->species_model, world->species_model_size);
		  if ((s == NULL) || (sindex != s->paramindex_mu))
		    {
		      continue;
		    }
		  else
		    {
		      mu = param0[s->paramindex_mu];
		      sigma = param0[s->paramindex_sigma];
		      kpop = k[s->to];
		      if (kpop>0)
			waitprob_spec += kpop * (*log_prob_wait_speciate)(t0,t1,mu,sigma,s);
		      assert(!isnan(waitprob_spec));
		    }
		}
	    }
	  long from;
	  long to;
	  if (tli->eventnode==NULL)
	    {
	      from = tli->from;
	      to = tli->to;
	    }
	  else
	    {
	    from = tli->eventnode->pop;
	    to = tli->eventnode->actualpop;
	    }
	  switch(type)
	    {
	    case 'i':
	      if (world->has_growth && growpops[from]!= 0 && (fabs(growth[growpops[from]-1])>EPSILON))
		{
		  x = mu_rate * param0[from];
		  g = growth[growpops[from]-1];
		  eventprob = LOG2 + g * (t1) - LOG(x);
		  assert(!isnan(eventprob));
		  //eventprob =  LOG2 + log((exp(g * t1)/(x * g));
		}
	      else
		{
		  x = mu_rate * param0[from];
		  eventprob = LOG2 - log(x);
		  assert(!isnan(eventprob));
		}
	      break;
	    case 'm':
	      if (usem)
		{
		  eventprob = log(param0[m2mmm(from,
					     to, (long) numpop)]);
		}
	      else
		{
		  eventprob = log(param0[m2mmm(from,
					     to, (long) numpop)]/param0[tli->eventnode->actualpop]);
		}
	      assert(!isnan(eventprob));
	      break;
	    case 'd':
	      s = get_fixed_species_model(from,to, world->species_model, world->species_model_size);
	      mu = param0[s->paramindex_mu];
	      sigma = param0[s->paramindex_sigma];
	      if (world->species_model_dist == NORMALSHORTCUT_DIST)
	      	tx = (t0+t1)/2.0;
	      else
		tx = t1;
	      eventprob = (*log_point_prob_speciate)(tx,mu,sigma,s);//-log(k[s->to]) no 1/ necessary;
	      assert(!isnan(eventprob));
	      break;
	    default:
	      error("point prob failed");
	    }
	}
      //sumprob += deltatime * (waitprob + waitprob_spec) + eventprob;
      if(has_mlalpha)
	{
	   //mittag-leffler
	  double pw = waitprobcoal + waitprobmig + waitprob_spec;
#ifdef WINDOWS
	  _Dcomplex pwc = {pw,0.0};
	  _Dcomplex mlfc = mittag_leffler(mlalpha, mlalpha, pwc);
#else
	  MYCOMPLEX mlfc = mittag_leffler(mlalpha, mlalpha, pw);
#endif
	  double mlfval = creal(mlfc);    
	  double mittag_result = LOG(deltatime2)*(mlalpha-1.0) + mlfval + eventprob;
	  sumprob += mittag_result;
	  assert(!isnan(sumprob));
	  //fprintf(stderr,"%i> locus %li sumprob=%f\n",myID,locus,sumprob);
	}
	  else
	{
	  sumprob += waitprobcoal + waitprobmig + waitprob_spec + eventprob;
	  assert(!isnan(sumprob));
	}
      //printf("%i> progtreetimes(): sumprob=%f waitcoal=%f + waitmig=%f + waitspec=%f + eventprob=%f\n",myID, sumprob,waitprobcoal,waitprobmig,waitprob_spec, eventprob);
    }
    assert(!isnan(sumprob));
    return sumprob;
}
//end new local

/// \brief Calculate Prob(g|param)Prob(D|G) from world->treetimes.
///
/// Calculate Prob(g|param)Prob(D|G) from world->treetimes
// uses as input    vtlist *tl = world->treetimes->tl;
// uses as input    world->treetimes->T
// to allow for custom p(g*|theta) calculation, but currently not needed
//MYREAL probg_treetimesX(world_fmt *world, vtlist *tl, long T);
MYREAL probg_treetimesCOMPLEX(world_fmt *world)
//{
//  return probg_treetimesX(world,world->treetimes->tl,world->treetimes->T);
//}
//MYREAL probg_treetimesX(world_fmt *world, vtlist *tl, long T)
{
  warning("probg_treetimesCOMPLEX needs mlalpha treatement [not done yet]");
    long numpop = world->numpop;
    long numpop2 = world->numpop2;
    long locus = world->locus;
    MYREAL like = world->likelihood[world->G];
    boolean notsky = (!world->options->skyline);
    // each locus is run with its native theta
    // we do not need to adjust for inheritance
    // before the combining over loci.
    // const MYREAL  inheritance = 1.0;
    // const MYREAL linheritance = 0.0;
    long T = world->treetimes->T;
    vtlist *tl = world->treetimes->tl; // perhaps comment this out to
    // making probgtreetimes more general
    // so that it can be used in the proposal structure in MCMC1.c
    vtlist *tli = &(tl[0]);
    vtlist *tli1= NULL;
    
    long i;
    long tlif = tli->from;
    long tlit = tli->to;
    
    MYREAL deltatime = tl[0].age;
    MYREAL sumprob = 0.;
    MYREAL eventprob=0.;
    //xcode MYREAL k;
    
    //const MYREAL *geo = world->data->geo;
    //const MYREAL *lgeo = world->data->lgeo;
    MYREAL *param0 = world->param0;
    MYREAL *locallparam;
    MYREAL *localparam;
    MYREAL *sm;
    MYREAL *skyparam;
    MYREAL skytime;
     long skyz;
    const MYREAL lmu_rate = world->options->lmu_rates[locus];
    MYREAL pw;
     long addition=0;
     long skypnum = world->timeelements*(numpop2+addition);
    MYREAL specw;
    locallparam = (MYREAL *) mycalloc ((numpop2 + 2 * numpop + skypnum), sizeof (MYREAL));
    
    // pointers into the localparam vector
    localparam = locallparam + numpop2;
    sm = localparam + numpop;
    skyparam = sm + numpop;
    // fill localparam vector with data
    memcpy (localparam, param0, sizeof (MYREAL) * (size_t) numpop);
    memcpy (skyparam, world->timek, sizeof (MYREAL) * (size_t) skypnum);
    skyz = 0;
    adjust_llocalparam(locallparam, world, locus, skyz);
    if(tl[0].eventnode->type == 't')
    {
        eventprob = 0.0;
        //deltatime = 0.0; //the first tip might be in the past,
        // this makes sure that if all sample have dual date > 0 that
        // we do not calculate strange probabilities,
        // this is a problem only for the first entry.
    }
    else
    {
        // k = tli->lineages[tlit];
        //eventprob = ((tlif==tlit) ? (locallparam[tlif]) : locallparam[mm2m(tlif,tlit, numpop)]);
        if(notsky)
        {
            if (tlif == tlit)
                eventprob = locallparam[tlif];
            else if (tl[0].eventnode->type == 'd')
                eventprob = 0.0;
            else
	      eventprob = locallparam[mm2m(tlif,tlit, (long) numpop)];
            sumprob = deltatime * probWait(tli->lineages, locallparam, (long) numpop, (long) numpop2, (long) skyz) + eventprob;
        }
        else
        {
            warning("check eventprob in bayes.c:466");
            eventprob = ((tlif==tlit) ? (locallparam[tlif]) : (tl[0].eventnode->type == 'd' ? 0.0 : locallparam[mm2m(tlif,tlit, (long) numpop)]));
            sumprob = deltatime * skyprobWait(world, tli->lineages, locallparam, (long) numpop, (long) numpop2, (long) skyz) + eventprob;
        }
        if (world->has_speciation)
        {
	  if (tl[0].eventnode->type == 'd')
	    {
	      specw = wait_event_species(world, tli, 0, tli->age, TRUE, &eventprob);
	      specw += eventprob;
	    }
	  else
	    {
	      specw = wait_event_species(world, tli, 0, tli->age, FALSE, &eventprob);
	    }
	  sumprob += specw; 
        }
    }
    skytime = world->times[1];
    // is is something
    for(i=1; i<T-1;i++)
    {
        tli = &tl[i];
        tli1 = &tl[i-1];
        //goes through all tipnodes if they are contemporary
        if(tli->age==0.0)
            continue;
        tlif = tli->from;
        tlit = tli->to;
        
        if(tli->eventnode ==NULL)
            continue;
       
        if(tli->eventnode->type == 't')
        {
            eventprob = 0.0;
        }
        else
        {
            // k = tli->lineages[tlit];
	  eventprob = ((tlif==tlit) ? (locallparam[tlif]) :  (tli->eventnode->type == 'd' ? 0.0 : locallparam[mm2m(tlif,tlit, (long) numpop)]));
        }
        if(skytime > tli->age || notsky)
        {
            deltatime = (tli->age - tli1->age);
            if(notsky)
            {
	      pw = probWait(tli->lineages, locallparam, (long) numpop, (long) numpop2, (long) skyz);
            }
            else
            {
	      pw = skyprobWait(world, tli->lineages, locallparam, (long) numpop, (long) numpop2, (long) skyz);
	      eventprob = ((tlif==tlit) ? (LOG2 - log(skyparam[skyz*(numpop2)+( long) tlit] * param0[tlit]) + lmu_rate) :  (tli->eventnode->type == 'd' ? 0.0 : locallparam[mm2m(tlif,tlit, (long) numpop)]));
            }
	    if(tli->eventnode->type != 'd')
	      {
		specw = wait_event_species(world, tli, tli1->age, tli->age, FALSE, &eventprob);
	      }
	    else
	      {
		specw = wait_event_species(world, tli, tli1->age, tli->age, TRUE, &eventprob);
	      }
            pw *= deltatime;
            pw += specw;
            pw += eventprob;
            sumprob += pw;
        }
        else
        {
            deltatime = (skytime - tli1->age);
            specw = wait_event_species(world, tli, tli1->age, skytime, FALSE, &eventprob);
            pw = skyprobWait(world, tli->lineages, locallparam, (long) numpop, (long) numpop2, (long) skyz);
            pw *= deltatime;
            pw += specw;
            sumprob += pw;
            //up to the next sky boundary
            deltatime = (tli->age - skytime);
            pw = skyprobWait(world, tli->lineages, locallparam, (long) numpop, (long) numpop2, (long) skyz);
            if (tli->eventnode->type != 'd')
	      {
                eventprob = ((tlif==tlit) ?
                             (LOG2 - log(skyparam[skyz*(numpop2)+ ( long) tlit] * param0[tlit]) + lmu_rate)
                             : locallparam[mm2m(tlif,tlit, (long) numpop)]);
		specw = wait_event_species(world, tli, skytime, tli->age, FALSE, &eventprob);
	      }
            else
	      {
		specw = wait_event_species(world, tli, skytime, tli->age, TRUE, &eventprob);
	      }
            pw *= deltatime;
            pw += specw;
            pw += eventprob;
            sumprob += pw;
            skyz++;
            skytime = world->times[skyz+1];
            adjust_llocalparam(locallparam, world, locus, skyz);
        }
    //printf(" bayes.c:505: i=%li sumprob=%f eventprob=%f specw=%f pw=%f\n",i, sumprob, eventprob, specw, pw);
}
myfree(locallparam);
//printf("ln p(d|m)=%f+%f=%f\n",sumprobWait, like ,sumprob + like);
return /*-LOG(world->treetimes->T-1) +*/ (sumprob + like + calculate_prior(world));
}

///
/// calculates the waiting time part for prob(g|param)
/// ignoring speciation nodes
MYINLINE
MYREAL probWait(long *lines, MYREAL *locallparam, long numpop, long numpop2, long skyz)
{
    MYREAL *invtheta = locallparam + numpop2;
    MYREAL *msum = invtheta + numpop;
    MYREAL *skyparam = msum + numpop + skyz * (numpop2);//why was this numpop2+1 ?? replace this
    long j;
    register double line0 = (double) lines[0];
    double line01 = (line0 * line0) - line0;
    register double line1;
    double line11;
    register double line2;
    double line21;
    register double line3;
    double line31;
    MYREAL probm = 0.;
    MYREAL probth = 0.;
    warning("this needs to be called with complex and is disable for now");
    switch(numpop)
    {
        case 1:
            probth = line01 * invtheta[0]/skyparam[0];
            //should be zero! probm = msum[0] * line0;
            return probth; // + probm;
        case 2:
            line1 = (double) lines[1];
            line11 = (line1 * line1) - line1;
            probm = msum[0] * line0 + msum[1] * line1;
            probth = invtheta[0]/skyparam[0] * line01 + invtheta[1]/skyparam[1] * line11;
            return probth + probm;
        case 3:
            line1 = (double) lines[1];
            line2 = (double) lines[2];
            line11 = line1 * (line1 - 1.0) * invtheta[1];
            probm = msum[0] * line0 + msum[1] * line1 + msum[2] * line2;
            line21 = line2 * (line2 - 1.0) * invtheta[2];
            //    probth = invtheta[0]/skyparam[0] * line01 + invtheta[1]/skyparam[1] * line11 + invtheta[2]/skyparam[2] * line21;
            probth = invtheta[0] * line01 + line11 + line21;
            return probth + probm;
        case 4:
            line1 = (double) lines[1];
            line11 = (line1 * line1) - line1;
            line2 = (double) lines[2];
            line21 = (line2 * line2) - line2;
            line3 = (double) lines[3];
            line31 = (line3 * line3) - line3;
            probm = msum[0] * line0 + msum[1] * line1 + msum[2] * line2 + msum[3] * line3;
            probth = invtheta[0]/skyparam[0] * line01 + invtheta[1]/skyparam[1] * line11 + invtheta[2]/skyparam[2] * line21 + invtheta[3]/skyparam[3] * line31;
            return probth + probm;
        default:
            line1 = (double) lines[1];
            line11 = (line1 * line1) - line1;
            line2 = (double) lines[2];
            line21 = (line2 * line2) - line2;
            line3 = (double) lines[3];
            line31 = (line3 * line3) - line3;
            probm = msum[0] * line0 + msum[1] * line1 + msum[2] * line2 + msum[3] * line3;
            probth = invtheta[0]/skyparam[0] * line01 + invtheta[1]/skyparam[1] * line11 + invtheta[2]/skyparam[2] * line21 + invtheta[3]/skyparam[3] * line31;
            for(j=4; j < numpop; j++)
            {
                const double linex = (double) lines[j];
                probth += ((linex * linex) - linex) * invtheta[j]/skyparam[j];
                probm += linex * msum[j];
            }
            return probth + probm;
    }
}

///
/// calculates the waiting time part for prob(g|param) with the skyline plot
MYINLINE
MYREAL skyprobWait(world_fmt *world, long *lines, MYREAL *locallparam, long numpop, long numpop2, long skyz)
{
    MYREAL *invtheta = locallparam + numpop2;
    MYREAL *msum = invtheta + numpop;// ignored and recalculated
    MYREAL *skyparam = msum + numpop + skyz * (numpop2);
    long j;
    MYREAL probm = 0.;
    MYREAL probth = 0.;
    MYREAL tmp=0.0;
    long msta, msto, i;
    const MYREAL mu_rate = world->options->mu_rates[world->locus];
    const MYREAL invmu_rate = 1.0/mu_rate;
    const MYREAL inheritance_scalar = world->options->inheritance_scalars[world->locus];
    const MYREAL iii = inheritance_scalar * mu_rate;
    for(j=0; j < numpop; j++)
    {
        const long linex = lines[j];
        MYREAL denom = iii *  world->param0[j] * skyparam[j];
        MYREAL invdenom = 1.0 / denom;
        probth -= ((linex * linex) - linex) * invdenom;
        msta = world->mstart[j];
        msto = world->mend[j];
        tmp = 0.0;
        for (i = msta; i < msto; i++)
        {
            //1229 if (!world->options->custm2[i]=='d')
            tmp += world->data->geo[i] * world->param0[i] * skyparam[i] * invmu_rate;
        }
        probm -= linex * tmp;
    }
    return probth + probm;
}

///
/// do we accept parameter update
boolean
bayes_accept (MYREAL newval, MYREAL oldval, MYREAL heat, MYREAL hastingsratio)
{
    MYREAL r;
    MYREAL diff = (newval - oldval + hastingsratio) * heat;
    // for thermodynamic integration we should not heat this
    // for the standard heating scheme we perhaps should but there is no visible improvement in the runs
    // when we heat this portion for swapping, if we will find some differences I will need to add a flag
    // that differentiates between thermodynamic integration and standard heating.
    if (diff >= 0.0)
        return TRUE;
    r = LOG (RANDUM ());
    if (r < diff)
        return TRUE;
    else
        return FALSE;
}

///
/// log prior ratios for all parameters
// MYREAL log_prior_ratio_all(world_fmt *world, MYREAL *newvals)
// {
//   long np = world->numpop;
//   long np2 = world->numpop2;
//   long npp = np2 + (world->bayes->mu * world->loci);
//   MYREAL logmax = -MYREAL_MAX;
//   long pop;
//   MYREAL ratio = 0;
//   MYREAL *vals;
//   vals = (MYREAL *) mycalloc(npp,sizeof(MYREAL));
//   for(pop = 0; pop < np; pop++)
//     {
//       vals[pop] = log_prior_ratiotheta(newvals[pop], -1., world->bayes,0L);
//       if(logmax < vals[pop])
// 	logmax = vals[pop];
//     }
//   for(pop = np; pop < np2; pop++)
//     {
//       vals[pop] = log_prior_ratiomig(newvals[pop], -1., world->bayes,0L);
//       if(logmax < vals[pop])
// 	logmax = vals[pop];
//     }
//   if(world->bayes->mu)
//     {
//       vals[pop] = log_prior_ratiorate(newvals[pop], -1., world->bayes,0L);
//       if(logmax < vals[pop])
// 	logmax = vals[pop];
//     }
//   ratio = 0;
//   for(pop = 0; pop < npp; pop++)
//     ratio += EXP(vals[pop]-logmax);
//   myfree(vals);
//   return log(ratio)+logmax;
// }

///
/// Log Prior distribution ratios between old and new parameter:
// UNIFORM
// $$p[x] = \frac{1}{b-a} $$
// in reality I use a window to propose new parameters
// but the distribution should be still the same $p[x]$
// the prior distribution for the new parameter and the old are the same
MYREAL log_prior_ratio_uni(MYREAL newparam,
                           MYREAL oldparam,
                           bayes_fmt * bayes,
                           long which)
{
  (void) oldparam;
    if((newparam > bayes->maxparam[which]) || (newparam < bayes->minparam[which]))
      return (double) -HUGE;
    else
        return 0.0;
}

//
// Uniform prior distribution for theta or migration rate used in heating acceptance
MYREAL log_prior_uni(world_fmt *world, long numparam)
{
  long numpop = (long) world->numpop;
  long start = (long) ((numparam <= numpop || numpop==1) ? 0 : numpop);
  long stop = (long) ((start == numpop) ? (long) world->numpop2 : numpop);
  long i;
    //  long frompop;
    //long topop;
  MYREAL * param0 = world->param0;
  bayes_fmt * bayes = world->bayes;
  MYREAL value=0.0;
  MYREAL mu;
  MYREAL p0;
  
  if(numparam>stop) // rate change
    {
      mu = world->options->mu_rates[world->locus];
      i = (long) (world->numpop2 + world->locus);
      if((mu > bayes->maxparam[i]) || (mu < bayes->minparam[i]))
	return (double) -HUGE;
      else
	return  -LOG(bayes->maxparam[i] - bayes->minparam[i]);
    }
  // migration or theta parameters
  for(i = start; i < stop; i++)
    {
      if(!strchr("0c", world->options->custm2[i]))
        {
            p0 = param0[i];
            if((p0 > bayes->maxparam[i]) || (p0 < bayes->minparam[i]))
	      return (double) -HUGE;
            else
                value += -LOG(bayes->maxparam[i] - bayes->minparam[i]);
        }
    }
    return value;
}

//
// uniform prior distribution, returns log(1/(b-a))
MYREAL log_prior_uni1(world_fmt *world, long numparam, MYREAL val)
{
    bayes_fmt * bayes = world->bayes;
    long i = numparam;
    //  long frompop;
    //long topop;
    MYREAL retval = (double) -HUGE;
    MYREAL p0;
    //@@if(i>=world->numpop && !world->options->usem)
    //@@  {
    //@@    m2mm(i,world->numpop,&frompop,&topop);
    //@@    p0 = val * world->param0[topop];
    //@@  }
    //@@else
    //@@  {
    p0 = val;
    //@@  }
    
    if((p0 <= bayes->maxparam[i]) && (p0 >= bayes->minparam[i]))
    {
        retval = bayes->maxparam[i] - bayes->minparam[i];
        if (retval > 0.0)
            return -log(retval);
    }
    return retval;
}

//
// Log Prior distribution ratios between old and new parameter:
// EXPONENTIAL
// $$p[x] = Integrate[Exp[-u/mean]/mean, {u, a, x}]/-Exp[-b/mean] + Exp[-a/mean] $$
// the ratio between old and new parameter prior will be then
// $$
//\frac{e^{-\frac{a}{\text{mean}}}-e^{-\frac{x}{\text{mean}}
//}}{e^{-\frac{a}{\text{mean}}}-e^{-\frac{\text{x0}}{\tex
//						     t{mean}}}}
// $$
//
// The log(hastingsratio) for this move is not zero but cancels with the log_prior_ratio_exp
// so instead of calculating stuff that goes away we set in both places the value to zero
MYREAL log_prior_ratio_exp(MYREAL newparam, MYREAL oldparam, bayes_fmt * bayes, long which)
{
  (void) oldparam;
    if((newparam > bayes->maxparam[which]) || (newparam < bayes->minparam[which]))
      return (double) -HUGE;
    else
        return 0.;
    //  MYREAL imean = 1./bayes->priormean[which];
    //MYREAL precalc = exp(-bayes->minparam[which] *imean);
    //MYREAL xo = EXP(-oldparam * imean);
    //MYREAL xn = EXP(-newparam * imean);
    //MYREAL idenom = 1./(precalc - xo );
    //MYREAL nom = (precalc - xn );
    //return log(nom * idenom);
}

//
// Exponential prior distribution for theta or migration rate used in heating acceptance
// Out[12]//TextForm=
//   b/mean                  (a - c)/mean
// E       (-a + c + (-1 + E            ) mean)
// --------------------------------------------
//                a/mean    b/mean
//              -E       + E
// a: lower bound
// b: upper bound
// c: parameter value
// (Power(E,b/m)*(-1 + Power(E,(a - c)/m)))/(Power(E,a/m) - Power(E,b/m))
MYREAL log_prior_exp(world_fmt *world, long numparam)
{
    return log_prior_gamma(world,numparam);
}

//
// return probability of prior_exp with value
MYREAL log_prior_exp1(world_fmt *world, long numparam, MYREAL value)
{
    return log_prior_gamma1(world,numparam,value);
}



//
// Log Prior distribution ratios between old and new parameter:
// EXPONENTIAL
// $$p[x] = Integrate[Exp[-u/mean]/mean, {u, a, x}]/-Exp[-b/mean] + Exp[-a/mean] $$
// the ratio between old and new parameter prior will be then
// $$
//\frac{e^{-\frac{a}{\text{mean}}}-e^{-\frac{x}{\text{mean}}
//}}{e^{-\frac{a}{\text{mean}}}-e^{-\frac{\text{x0}}{\tex
//						     t{mean}}}}
// $$
// see under log_prior_ratio_exp(), this needs more careful checking as here we
// depend on current theta and almost certainly the hastings ratio is different from the
// one with the exp proposal.
MYREAL log_prior_ratio_wexp(MYREAL newparam, MYREAL oldparam, bayes_fmt * bayes, long which)
{
  (void) oldparam;
    if((newparam > bayes->maxparam[which]) || (newparam < bayes->minparam[which]))
      return (double) -HUGE;
    else
        return 0.;
}

MYREAL log_prior_wexp(world_fmt *world, long numparam)
{
    return log_prior_exp(world, numparam);
}

MYREAL log_prior_wexp1(world_fmt *world, long numparam, MYREAL value)
{
    return log_prior_exp1(world, numparam, value);
}
//
// Log Prior distribution ratios between old and new parameter:
// MULTIPLIER
// $$p[x] = 1/(lambda m)
// the ratio between old and new parameter prior will be then
// $$
//
// $$
MYREAL log_prior_ratio_mult(MYREAL newparam, MYREAL oldparam, bayes_fmt * bayes, long which)
{
    if((newparam > bayes->maxparam[which]) || (newparam < bayes->minparam[which]))
      return (double) -HUGE;
    else
    {
        if(newparam > 0.0)
            return log(newparam/oldparam);
        else
            return -MYREAL_MAX;
    }
}


// incorrect
MYREAL log_prior_mult(world_fmt *world, long numparam)
{
  long numpop = (long) world->numpop;
  long numpop2 = (long) world->numpop2;
  long locus = (long) world->locus;
    long start = ((numparam <= numpop || numpop==1) ? 0 : numpop);
    long stop = ((start == numpop) ? numpop2 : numpop);
    MYREAL * param0 = world->param0;
    MYREAL minp;
    MYREAL maxp;
    MYREAL mean;
    MYREAL imean;
    bayes_fmt *bayes = world->bayes;
    MYREAL value = 0.0;
    MYREAL mu;
    long i;
    
    if(numparam>stop) // rate change
    {
        i = numpop2 + locus;
        mu = world->options->mu_rates[locus];
        if((mu > bayes->maxparam[i]) || (mu < bayes->minparam[i]))
	  return (double) -HUGE;
        else
        {
            maxp = bayes->maxparam[i];
            minp = bayes->minparam[i];
            mean = bayes->meanparam[i];
            imean = 1.0 / mean;
            value += log(EXP(maxp * imean) * (-minp + param0[i] + EXP((minp-param0[i]) * imean))/
                         (EXP(maxp * imean) - EXP(minp * imean)));
        }
        
    }
    for(i = start; i < stop; i++)
    {
        if(!strchr("0c", world->options->custm2[i]))
        {
            maxp = bayes->maxparam[i];
            minp = bayes->minparam[i];
            mean = bayes->meanparam[i];
            imean = 1.0 / mean;
            if((param0[i] > bayes->maxparam[i]) || (param0[i] < bayes->minparam[i]))
	      return (double) -HUGE;
            else
            {
                value += log(EXP(maxp * imean) * (-minp + param0[i] + EXP((minp-param0[i]) * imean))/
                             (EXP(maxp * imean) - EXP(minp * imean)));
            }
        }
    }
    return value;
}

//
//
MYREAL log_prior_mult1(world_fmt *world, long numparam, MYREAL val)
{
    bayes_fmt * bayes = world->bayes;
    long i = numparam;
    MYREAL retval = (double) -HUGE;
    MYREAL minp;
    MYREAL maxp;
    MYREAL mean;
    MYREAL imean;
    if((val <= bayes->maxparam[i]) && (val >= bayes->minparam[i]))
    {
        maxp = bayes->maxparam[i];
        minp = bayes->minparam[i];
        mean = bayes->meanparam[i];
        imean = 1.0 / mean;
        retval = log(EXP(maxp * imean) * (-minp + val + EXP((minp-val) * imean))/(EXP(maxp * imean) - EXP(minp * imean)));
    }
    return retval;
}


//
// scaling prior used for multiple loci posterior
MYREAL scaling_prior(world_fmt *world, long numparam, MYREAL val)
{
    return (*log_prior_1[numparam])(world,numparam,val);
}


MYREAL calculate_prior(world_fmt *world)
{
    const long np = world->numpop2 + world->bayes->mu + world->species_model_size * 2 + world->grownum;
    long i;
    MYREAL retval=0.0;
    
    for (i=0;i<np;i++)
    {
        if (world->bayes->map[i][1] != INVALID)
            retval += (*log_prior_1[i])(world,i,world->param0[i]);
    }
    return retval;
}

void calculate_plotter_priors(world_fmt *world)
{
    bayes_fmt *bayes = world->bayes;
    long np = world->numpop2 + world->bayes->mu + world->species_model_size * 2 + world->grownum;
    MYREAL *priors;
    MYREAL *maxpriors;
    MYREAL *totalspriors;
    boolean *visited;
    long *bins = world->bayes->histogram[0].bins;
    //MYREAL allpriortotal;
    long targetsumbin;
    long pa0, pa, bin;
    MYREAL h,v,w, logw;
    visited = (boolean *) mycalloc(np,sizeof(boolean));
#ifdef DEBUG
    // float * sums = (float *) mycalloc(np,sizeof(float));
#endif
    priors = (MYREAL *) mycalloc(bayes->histogram[0].binsum, sizeof(MYREAL));
    doublevec1d(&maxpriors, np);
    doublevec1d(&totalspriors, np);
    targetsumbin = 0;
    memset(visited,0,sizeof(boolean)*(size_t) np);
    for(pa0 = 0; pa0 < np; pa0++)
    {
        if(shortcut(pa0,world,&pa))
            continue;
        //1229if(world->options->custm2[pa] == 'd')
        //1229    continue;
        if(visited[pa]==TRUE)
            continue;
        visited[pa] = TRUE;
        w = world->bayes->deltahist[pa];
        logw = log(w);
        maxpriors[pa] = (double) -HUGE;
        for(bin=targetsumbin; bin < targetsumbin + bins[pa]; bin++)
        {
            h = world->bayes->deltahist[pa];
            v = world->bayes->histogram[0].minima[pa] + (bin-targetsumbin) * h + h/2;
            priors[bin] = scaling_prior(world,pa,v);
            if(priors[bin] > maxpriors[pa])
                maxpriors[pa] = (MYREAL) priors[bin];
        }
        totalspriors[pa]=0.0;
        for(bin=targetsumbin; bin < targetsumbin + bins[pa]; bin++)
        {
            totalspriors[pa] += exp(priors[bin] - maxpriors[pa]);
        }
        if(totalspriors[pa]>=0.0)
            totalspriors[pa] = log(totalspriors[pa]) + maxpriors[pa] + logw;
        else
            totalspriors[pa]= (double) -HUGE;
        for(bin=targetsumbin; bin < targetsumbin + bins[pa]; bin++)
        {
            priors[bin] = exp(priors[bin]);
#ifdef DEBUG
            //sums[pa] += priors[bin];
#endif
        }
        targetsumbin += bins[pa];
    }
#ifdef DEBUG
    //for(pa0 = 0; pa0 < np; pa0++)
    //    printf("sums[%li]=%f ==?== 1.0??\n",pa0,sums[pa0]);
    //myfree(sums);
#endif
    world->bayes->priors = priors;
    myfree(maxpriors);
    myfree(totalspriors);
    myfree(visited);
}


//
// verbose log to stdout: report all changes [too much output, this will slow down run a lot]
void print_bayes_verbose(long which, world_fmt *world, MYREAL newparam, boolean success)
{
    long i;
    //  long frompop;
    //long topop;
    //MYREAL theta;
    
    if(world->options->verbose)
    {
        fprintf(stdout,"%i> <%li> ", myID, which);
        for(i=0;i<world->numpop;i++)
            fprintf(stdout,"%f ", which==i ? newparam : world->param0[i]);
        //@@if(world->options->usem)
        //@@{
        for(i=world->numpop;i<world->numpop2;i++)
            fprintf(stdout,"%f ", which==i ? newparam : world->param0[i]);
        //@@}
        //@@else
        //@@	{
        //@@ for(i=world->numpop;i<world->numpop2;i++)
        //@@   {
        //@@m2mm(i,world->numpop,&frompop,&topop);
        //@@theta = world->param0[topop];
        //@@fprintf(stdout,"%f ", which==i ? newparam * theta :
        //@@      world->param0[i] * theta);
        //@@	  }
        //@@}
        if(world->bayes->mu)
            fprintf(stdout,"%f ", which==world->numpop2 ? newparam : world->options->mu_rates[world->locus]);
        fprintf(stdout,"%f %c\n",world->param_like,success ? '*': ' ');
    }
}

//
// standard metropolis proposal
MYREAL uniform_proposal(long which, world_fmt * world, MYREAL *oldparam, boolean *success)
{
#ifdef RANGEDEBUG
    int rangedebug=0;
#endif
    //long        topop;
    //long        frompop;
    //static long counter = 0;
    const long  numpop2 = world->numpop2;
    const long  npx = numpop2 + world->species_model_size * 2;
    //const long  npx2 = npx + world->grownum;
    const MYREAL themin = world->bayes->minparam[which];
    const MYREAL themax = world->bayes->maxparam[which];
    //  MYREAL      mean;
    MYREAL      delta;
    MYREAL      newparam= -999;
    MYREAL      oldval;
    MYREAL      r;
    MYREAL      newval;
    MYREAL      hastingsratio;
    //MYREAL      theta;
    //MYREAL      nm;
    long specid = -1;
    //MYREAL oldsigoma = 0.0;
    //MYREAL oldmu=0.0;
    //MYREAL newmu=0.0;
    //MYREAL newsigma=0.0;
    const MYREAL      murate = world->options->mu_rates[world->locus];
    MYREAL    * param0 = world->param0;
    bayes_fmt * bayes = world->bayes;
    boolean is_mu=FALSE;//this is the mean of the speciation time mu
    const boolean verbose = (world->heat > (1.0  - SMALLEPSILON)) && myID==MASTER && world->options->verbose;
    
    //boolean growthprop = FALSE;
    // calculate the probability for the old parameter set
    // this is done from scratch because the tree might have changed in the last step
    delta  = bayes->delta[which];
    oldval = probg_treetimes(world);
    r = UNIF_RANDUM();
    // draw a new parameter from the prior distribution
    // for migration parameters we need to distinguish whether the
    // prior is in terms of M or xNm
    if(which < world->numpop)
    {
        newparam = (*propose_new[which]) (param0[which],which, world, &r);
        check_min_max_param(&newparam,themin,themax);
        hastingsratio = (*hastings_ratio[which])(newparam, oldparam[which], delta, r, bayes, which);
        // add the prior ratio to the log of the hastings ratios
        hastingsratio += (*log_prior_ratio[which])(newparam,oldparam[which], bayes, which);
        bayes_set_param(world->param0,newparam,which,world->options->custm2, world->numpop);
        reprecalc_world(world, which);
        // calculate prob(G|params)prob(D|G)
        newval = probg_treetimes(world);
    }
    else if(which < numpop2)
        {
            newparam = (*propose_new[which]) (param0[which], which, world, &r);
            check_min_max_param(&newparam,themin,themax);
            hastingsratio = (*hastings_ratio[which])(newparam, oldparam[which], delta, r, bayes, which);
            // add the prior ratio to the log of the hastings ratios
            hastingsratio += (*log_prior_ratio[which])(newparam,oldparam[which], bayes, which);
            bayes_set_param(world->param0,newparam,which,world->options->custm2, world->numpop);
            reprecalc_world(world, which);
            // calculate prob(G|params)
            newval = probg_treetimes(world);
        }
    else if(bayes->mu && which == numpop2)
            {
                newparam = -1000.0;
                newparam = (*propose_new[which]) (murate, which, world, &r);
                check_min_max_param(&newparam,themin,themax);
                hastingsratio = (*hastings_ratio[which])(newparam, murate, delta, r, bayes, which);
                // add the prior ratio to the log of the hastings ratios
                hastingsratio += (*log_prior_ratio[which])(newparam, murate, bayes, which);
                world->options->mu_rates[world->locus] = newparam;
                world->options->lmu_rates[world->locus] = log(newparam);
                // change the tree length because of the rate
                reprecalc_world(world, which);
                recalc_timelist(world, world->options->mu_rates[world->locus] , oldparam[which]);
                // calculate prob(G|params)
                newval = probg_treetimes(world);
            }
    else if (which < npx)
            {
	//warning("GROWTH has issues with speciation");
                specid = propose_new_spec_mu(world, which, &is_mu, &newparam);
                hastingsratio = (*hastings_ratio[which])(newparam, oldparam[which], delta, r, bayes, which);
                hastingsratio += (*log_prior_ratio[which])(newparam, oldparam[which], bayes, which);
                bayes_set_param(world->param0,newparam,which,world->options->custm2, world->numpop);
                reprecalc_world(world, which);
                // calculate prob(G|params)
                newval = probg_treetimes(world);
            }
    else
      {
	long gw = which-npx;//world->options->growpops[which-npx]-1;
	if (gw<0)
	  {
	    error("proposal for new growth value tried to use a no-growth population\n");
	  }
	else
	  {
	    //growthprop = TRUE;
	    world->savegrowth[gw] = world->growth[gw];
	    world->growth[gw] = (*propose_new[which]) (world->growth[gw],which, world, &r);
	    newparam = world->growth[gw];
	    //fprintf(stdout,"%i> new growth[%li]=%f trial\n",myID, gw,newparam);
	    check_min_max_param(&world->growth[gw],themin,themax);
	    hastingsratio = (*hastings_ratio[which])(world->growth[gw], world->savegrowth[gw], delta, r, bayes, which);
	    // add the prior ratio to the log of the hastings ratios
	    hastingsratio += (*log_prior_ratio[which])(world->growth[gw], world->savegrowth[gw], bayes, which);
	    bayes_set_param(world->param0,newparam,which,world->options->custm2, world->numpop);
	    //reprecalc_world(world, which);
	  }
	// calculate prob(G|params)
	newval = probg_treetimes(world);
	//printf("%c Growth which=%li newval=%f oldval=%f diff=%f\n", newval-oldval>0 ? '+' : '-', which,newval,oldval,newval-oldval);
    }
    //Acceptance or rejection of the new value
    *success = bayes_accept(newval, oldval,world->heat, hastingsratio);
    if(*success)
      {
	//if (world->has_growth && which >= npx)
	  //  printf("+");fflush(stdout);
	  //{
	    //if(world->heat==1.0)
	    //  printf("+ %i> growth=%f theta0=%f\n", myID, world->growth[0], world->param0[0]);
	  //}
	//if(world->heat==1.0 && which >= world->numpop2)
	//  {
	//    printf("%i> + %li %f %f (%f,%f,%f)\n",myID,which,newval, oldval,world->param0[which],newparam,oldparam[which]);
	//  }
	if(verbose)
	  print_bayes_verbose(which,world, newparam, *success);
	return newval;
      }
    else
      {
        if(world->options->prioralone)
	  {
            *success = TRUE;
#ifdef DEBUG
            if(verbose)
	      print_bayes_verbose(which,world, newparam, *success);
#endif
            return newval;
	  }
#ifdef DEBUG
        if(verbose)
	  print_bayes_verbose(which,world, newparam, *success);
#endif
        if(which<numpop2)
	  {
	  bayes_set_param(world->param0,oldparam[which],which,world->options->custm2, world->numpop);
	    reprecalc_world(world, which);
	  }
        // mutation rate undo
        else if(bayes->mu && which==numpop2)
	  {
            recalc_timelist(world, oldparam[which], world->options->mu_rates[world->locus]);
            world->options->mu_rates[world->locus] = oldparam[which];
            world->options->lmu_rates[world->locus] = log(oldparam[which]);
            world->param0[which] = oldparam[which];
            reprecalc_world(world, which);
	  }
        // speciation stuff undo
        else if(which < npx)
	  {
            world->param0[which] = oldparam[which];
            if(is_mu)
	      world->species_model[specid].mu = (double) oldparam[which];
            else
	      world->species_model[specid].sigma = (double) oldparam[which];
	    if(world->species_model[specid].type == 't')
	      {
		long t;
		for (t=0;t<world->species_model_size;t++)
		  {
		    species_fmt *s     = &world->species_model[specid];
		    species_fmt *other = &world->species_model[t];
		    if(s->from == other->from && other->type == 't')
		      {
			other->mu = s->mu;
			other->sigma = s->sigma;
			long y = t*2 + world->bayes->mu + world->numpop2;
			world->param0[y] = (MYREAL) s->mu;
			world->param0[y+1] = (MYREAL) s->sigma;
		      }
		  }
	      }
		/*if(world->heat==1.0 && which >= world->numpop2)
	      {
	    	printf("%i>  -%li %f %f (%f,%f,%f)\n",myID,which,newval, oldval,world->param0[which],newparam,oldparam[which]);
		}	  */ 
	  }
	else
	  {
	    //  if (world->has_growth && which >= npx)
	    //{
	    //	if(world->heat==1.0)
	    //	  printf("- %i> growth=%f theta0=%f\n", myID, world->growth[0], world->param0[0]);
	    //  }
	    if (world->options->growpops[which-npx]!=0)
	      {
		long gw = world->options->growpops[which-npx] - 1;
		world->growth[gw] = world->savegrowth[gw];
		//printf(".");fflush(stdout);
	      }
	  }
        //printf("  %li %f %f\n",which,newval, oldval);
        //printf(".");fflush(stdout);
        return oldval;
      }
}

void traverse_adjust(node *theNode, MYREAL new_old_ratio)
{
    MYREAL age = 0.;
    node * right=NULL;
    if(theNode == NULL)
        return;
    
    if(theNode->type != 't')
    {
        //fprintf(stdout,"@traverse_adjust@(%i%c-%f/%f)\n",myID,theNode->type,theNode->tyme,theNode->tyme*new_old_ratio);
      right = theNode->next;
      while(!right->top)
	{
	  traverse_adjust(right->back, new_old_ratio);
	  right = right->next;
	}
      //if(right->top)
      //	    break;
      //	  if (theNode->type != 'm' && theNode->type != 'd' && right->back != NULL)
	  //traverse_adjust(theNode->next->back, new_old_ratio);
        //        if(theNode->type != 'm' || theNode->type != 'd')
        //if (theNode->type != 'm' && theNode->type != 'd' && theNode->next->next->back != NULL)
        //{
        //    traverse_adjust(theNode->next->next->back, new_old_ratio);
        //}
    }
#ifdef DEBUG
    //else
    //{
    //  fprintf(stdout,"(#traverse_adjust@%i%c-%f/%f)\n",myID,theNode->type,theNode->tyme,theNode->tyme*new_old_ratio);
    //}
#endif
    age = theNode->tyme * new_old_ratio;
    adjust_time_all (theNode, age);
}


//
// changes branchlength of all branches and the likelihood needs to be calculated completely
void
recalc_timelist (world_fmt * world, MYREAL new_ratio, MYREAL old_ratio)
{
    MYREAL new_old_ratio = new_ratio / old_ratio;
    //if (fabs(new_old_ratio - 1.0) < SMALLEPSILON)
    //    return;
#ifdef DEBUG
    //printf("+++++++++++++\n%f=%f/%f\n",new_old_ratio,new_ratio,old_ratio);
#endif
    first_smooth(world,world->locus);
    traverse_adjust(world->root->next->back, new_old_ratio);
    set_v (world->root->next->back);
    first_smooth(world,world->locus);
    construct_tymelist (world, &world->treetimes[0]);
#ifdef RATE_DEBUG
    //if(world->cold)
    printf("%i> HEAT: %f, RATE: newrate=%f oldrate=%f oldlike=%f ",myID, world->heat, new_ratio, old_ratio,world->likelihood[world->G]);
#endif
    world->likelihood[world->G] = treelikelihood(world);
#ifdef RATE_DEBUG
    //  if(world->cold)
    printf("newlike=%f\n",world->likelihood[world->G]);
#endif
}

void recalc_timelist_1 (world_fmt * world, MYREAL new_old_ratio)
{
    long k;
    MYREAL age;
    
    node *theNode;
    if(new_old_ratio < SMALLEPSILON)
    {
        new_old_ratio = SMALLEPSILON;
    }
    printf("\n\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    for (k = 0; k < world->treetimes[0].T; k++)
    {
        world->treetimes[0].tl[k].age *= new_old_ratio;
        age = world->treetimes[0].tl[k].age;
        theNode = world->treetimes[0].tl[k].eventnode;
        adjust_time_all (theNode, age);
    }
    for (k = 0; k < world->treetimes[0].T; k++)
    {
        theNode = world->treetimes[0].tl[k].eventnode;
        if (theNode->top == 0 )
            error("Node is not top\n");
        if(theNode->type=='i')
            printf("**p(%li):%f l(%li):%f r(%li):%f |||| %li %li\n",theNode->id, theNode->tyme, theNode->next->back->id, theNode->next->back->tyme,  theNode->next->next->back->id, theNode->next->next->back->tyme, world->treetimes[0].tl[k].lineages[0],world->treetimes[0].tl[k].lineages[1]);
        if(theNode->type=='m')
            printf("**p(%li):%f l(%li):%f r: ---- |||| %li %li\n",theNode->id, theNode->tyme,  theNode->next->back->id, theNode->next->back->tyme, world->treetimes[0].tl[k].lineages[0],world->treetimes[0].tl[k].lineages[1]);
        if(theNode->type=='t')
            printf("**p(%li):%f l:---- r: ---- |||| %li %li\n",theNode->id, theNode->tyme, world->treetimes[0].tl[k].lineages[0],world->treetimes[0].tl[k].lineages[1]);
        if(theNode->type !='t' && (theNode->next->back->tyme > theNode->tyme))
        {
            warning("***********************Time mess: above=%f this=%f\n,",theNode->next->back->tyme , theNode->tyme);
        }
    }
    treelikelihood(world);
    //    first_smooth(world,world->locus);
    //construct_tymelist (world, &world->treetimes[0]);
}

//
// Update for parameters using Bayesian inference
long
bayes_update (world_fmt * world)
{
    const boolean     mu = world->bayes->mu;
    const long        numpop = world->numpop;
    const long        numpop2 = world->numpop2;
    //const long        numpop2rate = world->numpop2 + mu;
    const long        npx2 = get_numparam(world);
    const long        npx = npx2 - world->grownum;
    const MYREAL      murate = world->options->mu_rates[world->locus];
    const boolean     writelog = (myID==MASTER &&
                                  world->options->progress &&
                                  world->cold);
    const long tt = world->options->lsteps * world->increment;
    const long ten = ((tt > 100000) ? (tt/10) : (tt/3));
    boolean           success=FALSE;
    long              ba=0;
    long              which;
    long              w; // w is a function of bayes->map[which]
    long              type;
    MYREAL            *oldparam=NULL;
    MYREAL            newval = (double) -HUGE;
    bayes_fmt         *bayes = world->bayes;
    MYREAL            paramval;
#ifdef DEBUG
    //double oldval = bayes->oldval;
#endif
    
    // progress reporter nested here so that we do need to ask whether we run in Bayes mode or not
    if(!world->in_burnin && writelog)
    {
        world->treesdone += 1;
        if(writelog && (world->treesdone % ten == 0))
        {
            bayes_progress(world, ten);
        }
    }
    // select parameter
    which /*= world->bayes->paramnum*/ = RANDINT(0L,npx2-1);
    //printf("select which: %li ",which);
    //    while(bayes->map[which][1] == INVALID)
    while(shortcut(which,world,&w))
      {
        which = RANDINT(0L,npx2-1);
        //    printf("%li ",which);
      }
    //  printf("\nselect which (chosen): %li\n",which);
    // savecopy the old parameters
    // we use a memcopy because the custom migration matrix can be filled with m and M and s and S
    // this will change multiple parameters at the same time
    doublevec1d(&oldparam,npx2);
    memcpy(oldparam, world->param0,sizeof(MYREAL) * (size_t) npx2);
    if(mu)
        oldparam[numpop2] = murate;
    //if(which<npx)
    //  w = world->bayes->map[which][1];
    //else
    //  w = which;   ## the while loop above should do the trick
#ifdef DEBUG
    //  printf("%i> which which: %li --> %li\n",myID,which,w);
#endif
    paramval = world->param0[w];
    if(which <numpop)
        type = THETAPRIOR;
    else  if (which < numpop2)
            type = MIGPRIOR;
    else if (mu && which==numpop2)
            {
                type = RATEPRIOR;
                paramval = world->options->mu_rates[world->locus];
            }
    else if (which<npx)
                type = SPECIESTIMEPRIOR;
    else
      type = GROWTHPRIOR;
    
    if(world->options->slice_sampling[type])
    {
      //warning("SLICE sampler does not work with MLF and GROWTH\n");
        newval = expslice(&paramval, which, world, log_prior_ratio[w]);
        success = TRUE;
    }
    else
    {
        newval = uniform_proposal(w, world, oldparam, &success);
    }
    world->logprior = calculate_prior(world);
    //if(which==0)
    //printf("%i> which=%li (%li) param=%f oldparam=%f\n", myID, which, w, world->param0[w],oldparam[w])
    double v = (newval + world->likelihood[world->G]) * world->heat + world->logprior;
    long ii = world->locus * world->options->heated_chains + world->heatid;
    if (world->steppingstone_counters[ii] < 1.0)
      {
      	world->steppingstone_scalars[ii] = v;
	world->steppingstone_counters[ii] += 1;
      }
    if (world->steppingstone_scalars[ii] < v)
      {
	world->steppingstones[ii] *= exp(world->steppingstone_scalars[ii] - v);
	world->steppingstone_scalars[ii] = v;
	world->steppingstones[ii] += 1.0;
      }
    else
      {
	world->steppingstones[ii] += exp(v-world->steppingstone_scalars[ii]);
      }
    
    //printf("[%li] @Steppingstone[%li]: %f %f\n",world->heatid,ii,world->steppingstones[ii], world->steppingstone_scalars[ii]);
    if(success)
      {
        bayes->oldval = newval;
        world->param_like = newval;
        ba = 1;
        bayes->accept[which] += ba;
	//if(which>=npx && world->cold)
	//  printf("%li/%li/%li ",which,bayes->accept[which],bayes->trials[which]);fflush(stdout);
	if(world->in_burnin && world->cold && !world->options->prioralone)
	  autotune_proposal(world,which);
      }
    else
      {
	if(world->in_burnin && world->cold && !world->options->prioralone)
	  autotune_proposal(world,which);
	memcpy(world->param0, oldparam, sizeof(MYREAL) * (size_t) npx2);
	if(type==RATEPRIOR)
	  {
            recalc_timelist(world, oldparam[which], world->options->mu_rates[world->locus]);
            //reprecalc_world(world, which);
	  }
        else
	  {
            reprecalc_world(world, which);
	  }
        if(mu && which==numpop2)
	  {
            world->options->mu_rates[world->locus] = oldparam[numpop2];
            world->options->lmu_rates[world->locus] = log(oldparam[numpop2]);
	  }
        bayes->oldval = newval; // uniform_proposal delivers either old or new value see success
        ba = 0;
	//if(which>=npx && world->cold)
	//  printf("%li/%li/%li ",which, bayes->accept[which],bayes->trials[which]);fflush(stdout);
      }
    bayes->trials[which] += 1;
#ifdef DEBUG
    // print a lot, left in here for other times
    //if (world->cold && which==3)//check migration
    //    printf("%i> which=%li param=%f oldval=%f logP=%f success=%c\n", myID, w, world->param0[w], oldval, bayes->oldval, success ? 'T' : 'F');
#endif
    myfree(oldparam);
    return ba;
}



//
// fill bayes record in array for histograms and bayesfile
// set parameter meaning according to option settings
void bayes_save_parameter(world_fmt *world, long pnum, long step)
{
  //needs addition of world->mlalpha
    long i;//,i0;
    //long frompop      = 0;
    //long topop        = 0;
    
    long n0           = world->numpop2 + (world->bayes->mu) + world->species_model_size * 2;
    long n            = n0 + world->grownum;
    long nn           = 2 + n;
    long numpop       = world->numpop;
    long numpop2      = world->numpop2;
    long nnpnum       = nn * pnum;
    boolean mu        = world->bayes->mu;
    //  MYREAL murate     = world->options->meanmu[world->locus] * world->options->mu_rates[world->locus];
    MYREAL murate     = world->options->mu_rates[world->locus];
    MYREAL inheritance_scalar = world->options->inheritance_scalars[world->locus];
    MYREAL * param0   = world->param0;
    worldoption_fmt *wopt = world->options;
    bayes_fmt * bayes = world->bayes;
    //MYREAL nm;
    (bayes->params+(nnpnum))[0] = bayes->oldval;
    (bayes->params+(nnpnum))[1] = world->likelihood[world->G];
    
    memcpy(bayes->params+(nnpnum+2), param0,sizeof(MYREAL) * (size_t) n);
    if(inheritance_scalar != 1.0)
    {
        for(i = 0; i < numpop; i++)
        {
            (bayes->params+(nnpnum + 2))[i] = param0[i] * inheritance_scalar;
        }
    }
    if(mu)
    {
        // we only write one rate per record because the locus is known
        (bayes->params+(nnpnum+2))[numpop2] = murate;
    }
    if(world->has_growth)
      {
	// check whether this is needed
	memcpy(bayes->params+(nnpnum+2+n0), world->growth, sizeof(double)* (size_t) world->grownum);
      }
    if(wopt->has_bayesmdimfile)
    {
      print_bayes_tofile(world->bayesmdimfile, bayes->params+nnpnum, bayes, world, numpop2, bayes->mdiminterval, world->locus, world->replicate, step);
    }
}

//
// print all parameters continously to a file. This may produce HUGE files
// and overrun your disk quota
#ifdef ZNZ
void	  print_bayes_tofile(znzFile mdimfile, MYREAL *params, bayes_fmt *bayes, world_fmt *world, long numpop2, long mdiminterval, long locus, long replicate, long step)
#else
void	  print_bayes_tofile(FILE *mdimfile, MYREAL *params, bayes_fmt *bayes, world_fmt *world, long numpop2, long mdiminterval, long locus, long replicate, long step)
#endif
{
#ifdef MPI
  (void) mdimfile;
#endif
  (void) numpop2;
    //  long i;
    long j, j0;

    //long grownum      = world->options->growpops_numalloc;
    long interval = mdiminterval > 0 ? mdiminterval : 1;
    boolean *visited;
    long nn = world->numpop2 + (world->bayes->mu) + world->species_model_size * 2 + world->grownum;
    // the recorder is dealing with the addition of species parameters and growth
    //#ifdef MPI
#if defined(MPI) && !defined(PARALIO)
    long addition = world->species_model_size*2 + world->grownum;
    float *temp;
    long z=0;
#else
#ifdef PARALIO
    MPI_Status status;
#endif
    int fmt = 5;
    long digits;
    long c=0;
    char *temp;
    double value = 0.0;
#endif
    bayes->mdimfilecount[locus] += 1;
    if (bayes->mdimfilecount[locus] % interval == 0)
    {
        visited = (boolean*) mycalloc(nn+2,sizeof(boolean));
        //#ifdef MPI
#if defined(MPI) && !defined(PARALIO)
        temp = (float*) mycalloc((size_t)(nn+9+3*world->options->heated_chains+2+(world->timeelements +  world->timeelements * (world->numpop2 + addition))), sizeof(float));
        temp[0] = (float) step;
        temp[1] = (float) locus;
        temp[2] = (float) replicate;
        temp[3] = (float) params[0];
        temp[4] = (float) params[1];
        temp[5] = (float) (params[0] - params[1]);
        temp[6] = (float) world->logprior;
        temp[7] = (float) (world->treetimes->T-1.0);
        temp[8] = (float) world->treelen;
        z=9;
#else
        temp = (char*) mycalloc(bayes->linesize, sizeof(char));
        c = sprintf(temp,"%li\t%li\t%li\t%f\t%f\t%f\t%f\t%li\t%f",step, locus+1, replicate+1, params[0], params[1],params[0] - params[1], world->logprior, world->treetimes->T-1, world->treelen);
#endif
        for (j0=0; j0 < nn; j0++)// the first value is the posterior Log(P(D|G)P(G|p)
            // second is log(P(D|G)
        {
	  if(shortcut(j0, world, &j))
	    continue;
	  else
            {
                j = j + 2;
                if(visited[j]==TRUE)
                    continue;
                else
                    visited[j]=TRUE;
            }
            //#ifdef MPI
#if defined(MPI) && !defined(PARALIO)
            temp[z] = (float) params[j];
            z++;
            //  fprintf (stdout, "%f (%li) ",temp[j+6],z);
#else
            if(fabs(params[j]) < SMALLEPSILON)
            {
                c += sprintf(temp+c,"\t0");
            }
            else
            {
	      if (params[j] >= 0.0)
		{
		  value = params[j];
		}
	      else
            {
		  value = -params[j];
		}
	      digits = (long) (floor(log10(value)));
                switch(digits)
                {
                    case -8:
                    case -6:
                    case -5:
                        fmt = 10;
                        break;
                    case -4:
                    case -3:
                        fmt = 8;
                        break;
                    case -2:
                    case -1:
                        fmt= 5;
                        break;
                    case 0:
                    case 1:
                        fmt = 4;
                        break;
                    case 2:
                        fmt = 2;
                        break;
                    case 3:
                        fmt = 1;
                        break;
                    case 4:
                    case 5:
                    case 6:
                    case 7:
                    case 8:
                        fmt = 0;
                        break;
                    default:
                        if(digits<-8)
                            fmt=20;
                        else
                            fmt = 5;
                        break;
                }
                c += sprintf(temp+c,"\t%.*f",fmt, params[j]);
            }
#endif
        }
#if defined(MPI) && !defined(PARALIO)
        //print_growth_record(temp, &z, world);
        print_marginal_like(temp, &z, world);
#else
        //print_growth_record(temp, &c, world);
        print_marginal_like(temp, &c, world);
#endif
#if defined(MPI) && !defined(PARALIO)
        //#ifdef MPI
        mpi_mdim_send(temp,z);
#else
        c += sprintf(temp+c,"\n");
        if(!bayes->has_linesize)
        {
            bayes->linesize = TWO * c;
            bayes->has_linesize = TRUE;
#ifdef ZNZ
	    long bytes = bayes->linesize > ONEMEGABYTE ? bayes->linesize : ONEMEGABYTE;
            znzbuffer(mdimfile,(size_t) bytes);
#endif
        }
#ifdef PARALIO
        //wrapper function,
        //      fprintf(stdout,"%i>%li %li\n@@%s@@\n\n",myID,c, strlen(temp),temp);
        if(c!=strlen(temp))
        {
            printf("%i> %li %li\n%s\n",myID,c,strlen(temp),temp);
            error("failed in bayes.c:2017");
        }
        mpi_mdim_send(&world->mpi_bayesmdimfile,temp,c);
#else
#ifdef ZNZ
        //znzprintf(mdimfile,"%s", temp);
        znzwrite(temp, (size_t) c, (size_t) sizeof(char), mdimfile);
#else
        /* this should do the trick, windows was printing an empty line*/
        FPRINTF(mdimfile,"%s", temp);
#endif /*znz*/
#endif /*paralio*/
#endif /*defined(mpi) and !defined(paralio)*/
        myfree(temp);
        myfree(visited);
    }
}



//
// Save the Bayesian results for printout into bayesfile
void bayes_save(world_fmt *world, long step)
{
  long np           = 2 + world->numpop2 + (world->bayes->mu) + world->species_model_size * 2 + world->grownum;// + growth
    long pnum         = world->bayes->numparams;
    long allocparams  = world->bayes->allocparams;
    
    if(world->options->has_bayesmdimfile)
    {
        bayes_save_parameter(world, 0L, step);
    }
    else
    {
        bayes_save_parameter(world, pnum, step);
        pnum++;
        if(pnum >= allocparams)
        {
            allocparams += 1000;
            world->bayes->params = (MYREAL *) myrealloc(world->bayes->params,sizeof(MYREAL)* (size_t) allocparams* (size_t) np);
        }
        world->bayes->numparams = pnum;
        world->bayes->allocparams = allocparams;
    }
}

 
void  set_map_groups(long numpop, long t, long size, longpair *map, long lastold)
{
  (void) size;
  long i;
  long j;
  long first = 0;

  for (j= -2; j >= lastold; j--)
    {
      for (i=0; i < numpop; i++)
	{      
	  if(map[i][1] == j)
	    {
	      first = map[i][0];
	      break;
	    }
	}
      for (i=0; i < numpop; i++)
	{      
	  if(map[i][1] == j)
	    {
	      map[i][1] = first;
	    }
	}
      for (i=numpop; i < t; i++)
	{      
	  if(map[i][1] == j)
	    {
	      first = map[i][0];
	      break;
	    }
	}
      for (i=numpop; i < t; i++)
	{      
	  if(map[i][1] == j)
	    {
	      map[i][1] = first;
	    }
	}
    }
  // do not adjust splits -- they do not work with the flexible
  // grouping used for migration parameters
  //for (i=t; i < size; i++)
  //  {
  //    map[i][0] = i;
  //    map[i][1] = i;
  //  }
}



long setup_bayes_map(longpair *map, world_fmt *world, long size)
{
  long specdistrib = world->species_model_dist;
    char *custm2 = world->options->custm2;
    const long numpop = world->numpop;
    const long numpop2 = world->numpop2;
    long invalid_count = 0;
    long j1, j2;
    long n = (long) strlen(custm2);
    long s = MIN(n,numpop2);
    long t = MIN(s,size);
    long i;
    long frompop;
    long topop;
    //long first = 0;
    boolean mu = world->bayes->mu;
    long species_start = mu + numpop2;
    long old = -2;
#ifdef DEBUG
    fprintf(stdout,"%i> @@@@@@@ heatid=%li custm=%s custm2=%s  (species: %li)\n",
	    myID, world->heatid, world->options->custm, world->options->custm2, specdistrib);
    
#endif
    for (i=0; i < t; i++)
    {
        map[i][0] = i;
        switch(custm2[i])
        {
	case 'D':
	case 'T':
	  //WARNING WARNING
	  //warning("This fractional coalescent version does not support DIVERGENCE ESTIMATION yet\n");
	  map[i][1] = i;
	  map[species_start][0]=species_start;
	  map[species_start][1]=species_start;
	  species_start++;
	  if(specdistrib == EXP_DIST)
	    {	      
	      map[species_start][0]= species_start;
	      map[species_start][1]= INVALID;
	    }
	  else
	    {
	      map[species_start][0]=species_start;
	      map[species_start][1]=species_start;
	    }
	  species_start++;
	  //custm2[i]='*';
	  break;
	case 't':
	case 'd':
	  //usererror("This fractional coalescent version does not support DIVERGENCE ESTIMATION yet");
	  map[i][1] = INVALID;
	  map[species_start][0]=species_start;
	  map[species_start][1]=species_start;
	  species_start++;
	  if(specdistrib == EXP_DIST)
	    {	      
	      map[species_start][0]=species_start;
	      map[species_start][1]= INVALID;
	    }
	  else
	    {
	      map[species_start][0]=species_start;
	      map[species_start][1]=species_start;
	    }
	  species_start++;
	  //custm2[i]='0';
	  break;
	case '*':
	  map[i][1] = i;
	  break;
	case '0':
	case 'c':
	  map[i][1] = INVALID;
	  invalid_count++;
	  break;
	case 's':
	case 'S':
	  m2mm(i,numpop,&frompop,&topop);
	  if(frompop != topop)
	    {
	      j1 = numpop + topop * (numpop-1) + ((frompop < topop) ? frompop : frompop-1);
	      j2 = numpop + frompop * (numpop-1) + ((topop < frompop) ? topop : topop-1);
	      //printf("j1=%li, j2=%li\n",j1,j2);
	      map[i][1] = j1 < j2 ? j1 : j2;
	    }
	  else
	    map[i][1] = -2; //solves miscoding of sizes with "symmetric" to "mean"
	  break;
	case 'm':
	case 'M':
	  map[i][1] = -2;
	  break;
	default:
	  if (i<numpop)
	    {
	      map[i][1] = ((int) 'a') - ((int) custm2[i]) - 3;
	      if (old > map[i][1])
		old = map[i][1];
	    }
	  else if (i<numpop*numpop)
	    {
	      map[i][1] = ((int) 'a') - ((int) custm2[i]) - 3;
	      if (old > map[i][1])
		old = map[i][1];
	    }
	  break;
        }
    }
    if (size>t)
      {
	// growth
	long npx2 = numpop2 + mu + 2 * world->species_model_size;
	for (i=npx2;i<world->grownum+npx2;i++)
	  {
	    map[i][0] = i;
	    map[i][1] = i;
	  }
      }
  set_map_groups(numpop, t, size, map, old);
#ifdef DEBUG
  fprintf(stdout,"%i> bayes map\n",myID);
  for (i=0;i<size;i++)
    {
      fprintf(stdout,"%i> %li map[0]=%li map[1]=%li ",myID,i, map[i][0],map[i][1]);
      if (i<t)
	{
	  fprintf(stdout,"%c\n",custm2[i]);
	}
      else
	  fprintf(stdout,"\n");
    }
#endif
    return invalid_count;
}

//
// Initialize the Bayesian framwork
void bayes_init(bayes_fmt *bayes, world_fmt *world, option_fmt *options)
{
  long sizeg = get_numparam(world);
  long size =  sizeg - world->grownum;
    long invalids;
    bayes->numpop2 = world->numpop2;
    bayes->mu = options->bayesmurates;
    bayes->linesize = MAXBUFSIZE;
    bayes->has_linesize=FALSE;
    bayes->progresslinesize = world->numpop2 * HUNDRED + 5 * HUNDRED;
    bayes->map = (longpair *) mycalloc(sizeg, sizeof(longpair));
    bayes->mapsize = sizeg;
    //printf("BAYES map size=%li\n",size);
    
    bayes->hyperprior = options->hyperprior;
    if (bayes->hyperprior)
    {
        bayes->hyperinterval=options->hyperinterval;
	bayes->hyperfactormean = options->hyperfactormean;
	bayes->hyperfactoralpha = options->hyperfactoralpha;
        init_hyperpriorrecord(&bayes->hyperp,size);
    }
    invalids = setup_bayes_map(bayes->map, world, sizeg);
    if(invalids == size)
    {
        // reset the choices for the updating to exlcude parameter changes
        //
        options->updateratio = 1.0;
        //world->options->updateratio = 1.0;
        set_updating_choices(world->options->choices, options, NOPARAMETER);
    }
    bayes->oldval = -MYREAL_MAX;
    bayes->allocparams = 1;
    bayes->numparams = 0;
    bayes->paramnum = 0;
    // datastore for several variables
    bayes->datastore = (MYREAL *) mycalloc((8 * sizeg + 2),sizeof(MYREAL));
    // pointers into datastore
    bayes->priormean = bayes->datastore;
    bayes->delta = bayes->priormean + sizeg;
    bayes->minparam = bayes->delta + sizeg;
    bayes->maxparam = bayes->minparam + sizeg;
    bayes->meanparam = bayes->maxparam + sizeg;
    bayes->deltahist = bayes->meanparam + sizeg;
    bayes->alphaparam = (MYREAL *) mycalloc((3*sizeg),sizeof(MYREAL));
    bayes->alphaorigparam = bayes->alphaparam + sizeg;
    bayes->betaparam =  bayes->alphaorigparam + sizeg;
    //old  bayes->deltahist = (MYREAL *) mycalloc(npp, sizeof(MYREAL));
    bayes->datastore2 = (long *) mycalloc((2 * (sizeg + 2)),sizeof(long));
    bayes->accept = bayes->datastore2;
    bayes->trials = bayes->accept + sizeg + 2;
    
    // records for all bayes derived values
    bayes->params = (MYREAL *) mycalloc((bayes->allocparams * (sizeg+2)),sizeof(MYREAL));
    //  bayes->params[0] = (MYREAL *) mycalloc(size+1,sizeof(MYREAL));
    
    // set counter for mdim file (typically called bayesallfile)
    if(world->cold && options->has_bayesmdimfile)
    {
        bayes->mdimfilecount = (long *) mycalloc(world->loci,sizeof(long));
    }
}


// initialize the Bayes histogram structure, adds an additional element for the
// summary over loci.
void bayes_init_histogram(world_fmt * world, option_fmt * options)
{
    long sumloc = world->loci > 1 ? 1 : 0;
    bayes_fmt *bayes = world->bayes;
    long loc;
    //long i;
    long np = get_numparam(world) - world->grownum;
    long npp = get_numparam(world); 
    //long npp = np + (bayes->mu) + 2 * world->species_model_size;
    //long nppt = np + (bayes->mu)  + 2 * world->species_model_size;
    long pa;
    long pa0;
    long bins;
    //long pa2;
    bayeshistogram_fmt *hist;
    bayes->scaling_factors = (MYREAL *) mycalloc(npp,sizeof(MYREAL));
    bayes->histtotal = (MYREAL *) mycalloc((world->loci * npp), sizeof(MYREAL));
    bayes->prettyhist = options->bayespretty;
    bayes->mdiminterval = options->bayesmdiminterval;
    bayes->histogram = (bayeshistogram_fmt *) mycalloc(world->loci + 1,sizeof(bayeshistogram_fmt));
    
    for(loc=0; loc < world->loci + sumloc; loc++)
    {
        hist = &(bayes->histogram[loc]);
        hist->bins = (long *) mycalloc(npp, sizeof(long));
        hist->binsum = 0; //DEBUG
        hist->results = NULL; //calloc(hist->binsum, sizeof(MYREAL));   // contains histogram, size is bins*numparam
        // on a per parameter basis
        // structure has a data storage vectors and the following are all pointers into it
        hist->numparam = npp;    // number of parameters: thetas + migrates + murate + growth + speciesnum
        hist->datastore = (MYREAL *) mycalloc((11*npp), sizeof(MYREAL)); // data storage, size is numparam*11
        // pointers into data storage
        hist->minima = hist->datastore;    // contains minimal values for each parameter
        hist->maxima = hist->datastore + npp;    // contains maximal values for each parameter
        hist->adjmaxima = hist->datastore + 2*npp;// holds maxima values used in histogram [are smaller than maxima]
        hist->cred50l  = hist->datastore + 3*npp;    // holds 50%-credibility margins (<all lower values>,
        hist->cred50u = hist->datastore + 4*npp;   //  <all high values>)
        hist->cred95l = hist->datastore + 5*npp;    // holds 95%-credibility margins (<all lower values>)
        hist->cred95u = hist->datastore + 6*npp;   //  <all high values>)
        hist->modes = hist->datastore + 7*npp;    // holds 95%-credibility margins (<all lower values>, <all high values>)
        hist->medians = hist->datastore + 8*npp;
        hist->means = hist->datastore + 9*npp;
        hist->stds = hist->datastore + 10*npp;
        //pa2 = 0;
        for(pa0=0; pa0 < npp; pa0++)
	  {
	    if(shortcut(pa0, world, &pa))
	      continue;
	    else
	      {
		if (pa>=np)
		  bins = options->bayes_priors[0].bins;
		else
                bins = options->bayes_priors[pa].bins;
                hist->bins[pa] = bins;
                hist->binsum += bins;
                hist->minima[pa] = (double) HUGE;
	      }
	  }
	/*
	// replace this with the above
        for(pa=0; pa < world->numpop; pa++)
        {
            if(!strchr("c", world->options->custm2[pa]))
            {
                bins = options->bayes_priors[pa].bins;
                hist->bins[pa] = bins;
                hist->binsum += bins;
                hist->minima[pa] = (double) HUGE;
            }
        }
        //pa2=0;
        for(pa=world->numpop; pa < world->numpop2; pa++)
        {
            if(!strchr("0cdt", world->options->custm2[pa]))
            {
                bins = options->bayes_priors[pa].bins;
                hist->bins[pa] = bins;
                hist->binsum += bins;
                hist->minima[pa] = (double) HUGE;
            }
            else
            {
                hist->bins[pa] = 0;
                hist->minima[pa] = (double) HUGE;
            }
        }
        if(world->bayes->mu)
        {
            bins = options->bayes_priors[pa].bins;
            hist->bins[pa] = bins;
            hist->binsum += bins;
            hist->minima[pa] = (double) HUGE;
        }
        if(world->has_speciation)
        {
            pa = world->bayes->mu + world->numpop2;
            for (i=0; i < 2 * world->species_model_size; i++)
            {
                bins = options->bayes_priors[pa+i].bins;
                hist->bins[pa+i] = bins;
                hist->binsum += bins;
                hist->minima[pa+i] = (double) HUGE;
            }
	    }*/
    }
}

//
// selects the specific set of prior parameters according to the prior options setting
// for each parameter with array_count i
MYINLINE  void select_prior_param(int selector, long i, bayes_fmt *bayes, prior_fmt *prior)
{
    bayes->minparam[i] = prior->min;
    bayes->maxparam[i] = prior->max;
    bayes->meanparam[i] = prior->mean;
    bayes->deltahist[i] = (prior->max - prior->min)/ (prior->bins);//bins-1
    switch(selector)
    {
        case MULTPRIOR:
        case WEXPPRIOR:
        case EXPPRIOR:
            bayes->priormean[i] = prior->mean; // fill mean for the call to bayes_epxb_newparam
            bayes->delta[i] =  prior->delta; //(prior->min + prior->max)/(20.); // 1/10 of the max span
            bayes->alphaparam[i] = 1.0;
            bayes->alphaorigparam[i] = 1.0;
            bayes->betaparam[i] = find_beta_truncgamma(prior->mean, 1.0, bayes->minparam[i],bayes->maxparam[i]);
            break;
        case UNIFORMPRIOR:
            bayes->priormean[i] = prior->mean; // fill mean for the call to bayes_epxb_newparam
            bayes->delta[i] =  prior->delta; //(prior->max - prior->min)/(10.); // 1/10 of the max span
            break;
        case GAMMAPRIOR:
            bayes->priormean[i] = prior->mean;
            bayes->alphaparam[i] = prior->alpha;
            bayes->alphaorigparam[i] = prior->alpha;
            bayes->delta[i] =  prior->delta; //(prior->min + prior->max)/(20.); // 1/10 of the max span
            bayes->betaparam[i] = find_beta_truncgamma(prior->mean, prior->alpha, bayes->minparam[i],bayes->maxparam[i]);
	    break;
        case NORMALPRIOR:
            bayes->priormean[i] = prior->mean;
            bayes->alphaparam[i] = prior->alpha; //check why wrote this: prior->std;
            bayes->alphaorigparam[i] = prior->alpha;
            bayes->delta[i] =  prior->delta; //(prior->min + prior->max)/(20.); // 1/10 of the max span
	    //bayes->betaparam[i] = find_beta_truncgamma(prior->mean, prior->alpha, bayes->minparam[i],bayes->maxparam[i]);
            break;
        default:
            error("Problems with the specification of the prior distribution");
            //break;
    }
}

// fill the Bayesian framework with values
void bayes_fill(world_fmt *world, option_fmt *options)
{
    long i;
    long locus;
    bayes_fmt * bayes = world->bayes;
    long numpop2 = world->numpop2;
    long np = numpop2 + world->bayes->mu + world->species_model_size * 2 + world->grownum;
    if (world->cold)
    which_prior(options->bayes_priors,np);
    
    for(i=0; i< np;i++)
    {
        select_prior_param(options->bayes_priors[i].kind, i, bayes, &options->bayes_priors[i]);
    }
    if(bayes->mu)
    {
        for(locus=0; locus < options->muloci;locus++)
        {
            world->options->mu_rates[locus] = options->bayes_priors[numpop2].mean;
            world->options->lmu_rates[locus] = log(options->bayes_priors[numpop2].mean);
        }
        for(locus=options->muloci;locus < world->loci; locus++)
        {
            world->options->mu_rates[locus] = world->options->mu_rates[options->muloci-1];
            world->options->lmu_rates[locus] = world->options->lmu_rates[options->muloci-1];
        }
        
    }
    // attach to the custm migration matrix
    bayes->custm2 = world->options->custm2;

    set_speciate_functions(world);
}

// resetting the bayes storage machinery
void bayes_reset(world_fmt * world)
{
    bayes_fmt *bayes = world->bayes;
    if(world->options->bayes_infer)
    {
        //bayes->allocparams = 1;
        //bayes->params = (MYREAL *) myrealloc(bayes->params, bayes->allocparams * size * sizeof(MYREAL));
        bayes->numparams = 0; // each locus start a new set overwriting the old, allocparam is not reset
    }
}

// free the Bayesian framework
void bayes_free(world_fmt *world)
{
    long i;
    long sumloc = world->loci > 1 ? 1 : 0;
    if(world->bayes->histogram != NULL)
    {
        for(i=0;i < world->loci + sumloc; i++)
        {
            myfree(world->bayes->histogram[i].bins);
            myfree(world->bayes->histogram[i].datastore);
            myfree(world->bayes->histogram[i].covariance);
            if(myID==MASTER && !world->data->skiploci[i])
            {
                myfree(world->bayes->histogram[i].results);
		if (world->bayes->histogram[i].results2 != NULL)
		  myfree(world->bayes->histogram[i].results2);
                //	      printf("bayes results freed for locus %li\n",i);
                myfree(world->bayes->histogram[i].set95);
                //printf("bayes set95 freed for locus %li\n",i);
                
            }
        }
        myfree(world->bayes->histogram);
    }
    myfree(world->bayes->datastore);
    myfree(world->bayes->datastore2);
    myfree(world->bayes->alphaparam);
    //betaparam and alpharorigparam are freed with alphaparam
    myfree(world->bayes->params);
    myfree(world->bayes->map);
    if(world->bayes->priors != NULL)
      myfree(world->bayes->priors);
    myfree(world->bayes);
}


// calculate the Bayesian statistics and prints the statistics
void bayes_stat(world_fmt *world, data_fmt *data)
{
  (void) data;
    MYREAL threshold;
    MYREAL mutot=0;
    long n=0;
    double mu;
    long lmu;
    double meanmu;
    bayes_fmt * bayes = world->bayes;
    bayeshistogram_fmt *hist;
    long locus;
    long l;
    long frompop=0;
    long topop=0;
    long j=0, j0, j1;
    long numpop2 = world->numpop2;
    long numparam1 = numpop2 + bayes->mu+ world->species_model_size * 2;
    long numparam = numparam1 + world->grownum;
    //long np = numparam - world->grownum;
    int fmt;
#ifdef DEBUG
    //print_bf_values(world);
#endif
    // for single locus data one is not calculating the overall values
    long lozi = world->loci > 1 ? world->loci : 0;
#ifdef LONGSUM
    long addon = (world->options->fluctuate ? (world->numpop * 3) : 0);
#else
    //xcode  long addon = 0;
#endif
    char st[7];
    char *stemp;
    
    //MYREAL ***cov;
    //long nsamp=0;
    //long n1=0;
    long z;
    stemp = (char *) mycalloc(LINESIZE,sizeof(char));
    // if hyperpriors were used then adjust the priors so that we
    // use the averaged parameters
    if(world->bayes->hyperprior)
    {
        for (z=0;z<numparam;z++)
        {
            if(bayes->hyperp[z].alphan > 1)
            {
                bayes->alphaparam[z] = bayes->hyperp[z].alpha;
                bayes->meanparam[z] = bayes->hyperp[z].mean;
                bayes->betaparam[z] = find_beta_truncgamma(bayes->meanparam[z],
                                                           bayes->alphaparam[z],
                                                           bayes->minparam[z],
                                                           bayes->maxparam[z]);
#ifdef DEBUG
                fprintf(stdout,"new prior parameter: %f %f %f\n",bayes->meanparam[z],bayes->alphaparam[z],bayes->betaparam[z]);
#endif
            }
        }
    }
    if(world->loci>1)
    {
        bayes_combine_loci(world);
    }
    else
    {
        calculate_plotter_priors(world);
    }
    // print raw histogram data into the bayesfile
    if(world->options->has_bayesfile)
    {
        print_locus_histogram_header(world->bayesfile, bayes->deltahist, bayes->custm2, world->numpop,
                                     numparam, world->options->usem, bayes->mu, world);
        for(locus=0; locus < world->loci; locus++)
        {
            if(world->data->skiploci[locus])
                continue;
            
            print_locus_histogram(world->bayesfile, world, locus, numparam);
        }
        if(world->loci>1)
        {
            print_loci_histogram(world->bayesfile, world, world->loci, numparam);
        }
    }
#ifdef PRETTY
    pdf_print_bayestable(world);
    if(world->options->allposteriors==TRUE)
    {
        for(locus=0;locus<world->loci;locus++)
            page_height = pdf_locus_histogram(world,locus);
    }
    page_height = pdf_loci_histogram(world);
#endif
    FPRINTF(world->outfile,"\n\n\nBayesian estimates\n");
    FPRINTF(world->outfile,"==================\n\n");
    FPRINTF(world->outfile,"Locus Parameter        2.5%%      25.0%%    mode     75.0%%   97.5%%     median   mean\n");
    FPRINTF(world->outfile,"-----------------------------------------------------------------------------------\n");
    for(locus=0; locus <= lozi; locus++)
    {
        if(locus<world->loci)
        {
            if(world->data->skiploci[locus])
                continue;
        }
        hist = &bayes->histogram[locus];
        if(locus == world->loci)
            strcpy(st,"  All ");
        else
            sprintf(st,"%5li ",locus + 1);
        
        for(j0=0; j0< numparam; j0++)
        {
	  if(shortcut(j0,world,&j))
	    {
                continue;
	    }
            
            if(j < world->numpop)
            {
                FPRINTF(world->outfile,"%5s ", st);
                FPRINTF(world->outfile,"Theta_%-3li      ",j0+1);
                FPRINTF(world->outfile, "%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",
                        hist->cred95l[j], hist->cred50l[j], hist->modes[j],
                        hist->cred50u[j], hist->cred95u[j],hist->medians[j], hist->means[j]);
            }
            else
            {
                if (j < numpop2)
                {
                    //1229if(world->options->custm2[j]!='d')
                    //1229  {
                    m2mm(j0,world->numpop,&frompop, &topop);
                    if(world->options->usem)
                        sprintf(stemp,"M_%li->%li", frompop+1, topop+1);
                    else
                        sprintf(stemp,"Theta_%li*M_%li->%li", topop+1, frompop+1, topop+1);
                    //1229 }
                }
                else
                {
                    if(j0==numpop2 && bayes->mu)
                        continue;
		    if(world->has_speciation && j0 < numparam1)
		      {
                    j1 = (j0 - world->numpop2 - bayes->mu);
                    z = j1 % 2;
                    j1 = j1/2;
                    species_fmt *s = &world->species_model[j1];
                    if (z==0)
                    {
                        sprintf(stemp,"D_%li->%li", s->from+1, s->to+1);
                    }
                    else
                    {
                        sprintf(stemp,"S_%li->%li", s->from+1, s->to+1);
                        //j1 = j1 + 1;
                    }
                }
		    if(world->has_growth && j0 >= numparam1)
		      {
			long d=j0-numparam1;
			while (world->options->growpops[d]==0 && d < world->options->growpops_numalloc)
			  d++;
			
			sprintf(stemp,"Growth_%-3li",d+1);
		      }
		}
                //fmt = 2;
                int dig = finddigits((long) hist->modes[j0]);
                switch (dig)
                {
                    case 1:
                    case 2:
                        fmt= 5;
                        break;
                    case 3:
                        fmt = 4;
                        break;
                    default:
                        fmt = 2;
                        break;
                }
                FPRINTF(world->outfile,"%5s ", st);
                FPRINTF(world->outfile, "%-15.15s",stemp);
                FPRINTF(world->outfile,"%8.*f %8.*f %8.*f %8.*f %8.*f %8.*f %8.*f\n",
                        fmt, hist->cred95l[j], fmt, hist->cred50l[j], fmt, hist->modes[j],
                        fmt, hist->cred50u[j], fmt, hist->cred95u[j], fmt, hist->medians[j], fmt,
                        hist->means[j]);
                // warning block: issuing a warning when the 75% and the95% percentiles are within 10% of the
                // upper prior boundary
                threshold =  0.9 * bayes->maxparam[j];
                if(hist->cred95u[j] > threshold && hist->cred50u[j] > threshold)
                {
                    if(locus == lozi && locus>0)
                        record_warnings(world,"Param %li (all loci): Upper prior boundary seems too low! ",j+1);
                    else
                        record_warnings(world,"Param %li (Locus %i): Upper prior boundary seems too low! ", j+1,locus+1);
                }
            }
        }
        if(bayes->mu)
        {
            FPRINTF(world->outfile,"%5.5s ", st);
            j0=world->numpop2;
            if(locus==lozi && lozi>1)
            {
                meanmu = 0.;
                for(l=0;l<world->loci;l++)
                {
                    meanmu += world->options->meanmu[l];
                }
                meanmu /= world->loci;
                if (meanmu < SMALL_VALUE)
                    meanmu = 10e-8;
                lmu = (long) (floor( log10(hist->modes[j0] * meanmu)));
                mu = meanmu * pow(10. , -lmu);
            }
            else
            {
	      lmu = (long) (floor( log10(hist->modes[j0] * world->options->meanmu[locus])));
                mu = world->options->meanmu[locus] * pow(10., -lmu);
            }
            FPRINTF(world->outfile," mu [10^%li]     ",lmu);
            FPRINTF(world->outfile, "%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",
                    mu*hist->cred95l[j0], mu*hist->cred50l[j0] , mu*hist->modes[j0],
                    mu*hist->cred50u[j0], mu*hist->cred95u[j0] , mu*hist->medians[j0],
                    mu*hist->means[j0]);
        }
    } // over all +1 loci
    FPRINTF(world->outfile,"-----------------------------------------------------------------------------------\n");
    if(bayes->mu)
    {
        mutot = 0;
        n = 1;
        for(locus=0; locus < world->loci; locus++)
        {
            if(world->data->skiploci[locus])
                continue; //invariantloci: not resolved yet
            mu = (double) world->options->meanmu[locus];
            mutot += (mu - mutot) / n;
            n++;
            FPRINTF(world->outfile,"(*) Mutation rate for locus %li was set to %8.5e", locus, mu);
        }
        if(world->loci>1)
        {
            FPRINTF(world->outfile,"(*) Average mutation rate for all loci is  %8.5e", mutot);
        }
    }
    //if (world->bayes->histogram[0].covariance != NULL)
    //{
    //    cov = (MYREAL ***) mycalloc(world->loci+1,sizeof(MYREAL **));
    //    for (locus=0;locus<world->loci;locus++)
    //    {
    //        cov[locus] = world->bayes->histogram[locus].covariance;
    //        n1 = world->bayes->histogram[locus].n;
    //       nsamp += n1/world->loci;
    //adjust_covariance(cov[locus],world->numpop2, n1);
    //    }
    //if(world->loci>1)
    //{
    //    cov[locus] = world->bayes->histogram[locus].covariance;
    //    //adjust_covariance(cov[locus],world->numpop2, nsamp);
    // }
    //print_cov2 (world, world->numpop, world->loci, cov);
    //myfree(cov);
    // print out the acceptance ratios for every parameter and the tree
    //}
    if(world->options->progress)
    {
        // final acceptance ratios for the Bayesian run
        bayes_print_accept(stdout,world);
    }
    myfree(stemp);
}


//
// print out the acceptance ratios for all the different Bayesian updates
void
bayes_print_accept(FILE * file,  world_fmt *world)
{
  long j0, j;          
    long topop    = 0;  
    long frompop  = 0;  
    char *stemp;       
    long trials   = 0;
    long accept   = 0;
    long tc = world->numpop2 + world->bayes->mu + world->species_model_size * 2;
    long tcg = tc+world->grownum;
    bayes_fmt *bayes = world->bayes;
    //long estimated_trials=0;
    //species_fmt *s = NULL;
    //long z=0;
    FPRINTF(file,"\n\n\nAcceptance ratios for all parameters and the genealogies\n");
    FPRINTF(file,"---------------------------------------------------------------------\n\n");
    // This needs more attention but will need more stuff to safe
    //if(world->options->checkpointing)
    //{
    //    FPRINTF(file,"\n\nMay not work correctly\n\n");
    //}
    
    FPRINTF(file,"Parameter           Accepted changes               Ratio\n");
    // population sizes
    for(j0=0; j0 < world->numpop; j0++)
    {
      //        if(!strchr("0c", bayes->custm2[j]))
      if(shortcut(j0,world,&j))
        {
	  continue;
	}
      else
	{
	  accept = world->accept_archive[j];
	  trials = world->trials_archive[j];
	  if(trials > 0)
            {
	      FPRINTF(file,"Theta_%-3li             %8li/%-8li         %8.5f\n", j+1,
		      accept, trials, (MYREAL) accept/trials);
	      //estimated_trials += trials;
            }
        }
    }
    // migration rates
    stemp = (char *) mycalloc(LINESIZE,sizeof(char));
    for(j0=world->numpop; j0 < world->numpop2; j0++)
    {
      //if(!strchr("0cd", bayes->custm2[j]))
      if(shortcut(j0,world,&j))
        {
	  continue;
	}
      else
	{
	  accept = world->accept_archive[j];
            trials=world->trials_archive[j];
            if(trials>0)
            {
                m2mm (j, world->numpop, &frompop, &topop);
                if(world->options->usem)
                {
                    sprintf(stemp, "M_%li->%li", frompop+1, topop+1);
                }
                else
                {
                    sprintf(stemp, "xN_%lim_%li->%li", topop+1, frompop+1, topop+1);
                }
                FPRINTF(file, "%-12.12s          %8li/%-8li         %8.5f\n", stemp,
			accept, trials, (MYREAL) accept/trials);
                //estimated_trials += trials;
            }
        }
    }
    // accepted rate of mutation rate changes for each locus mutation rate
    if(bayes->mu)
    {
        //      for (j=world->numpop2;j<world->numpop2 + bayes->mu*world->loci;j++)
        for (j=world->numpop2;j<world->numpop2 + bayes->mu;j++)
        {
            FPRINTF(file, "Rate of mutation rate (%li) %8li/%-8li         %8.5f\n", j+1,
                    world->accept_archive[j],
                    world->trials_archive[j],
		    (MYREAL) world->accept_archive[j]/world->trials_archive[j]);
            //estimated_trials += world->trials_archive[j];
        }
    }
    // accepted species events
    if (world->has_speciation)
    {
      //z=0;
        //	for(j=world->numpop2+world->loci * bayes->mu;j <= tc;j+=2)
        for(j0=world->numpop2+bayes->mu;j0 < tc;j0++)
        {
	  if(shortcut(j0,world,&j))
	    {
	      continue;
	    }
	  else
	    {
	      species_fmt * s = get_which_species_model(j, world->species_model, world->species_model_size);
	      long from = s->from;
	      long to = s->to;
	      accept = world->accept_archive[j];
	      trials=world->trials_archive[j];
	      sprintf(stemp,"%c_%li->%li",j==s->paramindex_mu ? 'D' : 'S',1+from,1+to);
	      FPRINTF(file, "D%-12.12s          %8li/%-8li         %8.5f\n", stemp, accept,
		      trials, (MYREAL) accept/trials);
	      //estimated_trials += trials;
	    }
	}
    }
    // accepted growth
    if(world->has_growth)
    {

        for (j=tc; j<tcg; j++)
        {
	  long d=j-tc;
	  while (world->options->growpops[d]==0 && d < world->options->growpops_numalloc)
	    d++;
	  FPRINTF(file, "Growth_%li:               %8li/%-8li         %8.5f\n", d+1,
		  world->accept_archive[j],
		  world->trials_archive[j],
		  (MYREAL) world->accept_archive[j]/world->trials_archive[j]);
	  //estimated_trials += world->trials_archive[j];
        }
    }

    // accepted trees
    trials=world->trials_archive[tcg];
    if(trials>0)
    {
        FPRINTF(file, "Genealogies           %8li/%-8li          %8.5f\n", world->accept_archive[tcg], (long)
                trials, (MYREAL) world->accept_archive[tcg]/trials);
        //estimated_trials += trials;
    }
    //    FPRINTF(file, "Sum of all acceptances    %8li/%-8li\n", estimated_trials, world->options->lsteps*world->maxreplicate*world->loci*world->options->lincr);
    
    myfree(stemp);
}


//
// print out the hyperprior information if any
void
bayes_print_hyperprior(FILE * file,  world_fmt *world)
{
  //WARNING not adapted to growth yet
    long j;             //used to loop over all parameters
    long topop    =0;   // indicator into the parameter vector, specifying originating population
    long frompop  =0;   // receiving population
    char *stemp;       // string variable holding print-string
    long trials   =0;   //
    long tc = world->numpop2 + world->bayes->mu + world->species_model_size * 2;
    bayes_fmt *bayes = world->bayes;
    hyper_fmt *hyperp = bayes->hyperp;
    //long estimated_trials=0;
    long z=0;
    
    if (!world->bayes->hyperprior)
        return; //no hyperprior used
    
#ifndef MPI
    for (j=0;j<world->bayes->mapsize;j++)
    {
        if(hyperp[j].meann > 1)
        {
            onepass_mean_std_end(&hyperp[j].mean,&hyperp[j].meanstd,&hyperp[j].meann);
            onepass_mean_std_end(&hyperp[j].alpha,&hyperp[j].alphastd,&hyperp[j].alphan);
        }
    }
#endif
    
    FPRINTF(file,"\n\n\nHyperprior information for all parameters\n");
    FPRINTF(file,"---------------------------------------------------------------------\n\n");
    
    FPRINTF(file,"Parameter        Priormean/std     PriorAlpha/std       N\n");
    for(j=0; j < world->numpop; j++)
    {
        trials=hyperp[j].meann;
        if(trials>1)
        {
            FPRINTF(file,"Theta_%-3li       %8.5f/%-8.5f    %8.5f/%-8.5f     %li\n", j+1, hyperp[j].mean,hyperp[j].meanstd,hyperp[j].alpha,hyperp[j].alphastd,trials);
        }
    }
    // migration rates
    stemp = (char *) mycalloc(LINESIZE,sizeof(char));
    for(j=world->numpop; j < world->numpop2; j++)
    {
        trials=hyperp[j].meann;
        if(trials>1)
        {
            m2mm (j, world->numpop, &frompop, &topop);
            if(world->options->usem)
            {
                sprintf(stemp, "M_%li->%li", frompop+1, topop+1);
            }
            else
            {
                sprintf(stemp, "xN_%lim_%li->%li", topop+1, frompop+1, topop+1);
            }
            FPRINTF(file,"%s          %8.5f/%-8.5f %8.5f/%-8.5f     %li\n", stemp, hyperp[j].mean,hyperp[j].meanstd,hyperp[j].alpha,hyperp[j].alphastd,trials);
        }
    }
    // accepted rate of mutation rate changes for each locus mutation rate
    if(bayes->mu)
    {
        trials=hyperp[world->numpop2].meann;
        //      for (j=world->numpop2;j<world->numpop2 + bayes->mu*world->loci;j++)
        for (j=world->numpop2;j<world->numpop2 + bayes->mu;j++)
        {
            FPRINTF(file,"Rate of mutation rate (%li)  %8.5f/%-8.5f         %8.5f/%-8.5f     %li\n", j+1, hyperp[j].mean,hyperp[j].meanstd,hyperp[j].alpha,hyperp[j].alphastd,trials);
        }
    }
    // accepted species events
    if (world->has_speciation)
    {
        z=0;
        //	for(j=world->numpop2+world->loci * bayes->mu;j <= tc;j+=2)
        for(j=world->numpop2+bayes->mu;j < tc;j+=2)
        {
            trials=hyperp[j].meann;
            long from = world->species_model[z].from;
            long to = world->species_model[z].to;
            z++;
            if(trials>1)
            {
                sprintf(stemp,"_%li->%li",1+from,1+to);
                FPRINTF(file,"D%s          %8.5f/%8.5f %8.5f/%8.5f     %li\n",
                        stemp, hyperp[j].mean, hyperp[j].meanstd,
                        hyperp[j].alpha, hyperp[j].alphastd, trials);
                
            }
            trials=hyperp[j+1].meann;
            if(trials>1)
            {
                FPRINTF(file,"S%s          %8.5f/%8.5f %8.5f/%8.5f     %li\n",
                        stemp, hyperp[j].mean,hyperp[j].meanstd,
                        hyperp[j].alpha,hyperp[j].alphastd,trials);
            }
        }
    }
    myfree(stemp);
}


//
// print out the acceptance ratios for all the different Bayesian updates
void
bayes_progress(world_fmt *world, long ten)
{
  (void) ten;
  const boolean writelog = world->options->writelog; 
    const boolean progress = world->options->progress;
    long np = world->numpop2 + world->bayes->mu + 2 * world->species_model_size;// growth is dealt independently
    char *buffer;
    long bufsize=0;
    char spacer[]="";
    //
    char nowstr[STRSIZE]; // will hold time of day
    //    char temp[LINESIZE];  // for locus printing with rates
    long j0;              // used to loop over all parameters
    long j;
    char *paramstr;       // string variable holding print-string
    bayes_fmt *bayes = world->bayes;
    MYREAL * autocorr = world->autocorrelation;
    MYREAL * effsample = world->effective_sample;
    MYREAL param0;
    MYREAL accrat=0.0;
    buffer = (char *) mycalloc(bayes->progresslinesize,sizeof(char));
    paramstr = (char *) mycalloc(LINESIZE,sizeof(char));
    
    get_time (nowstr, "%H:%M:%S");
    
    //if(world->options->heating)
    //{
    //    bufsize += sprintf(buffer+bufsize,
    //                       "\n%s   [NODE:%i, Locus: %li, Heating:ON (Swaps: %li)] ",
    //                       nowstr,myID, world->locus + 1, world->swapped);
    //}
    //else
    //{
    //    bufsize += sprintf(buffer+bufsize,"\n%s   [NODE:%i, Locus: %li] ", nowstr,myID, world->locus + 1);
    //}
    bufsize += sprintf(buffer+bufsize,"%s   ",nowstr);
    prognose_time (nowstr, world, 1, world->treesdone, spacer, TRUE);
    bufsize += sprintf(buffer+bufsize, "Sampling: prognosed end of run is %s [%f done])\n", nowstr,
                       ((MYREAL) world->treesdone / (MYREAL) world->treestotal));
    
    if(world->options->has_autotune)
    {
        bufsize += sprintf(buffer+bufsize,"           Parameter     Acceptance Current      PropWindow AutoCorr ESS\n");
        bufsize += sprintf(buffer+bufsize,"           ------------  ---------- ------------ ---------- -------- --------\n");
        //                                            123456789012  1234567890 123456789012 1234567890 12345678 12345678
    }
    else
    {
        bufsize += sprintf(buffer+bufsize,"           Parameter     Acceptance Current      AutoCorr ESS\n");
        bufsize += sprintf(buffer+bufsize,"           ------------  ---------- ------------ -------- --------\n");
    }
    for(j0=0; j0 < np; j0++)
    {
        if(shortcut(j0,world,&j))
        {
            continue;
        }
        else
        {
            set_paramstr(paramstr, j,world);
            param0 = world->param0[j];
            
            //if(world->options->has_autotune)
            //{
            if(bayes->trials[j]>0)
                accrat = (MYREAL) bayes->accept[j]/bayes->trials[j];
            else
                accrat = 0.0;
            double propwindow = bayes->delta[j0];
	    if(world->options->has_autotune)
	      {
		bufsize += sprintf(buffer+bufsize, "           %-12.12s     % 5.3f    % 10.5f  %9.3f % 8.3f % 8.2f\n",
                               paramstr, accrat, param0, propwindow, autocorr[j],effsample[j]);
	      }
	    else
	      {
		bufsize += sprintf(buffer+bufsize, "           %-12.12s     % 3.2f % 10.5f % 8.3f % 8.2f\n",
				   paramstr, accrat, param0, autocorr[j],effsample[j]);
	      }
        }
    }
    if (world->has_growth)
      {
	long g;
	for(g=0;g<world->grownum;g++)
	  {
	    double propwindow = bayes->delta[np+g];
	    long d=g;
	    while (world->options->growpops[d]==0 && d < world->options->growpops_numalloc)
	      d++;			
	    sprintf(paramstr,"Growth_%li",d+1);
	    j = np + g;
	    if(bayes->trials[j]>0)
	      accrat = (MYREAL) bayes->accept[j]/bayes->trials[j];
            else
	      accrat = 0.0;
	    if(world->options->has_autotune)
	      {
		bufsize += sprintf(buffer+bufsize, "           %-12.12s     % 5.3f    % 10.5f  %9.3f % 8.3f % 8.2f\n",
				   paramstr, accrat, world->growth[g], propwindow, autocorr[np+g],effsample[np+g]);
	      }
	    else
	      {
		bufsize += sprintf(buffer+bufsize, "           %-12.12s     % 3.2f % 10.5f % 8.3f % 8.2f\n",
				   paramstr, accrat, world->growth[g], autocorr[np+g], effsample[np+g]);
	      }
	  }
	j0 = np+world->grownum; //the genealogy is the last element in autocorr and effsample
      }
    if(world->options->has_autotune)
    {
        sprintf(buffer+bufsize, "           Genealogies      % 5.3f % 13.5f      ---   % 8.3f % 8.2f\n",
                (MYREAL) bayes->accept[j0] / (MYREAL) bayes->trials[j0],world->likelihood[world->numlike-1], autocorr[j0],effsample[j0]);
    }
    else
    {
        sprintf(buffer+bufsize, "           Genealogies      % 5.3f % 13.5f % 8.3f % 8.2f\n",
                (MYREAL) bayes->accept[j0]/bayes->trials[j0],world->likelihood[world->numlike-1], autocorr[j0],effsample[j0]);
    }
    if(progress)
    {
        FPRINTF(stdout,"%s",buffer);
    }
    
    if(writelog)
    {
        FPRINTF(world->options->logfile,"%s",buffer);
    }
    myfree(buffer);
    myfree(paramstr);
    //print_bayes_ess(stdout, world, world->autocorrelation, world->effective_sample); //offset was 11
}




// adjusts the allocations for the histogram bins for the current locus
void adjust_bayes_bins(world_fmt * world, long locus)
{
    bayeshistogram_fmt *hist = &(world->bayes->histogram[locus]);
    if(world->bayes->histogram[locus].results == NULL)
    {
        world->bayes->histogram[locus].results = (double *) mycalloc(hist->binsum, sizeof(double));
        world->bayes->histogram[locus].set95 = (char *) mycalloc(hist->binsum * 2, sizeof(char));
        world->bayes->histogram[locus].set50 = world->bayes->histogram[locus].set95 + hist->binsum;
    }
}


// assembles the histogram from the sampled parameter values in params
// there is a similar function that reads from the bayesallfile mdimfile
// synonym: construct_locus_histogram
void construct_param_hist(world_fmt *world, long locus, long npa, long pa, long numbin, MYREAL *mini,
                          MYREAL *maxi, double **results,
                          long *total, MYREAL *themean, MYREAL *thestd)
{
  (void) numbin;
    bayes_fmt *bayes = world->bayes;
    long      j, j0;
    long      i;
    long      np = npa;
    long      npx = np + 2;
    long      bin;
    long      nb;
    
    MYREAL    delta; 
    MYREAL    value;
    MYREAL    *params = bayes->params;
    long    floorindex;
    float *p = (float *) mycalloc(bayes->numparams,sizeof(float));
    //float *q = (float *) mycalloc(bayes->numparams,sizeof(float));
    for(i=0;i < bayes->numparams; i++)
    {
        floorindex = npx * i + 2;
        value = params[floorindex+pa];
        p[i] = (float) value;
    }
    qsort(p, (size_t) bayes->numparams, sizeof(float), floatcmp);
    delta = bayes->deltahist[pa];
    for(i=0;i < bayes->numparams; i++)
    {
      value = (MYREAL) p[i]; //params[floorindex+pa];
      *themean += value;
      *thestd += value * value;
      if(value > maxi[pa])
        {
            bin = bayes->histogram[locus].bins[pa] - 1;
        }
        else
        {
            bin = (long) ((value - mini[pa])/delta);
            if(bin<0)
                bin=0;
        }
        if((bin) > bayes->histogram[locus].bins[pa])
        {
            warning("%i> value not counted for histogram: bin=%li > histbins=%li\n", myID,bin, bayes->histogram[locus].bins[pa]);
            continue;
        }
        nb = 0;
        for(j0=0;j0<pa;j0++)
        {
	  if(shortcut(j0,world,&j))
            {
                continue;
            }
            nb += bayes->histogram[locus].bins[j];
        }
        (*results)[nb + bin] += 1.;
        *total += 1;
    }
    //
    // adjusting the integral =1
    nb=0;
    for(j0=0;j0<pa;j0++)
      {
	if(shortcut(j0,world,&j))
	  {
	    continue;
	  }
	nb += bayes->histogram[locus].bins[j];
      }    
    for(bin=0;bin < bayes->histogram[locus].bins[pa]; bin++)
    {
        (*results)[nb + bin] /= (double)(*total)  * delta;
    }
    myfree(p);
    //myfree(q);
}

// construct bayes histogram: make a at least BAYESNUMBIN slices through the min-max range for each locus
// while calculating the histogram calculate also the mean and standard deviation.
// This is done in the same loop, but is somehwat opaque. This is used from bayesallfile!
// Synonym:  construct_param_hist() for ram based parameters
void construct_locus_histogram(world_fmt *world, long locus, MYREAL *mini, MYREAL *maxi, double **results)
{
    long pa;
    long pa0;
    long rpa;
    long total=0;
    long numbin=0;
    boolean *visited;
    bayes_fmt *bayes = world->bayes;
    long npa = world->numpop2 + bayes->mu;
    long np = npa + 2 * world->species_model_size;
    long npg = np + world->grownum;
    MYREAL themean=0.0;
    MYREAL thestd=0.0;
    visited = (boolean *) mycalloc(npg, sizeof(boolean));
    for(pa0=0; pa0 < npg; pa0++)
    {
      if(shortcut(pa0,world,&pa))
        {
            continue;
        }
      rpa = pa;
        // handles cases with symmetric migration rates and mean of migration rates
        // TODO needs revision for new scheme that allows more freedom in setting migration rate groups
      if(!visited[rpa])
        {
            themean = 0.0;
            thestd = 0.0;
            total = 0;
            construct_param_hist(world,locus,npg, rpa,numbin, mini, maxi, results, &total,&themean,&thestd);
            world->bayes->histogram[locus].means[pa] = themean/total;
            world->bayes->histogram[locus].stds[pa]  = thestd / total;
            world->bayes->histtotal[locus*npg+pa] = (MYREAL) total;
            //printf("total=%li\n",total);
	    //printf("@@1@@@ locus=%li parameter=%li total=%li mean=%f \"std\"=%f\n",locus, rpa, total,themean/total,thestd/total);
            //	  printf("histtotal=%li locus=%li pa=%li\n", total, locus, pa);
            visited[rpa] = TRUE;
        }
    }
    //if (world->has_growth) // dealth with above?
    //  {
    //    construct_locusgrowth_histogram(world, locus, mini, maxi, results);
    //  }
    covariance_bayes(world,locus);
    myfree(visited);
}

//
// find the largest interval for all loci looking at the stored minima and maxima of the histogram
// and reports these into adjmini and adjmaxi
// this allows combination of different histograms without storing all raw data
void adjust_bayes_min_max(world_fmt* world, MYREAL **mini, MYREAL **maxi, MYREAL **adjmini, MYREAL **adjmaxi)
{
    //MYREAL delta;
    long pa0, pa;
    const long numpop2 = world->numpop2;
    const long np = numpop2 + (world->bayes->mu) + 2 * world->species_model_size + world->grownum;
    
    for(pa0=0; pa0 < np; pa0++)
    {
      // if custom migration matrix is set to zero
      // continue
      if(shortcut(pa0,world,&pa))
        {
            continue;
        }
        
        if((*mini)[pa] < (*adjmini)[pa0])
        {
            (*adjmini)[pa0] = (*mini)[pa];
        }
        if((*maxi)[pa] > (*adjmaxi)[pa0])
        {
            (*adjmaxi)[pa0] = (*maxi)[pa];
        }
    }
}


// find minimum and maximum values for each parameters
void find_bayes_min_max(world_fmt* world, MYREAL **mini, MYREAL **maxi, MYREAL **adjmaxi)
{
  (void) adjmaxi;
    long pa0, pa;
    long np = world->numpop2 + (world->bayes->mu) + 2 * world->species_model_size + world->grownum;
    
    for(pa0=0; pa0 < np; pa0++)
    {
        if(shortcut(pa0,world,&pa))
        {
            continue;
        }
        (*mini)[pa0] =  world->bayes->minparam[pa];
        (*maxi)[pa0] =  world->bayes->maxparam[pa];
    }
}

//
// prints order of parameters for header files in bayesfile and bayesallfile
void print_param_order(char **buf, long *bufsize, long *allocbufsize, world_fmt *world, long numparam)
{
  long i;
    long pa;
    long mypa;
    long pa0;
    long numpop = world->numpop;
    long numpop2 = world->numpop2;
    long frompop, topop;
    long frompop2, topop2;
    char *custm2 = world->options->custm2;
    boolean usem = world->options->usem;
    bayes_fmt *bayes = world->bayes;
    *bufsize += sprintf(*buf + *bufsize,"# Order of the parameters:\n");
    *bufsize += sprintf(*buf + *bufsize,"# Parameter-number Parameter\n");
    for(pa0=0;pa0<numparam;pa0++)
    {
        if(*bufsize > (*allocbufsize - 150))
        {
            *allocbufsize += LINESIZE;
            *buf = (char *) realloc(*buf, sizeof(char) * (size_t) (*allocbufsize));
        }
        if(shortcut(pa0,world,&pa))
        {
            continue;
        }
        if(pa0 < numpop)
        {
            if((pa0 == pa) && (custm2[pa0] == '*'))
                *bufsize += sprintf(*buf + *bufsize,"#@ %6li    %s_%li\n", pa0+1, "Theta",pa0+1);
            else
                *bufsize += sprintf(*buf + *bufsize,"#@ %6li    %s_%li = %s_%li   [%c]\n",
                                    pa0+1, "Theta",pa0+1, "Theta",pa+1, custm2[pa0]);
        }
        else
        {
            // do we estimate mutation rate changes?
            if(bayes->mu && pa0 == numpop2)
            {
                *bufsize += sprintf( *buf + *bufsize,"#@ %6li    %s\n", pa0+1, "Rate");
            }
            else
            {
	      	  if(pa0<numpop2)
		    {
		      m2mm(pa0,numpop,&frompop,&topop);
		      if((pa0==pa) && (custm2[pa0]=='*'))
			{
			  if(usem)
			    *bufsize += sprintf( *buf + *bufsize,"#@ %6li    %s_(%li,%li)\n", pa+1, "M", frompop+1, topop+1);
			  else
			    *bufsize += sprintf( *buf + *bufsize,"#@ %6li    %s_(%li,%li) = %s_(%li,%li)*%s_%li\n", pa+1, "xNm", frompop+1, topop+1, "M", frompop+1, topop+1, "Theta", topop+1);
			}
		      else
			{
			  m2mm(pa,numpop,&frompop2,&topop2);
			  if(usem)
			    *bufsize += sprintf(*buf + *bufsize,"#@ %6li    %s_(%li,%li) =  %s_(%li,%li) [%c]\n", pa+1, "M", frompop+1, topop+1, "M", frompop2+1, topop2+1, custm2[pa0]);
			  else
			    *bufsize += sprintf(*buf + *bufsize,"#@ %6li    %s_(%li,%li) = %s_(%li,%li)*%s_%li  = %s_(%li,%li) = %s_(%li,%li)*%s_%li [%c]\n",
						pa+1, "xNm", frompop+1, topop+1, "M", frompop+1, topop+1, "Theta", topop+1,
						"xNm", frompop2+1, topop2+1, "M", frompop2+1, topop2+1, "Theta", topop2+1,
						custm2[pa0]);
			}
		    }
	    }
	}
    }
    if (world->has_speciation)
    {
      mypa = world->numpop2 + world->bayes->mu;
      for(pa0=mypa;pa0<mypa + 2 * world->species_model_size;pa0++)
        {
	  if(shortcut(pa0,world,&pa))
	    {
	      continue;
	    }
	  else
	    {
	      species_fmt * s = get_which_species_model(pa,world->species_model,world->species_model_size);
	      frompop = s->from;
	      topop = s->to;
	      if(pa == s->paramindex_mu)
		*bufsize += sprintf( *buf + *bufsize,"#@ %6li    %s_(%li,%li)\n", pa+1, "D", frompop+1, topop+1);
	      else
		*bufsize += sprintf( *buf + *bufsize,"#@ %6li    %s_(%li,%li)\n", pa+1, "S", frompop+1, topop+1);
	    }
	}
    }
    if (world->has_growth)
      {
	*bufsize += sprintf(*buf+ *bufsize,"#@ Growpopulations:%li: ",world->options->growpops_numalloc);
	for(i=0;i<world->options->growpops_numalloc;i++)
	  {
	    *bufsize += sprintf(*buf+ *bufsize," %li", world->options->growpops[i]);
	  }
	*bufsize += sprintf(*buf+ *bufsize,"\n");
	for(i=0;i<world->grownum;i++)
	  {
	    *bufsize += sprintf(*buf+ *bufsize,"#@ %s_%li\n","Growth", i);
	  }
      }
}
//
// prints a comment header, using shell script comments for the output of the raw histogram data for bayesdata
//
// # Raw data for the histogram of the posterior probabilities for all parameters and loci\n
// # produced by the program migrate-n VERSIONNUMBER (popgen.csit.fsu.edu/migrate.hml)\n
// # written by Peter Beerli 2004, Tallahassee, if you have problems email to beerli@fsu.edu\n
// #
// # The HPC values are indicators whether the parametervalue is in the highest-posterior credibility set,\n
// # a 0 means it is outside and a 1 means the value is inside the credibility set.\n
// #
// # --------------------------------------------------------------\n
// # Locus Parameter 50%HPC  95%HPC parameter-value count frequency\n
// # --------------------------------------------------------------\n
// #
void print_locus_histogram_header(FILE *bayesfile, MYREAL *deltas, char *custm2, long numpop, long numparam, boolean usem, boolean mu, world_fmt *world)
{
  (void) numpop;
  (void) usem;
    long pa;
    long bufsize = 0;
    long allocbufsize = LONGLINESIZE;
    char *buf = (char *) mycalloc(allocbufsize, sizeof(char));
    fprintf(bayesfile, "# Raw data for the histogram of the posterior probabilities for all parameters and loci\n");
    fprintf(bayesfile, "# produced by the program migrate-n %s (http://popgen.sc.fsu.edu/migrate.hml)\n",MIGRATEVERSION);
    fprintf(bayesfile, "# written by Peter Beerli 2004-2013, Tallahassee, if you have problems email to beerli@fsu.edu\n");
    fprintf(bayesfile, "#\n");
    fprintf(bayesfile, "# The HPC values are indicators whether the parametervalue is in the highest-posterior credibility set,\n");
    fprintf(bayesfile, "# a 0 means it is outside and a 1 means the value is inside the credibility set.\n");
    fprintf(bayesfile, "#\n");
    print_param_order(&buf, &bufsize,&allocbufsize, world,numparam);
    fprintf(bayesfile, "#%s\n",buf);
    fprintf(bayesfile, "# Delta for Theta and M ");
    for(pa=0;pa<numparam-1; pa++)
    {
        fprintf(bayesfile,"%f ", custm2[pa]!='0' ? deltas[pa] : -99);
    }
    if(!mu)
        fprintf(bayesfile,"%f ", custm2[pa]!='0' ? deltas[pa] : -99);
    
    fprintf(bayesfile, "\n# -----------------------------------------------------------------------------------------\n");
    fprintf(bayesfile, "# Locus Parameter 50%%HPC  95%%HPC parameter-value count frequency cummulative_freq priorfreq\n");
    fprintf(bayesfile, "# -------------------------------------------------------------------------------------------\n");
    myfree(buf);
}

//
// print the locus data for a histogram to the file bayesfile the data has a header starting with # so that
// other programs easiliy can remove ot for processing. the function calculates the total of all bins and
// then is printing locusnumber parameternumber 95%HPC 50%HPC bincounts frequency for all bins and parameters
// the HPC columns are 0 and 1, where 0 means no in the credibiliity set.
// The header is printed in print_locus_histogram_header(bayesfile)
void print_locus_histogram(FILE *bayesfile, world_fmt *world, long locus, long numparam)
{
    //  char indicator;
    bayes_fmt *bayes = world->bayes;
    long bin;
    long pa0, pa;
    long numbins = 0;
    long numbinsall = 0;
    long counterup;
    long counterlow=0;
    //long j;
    MYREAL delta;
    MYREAL value;
    //MYREAL total=0.;
    MYREAL freq=0.;
    MYREAL sumfreq=0.;
    long *bins = bayes->histogram[locus].bins;
    double *results = bayes->histogram[locus].results;
    MYREAL *mini = bayes->histogram[locus].minima;
    char *set50 = bayes->histogram[locus].set50;
    char *set95 = bayes->histogram[locus].set95;
    if(results == NULL)
        return;
    for(pa0=0; pa0 < numparam; pa0++)
    {
        if(shortcut(pa0,world,&pa))
	  {
            continue;
	  }
        delta = bayes->deltahist[pa];
        //(maxi[pa] - mini[pa])/bins[pa]@@@@;
        value = mini[pa] + delta/2.;
        //total = 0.;
        counterup=0;
        counterlow = 0;
        //	numbins = 0;
        //for(j=0; j < pa; j++)
        //  numbins += bins[j];
        if(pa<pa0)
            continue;
        numbinsall += bins[pa];
        numbins = numbinsall - bins[pa];
        for(bin=0;bin<bins[pa];bin++)
        {
	  //total += results[numbins + bin];
            
            if(set95[numbins+bin]=='1')
            {
                counterlow++;
                break;
            }
        }
        for(bin=counterlow;bin<bins[pa];bin++)
        {
	  //total += results[numbins + bin];
            if(set95[numbins+bin]=='1')
                counterup++; //upper 97.5
        }
        sumfreq = 0.;
        for(bin=0;bin<bins[pa];bin++)
        {
            freq = results[numbins + bin];// /total;
            sumfreq += freq;
            fprintf(bayesfile,"%li %li %c %c %f %f %f %f %f\n", locus+1, pa+1, set50[numbins+bin], set95[numbins+bin], value,
                    results[numbins + bin],freq,sumfreq, delta * exp(scaling_prior(world,pa,value)));
            value += delta;
        }
        fprintf(bayesfile,"\n");
    }
}


//
// print the loci data for a histogram to the file bayesfile the data has a header starting with # so that
// other programs easiliy can remove ot for processing. the function calculates the total of all bins and
// then is printing locusnumber parameternumber 95%HPC 50%HPC bincounts frequency for all bins and parameters
// the HPC columns are 0 and 1, where 0 means no in the credibiliity set.
// The header is printed in print_locus_histogram_header(bayesfile)
void print_loci_histogram(FILE *bayesfile, world_fmt * world, long locus, long numparam)
{
    bayes_fmt *bayes=world->bayes;
    long bin;
    //long j;
    long pa0, pa;
    long numbins = 0;
    long numbinsall = 0;
    MYREAL delta;
    MYREAL value;
    //MYREAL total=0.;
    MYREAL freq=0.;
    MYREAL sumfreq=0.;
    long *bins = bayes->histogram[locus].bins;
    double *results = bayes->histogram[locus].results;
    MYREAL *mini = bayes->histogram[locus].minima;
    //MYREAL *maxi = bayes->histogram[locus].maxima;
    char *set50 = bayes->histogram[locus].set50;
    char *set95 = bayes->histogram[locus].set95;
    
    for(pa0=0; pa0 < numparam; pa0++)
    {
        if(shortcut(pa0,world,&pa))
        {
            continue;
        }
        
        delta = bayes->deltahist[pa];//(maxi[pa] - mini[pa])/bins[pa];
        value = mini[pa] + delta/2.;
        sumfreq =0.;
        
        //numbins = 0;
        //for(j=0; j < pa; j++)
        //	numbins += bins[j];
        if(pa<pa0)
            continue;
        numbinsall += bins[pa];
        numbins = numbinsall - bins[pa];
        
        for(bin=0;bin<bins[pa];bin++)
        {
            freq = (MYREAL) results[numbins + bin];
            sumfreq += freq;
            //total = (MYREAL) bayes->trials[pa];
            fprintf(bayesfile,"%li %li %c %c %f %f %f %f %f\n", locus+1, pa+1, set50[numbins+bin], set95[numbins+bin], value,
                    results[numbins+bin] , freq, sumfreq, bayes->priors[numbins+bin]);
            value += delta;
        }
        fprintf(bayesfile,"\n");
    }
}


void print_bayes_credibility(FILE *file, MYREAL *cred, long numpop)
{
    long pa;
    long numpop2 = numpop * numpop;
    for(pa=0; pa < numpop; pa++)
        fprintf(file,"#Theta lower=%f upper=%f\n",cred[pa], cred[pa+ numpop2]);
    for(pa=numpop; pa < numpop2; pa++)
        fprintf(file,"#Migration lower=%f upper=%f\n",cred[pa], cred[pa+ numpop2]);
}

long number_replicates(world_fmt *world)
{
    long repmax = 1;
    if (world->options->replicate)
    {
        if (world->options->replicatenum == 0)
            repmax = world->options->lchains;
        else
            repmax = world->options->replicatenum;
    }
    else
        repmax = 1;
    return repmax;
}
long number_replicates2(option_fmt *options)
{
    long repmax = 1;
    if (options->replicate)
    {
        if (options->replicatenum == 0)
            repmax = options->lchains;
        else
            repmax = options->replicatenum;
    }
    else
        repmax = 1;
    return repmax;
}

//
// print header for complete (raw) parameter outfile (options->bayesmdimfilename)
// Header uses shell script escapes
// # Migrate MIGRATEVERSION (Peter Beerli, (c) 2008)
// # Raw results from Bayesian inference: these values can be used to generate
// # joint posterior distribution of any parameter combination
// # and also used for TRACER results [splitting utility may be needed]
#ifdef ZNZ
void print_bayes_mdimfileheader(znzFile file, long interval, world_fmt* world, data_fmt *data)
#else
void print_bayes_mdimfileheader(FILE *file, long interval, world_fmt* world, data_fmt *data)
#endif
{
#ifdef MPI
#ifdef PARALIO
    MPI_Status status;
#endif
#endif
    long i;
    long pa;
    long pa0;
    long pop;
    long frompop;
    long topop;
    long repmax;
    long numpop = world->numpop;
    long numpop2 = world->numpop2;
    long specstart = numpop2 + world->bayes->mu;
    long numparam = specstart + 2 * world->species_model_size;
    boolean *visited;
    bayes_fmt *bayes = world->bayes;
    long bufsize = 0;
    long allocbufsize = LONGLINESIZE;
    char *buf = (char *) mycalloc(allocbufsize,sizeof(char)); //this should be plenty for the header (currently 10^4 characters)
    bufsize = sprintf(buf,"# Migrate %s (Peter Beerli, (c) 2014)\n", MIGRATEVERSION);
    bufsize += sprintf(buf+bufsize,"# Raw results from Bayesian inference: these values can be used to generate\n");
    bufsize += sprintf(buf+bufsize,"# joint posterior distribution of any parameter combination\n");
    bufsize += sprintf(buf+bufsize,"# Writing information on parameters (Thetas, M or xNm)\n");
    bufsize += sprintf(buf+bufsize,"# every %li parameter-steps\n", interval+1);
    bufsize += sprintf(buf+bufsize,"# \n");
    bufsize += sprintf(buf+bufsize,"# --  %s\n", "Steps");
    bufsize += sprintf(buf+bufsize,"# --  %s\n", "Locus");
    bufsize += sprintf(buf+bufsize,"# --  %s\n", "Replicates");
    bufsize += sprintf(buf+bufsize,"# --  %s\n", "log(Posterior)");
    bufsize += sprintf(buf+bufsize,"# --  %s\n", "log(prob(D|G))");
    bufsize += sprintf(buf+bufsize,"# --  %s\n", "log(prob(G|Model))");
    bufsize += sprintf(buf+bufsize,"# --  %s\n", "log(prob(Model))");
    bufsize += sprintf(buf+bufsize,"# --  %s\n", "Sum of time intervals on G");
    bufsize += sprintf(buf+bufsize,"# --  %s\n", "Total tree length of G");
    print_param_order(&buf, &bufsize, &allocbufsize, world, world->numpop2);
    
    //if(world->bayes->mu)
    //    bufsize += sprintf(buf+bufsize,"# %li  %s\n",  world->numpop2+1, "Rate");
    bufsize += sprintf(buf+bufsize,"# \n");
    print_marginal_order(buf, &bufsize, world);
    repmax = number_replicates(world);
    bufsize += sprintf(buf+bufsize,"# \n");
    bufsize += sprintf(buf+bufsize,"#$ ------------------------------------------------------------------ \n");
    bufsize += sprintf(buf+bufsize,"#$ begin  [do not change this content]\n");
    bufsize += sprintf(buf+bufsize,"#$ Model=%s\n",world->options->custm);
    bufsize += sprintf(buf+bufsize,"#$ Mode2=%s\n",world->options->custm2);
    bufsize += sprintf(buf+bufsize,"#$ %li %li %li %li %li %i\n",world->loci, world->numpop,
                       world->numpop2, (long) world->options->replicate, repmax, (int) world->options->usem);
    for (pop = 0; pop < numpop; pop++)
    {
        bufsize += sprintf(buf+bufsize,"#$ %s\n",data->popnames[pop]);
    }
    bufsize += sprintf(buf+bufsize,"#$ end\n");
    bufsize += sprintf(buf+bufsize,"#$ ------------------------------------------------------------------ \n");
    bufsize += sprintf(buf+bufsize,"# \n");
    bufsize += sprintf(buf+bufsize,"# remove the lines above and including @@@@@, this allows to use\n");
    bufsize += sprintf(buf+bufsize,"# Tracer (http://tree.bio.ed.ac.uk/software/tracer/) to inspect\n");
    bufsize += sprintf(buf+bufsize,"# this file. But be aware that the Tracer program (October 2006)\n");
    bufsize += sprintf(buf+bufsize,"# only works with single-locus, single-replicate files\n");
    bufsize += sprintf(buf+bufsize,"# The migrate contribution folder contains a command line utility written\n");
    bufsize += sprintf(buf+bufsize,"# in PERL to split the file for Tracer, it's name is mba\n");
    bufsize += sprintf(buf+bufsize,"# @@@@@@@@\n");
    bufsize += sprintf(buf+bufsize,"Steps\tLocus\tReplicate\tlnPost\tlnDataL\tlnPrbGParam\tlnPrior\ttreeintervals\ttreelength");
    
    visited = (boolean*) mycalloc(numpop2,sizeof(boolean));
    for(pa0=0;pa0<numpop2;pa0++)
    {
        if(shortcut(pa0,world,&pa))
        {
            continue;
        }
        else
        {
            if(visited[pa]==TRUE)
                continue;
            else
                visited[pa]=TRUE;
        }
        if(pa < numpop)
        {
            bufsize += sprintf(buf+bufsize,"\t%s_%li", "Theta",pa+1);
        }
        else
        {
            m2mm(pa,numpop,&frompop,&topop);
            bufsize += sprintf(buf+bufsize,"\t%s_%li_%li",
                               (world->options->usem ? "M" : "xNm"), frompop+1, topop+1);
        }
    }
    if(bayes->mu)
    {
        bufsize += sprintf(buf+bufsize,"\t%s", "murate");
    }
    for(pa0=specstart; pa0 < numparam; pa0++)
    {
        if(shortcut(pa0,world,&pa))
        {
            continue;
        }
	species_fmt * s = get_which_species_model(pa,world->species_model,world->species_model_size);
        frompop = s->from;
        topop = s->to;
	if (pa == s->paramindex_mu)
	  bufsize += sprintf(buf+bufsize,"\t%s_%li_%li","D", frompop+1, topop+1);
	else
	  bufsize += sprintf(buf+bufsize,"\t%s_%li_%li","S", frompop+1, topop+1);
    }
    if (world->has_growth)
      {
	for(i=0;i<world->grownum;i++)
	  {
	    bufsize += sprintf(buf+bufsize,"\t%s_%li","Growth", i);
	  }
    }
    for(i=0;i<world->options->heated_chains;i++)
        bufsize += sprintf(buf+bufsize,"\tmL_%f", world->options->heat[i]);
    if(world->options->heated_chains>1)
        bufsize += sprintf(buf+bufsize,"\tmL_thermo");
    bufsize += sprintf(buf+bufsize,"\tmL_harmonic");
    bufsize += sprintf(buf+bufsize,"\n");
    
#ifdef MPI
#ifdef PARALIO
    //wrapper function [see migrate_mpi.c]
    //fprintf(stdout,"%i>%li %li\n@@%s@@\n\n",myID,bufsize, strlen(buf),buf);
    mpi_mdim_send(&world->mpi_bayesmdimfile,buf,bufsize);
#else
#ifdef ZNZ
    znzprintf(file,"%s",buf);
#else
    fprintf(file,"%s",buf);
#endif
#endif
#else
#ifdef ZNZ
    znzprintf(file,"%s",buf);
#else
    fprintf(file,"%s",buf);
#endif
#endif
    
    myfree(visited);
    myfree(buf);
}

//
// finds the mode of the results histogram
// fills the bayes->histogram[locus].modes
void     find_bayes_mode(world_fmt *world, long locus, long numparam)
{
  bayes_fmt *bayes = world->bayes;
    long bin;
    long j;
    long jj;
    long numbin=0;
    long pa0, pa;
    long *bins = bayes->histogram[locus].bins;
    
    MYREAL tmp;
    MYREAL *modes = bayes->histogram[locus].modes;
    MYREAL *mini = bayes->histogram[locus].minima;
    double *results = bayes->histogram[locus].results;
    
    for(pa0=0; pa0 < numparam; pa0++)
    {
        if(shortcut(pa0,world,&pa))
        {
            continue;
        }
	numbin = 0;// DEBUG/CHECK  THIS
        for(j=0; j < pa; j++)
        {
            if(shortcut(j,world,&jj))
            {
                continue;
            }
            numbin += bins[jj];
        }
        tmp = 0.;
        for(bin=0;bin<bins[pa]; bin++)
        {
            if (tmp < (double) results[numbin + bin])
            {
                tmp = (double) results[numbin+bin];
                modes[pa] = (MYREAL) bin * bayes->deltahist[pa] + mini[pa];
            }
        }
        //      numbin += bins[pa];
    }
}


//
// moving average smoothing
// a window of size 2*el + 1 is averaged and replaces the central value
// the array with the values is prepended with el zeros and also appended with el zeros,
// and the output replaces the input array,
// this methods uses an Epanechnikov kernel [k(x) = 3/4 (1-x^2) for -1<=x<=1
// my method assumes a boundary with boolean 'boundary' on the left,
// if this is set to TRUE then the data will be duplicated, reversed and attached to the
// left before run, return will be only the second half
void bayes_smooth(double *x, long xelem, long el, boolean lastfirst, boolean boundary)
{
    double d;
    double dd;
    double xsum;
    double *xn;
    double *xx;
    double *weight;
    double total;
    double totalx=0.;
    double totalxsum = 0.;
    long nxelem;
    long i, j, jj, w;
    long el2 = 2 * el + 1;
    // boundary code
    if (boundary)
      {
	nxelem = 2 * xelem;
	xn =  (double *) mycalloc(nxelem, sizeof(double));
	memcpy(xn+xelem, x, sizeof(double) * (size_t) xelem);
	for(i=1;i<xelem;i++)
	  {
	    xn[xelem-i] = x[i];
	  }
      }
    else
      {
	nxelem = xelem;
	xn = x;
      }
    // end boundary preparation
    weight = (double *) mycalloc(el2, sizeof(double));
    if(xelem==0)
        return;
    xx = (double *) mycalloc((1 + el2 + nxelem),sizeof(double));
    
    memcpy(xx+el, xn, sizeof(double) * (size_t) nxelem);
    if(lastfirst)
    {
        for(i=0; i < el; i++)
        {
            xx[i] = xx[el];
            xx[i+nxelem] = xx[nxelem-1];
        }
    }
    weight[el] = 0.75;
    total = 1.;
    d =  (2.0/(el2-1.));
    for(j = el+1; j < el2; j++)
    {
        dd = (-1.0 + d*j);
        dd *= dd;
        weight[j] =  0.75 * (1.0 - dd);
        total += 2 * weight[j];
    }
    //    printf("weight=%f\n",total);
    weight[el] /= total;
    
    for(j=el+1, jj = el-1; j < el2; j++, jj--)
    {
        weight[j] /= total;
        weight[jj] = weight[j];
    }
    
    for(i=el; i < nxelem + el; i++)
    {
        xsum = 0.;
        for(j=i - el, w=0; j < i + el + 1; j++, w++)
            xsum += xx[j] * weight[w];
        totalxsum += xsum;
        totalx += xn[i-el];
        //fprintf(stdout,"%20.20f %20.20f\n",x[i-el], xsum);
        xn[i-el] = xsum;
        //nx[i-el] = xsum/el2;
    }
    //    fprintf(stdout,"%20.20f %20.20f\n",totalx, totalxsum);
    if(totalxsum>0.0)
        totalx /= totalxsum;
    else
    {
        myfree(xx);
        myfree(weight);
        return;
    }
    for(i=0; i < nxelem; i++)
    {
        xn[i] *= totalx;
        //fprintf(stdout,"%20.20f\n",x[i]);
    }
    if(boundary)
      {
	memcpy(x, xn+xelem, sizeof(double) * (size_t) xelem);
	myfree(xn);
      }
    myfree(xx);
    myfree(weight);
}

//
// find the credibility set by using the Highest Posterior Probability Density(HPD) and the standard point
// descriptors. the statistics are filled into the statistics part of the bayes structure (datastore)
void calc_hpd_credibility(world_fmt *world,long locus, long numpop2, long numparam)
{
  (void) numpop2;
    const MYREAL alpha95 = 0.95;
    const MYREAL alpha50 = 0.50;
    bayes_fmt *bayes = world->bayes;
    //long j;
    long li;
    long pa0, pa;
    long rpa; // the many locus-rates are using only a single rate prior and recording (per locus)
    // and we need to distinguish here
    
    long numbins=0;
    long numbinsall = 0;
    long locmedian;
    
    pair *parts; // pairs of doubles
    long csum = 0;
    boolean *smoothy=NULL;
    
    MYREAL biggestcolumn;
    MYREAL total;
    MYREAL total1;
    MYREAL total2;
    MYREAL cdf;
    MYREAL cutoff95;
    MYREAL cutoff50;
    MYREAL delta;
    MYREAL tmp;
    
    MYREAL * mini;
    MYREAL * maxi;
    double *results;
    double *results2;
    long *bins;
    MYREAL *modes;
    MYREAL *medians;
    MYREAL *cred50;
    MYREAL *cred95;
    char *set50;
    char *set95;
    //double xvals[3];
    bins = bayes->histogram[locus].bins;
    modes = bayes->histogram[locus].modes;
    medians = bayes->histogram[locus].medians;
    mini =  bayes->histogram[locus].minima;
    maxi =  bayes->histogram[locus].maxima;
    cred50 =  bayes->histogram[locus].cred50l;
    cred95 =  bayes->histogram[locus].cred95l;
    set50 =  bayes->histogram[locus].set50;
    set95 =  bayes->histogram[locus].set95;
    results =  bayes->histogram[locus].results;
    if (bayes->histogram[locus].results2 == NULL)
    {
      bayes->histogram[locus].results2 = (double *) mymalloc(sizeof(double) * (size_t) bayes->histogram[locus].binsum);
      memcpy(bayes->histogram[locus].results2,bayes->histogram[locus].results,sizeof(double)* (size_t) bayes->histogram[locus].binsum);
    }
    else
    {
      myfree(bayes->histogram[locus].results2);
      bayes->histogram[locus].results2 = (double *) mymalloc(sizeof(double) * (size_t) bayes->histogram[locus].binsum);
      memcpy(bayes->histogram[locus].results2,bayes->histogram[locus].results,sizeof(double)* (size_t) bayes->histogram[locus].binsum);
    }
    results2 =  bayes->histogram[locus].results2;
    //smoothy = (boolean *) mycalloc(numparam, sizeof(boolean));
    if(bayes->histogram[locus].smoothed == NULL)
        bayes->histogram[locus].smoothed = (boolean *) mycalloc(numparam, sizeof(boolean));
    smoothy = bayes->histogram[locus].smoothed;
    for(pa0=0; pa0 < numparam; pa0++)
    {
        if(shortcut(pa0,world,&pa))
        {
            continue;
        }
        
        // so that we can check whether the bins got smoothed already or not
        smoothy[pa0] = FALSE;
        
        rpa = pa;
        
        if(csum<bins[rpa])
            csum = bins[rpa];
    }
    parts = (pair *) mycalloc(csum, sizeof(pair));
    
    for(pa0=0; pa0 < numparam; pa0++)
    {
        if(shortcut(pa0,world,&pa)) //1229|| (pa0<numpop2 && bayes->custm2[pa0]=='d'))
        {
            continue;
        }
        if(pa<pa0)
            continue;
        
        rpa = pa;
        
        numbinsall += bins[rpa];
        numbins = numbinsall - bins[rpa];
        
        delta = bayes->deltahist[rpa];//(maxi[rpa] - mini[rpa]) / bins[rpa];
        //printf("%i> pa=%li delta = %f",myID, pa, delta);
        parts = (pair *) memset(parts,0,(size_t) csum * sizeof(pair));
        
        //
        // calculate the total and then the average
        //total = 0.;
        biggestcolumn = 0.;
        locmedian = 0;
        //
        // smooth the results before we calculate means etc
        if(!smoothy[rpa])
        {
            //	  bayes_smooth(results+numbins,bins[pa], bins[pa]/50);
            //bayes_smooth(results+numbins,bins[rpa], MAX(5.0,bins[rpa]/50.0),FALSE);
	    //if(locus!=world->loci)
	    //kernel_smooth(results+numbins, bins[rpa], results+numbins, bins[rpa], (long) (sqrt(world->options->lsteps)), mini[rpa], maxi[rpa]);
	    bayes_smooth(results+numbins,bins[rpa], MAX(5.0,(long) (sqrt((double) bins[rpa]))),FALSE,FALSE);
#ifdef DEBUG
	  printf("%i> after bayes_smooth: bins[%li]=%li (%li)\n",myID, rpa,bins[rpa], (long) (sqrt(world->options->lsteps)));
#endif
            smoothy[rpa] = TRUE;
        }
        total1 = 0.0;
	total2 = 0.0;
        for(li=0;li<bins[rpa]; li++)
        {
            tmp = (double) results[numbins + li];
            total1 += tmp;//* delta;
            total2 += (double) results2[numbins + li];//* delta;
        }
        total=0.0;
	double tmp2=0.0;
	double total22=0.0;
        for(li=0;li<bins[rpa]; li++)
        {
            //@@@@@@@@@@
            results[numbins + li] /= (double) total1;
            results2[numbins + li] /= (double) total2;
            tmp = (double) results[numbins + li];
            tmp2 = (double) results2[numbins + li];
            parts[li][1] = tmp;
            total += tmp;
	    total22 += tmp2;
            parts[li][0] = delta/2.0 + mini[rpa] + li * delta;
            if(tmp > biggestcolumn && li>1)
            {
                biggestcolumn = tmp;
                locmedian = li;
            }
        }
#ifdef DEBUG 
	printf("%i> locus=%li pa=%li sum(result,result2)=%f,%f\n",myID, locus, rpa, total, total22);
#endif
        //
        // mode is the value with largest column
        modes[pa] = parts[locmedian][0];
        //
        // median is at 0.5 of the cumulative probability distribution
        li = 0;
        tmp = 0;
        while(tmp < total/2  && li < bins[rpa]-1)
        {
            tmp += parts[li][1];
            li++;
        }
        medians[pa] = parts[li][0];
        //
        // sort parts using the histogram y-axis for easy finding of quantiles
        paired_qsort2(parts, bins[rpa]);
        
        // find HPD crediblity intervals for 50% and 95% credibility intervals
        // for 50% intervals, starting from the highest value (mode) and moving down
        // until the cumulative sum / total is >=0.5
        cdf = 0.;
        li = bins[rpa];
        while(cdf < alpha50 && li>0)
        {
            li--;
            cdf += parts[li][1] /total;
        }
        cutoff50 = parts[li][1]; // anything lower than this will be in the 50% credibility set
#ifdef DEBUG
	printf("pa=%li c50=%f\n",rpa,cutoff50);
#endif
        // or 95% interval
        while(cdf < alpha95 && li>0)
        {
            li--;
            cdf += parts[li][1] /total;
        }
        cutoff95 = parts[li][1]; // anything lower than this will be in the 95% credibility set
#ifdef DEBUG
	printf("pa=%li c95=%f\n",rpa,cutoff95);
#endif
        
        // fill the credibility sets (are printed into Bayesfile
        for(li=numbins;li<numbins + bins[rpa]; li++)
        {
            if(results[li] < cutoff95)
            {
                set50[li] = '0';
                set95[li] = '0';
            }
            else
            {
                set95[li] = '1';
                if(results[li] < cutoff50)
                    set50[li] = '0';
                else
                    set50[li] = '1';
            }
        }
        // fill the innermost 50% levels, smooth over adjacent bins
        //
        // start at the modus and go left
        li = numbins + locmedian;
        while (set50[li] == '1' && li > numbins)
            --li;
        cred50[pa] = mini[rpa] + (li-numbins) * delta;// + delta/2.;
        while (set95[li] == '1' && li > numbins)
            --li;
        cred95[pa] = mini[rpa] + (li-numbins) * delta;// + delta/2.;

        // start at the modus and go right
        li = numbins + locmedian;
        while (set50[li] == '1' && li < numbins + bins[rpa]-1)
            ++li;
        cred50[pa + numparam] = maxi[rpa] - (bins[rpa] - li + numbins) * delta;// - delta/2.;
        while (set95[li] == '1' && li < numbins + bins[rpa]-1)
            ++li;
        cred95[pa + numparam] = maxi[rpa] - (bins[rpa] - li + numbins) * delta;// - delta/2.;
#ifdef DEBUG
	printf("pa=%li cred50=%f..%f [%li]\n",rpa,cred50[pa],cred50[pa+numparam],numparam);
	printf("pa=%li cred95=%f..%f [%li]\n",rpa,cred95[pa],cred95[pa+numparam],numparam);
#endif

    }
    myfree(parts);
    myfree(smoothy);
}
//
// combines over loci
void bayes_combine_loci(world_fmt * world)
{
  const long loci    = world->loci;
  const long numpop2 = world->numpop2;
  //const long np      = numpop2 + world->bayes->mu + 2 * world->species_model_size;
  const long np2     = numpop2 + world->bayes->mu + 2 * world->species_model_size + world->grownum;
  long    targetstartbin;
  long    locus;
  long    pa0;
  long    pa=0;
  long    sourcebin;
  long sourcesumbin;
  long targetsumbin;
  long    bin;
  long    *bins;
  MYREAL  v;
  MYREAL  count;
  MYREAL  sumprob;
  MYREAL  *maxvala;
  MYREAL  *maxvalaloc;
  MYREAL  total;
  MYREAL  *priors;
  bayes_fmt * bayes = world->bayes;
  // records whether a cell in the all-loci histogram was filled or not
  char    *touched;
  
  // source is the locus histogram, already filled in by the locus workers or then by the locus loop
  bayeshistogram_fmt * source = NULL;
  // target is the summary over all loci
  bayeshistogram_fmt *target = &bayes->histogram[world->loci];
  //@Dec16/15 float    *results;
  double    *results;
  MYREAL   *minima;
  boolean  *visited;
  MYREAL   *mini;
  MYREAL   *maxi;
  MYREAL   *adjmaxi;
  MYREAL   w;// width of bin
  MYREAL   logw; //log width of bins
  long tts;
  //MYREAL   countsum;
  long     *counts;
  MYREAL   *locusweight = world->data->locusweight;
  MYREAL   allpriortotal;
  // allocation of space for maximal values for parameters
  doublevec1d(&maxvala, np2);
  doublevec1d(&maxvalaloc, np2);
  visited = (boolean *) mycalloc(np2,sizeof(boolean));
  // Not every locus is filling all bins in the posterior
  // counts keeps track of the number of loci that are used per bin
  counts = (long *) mycalloc(bayes->histogram[world->loci].binsum,sizeof(long));
  //
  // set the target ('sum' over all loci) minima and maxima
  for(pa0=0;pa0<np2; pa0++)
    {
      bayes->histogram[loci].minima[pa0] = MYREAL_MAX;
      bayes->histogram[loci].maxima[pa0] = -MYREAL_MAX;
    }
  //
  // calculate the number of used loci (under some option invariants are skipped)
  // and adjust the min/max for the working loci
  for(locus=0;locus<world->loci; locus++)
    {
      if(world->data->skiploci[locus])
        {
	  continue;
        }
      mini =  world->bayes->histogram[locus].minima;
      maxi =  world->bayes->histogram[locus].maxima;
      adjmaxi =  world->bayes->histogram[locus].adjmaxima;
      
      find_bayes_min_max(world, &mini, &maxi, &adjmaxi);
      adjust_bayes_min_max(world,&mini,
			   &maxi,
			   &bayes->histogram[loci].minima,
			   &bayes->histogram[loci].maxima);
    }
  // allocation of all-loci-histogram table into results, and for the 50% and 95% credibility sets
  adjust_bayes_bins(world, world->loci);
  bins = target->bins;
  results = (double *) mycalloc((size_t) (1+bayes->histogram[loci].binsum), sizeof(double));
  minima = target->minima;
  // records whether a cell in the all-loci histogram was filled or not
  touched = (char *) mycalloc((size_t) (1 + bayes->histogram[loci].binsum), sizeof(char));
  // record prior
  priors = (MYREAL *) mycalloc((size_t) (1+bayes->histogram[loci].binsum), sizeof(MYREAL));
  memset(touched,' ', (size_t) bayes->histogram[world->loci].binsum * sizeof(char));
  // combine previous results into target->results
  // over all loci
  // ln mL = ln K * Sum ln mL_i
  // ln K = Prod_param Sum_bins Exp[ln(dx) + ln p(param,bin)*(1-loci) + Sum_loci ln posterior(param,bin, locus)]
  //             1.         2.   3.    7.        4.                          5.             6.
  // Kahan sum error compensation
  // [wikipedia 01-05-2016]
  // sum is results[i]
  // temp variables: y,t
  // compensation for lost low-order bits
  //@ double * kahan_c = (double *) mycalloc(bayes->histogram[world->loci].binsum, sizeof(double));
  //FILE * gugu = fopen("gugus","w");
  for(locus=0;locus<loci; locus++)
    {
      if(world->data->skiploci[locus])
	continue;
      sourcesumbin = 0;
      targetsumbin = 0; //records the already worked on samples
      source = &bayes->histogram[locus];
      memset(visited,0,sizeof(boolean)*(size_t) np2);
      // over all parameters
      for(pa0=0;pa0<np2;pa0++)
	{
	  if(shortcut(pa0,world,&pa))
	    {
	    continue;
	    }
	  if(visited[pa]==TRUE)
	    continue;
	  visited[pa] = TRUE;
	  //w = world->bayes->deltahist[pa];
	  maxvala[pa] = (double) -HUGE;
	  // set the start bin for the target
	  // this is safeguarding against different min max per locus
	  targetstartbin = 0;//(long) ((source->minima[pa] - minima[pa]) / w);
	  if(source->results2 != NULL)
            {
	      //      countsum=0.0;
	      tts = targetsumbin + targetstartbin;
	      //const long otts = tts;
	      const long sourcebins = source->bins[pa] + sourcesumbin;
#ifdef DEBUG
	      printf("@#$ %li $$ %li %li $$ ",locus, source->bins[pa], tts);
#endif
	      for(sourcebin=sourcesumbin; sourcebin < sourcebins; tts++, sourcebin++)
                {
		  //if (locus>20)
		  //  continue;
		  if(source->results2[sourcebin] > 0.0)
                    {
		      // frequencies are normalized per locus
		      //count = (double) source->results2[sourcebin];
		      count = source->results2[sourcebin];
		      //v = tts * w + w/2;
		      //priors
		      //double pbin = scaling_prior(world,pa,v); //calculate log prior for parameter
		      //if (log(count) < -9 || locus > 20)
		      //continue;
		      double logcount = log(count); // - kahan_c[tts];
		      touched[tts] = '1';
		      //results[tts] += /*(float)*/ log(count);//- pbin ; //5.
		      double t = results[tts] + logcount;
		      //kahan_c[tts] = (t - results[tts]) - logcount;
		      results[tts] = t;
		      
		      if (results[tts] > maxvala[pa])
			{
			  maxvala[pa] = results[tts];
			  maxvalaloc[pa] = sourcebin;
			}
		      counts[tts] += (long) locusweight[locus];//1.
		      //	      countsum += count;
                    } /**if sourcebin is not zero**/
                }/**over all sourcebins**/
                sourcesumbin += source->bins[pa];
                targetsumbin += target->bins[pa];
            }/**source!=null**/
        }/**parameter**/
    }/**loci**/
    // then
    // adjust the total and use the overflow safeguard
    targetsumbin = 0;
  memset(visited,0,sizeof(boolean)* (size_t) np2);
  for(pa0 = 0; pa0 < np2; pa0++)
    {
      if(shortcut(pa0,world,&pa))
	{
	  continue;
	}
        if(visited[pa]==TRUE)
	  continue;
        visited[pa] = TRUE;
        w = world->bayes->deltahist[pa];
      logw = log(w); //7. ln(dx)
      for(bin=targetsumbin; bin < targetsumbin + bins[pa]; bin++)
	{
	  v = source->minima[pa] + (bin-targetsumbin) * w + w/2;
	  //priors
	  priors[bin] = logw + scaling_prior(world,pa,v); //calculate log prior for parameter
            if(touched[bin] == '1')
	      {
                // posteriors
	      if(counts[bin]==0)
		  {
		    results[bin] = (double) -HUGE;
		  }
	      if(results[bin] > maxvala[pa])
		maxvala[pa] = (MYREAL) results[bin];
	    }
	}
      allpriortotal=0.0;
	maxvala[pa]= (double) -HUGE;
	for(bin=targetsumbin; bin < targetsumbin + bins[pa]; bin++)
	  {
	    if(touched[bin]=='1')
	    {
	      if (priors[bin] > -11.0)  //this arbitrary cutoff seems to work, but
		// does not correct some results[bin], new experiment with setting
		{
		  //if (counts[bin] != loci)
		  ///if (counts[bin] < loci/10)
		  ///  {
		  ///    touched[bin] = ' ';
		  ///  }
		  ///else
	      {
		    results[bin] += priors[bin] * (1.0 - counts[bin]);
		if(results[bin] > maxvala[pa])
		  maxvala[pa] = (MYREAL) results[bin];
		  }
	      }
	    else
	      {
		  touched[bin] = ' ';
		}
	    }
	  priors[bin] = exp(priors[bin]); //@Dec9/15//non-log the priors for print/plotting later
	  allpriortotal += priors[bin]; //@Dec9/15  checksum must total to 1.0
	}
	target->means[pa] = 0.;
	sumprob = 0.;
	for(bin=targetsumbin; bin < targetsumbin + bins[pa]; bin++)
	  {
	    if(touched[bin]=='1')
	      {
		double x = exp(results[bin] - maxvala[pa]);
	      results[bin] = /*(float)*/ x;
		sumprob += x;
	      }
	    else
	      {
		results[bin] = 0.;
	      }
	  }
	if(sumprob > 0.0)
	  {
	    //log the integral(SUMPROB) and adjust with the biggest value
	  bayes->scaling_factors[pa] = logw + log(sumprob) + maxvala[pa];//logw was commented out
	    
#ifdef DEBUG
	    printf("%i> scalingfactor[%li]=%f\n",myID, pa, bayes->scaling_factors[pa]);
#endif
	  }
	else
	  bayes->scaling_factors[pa] = (double) -HUGE;
	total = 0.0;
	for(bin=targetsumbin; bin < targetsumbin + target->bins[pa]; bin++)
	  {
	    if(touched[bin]=='1')
	      {
	      results[bin] /= /*(float)*/ (sumprob);//this takes care of the maxvala scaling
		total +=  results[bin];
	      target->means[pa] += results[bin] * (minima[pa] +  w/2. + w * (bin-targetsumbin));
	      }
	  }
	target->means[pa] /= total;
#ifdef DEBUG
	fprintf(stdout,"allpriortotal=%f total=%f sumprob=%f w=%f\n== end summary \n",allpriortotal,total, sumprob, w);
#endif
	targetsumbin += bins[pa];
    }
  long i;
  for (i=0; i<bayes->histogram[loci].binsum;i++)
    target->results[i] = results[i];
    //covariance_summary(world);
    /// calculate the credibility intervals, histograms, means etc
  calc_hpd_credibility(world, loci, world->numpop2, np2);
    bayes->priors = priors; //save the priors
    myfree(touched);
  //myfree(kahan_c);
    myfree(maxvala);
    myfree(visited);
    myfree(counts);
}
////////////////////////////////////////////////////////////////////////////////

// calculates the HPC credibility set and also calculates the posterior singlelocus-histogram
void calculate_credibility_interval(world_fmt * world, long locus)
{
    MYREAL *mini;
    MYREAL *maxi;
    MYREAL *adjmaxi;
    //    long i;
    long np = world->numpop2 + world->bayes->mu + 2 * world->species_model_size + world->grownum;
    
    mini =  world->bayes->histogram[locus].minima;
    maxi =  world->bayes->histogram[locus].maxima;
    adjmaxi =  world->bayes->histogram[locus].adjmaxima;
    
    find_bayes_min_max(world, &mini, &maxi, &adjmaxi);
    adjust_bayes_bins(world, locus);
    
    // construct histogram
    construct_locus_histogram(world, locus, mini, maxi, &world->bayes->histogram[locus].results);
    
    // calc_credibility
    calc_hpd_credibility(world, locus, world->numpop2, np);
}


/// 
/// sets parameters so that they represent an average, this allows to have multiple groups
/// with different 'averages'
void handle_averages(char representant, MYREAL newparam, long which, long numpop, char *custm2,  MYREAL * param)
{
  long i;
  if(which<numpop)
    {
      for(i=0;i<numpop; ++i)
	{
	  if(custm2[i]==representant)
	    param[i] = newparam;
	}
    }
  else
    {
      for(i=numpop;i<numpop*numpop; ++i)
	{
	  if(custm2[i]==representant)
	    param[i] = newparam;
	}
    }
}
//
// adjust the parameters so that the new set is consistent with the custom migration matrix
void
bayes_set_param (MYREAL *param, MYREAL newparam, long which, char *custm2, long numpop)
{
    char    migtype;
    
    long    frompop = 0;
    long    topop   = 0;
    long    i;
    
    MYREAL  nmig;
    
    if (which >= numpop*numpop)
    {
        param[which] = newparam;
        return;
    }
    migtype = custm2[which];
    //  check custm matrix and then decide
    switch(migtype)
      {
      case 'C':
      case 'c':
	break;
      case 's':
	m2mm (which, numpop, &frompop, &topop);
	param[mm2m(frompop,topop,numpop)] = newparam; // makes the
	param[mm2m(topop,frompop,numpop)] = newparam; // two parameter equal
	break;
      case 'S':
	m2mm (which, numpop, &frompop, &topop);
	param[mm2m(frompop,topop,numpop)] = newparam; // makes the
	param[mm2m(topop,frompop,numpop)] = param[mm2m(frompop,topop,numpop)] * param[topop] / param[frompop]; // two parameter equal
	break;
      case 'm':
	handle_averages('m', newparam, which, numpop, custm2,  param);
	break;
      case 'M':
	if(which<numpop)
	  {
	    for(i=0;i<numpop; ++i)
	      {
		if(strchr("Mm",custm2[i]))
		  param[i] = newparam;
	      }
	  }
	else
	  {
	    m2mm (which, numpop, &frompop, &topop);
	    nmig = param[topop] * newparam;
	    for(i=numpop;i<numpop*numpop; ++i)
	      {
		if(custm2[i]=='M')
		  {
		    m2mm (i, numpop, &frompop, &topop);
		    param[i] = nmig / param[topop]; // makes the 
		  }
	      }
	  }
	break;
      case 'd':
      case 't':
	param[which] = 0.0;
	warning("custm2 should be set to migration things already!");
	break;
      case 'x':
      case '*':
	param[which] = newparam;
	break;
      default:
	// case allows to specify different groups beyond mM	
	handle_averages(custm2[which], newparam, which, numpop, custm2,  param);	
	break;
    }
}
