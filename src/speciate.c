 /*------------------------------------------------------
 inference of population parameters
 using a Metropolis-Hastings Monte Carlo algorithm
 -------------------------------------------------------
 speciation (fusion [and later fission] routines)
 
 Peter Beerli 2013, Tallahassee
 beerli@fsu.edu
 
 Copyright 2013-2015 Peter Beerli
 
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
 
$Id:$
*/
/* \file speciate.c 
speciation tools
*/
#include "migration.h"
#include "bayes.h"
#include "random.h"
#include "tools.h"
#include "tree.h"
#include "random.h"
#include "uep.h"
#include "random.h"
#include "sighandler.h"
#include "migrate_mpi.h"
#include "assignment.h"
#include "mcmc.h"
#include "mcmc2.h"
#include "priors.h"
#include "pretty.h"
#include "speciate.h"
#include "mittag_leffler.h"
#include <assert.h>

extern profuncptr *propose_new;
log_prob_wait_speciate_func log_prob_wait_speciate;
log_point_prob_speciate_func log_point_prob_speciate;
time_to_speciate_func time_to_speciate;



MYREAL eventtime_single(proposal_fmt *proposal, world_fmt *world, long pop, long timeslice, long *lineages, double age, char * event, long *to, long *from);
long init_speciesvector(world_fmt * world, option_fmt * options);
int beyond_last_node(proposal_fmt* proposal, vtlist *tentry, long gte, long *slider);
void fill_speciesvector(world_fmt*world, option_fmt * options);
void species_datarecorder(world_fmt *world);
#if defined(MPI) && !defined(PARALIO) /* */
void print_species_record(float *temp, long *z, world_fmt *world);
#else
void  print_species_record(char *temp, long *c, world_fmt * world);
#endif
void read_species_record(long masternumbinall, long locus, MYREAL *params, long *n, MYREAL *oldmeans, MYREAL *lowerbound, MYREAL *upperbound, MYREAL *delta, world_fmt *world, char **inptr);
void check_speciateparameters(species_fmt *s);
species_fmt * get_which_species_model( long which, species_fmt * s ,  long ssize);
species_fmt * get_fixed_species_model( long from, long to, species_fmt * s ,  long ssize);
species_fmt * get_species_model( long to, species_fmt * s ,  long ssize);




MYREAL log_prob_wait_speciate_exp(MYREAL t0, MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt *s);
MYREAL log_point_prob_speciate_exp(MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt * s);
MYREAL time_to_speciate_exp(world_fmt *world, long pop, MYREAL t0, char *event, long *to, long *from);

MYREAL log_prob_wait_speciate_weibull(MYREAL t0, MYREAL t1, MYREAL lambda, MYREAL k, species_fmt *s);
MYREAL log_point_prob_speciate_weibull(MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt * s);
MYREAL time_to_speciate_weibull(world_fmt *world, long pop, MYREAL t0, char *event, long *to, long *from);

MYREAL log_prob_wait_speciate_normal(MYREAL t0, MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt *s);
MYREAL log_point_prob_speciate_normal(MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt * s);
MYREAL time_to_speciate_normal(world_fmt *world, long pop, MYREAL t0, char *event, long *to, long *from);

MYREAL time_to_speciate_normal_erf(world_fmt *world, long pop, MYREAL t0, char *event, long *to, long *from);

MYREAL log_prob_wait_speciate_normalorig(MYREAL t0, MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt *s);
MYREAL log_point_prob_speciate_normalorig(MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt * s);
MYREAL time_to_speciate_normalorig(world_fmt *world, long pop, MYREAL t0, char *event, long *to, long *from);

MYREAL log_prob_wait_speciate_normalshortcut(MYREAL t0, MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt *s);
MYREAL log_point_prob_speciate_normalshortcut(MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt * s);
MYREAL time_to_speciate_normalshortcut(world_fmt *world, long pop, MYREAL t0, char *event, long *to, long *from);

char   set_type(world_fmt *world, long topop, long frompop, char *custm2, long numpop);
node * set_type2(world_fmt *world, node *p, node *q, char *custm2);

MYREAL time_to_coalmig (world_fmt * world, long pop, long timeslice, long *lineages, char * event,  long * to, long * from);
MYREAL time_to_coalescence(world_fmt * world, long pop, double age, long timeslice, long *lineages, char * event, long * to, long *from);
MYREAL time_to_migration(proposal_fmt *proposal, world_fmt *world, long pop, long timeslice, long * lineages, char * event, long * to, long * from);
MYREAL time_to_speciation(world_fmt * world, long pop, double age, char * event, long * to, long * from);
void keep_min_eventtime(MYREAL *the_eventtime, MYREAL eventtime, char *the_event, char myevent,
			long *to, long tox, long *from, long fromx);
long speciation_from(long to, proposal_fmt * proposal);
void loopcleanup(boolean assign, world_fmt * world, proposal_fmt *proposal, long oldpop, long newpop, timelist_fmt *timevector);
long newtree_update (world_fmt * world, long g, boolean assign);
void set_things(MYREAL t1, MYREAL t2, char e1, char e2, long to1, long from1, long to2, long from2, 
		MYREAL *time, char *event, long *to, long *from);
void calc_time_per_line(proposal_fmt * proposal, char type, long popto, boolean same,MYREAL *time, char *event, long *to, long *from);
int beyond_last_node(proposal_fmt* proposal, vtlist *tentry, long gte, long *slider);
long get_species_record(world_fmt * world,long which);
long propose_new_spec_mu(world_fmt * world, long which, boolean *is_mu, MYREAL *newparam);
void adjust_averagediv(world_fmt * world, long which, MYREAL *newparam);
MYREAL wait_event_species(world_fmt *world, vtlist *tli, MYREAL t0, MYREAL t1, boolean waitonly, MYREAL *eventprob);
MYREAL wait_D(long pop, MYREAL t0, MYREAL t1, long *lineages, world_fmt *world);
void set_first_speciestree(node *mrca, world_fmt *world);
void construct_species_hist(species_fmt *db, world_fmt *world, long locus,
			    long npa, long pai, long offset, long numbin, 
			    MYREAL *mini, MYREAL *maxi, float **results,
			    float *total, float *themean, float *thestd);
void construct_locusspecies_histogram(world_fmt *world, long locus, MYREAL *mini, MYREAL *maxi, float **results);
void get_div_from_tree(node *root, divtime_fmt **div, long *divtime_alloc, long *divtime_num);
void record_parameters(world_fmt *world);
void set_speciate_functions(world_fmt *world);

double mylogerf(double x);

#ifndef SPECIATETESTER
//##

void set_speciate_functions(world_fmt *world)
{
  switch(world->species_model_dist)
    {
    case WEIBULL_DIST:
      log_prob_wait_speciate = (double (*) (double, double, double, double, species_fmt *)) log_prob_wait_speciate_weibull;
      log_point_prob_speciate = (double (*) (double, double, double, species_fmt *)) log_point_prob_speciate_weibull;
      time_to_speciate = (double (*) (world_fmt *, long, double, char *, long *, long *)) time_to_speciate_weibull;
      //printf("%i> divergence method: Weibull\n",myID);
      break;
    case 9:
      error("what is divergence distribution case 9?");
      //break;
    case NORMAL_DIST:
      log_prob_wait_speciate = (double (*) (double, double, double, double, species_fmt *)) log_prob_wait_speciate_normalorig;
      log_point_prob_speciate = (double (*) (double, double, double, species_fmt *)) log_point_prob_speciate_normalorig;
      time_to_speciate = (double (*) (world_fmt *, long, double, char *, long *, long *)) time_to_speciate_normalorig;
      //printf("%i> divergence method: Normal\n",myID);
      break;
    case NORMALSHORTCUT_DIST:
      log_prob_wait_speciate = (double (*) (double, double, double, double, species_fmt *)) log_prob_wait_speciate_normalshortcut;
      log_point_prob_speciate = (double (*) (double, double, double, species_fmt *)) log_point_prob_speciate_normalshortcut;
      time_to_speciate = (double (*) (world_fmt *, long, double, char *, long *, long *)) time_to_speciate_normalshortcut;
      //printf("%i> divergence method: Normal (Shortcut)\n",myID);
      break;
    case EXP_DIST:
    default:
      log_prob_wait_speciate = (double (*) (double, double, double, double, species_fmt *)) log_prob_wait_speciate_exp;
      log_point_prob_speciate = (double (*) (double, double, double, species_fmt *)) log_point_prob_speciate_exp;
      time_to_speciate = (double (*) (world_fmt *, long, double, char *, long *, long *)) time_to_speciate_exp;
      //printf("%i> divergence method: Exponential\n",myID);
      break;
    }

}
  


long init_speciesvector(world_fmt * world, option_fmt *options)
{
  const long numpop = world->numpop;
  const long numpop2 = world->numpop2;
  long i;
  long ssize = 0;
  char * custm2 = world->options->custm2;
  if (world->has_speciation)
    {
      for (i=numpop; i<numpop2; i++)
	{
	  if(strchr("dDtT",custm2[i]))
	    {
	      ssize++;
	    }
	}
#ifdef DEBUG
      printf("%i> speciation events: %li %s\n",myID, ssize, custm2);
#endif
      world->species_model = (species_fmt *) mycalloc((ssize+1),sizeof(species_fmt));
      world->species_model_size = ssize;
      world->species_model_dist = options->species_model_dist;
    }
  return ssize;
}

// init speciation vector
//
void fill_speciesvector(world_fmt*world, option_fmt * options)
{
  const long numpop = world->numpop;
  const long numpop2 = world->numpop2;
  long i;
  long from;
  long to;
  long ssize = 0;
  char * custm2 = world->options->custm2;
  species_fmt *s = world->species_model;
  bayeshistogram_fmt *hist;
  long z;
  long zz;
  long z3;
  if (world->has_speciation)
    {
      z = numpop2+world->bayes->mu;
      zz = numpop2+world->bayes->mu;
      z3 = numpop2+world->bayes->mu;
      for (i=numpop; i<numpop2; i++)
	{
	  if(strchr("dDtT", custm2[i]))
	    {
	      s[ssize].type = (char) tolower(custm2[i]);//this contains either 'd' or 't'
	      custm2[i] = (strchr("dt", custm2[i]) ? '0' : '*');
	      m2mm(i,numpop,&from,&to);
	      s[ssize].id = ssize;
	      s[ssize].from = from;
	      s[ssize].to   = to;
	      s[ssize].mu   = (double) options->bayes_priors[z3].mean;
	      s[ssize].min  = (double) options->bayes_priors[z3].min;
	      s[ssize].max  = (double) options->bayes_priors[z3].max;
	      s[ssize].paramindex_mu = z3++;
	      s[ssize].sigma= (double) options->bayes_priors[z3].mean;
	      s[ssize].sigmamin  = (double) options->bayes_priors[z3].min;
	      s[ssize].sigmamax  =  (double) options->bayes_priors[z3].max;
	      s[ssize].paramindex_sigma = z3++;
	      s[ssize].allocsize = HUNDRED;
	      s[ssize].data = (double*) mycalloc(s[ssize].allocsize,sizeof(double));
	      world->bayes->minparam[zz]=(double) s[ssize].min;
	      world->bayes->maxparam[zz++]=(double) s[ssize].max;
	      world->bayes->minparam[zz]=(double) s[ssize].sigmamin;
	      world->bayes->maxparam[zz++]=(double) s[ssize].sigmamax;
	      if (world->cold)
		{
		  hist = &(world->bayes->histogram[world->locus]);
		  hist->bins[z] = options->bayes_priors[z].bins; 
		  hist->minima[z] = (double) s[ssize].min;
		  hist->maxima[z++] = (double) s[ssize].max;
		  //sigma
		  hist->bins[z] = options->bayes_priors[z].bins;  
		  hist->minima[z] = (double) s[ssize].sigmamin;
		  hist->maxima[z++] = (double) s[ssize].sigmamax;	       
		}
	      ssize += 1;
	    }
	}
    }
}

void species_datarecorder(world_fmt *world)
{
  species_fmt * adb = world->species_model;
  species_fmt * db;
  long s;
  long i;
  // records the mu and sigma
  // for posterior distributions 
  for (i=0;i<world->species_model_size;i++)
    {
      db = &(adb[i]);
      s = db->size;
      if (s+3 >= db->allocsize)
	{
	  db->allocsize += HUNDRED;
	  db->data = (double *) myrealloc(db->data,sizeof(double)*(size_t) db->allocsize);
	}
      // this is always one single locus
      //
      db->data[s++] = db->mu;
      db->data[s++] = db->sigma;
      //printf("#EXPERIMENT: %f %f ",db->mu,db->sigma);
      //for(i=0;i<world->numpop2+ world->bayes->mu + 2* world->species_model_size;i++)
      //	printf("%f ",world->param0[i]);
      //printf("\n");
      db->size = s;
    }
}

#if defined(MPI) && !defined(PARALIO) /* */
void print_species_record_old(float *temp, long *z, world_fmt *world)
{
  long i;
  long j;
  //long locus = world->locus;
  if (world->has_speciation)
    {
      j = world->numpop2 + world->bayes->mu;
      for (i=0; i<world->species_model_size;i++,j+=2)
	{
	  temp[(*z)++] = world->species_model[i].mu; 
	  temp[(*z)++] = world->species_model[i].sigma; 
	  //temp[j] = world->param0[j]; 
	  //temp[j+1] = world->param0[j+1];
	  //(*z) += 2;
	}
    }
}

void print_species_record(float *temp, long *z, world_fmt *world)
{
  long i;
  long j;
  //long locus = world->locus;
  if (world->has_speciation)
    {
      j = world->numpop2 + world->bayes->mu;
      for (i=0; i<world->species_model_size;i++,j+=2)
	{
	  temp[(*z)++] = world->param0[world->species_model[i].paramindex_mu];
	  temp[(*z)++] = world->param0[world->species_model[i].paramindex_sigma]; 
	}
    }
}
#else /*not MPI or MPI & PARALIO*/
void  print_species_record(char *temp, long *c, world_fmt * world)
{
  //long locus = world->locus;
  long i;
  if (world->has_speciation)
    {
      for (i=0; i<world->species_model_size;i++)
	{
	  *c += sprintf(temp+ *c,"\t%f", world->param0[world->species_model[i].paramindex_mu]); 
	  *c += sprintf(temp+ *c,"\t%f",world->param0[world->species_model[i].paramindex_sigma]); 
	} 
    }
}
#endif


void read_species_record(long masternumbinall, long locus, MYREAL *params, long *n, MYREAL *oldmeans, MYREAL *lowerbound, MYREAL *upperbound, MYREAL *delta, world_fmt *world, char **inptr)
{
  const long numpop2 = world->numpop2;
  bayes_fmt *bayes = world->bayes;
  long j0;
  long j;
  long np = numpop2 + bayes->mu + world->species_model_size * 2;
  bayeshistogram_fmt *hist = &(world->bayes->histogram[locus]);
  long numbinsall = masternumbinall;
  long numbins;
  long bin;

  for(j0=numpop2+bayes->mu;j0 < numpop2+world->species_model_size*2; j0++)
    {
      if(shortcut(j0,world,&j))
	{
	  continue;
	}
      params[j+2] =  atof(strsep(inptr,"\t"));
      n[j] += 1;
      oldmeans[j] = hist->means[j];
      hist->means[j] += (params[j+2] - hist->means[j]) / n[j];
      numbinsall += hist->bins[j];
      numbins = numbinsall - hist->bins[j];
      
      if (params[j+2]>upperbound[j])
	{
	  warning("above upper bound: %f\n",params[j+2]);
	  continue;
	}
      bin = (long) ((params[j+2]-lowerbound[j]) / delta[j]);
      hist->minima[j0] = lowerbound[j0];
      hist->maxima[j0] = upperbound[j0];
      hist->results[numbins + bin] += 1.;
      bayes->histtotal[locus * np + j0] += 1;
    } 
}

void check_speciateparameters(species_fmt *s)
{
  if (s->mu < 0.0)
    {
      s->mu = 0.1;
    }
  if (s->sigma < 0.0)
    {
      s->sigma= 1.0;
    }
}
#endif /* #ifndef SPECIATETESTER */
species_fmt * get_which_species_model( long which, species_fmt * s ,  long ssize)
{
  long i;
  for(i=0;i<ssize;i++)
    {
      if (which == s[i].paramindex_mu || which == s[i].paramindex_sigma)
	return &(s[i]);
    }
  return NULL;
}

species_fmt * get_fixed_species_model( long from, long to, species_fmt * s ,  long ssize)
{
  const long siz = (long) ssize;
  long i;
  long z = 0;
  long *store = (long *) mycalloc(siz,sizeof(long));
  for ( i=0; i<siz; i++ )
    {
      if (from == s[i].from && to == s[i].to)
	store[z++] = i;
    }
    switch(z)
    {
    case 0:
      myfree(store);
      return NULL;
    case 1:
      i = store[0];
      myfree(store);
      return &(s[i]);
    default:
      //warning("more than one candidate speciation time");
      i = store[RANDINT(0,z-1)];
      myfree(store);
      return &(s[i]);
    }
}

species_fmt * get_species_model( long to, species_fmt * s ,  long ssize)
{
  const long siz = (long) ssize;
  long i;
  long z = 0;
  long *store = (long *) mycalloc(siz,sizeof(long));
  for ( i=0; i<siz; i++ )
    {
      if (to == s[i].to)
	store[z++] = i;
    }
  switch(z)
    {
    case 0:
      myfree(store);
      return NULL;
    case 1:
      i = store[0];
      myfree(store);
      return &(s[i]);
    default:
      //warning("more than one candidate speciation time");
      i = store[RANDINT(0,z-1)];
      myfree(store);
      return &(s[i]);
    }
}
#ifndef SPECIATETESTER
// using normal
MYREAL log_prob_wait_speciate_normal(MYREAL t0, MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt *s)
{
  (void) s;
  //12/3/16 normal_hazard_new.nb
  //-k Log[(1 + Erf[(mu - t0)/(Sqrt[2] sigma)])/(
  //1 + Erf[(mu - t1)/(Sqrt[2] sigma)])]
  // the k will be introduced in the calling function, ommitted here
  const MYREAL qsd = 1./(SQRT2 * sigma);
  if (0.0 > t0)
    error("negative time in speciate_normal");
  double e1 = 1 + ERF((mu - t0)*qsd);
  double e0 = 1 + ERF((mu - t1)*qsd);
  if (e0 == 0.0)
    return (double) -HUGE;
  else
    return -log(e1/e0);
}

MYREAL lower_bound_erf_approx(double x);

//MatrixForm[
// Table[{i, Log[1 - N[Erf[-i], 100]], Log[Erf[i] - -1], 
//    Log[-Erfc[i] + Erfc[-600]]}, {i, 0, -20, -1}] // N]
// this delivers the log of erf(x) - -1 which is equivalent to Erfc[-large] - Erfc[x]
// 1-erf(-large) - (1-erf(x)) = 1+1-1+erf(x) = 1 + erf(x)
MYREAL lower_bound_erf_approx(double x)
{
  return -0.3926792277130407 + 0.4220596907084323 * x - 0.9839405393215818 *  x * x;
}

// same as log_prob_wait_speciate_normal() but with b1 < Infinity
MYREAL log_prob_wait_speciate_normalorig(MYREAL t0, MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt *s)
{
  //(Erf[(b - mu)/(Sqrt[2] sigma)] + Erf[(mu - tx)/(Sqrt[2] sigma)])/
  //(Erf[(b - mu)/(Sqrt[2] sigma)] + Erf[(mu - t0)/(Sqrt[2] sigma)])
  const double b1 = (double) s->max;
  const double qsd = (double) 1./((double) SQRT2 * (double) sigma);
  if(t1 > (double) b1)
    return (double) -HUGE;
  //double lmu = (double) mu;
  double mub = mu - b1;
  //double eb1 = ERF((lmu-b1) * qsd);
  double eb1 = ERF(mub * qsd);
  double loge0;
  double loge1;
  double x0 = (mu - (double) t0) * qsd;
  double x1 = (mu - (double) t1) * qsd;
  //double e0 = ERF(x0) - eb1;
  //double e1 = ERF(x1) - eb1;
  double e0 = ERF(x0) - eb1;
  double e1 = ERF(x1) - eb1;
#ifdef DEBUG
  //printf ("@wait t1=%g x0=%g x1=%g eb1=%g e0=%g e1=%g r=%g | ", t1, x0, x1, eb1, e0,e1, e1/e0);
#endif
  if (isnan(e0) || isnan(e1) || isnan(eb1))
    error("found nan in wait_speciate_normalorig");

  //if (MYCLOSE(e1,0.0,DBL_EPSILON))
  //  {
  //    if (e0>0.0)
  //	{
  //     //Log[1 - N[Erf[-i] using a quadratic approximatioin
  //	  double  e1a = lower_bound_erf_approx(x1);
  //	  if (isnan(e1a) || isinf(e1a))
  //	    e1a = -(x1 * x1);
  //	  printf(" [%g] \n",e1a - log(e0));
  //	  return e1a - log(e0);
  //	}
  //}
  if (MYCLOSE(e0,0.0,DBL_EPSILON))
    {
      //for negative x
      loge0 = lower_bound_erf_approx(x0);
    }
  else
    loge0 = log(e0);
  if (MYCLOSE(e1,0.0,DBL_EPSILON))
    {
      //for negative x
      loge1 = lower_bound_erf_approx(x1);
    }
  else
    loge1 = log(e1);
  return loge1 - loge0;
}


// using weibull
MYREAL log_prob_wait_speciate_weibull(MYREAL t0, MYREAL t1, MYREAL lambda, MYREAL k, species_fmt *s)
{
  (void) s;
  // using mu but std is equal to k
  double mylambda = exp(log(lambda) - lgamma(1.0 + 1.0/k));
  // lambda^-k (t0^k - t1^k)
  return (pow(mylambda,-k)) * (pow(t0,k) - pow(t1,k));
}


// using an exponential (instead of normal)
// 
MYREAL log_prob_wait_speciate_exp(MYREAL t0, MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt *s)
{
  (void) s;
  (void) sigma;
  return (t0-t1) /  mu;
}



// weibull distribution instead of normal
// sigma = cv!
MYREAL log_point_prob_speciate_weibull(MYREAL t1, MYREAL lambda, MYREAL k, species_fmt * s)
{
  (void) s;
  // using mu but std is equal to k which the parameter not a #lineage
  double mylambda = exp(log(lambda) - lgamma(1.0 + 1.0/k));
  return log(k) - log(t1) + k * (-log(mylambda) + log(t1));
}


MYREAL log_point_prob_speciate_normal(MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt * s)
{
  (void) s;
  // 12/3/16
  //(Log[E^(-((mu - t1)^2/(
  // 2 sigma^2))) Sqrt[2/\[Pi]])/(k sigma Erfc[(-mu + t1)/(
  // Sqrt[2] sigma)]])
  // 1/k is used in the calling function and ommitted here
  const MYREAL qsd = 1./(SQRT2 * sigma);
  const MYREAL mut = (mu-t1);
  const MYREAL mut2 = mut*mut;
  const MYREAL var = sigma*sigma;
  MYREAL prob =  -(mut2/(2.0 * var)) + HALFLOG2DIVPI - log(sigma) - log(ERFC(-mut*qsd));
  return prob;
}

MYREAL log_point_prob_speciate_normalorig(MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt * s)
{
  const MYREAL b1 = (double) s->max;
  if (b1<t1)
    return (double) -HUGE;
  const MYREAL qsd = 1./(sqrt(2.0) * sigma);
  const MYREAL mut = (mu-t1);
  const MYREAL mub = (mu - b1);
  const MYREAL mut2 = mut*mut;
  const MYREAL var = sigma*sigma;
  // using a truncated normal (0..b1)
  //-((mu - t)^2/(2 sigma^2)) + 1/2 (-Log[2] - Log[\[Pi]]) - Log[sigma] - 
  //Log[Erfc[(-b1 + mu)/(Sqrt[2] sigma)] - 
  //    Erfc[(mu - t)/(Sqrt[2] sigma)]]
  double m = ERFC(mub*qsd)-ERFC(mut*qsd);// -erf(mub*qsd)+erf(mut*qsd)
  double logm;
  if (m < 0.0)
    error("log_point_prob_speciate_normalorig had a negative denominator!");
  
  if (m<DBL_EPSILON)
    {
      logm = lower_bound_erf_approx(mut*qsd);
      // this is an approximation using a quadratic polynomial
      //printf("hack: logm=%g\n",logm);
    }
  else
    {
      logm = log(m);
      //printf("----: logm=%g\n",logm);
    }
#ifdef DEBUG
  //printf("wait: %g %g %g \n", m,  ERFC(mub*qsd),ERFC(mut*qsd));
#endif
  MYREAL prob =  -(mut2/(2.0 * var)) + HALFLOG2DIVPI - log(sigma) - logm;//+ ERF(-mub*qsd));
  if (isnan(prob) || isinf(prob))
    error("failed in log_point_prob_normalorig\n");
  return prob;
}


MYREAL log_prob_wait_speciate_normalshortcut(MYREAL t0, MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt *s)
{
  warning("shortuct looks wrong! and also does not deliver correct results");
  const MYREAL b1 = (double) s->max;
  double tx;
  double deltatime;
  if (t1 >= 0)
    {
      tx = t1;//(t1+t0)/2.0;
      deltatime = t1 - t0;
    }
  else
    {
      tx = t0;
      deltatime = 1.0;
    }
  if (tx > b1)
    return (double) -HUGE;

  const double qsd = 1./(sqrt(2.0) * sigma);
  const double mut = (mu-tx);
  const double mub = (mu - b1);
  const double mut2 = mut*mut;
  const double var = sigma*sigma;
  double m = ERFC(mub*qsd)-ERFC(mut*qsd);
  double logm;
  if (m < 0.0)
    {
      error("log_prob_wait_speciate_normalshortcut() had a negative denominator!");
    }
  else if (m<DBL_EPSILON)
    {
      logm = lower_bound_erf_approx(mut*qsd);
    }
  else
    {
      logm = log(m);
    }
  double prob =  -(mut2/(2.0 * var)) + HALFLOG2DIVPI - log(sigma) - logm;//+ ERF(-mub*qsd));
  double expprob;
  if (isnan(prob) || isinf(prob))
    error("failed in log_prob_wait_normalorig\n");
  else
    expprob = -exp(log(deltatime) + prob);
  //printf("prob: t0=%f %f %f %f %f %f\n",t0, -deltatime, expprob, prob, -exp(log(deltatime)+prob),logm);
  return expprob;
}



MYREAL log_point_prob_speciate_normalshortcut(MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt * s)
{
  const MYREAL b1 = (double) s->max;
  if (b1<t1)
    return (double) -HUGE;
  const MYREAL qsd = 1./(sqrt(2.0) * sigma);
  const MYREAL mut = (mu-t1);
  const MYREAL mub = (mu - b1);
  const MYREAL mut2 = mut*mut;
  const MYREAL var = sigma*sigma;
  double m = ERFC(mub*qsd)-ERFC(mut*qsd);
  double logm;
  if (m < 0.0)
    {
      error("log_point_prob_speciate_normalorig had a negative denominator!");
    }
  else if (m<DBL_EPSILON)
    {
      logm = lower_bound_erf_approx(mut*qsd);
    }
  else
    {
      logm = log(m);
    }
  double prob =  -(mut2/(2.0 * var)) + HALFLOG2DIVPI - log(sigma) - logm;
  if (isnan(prob) || isinf(prob))
    error("failed in log_point_prob_normalorig\n");
  return prob;
}


// exp distribution instead of normal
MYREAL log_point_prob_speciate_exp(MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt * s)
{
  (void) s;
  (void) sigma;
  (void) t1;
  return -log(mu);
}


char set_type(world_fmt *world, long topop, long frompop, char *custm2, long numpop)
{
  (void) frompop;
  (void) custm2;
  (void) numpop;
  char type = 'm';
  //long i = mm2m(frompop, topop, numpop);
  species_fmt *s = NULL; 
  s = get_fixed_species_model(frompop, topop, world->species_model, world->species_model_size); 
  if (s!=NULL)
    {
      type = 'd';
    }
  //1229if (custm2[i]=='d' || custm2[i]=='D')
  //1229   type = 'd';
  return type;
}

node * set_type2(world_fmt *world, node *p, node *q, char *custm2)
{
  (void) custm2;
  
  long topop = p->actualpop;
  long frompop = q->actualpop;
  //long numpop = world->numpop;
  node *tmp=NULL;
  //char *custm2 = world->options->custm2;
  char type = 'm';
  species_fmt *s = NULL; 
  s = get_fixed_species_model(frompop, topop, world->species_model, world->species_model_size); 
  if (s!=NULL)
    {
      type = 'd';
      tmp = add_migration (world, p, type, q->actualpop, p->actualpop,
			   (MYREAL) RANDDOUBLE(0.0, q->tyme - p->tyme));
      //while (frompop != s->from)
      //	{  
      //	  s = get_species_model(s->from, world->species_model, world->species_model_size);
      //	  if (s==NULL)
      //	    break;
      //	  tmp = add_migration (world, tmp, type, q->actualpop, tmp->actualpop,
      //		       (MYREAL) RANDDOUBLE(0.0, q->tyme - tmp->tyme));
      //	}
    }
  return tmp;
}


//invert the speciation probability to get the time of event of split
// getprob[mu_, std_, b1_, t0_, t1_] := 
//                   Log[Erf[(-b1 + mu)/(Sqrt[2] std)] - Erf[(mu - t1)/(Sqrt[2] std)] ] - 
//                   Log[Erf[(-b1 + mu)/(Sqrt[2] std)] - Erf[(mu - t0)/(Sqrt[2] std)]]
// In: Solve[Log[r] == getprob[mu, std, b1, t0, t1], t1]
// Out: t1 <- mu - Sqrt[2] std InverseErf[(-1 + r) Erf[(b1 - mu)/(Sqrt[2] std)] + 
//		       r Erf[(mu - t0)/(Sqrt[2] std)]]

// exponential distribution
MYREAL time_to_speciate_exp(world_fmt *world, long pop, MYREAL t0, char *event, long *to, long *from)
{
  (void) t0;
  long fromi;
  long ssize = world->species_model_size;
  species_fmt *s;
  MYREAL mu;
  double interval= (double) HUGE;
  double intervalnew;
  boolean touched=FALSE;
  double mlalpha = world->mlalpha;
  //double mlinheritance = world->mlinheritance;
  boolean has_mlalpha = mlalpha < 1.0;

  *from = -1;
  *to = -1;
  for (fromi=0; fromi < world->numpop; fromi++)
    {
      if (fromi==pop)
	continue;
      s= get_fixed_species_model(fromi, pop, world->species_model, ssize);
      if(s==NULL)
	continue;
      mu = world->param0[s->paramindex_mu];
      //double priormin = s->min;
      //double priormax = s->max; 
      if(has_mlalpha)
	{
	  intervalnew = propose_new_mlftime(1.0/mu, mlalpha, UNIF_RANDUM(), UNIF_RANDUM());
	  //intervalnew = interval_mittag_leffler(UNIF_RANDUM(), mlalpha, 1.0/mu,tmin,tmax);
	}
      else
	{
      intervalnew =  (-(LOG (UNIF_RANDUM ())) * mu);
	}
      if (intervalnew < interval)
	{
	  interval = intervalnew;
	  *event = 'd';
	  *from = s->from; //shoudl be the same as from
	  *to = s->to;//should be the same as pop
	  touched=TRUE;
	}
    }
  //#if DEBUG
  //fprintf(proptimes,"proposedtimes> %li %li %f\n", *from, *to, interval);
  //#endif
  if(!touched)
    return (double) HUGE;
  else
    return interval;
}


// weibull
MYREAL time_to_speciate_weibull(world_fmt *world, long pop, MYREAL t0, char *event, long *to, long *from)
{
  //  f[a_, b_, t0_] := 
  //        Module[{value = 0, count = 0}, 
  //           While[value < t0 && count++ < 100, 
  //                 value = b (-Log[RandomReal[]])^(1/a)]; Return[value]]
  double r;
  long ssize = world->species_model_size;
  species_fmt *s;
  double t1 = (double) HUGE;
  double tnew;
  MYREAL lambda;
  MYREAL k;
  double mu;
  //double tmax;
  MYREAL interval;
  double lg;
  long sfrom = -1;
  long fromi;
  boolean touched=FALSE;
  double mlalpha = world->mlalpha;
  boolean has_mlalpha = mlalpha < 1.0;

  *from = -1;
  *to = -1;
  for (fromi=0; fromi < world->numpop; fromi++)
    {
      if (fromi==pop)
	continue;
      s= get_fixed_species_model(fromi, pop, world->species_model, ssize);
      if(s==NULL)
	continue;
      mu = world->param0[s->paramindex_mu];
      k = world->param0[s->paramindex_sigma];
      //double priormin = s->min;
      //double priormax = s->max;
      if (has_mlalpha)
	{
	  //r = UNIF_RANDUM();
	  interval = propose_new_mlftime(1.0/mu, mlalpha, UNIF_RANDUM(), UNIF_RANDUM());
	    //interval_mittag_leffler_func(r, mlalpha, t0, mu, k, s,priormin,priormax);
	  tnew = t0 + interval;
	}
      else
	{
      lg = lgamma(1.0 + 1.0/k);
      lambda = exp(log(mu) - lg);
      //tmax = (double) s->max;
      r = log(UNIF_RANDUM());
      //{{t1 -> (t0^k - lambda^k Log[rr])^(1/k)}}
      tnew = pow(pow(t0,k) - pow(lambda,k) * r,(1.0/k));
	}
      if(tnew < t1)
	{
	  t1= tnew;
	  sfrom = s->from;
	  touched=TRUE;
	}
    }
  if (!touched)
    return (double) HUGE;

  if (t1<t0)
    t1 = t0 + EPSILON;

  *event = 'd';
  *from = sfrom;
  *to = pop;

  interval =  t1-t0;  
  return interval;
}

// version of this function are in the section SPECIESTESTER(with main)  see below!!!
MYREAL time_to_speciate_normal_erf(world_fmt *world, long pop, MYREAL t0, char *event, long *to, long *from)
{
  long ssize = world->species_model_size;
  species_fmt *s;
  MYREAL r;
  MYREAL tmax;
  MYREAL t1= (double) HUGE;
  double tnew;
  double newtmax = 0.0;
  MYREAL mu;
  double newmu = (double) HUGE;
  MYREAL sigma;
  double newsigma= (double) HUGE;
  double qsd;
  MYREAL interval;
 long sfrom = -1;
 long fromi;
 boolean touched=FALSE;
 *from = -1;
 *to = -1;
 for (fromi=0; fromi < world->numpop; fromi++)
    {
      if (fromi==pop)
	continue;
      s= get_fixed_species_model(fromi, pop, world->species_model, ssize);
      if(s==NULL)
	continue;
      mu = world->param0[s->paramindex_mu];
      sigma = world->param0[s->paramindex_sigma];
      tmax = (double) s->max;
      qsd = 1./(SQRT2 * sigma);
      r = UNIF_RANDUM();
      //double eb1 = ERF((tmax-mu) * qsd); // is always +
      double eb0 = ERF((mu - t0) * qsd); // if t0 < mu + otherwise -
      double cdfval = -1.0 + r + r * eb0;
      if (cdfval <= -1.0)
	{
	  tnew = (double) HUGE; //4.2.9
	}
      else if (cdfval >= 1.0)
	{
	  tnew = t0>mu ? t0+EPSILON : mu+EPSILON; //4.2.9
	}
      else
	{
	  tnew = mu - SQRT2 * sigma * erfinv(cdfval);
	}
      if (tnew < t1)
	{
	  t1 = tnew;
	  sfrom = s->from;
	  newmu = mu;
	  newsigma=sigma;
	  newtmax= tmax;
	  touched=TRUE;
	}
    }
 if(!touched)
   return (double) HUGE;
 *event = 'd';
 *from = sfrom;
 *to = pop;
  if (t0 <= t1 && t1 < newtmax)
    {
      interval =  t1-t0;  
      return interval;
    }
  else
    {  //4.2.9
      if (newmu + 2.0 * newsigma < t0)
	return EPSILON;
      if (t1 - t0 > -10e-8)
	return EPSILON;
      else
	return EPSILON;
    }
}


// normal distribution using rannor
MYREAL time_to_speciate_normal(world_fmt *world, long pop, MYREAL t0, char *event, long *to, long *from)
{
  //12/3/16
  //mu - Sqrt[2] sigma InverseErf[-1 + r + r Erf[(mu - t0)/(Sqrt[2] sigma)]]
  long ssize = world->species_model_size;
  species_fmt *s;// = get_species_model(pop, world->species_model, ssize);
  MYREAL tmax;
  MYREAL t1 =(double) HUGE;
  MYREAL newtmax = 0.0;
  double tnew;
  MYREAL mu;
  MYREAL sigma;
  MYREAL interval;
  long sfrom = -1;
  long fromi;
  boolean touched=FALSE;
  *from = -1;
  *to = -1;
  for (fromi=0; fromi < world->numpop; fromi++)
    {
      if (fromi==pop)
	continue;
      s= get_fixed_species_model(fromi, pop, world->species_model, ssize);
      if(s==NULL)
	continue;
      mu = world->param0[s->paramindex_mu];
      sigma = world->param0[s->paramindex_sigma];
      tmax = (double) s->max;
      if (t0<tmax)
	tnew = trunc_random_normal2(t0, tmax, mu, sigma);
      else
	tnew = EPSILON;
      if(tnew < t1)
	{
	  t1 = tnew;
	  touched=TRUE;
	  newtmax = tmax;
	  sfrom = s->from;
	}
    }
  if(!touched)
    return (double) HUGE;
  *event = 'd';
  *from = sfrom;
  *to = pop;
  if (t0 <= t1 && t1 < newtmax)
    {
      interval =  t1-t0;  
      return interval;
    }
  else
    {
      if ((t1>newtmax) || t1 < 0.0)
	return EPSILON;
      else
	return (double) HUGE;
    }
}

//invert the speciation probability to get the time of event of split
// getprob[mu_, std_, b1_, t0_, t1_] := 
//                   Log[Erf[(-b1 + mu)/(Sqrt[2] std)] - Erf[(mu - t1)/(Sqrt[2] std)] ] - 
//                   Log[Erf[(-b1 + mu)/(Sqrt[2] std)] - Erf[(mu - t0)/(Sqrt[2] std)]]
// In: Solve[Log[r] == getprob[mu, std, b1, t0, t1], t1]
// Out: t1 <- mu - Sqrt[2] std InverseErf[(-1 + r) Erf[(b1 - mu)/(Sqrt[2] std)] + 
//		       r Erf[(mu - t0)/(Sqrt[2] std)]]
MYREAL time_to_speciate_normalorig(world_fmt *world, long pop, MYREAL t0, char *event, long *to, long *from)
{
  long ssize = world->species_model_size;
  species_fmt *s;
  MYREAL r;
  MYREAL tmax;
  double newtmax = 0.0;;
  MYREAL t1 = (double) HUGE;
  MYREAL tnew;
  MYREAL mu;
  //double newmu = (double) HUGE;
  MYREAL sigma;
  double qsd;
  MYREAL interval;
  long sfrom = -1;
  long fromi;
  boolean touched=FALSE;
  *from = -1;
  *to = -1;
  for (fromi=0; fromi < world->numpop; fromi++)
    {
      if (fromi==pop)
	continue;
      s= get_fixed_species_model(fromi, pop, world->species_model, ssize);
      if(s==NULL)
	continue;
      mu = world->param0[s->paramindex_mu];
      sigma = world->param0[s->paramindex_sigma];
      tmax = (double) s->max;
      qsd = 1./(SQRT2 * sigma);
      //double myhalftheta = world->param0[s->to]/2.;
      r = UNIF_RANDUM(); //
      //r = pow(UNIF_RANDUM(),myhalftheta);
      double eb1 = ERF((mu-tmax) * qsd); 
      double eb0 = ERF((mu - t0) * qsd);
      double cdfval = r*eb0 + (1.0-r)*eb1;
      if (cdfval <= -1.0)
	{
	  tnew = (double) HUGE;
	}
      else if (cdfval >= 1.0)
	{
	  //tnew = t0>mu ? t0+EPSILON : mu+EPSILON;
	  tnew = t0 + EPSILON * UNIF_RANDUM();
	}
      else
	{
	  tnew =  mu - SQRT2 * sigma * (erfinv(cdfval));
	}
      if(tnew < t1)
	{
	  t1 = tnew;
	  //  newmu = mu;
	  newtmax = tmax;
	  touched=TRUE;
	  sfrom = s->from;
	}
    }
  if(!touched)
    return (double) HUGE;
  *event = 'd';
  *from = sfrom;
  *to = pop;
  if (t0 <= t1 && t1 < newtmax)
    {
      interval =  t1-t0;  
      return interval;
    }
  else 
    {
      return EPSILON;
    }
}

// normal shortcut distribution
MYREAL time_to_speciate_normalshortcut(world_fmt *world, long pop, MYREAL t0, char *event, long *to, long *from)
{
  long ssize = world->species_model_size;
  species_fmt *s;
  MYREAL r;
  MYREAL tmax;
  double newtmax = 0.0;
  MYREAL t1 = (double) HUGE;
  MYREAL tnew;
  MYREAL mu;
  double newmu = (double) HUGE;
  MYREAL sigma;
  double qsd;
  MYREAL interval;
  long sfrom = -1;
  long fromi;
  boolean touched=FALSE;
  double mlalpha = world->mlalpha;
  boolean has_mlalpha = mlalpha < 1.0;

  *from = -1;
  *to = -1;
  for (fromi=0; fromi < world->numpop; fromi++)
    {
      if (fromi==pop)
	continue;
      s= get_fixed_species_model(fromi, pop, world->species_model, ssize);
      if(s==NULL)
	continue;
      mu = world->param0[s->paramindex_mu];
      sigma = world->param0[s->paramindex_sigma];
      tmax = (double) s->max;
      qsd = 1./(SQRT2 * sigma);
      r = UNIF_RANDUM();
      //double priormin = s->min;
      //double priormax = s->max;
      if (has_mlalpha)
	{
	  interval = propose_new_mlftime(mu, mlalpha, UNIF_RANDUM(), UNIF_RANDUM());
	  //interval = interval_mittag_leffler_func(r, mlalpha, t0, mu, sigma, s, priormin, priormax);
	  tnew = t0 + interval;
	}
      else
	{
	  double eb1 = ERF((tmax-mu) * qsd); 
	  double eb0 = ERF((mu - t0) * qsd);
	  double cdfval = r*eb0 + (-1.0+r)*eb1;
	  if (cdfval <= -1.0)
	    {
	      //warning("cdval is <= -1 ==> infinite waiting");
	      tnew = (double) HUGE;
	    }
	  else if (cdfval >= 1.0)
	    {
	      if(t0>mu)
		tnew = t0 + EPSILON * UNIF_RANDUM();
	      else
		tnew = mu + EPSILON * UNIF_RANDUM();
	    }
	  else
	    {
	      tnew =  mu - SQRT2 * sigma * (erfinv(cdfval));
	      if (tnew < t0)
		tnew = t0 + EPSILON * UNIF_RANDUM();
	    }
	}
      if(tnew < t1)
	{
	  t1 = tnew;
	  newmu = mu;
	  newtmax = tmax;
	  touched=TRUE;
	  sfrom = s->from;
	}
    }
  if(!touched)
    return (double) HUGE;
  *event = 'd';
  *from = sfrom;
  *to = pop;
  if (t0 <= t1 && t1 < newtmax)
    {
      interval =  t1-t0;  
      return interval;
    }
  else 
    {
      if (100.0 * newmu < t0)
	return EPSILON;
      if (t1 - t0 > -10e-8)
	return EPSILON;
      else
	return EPSILON;
    }
}



MYREAL time_to_coalescence(world_fmt * world, long pop, double age, long timeslice, long *lineages, char * event, long * to, long *from)
{
  
  double * growth = NULL;
  long * growpops = world->options->growpops;
  if (world->has_growth)
    {
      growth = world->growth;
    }
  world_fmt *w = world;
  long  lines;
  MYREAL  rate = w->options->mu_rates[w->locus];
  //long    timeslice;
  long    addition=0;
  long    timepop;
  MYREAL  inheritance = w->options->inheritance_scalars[w->locus];
  MYREAL  timethetarate;
  MYREAL denom=0.0;
  MYREAL invdenom;
  MYREAL interval= (double) HUGE;
  MYREAL r=0.0;
  MYREAL logr = (double) -HUGE;
  double mlalpha = world->mlalpha;
  //double mlinheritance = world->mlinheritance;
  boolean has_mlalpha = mlalpha < 1.0;
  *from = -1;
  *to = -1;
  //timeslice=tentry->timeslice;
  timepop = (w->numpop2+addition)*timeslice + pop;
  timethetarate = w->timek[timepop] * w->param0[pop];
  //double priormin = -10.;//w->bayes->minparam[pop];
  //double priormax = 10;//w->bayes->maxparam[pop];
  lines    = 2 * lineages[pop];
  if (lines == 0)
    {
      interval= (double) HUGE;
    }
  else
    {
      if(has_mlalpha)
	{	  
	  //DEBUG N^a discussion: denom = lines /(rate*timethetarate);
	  //-(1/(k (k - 1)/(theta^a)))^(1/
	  //a) ((Sin[Pi a]/(Tan[Pi a (1 - RandomReal[])])) - Cos[Pi a])^(1/
	  //a) Log[RandomReal[]]
	  denom = lines /(rate*timethetarate);
	  //interval = interval_mittag_leffler(r, mlalpha, denom,priormin,priormax);
	  interval = propose_new_mlftime(denom, mlalpha, UNIF_RANDUM(), UNIF_RANDUM());
	}
      else
	{
	  if(world->has_growth && growpops[pop]!=0 && fabs(growth[growpops[pop]-1])>EPSILON)
	    {
	      //invdenom    = (rate*timethetarate) / lines;		
	      //interval =  -(logr * invdenom);	     
	      //interval = interval_growth(r, age, world->param0[pop], growth[growpops[pop]-1], lines, 0.0, 1.0);
	      interval = get_time_for_growth(rate*timethetarate, growth[growpops[pop]-1], lines, age);
	    }
	  else
	    {
	      r = UNIF_RANDUM ();
	      logr = log(r);
	      invdenom    = (rate*timethetarate) / lines;
	      interval =  -(logr * invdenom);
	    }
	}
    }
  if (interval<0.0)
    warning("%i> x=%f in timecoal [denom=%f lines=%li r=%f inh=%f rate=%f timethetarate=%f]\n",myID, 
	    interval,denom,lines,r,inheritance,rate,timethetarate);
  *event = 'c';
  *to = pop;
  *from = pop;
  return interval;
}


MYREAL time_to_migration(proposal_fmt *proposal, world_fmt *world, long pop, long timeslice, long * lineages, char * event, long * to, long * from)
{
  (void) lineages;
  (void) proposal;
  long    addition=0;
  long    numpop = world->numpop;
  MYREAL  *skyparam = world->timek+(world->numpop2+addition)*timeslice;
  char    *custm2 = world->options->custm2;
  MYREAL  rate = world->options->mu_rates[world->locus];
  //MYREAL  invrate = 1./ rate;
  //long    timepop = (proposal->world->numpop2+addition)*timeslice + pop;
  //MYREAL  inheritance = proposal->world->options->inheritance_scalars[proposal->world->locus];
  //MYREAL  timethetarate = proposal->world->timek[timepop] * proposal->param0[pop];
  //MYREAL denom;
  MYREAL invdenom;
  long msta = world->mstart[pop];
  long msto = world->mend[pop];
  MYREAL mm = 0.0;
  long i;
  MYREAL interval=0.0;
  MYREAL the_eventtime= (double) HUGE;
  //MYREAL eventtime;
  char the_event=' ';
  //char myevent=' ';
  long tox = -1;
  long fromx = -1;
  long the_to = -1;
  long the_from = -1;
  long offlimit_from = -1; //nothing is offlimit
  double mlalpha = world->mlalpha;
  boolean has_mlalpha = mlalpha < 1.0;
  
  //if (proposal!=NULL && proposal->divlist->div_elem > 0) 
  //  offlimit_from = proposal->divlist->divlist[0][1]; //4.2.9
#ifdef OFFLIMIT_ON  
  if(proposal!=NULL 
     && proposal->migr_table_counter>0 
     && proposal->migr_table[proposal->migr_table_counter-1].event=='d')
    {
      // if we have encountered a divergence then we will not get a migration
      // back to the other population: e.g. D: from=1 to=2 --> no M from 2-->1
      // for the lineages proposed: custom={*D**} eseentially results for this 
      // lineage into {*d**}
      offlimit_from = proposal->migr_table[proposal->migr_table_counter-1].to; 
    }
#endif
  for (i = msta; i < msto; i++)
    { 
      //double priormin = world->bayes->minparam[i];
      //double priormax = world->bayes->maxparam[i];
      switch (custm2[i])
	{
	case 'd':
	  error("should not go here");
	  // break;
	case '0':
	  continue;
	default: /*migration with *, m, M, s, S, or c */
	  if(skyparam[i]<=0.0)
	    {
	      warning("skyparam[%li]=%f rate=%f\n",i, skyparam[i],rate);
	      error("mismatch with skyparams\n");
	    }
	  m2mm(i,numpop,&fromx,&tox);
	  if(offlimit_from == fromx)
	    {
	      interval= (double) HUGE;
	    }
	  else
	    {
	      mm =  world->data->geo[i] * world->param0[i] * skyparam[i]/rate;
	      if(has_mlalpha)
		{
		  //DEBUG N^a 		  
		  interval = propose_new_mlftime(mm, mlalpha, UNIF_RANDUM(), UNIF_RANDUM());
		  //interval = interval_mittag_leffler(UNIF_RANDUM(), mlalpha, mm, priormin, priormax);
		}
	      else
		{
	      invdenom = 1.0 / (mm);
	      interval =  (-(LOG (UNIF_RANDUM ())) * invdenom);
	    }
	    }
	  if (the_eventtime > interval)
	    {
	      the_eventtime = interval;
	      the_event = 'm';
	      the_to = tox;
	      the_from = fromx;
	    }
	}
    }
  *to = the_to;
  *from = the_from;
  *event = the_event;
#ifdef DEBUG
  //if (the_event == 'd')
  //  printf("%i@ %li <- %li: @time:%f (%f)\n",myID, the_to, the_from, tentry->age + interval, interval);
#endif
  return the_eventtime;
}

///
/// standard time evaluator for a single line, the time_to_speciate() is also called
/// when the last two lines in the tree need to be worked on jointly
MYREAL time_to_speciation(world_fmt * world, long pop, double age, char * event, long * to, long * from)
{
  MYREAL the_eventtime= (double) HUGE;
  //char the_event;
  char myevent=' ';
  long tox = -1;
  long fromx = -1;
  *from = -1;
  *to = -1;

  MYREAL interval;
  
  interval = (*time_to_speciate)(world,pop,age, &myevent, &tox, &fromx);

  if (interval < the_eventtime)
    {
      the_eventtime = interval;
      *event = myevent;
      *to = tox;
      *from = fromx;
    }
  return the_eventtime;
}


void keep_min_eventtime(MYREAL *the_eventtime, MYREAL eventtime, char *the_event, char myevent,
			long *to, long tox, long *from, long fromx)
{
  if (eventtime < *the_eventtime)
    {
      *the_eventtime = eventtime;
      *the_event = myevent;
      *to = tox;
      *from = fromx;
    }
}

MYREAL eventtime_single(proposal_fmt *proposal, world_fmt *world, long pop, long timeslice, long *lineages, double age, char * event, long *to, long *from)
{
  (void) proposal;
  char the_event = ' ';
  double the_eventtime= (double) HUGE;
  double ceventtime;
  double meventtime;
  double deventtime;
  char myevent=' ';
  long fromx = *from;
  long tox = *to;
  ceventtime = time_to_coalescence(world, pop, age, timeslice, lineages, &myevent,&tox, &fromx);
  keep_min_eventtime(&the_eventtime, ceventtime, &the_event, myevent, to, tox, from, fromx);
  meventtime = time_to_migration(proposal, world, pop, timeslice, lineages, &myevent,&tox, &fromx);
  keep_min_eventtime(&the_eventtime, meventtime, &the_event, myevent, to, tox, from, fromx);
  deventtime = time_to_speciation(world, pop, age, &myevent, &tox, &fromx);
  keep_min_eventtime(&the_eventtime, deventtime, &the_event, myevent, to, tox, from, fromx);
  //  }
  *event = the_event;
  if (*event == 'm' && *from == *to)
    error("FAIL");
  return the_eventtime;
}

long speciation_from(long to, proposal_fmt * proposal)
{
  species_fmt * s = proposal->world->species_model;
  const long siz = proposal->world->species_model_size;
  long i;
  for ( i=0; i<siz; i++ )
    {
      if (to == s[i].to)
	return s[i].from;
    }
  warning("Inclusion of a population splitting event failed\n");
  warning("because target was not found, this sample was ignored, but run continues\n");
  return -1;
}


void loopcleanup(boolean assign, world_fmt * world, proposal_fmt *proposal, long oldpop, long newpop, timelist_fmt *timevector)
{
  (void) world;
  if(assign)
    {
      reassign_individual(proposal->origin,oldpop,newpop);
      if(world->cold)
	{
	  record_assignment(world->locus, world); // records failures
	}
    }
  free_timevector (timevector);
#ifndef TESTING2
  free_masterproposal (proposal);
#endif
}

/*
return 1 if tree was accepted, 0 otherwise 
assign is set for assigning individuals, caller is responsible to reset the 
origin back to the original state when fail
*/
long newtree_update (world_fmt * world, long g, boolean assign)
{
  long assignindex = -1;
  boolean coalesced;
  char event;
  long slider;
  long bordernum;
  long actualpop = -99, zz;
  long skyz=0;
  MYREAL endtime, nexttime, age;
  long newpop = -1;;
  long oldpop;
  long from = -1;
  long to = -1;
  proposal_fmt *proposal=NULL; 
  timelist_fmt *timevector = NULL; /* local timelist */
  vtlist *tentry = NULL; /*pointer into timeslice */
  MYREAL x;
  //long i,ii;
  //long np = world->numpop2 + world->bayes->mu + world->species_model_size * 2;
  
  //if (assign && world->unassignednum<2)
  //  return 0;
  
  new_localtimelist (&timevector, &world->treetimes[0], world->numpop);
  new_proposal (&proposal, &world->treetimes[0], world);
  
    if(assign)
      {
	assignindex = chooseUnassigned (proposal);
	if (proposal->origin == NULL)
	  {
	    free_timevector(timevector);
#ifndef TESTING2
	    free_masterproposal (proposal);
#endif
	    //warning("No origin specified\n");
	    return 0;
	  }
	oldpop = proposal->origin->actualpop;
	newpop = oldpop;
	if (world->numpop>1)
	  {
	    while (newpop == oldpop)
	      {
		if (world->has_unassignedfreq)
		  newpop = assign_bypastfreq(proposal->origin,world->numpop);
		else
		  {
		    if (world->has_unassignedpoplist)
		      newpop = ( long) world->unassignedpoplist[RANDINT(0,(long) world->unassignedpoplistnum-1)];
		    else
		      newpop = RANDINT(0, (long) world->numpop-1);
		  }
	      }
	  }
	reassign_individual(proposal->origin,newpop, -1);
      }
    else
      {
	chooseOrigin (proposal);
	oldpop = proposal->origin->actualpop;
      }

    if(world->options->bayes_infer)
    {
      //negative value to flag for no-change later
      world->bayes->starttime = -proposal->origin->tyme;    
    }

    construct_localtimelist (timevector, proposal);
    //if(proposal->divlist->div_elem>0) //4.2.9
	// {
    //	memcpy(proposal->param0save, proposal->world->param0, sizeof(MYREAL)* (size_t) np);
    // 	for (i=0;i<proposal->divlist->div_elem;i++)
    //	  {
    //	    ii = m2mmm(proposal->divlist->divlist[i][0],
    //		       proposal->divlist->divlist[i][1],proposal->world->numpop);
    //	    proposal->world->param0[ii]=0.0;
    //	  }
    //}
    tentry = &(*timevector).tl[0];
    age = proposal->origin->tyme;
    if(age<0.0)
      {
	error("Abort time was negative");
      }
    zz = 0;
    // finding timeslice int the timelist with age
    while ((tentry->age < age || tentry->age - age < SMALLEPSILON)&& zz < (*timevector).T)
    {
        tentry = &(*timevector).tl[zz];
        zz++;
    }
    zz--;
 
    // adjusting the timeinterval for the skylineplots
    skyz=world->timeelements-1;
    while (skyz > 0 && tentry->age <  proposal->world->times[skyz])
      skyz--;
    skyz++;
    nexttime =  proposal->world->times[skyz];
    tentry->timeslice=skyz-1;
   
    if (tentry->age < nexttime)
      nexttime = tentry->age;

    if ((*timevector).T > 1)
        endtime = (*timevector).tl[(*timevector).T - 2].age;
    else
        endtime = 0.0;

    proposal->time = age;
    coalesced = FALSE;
    /*------------------------------------
     main loop: sliding down the tree  */
    slider = 0;
    //printf("-@---------------------------------------------------\n");
    while (nexttime <= endtime)
    {
        actualpop =
        (proposal->migr_table_counter >
         0) ? proposal->migr_table[proposal->migr_table_counter -
                                   1].from : proposal->origin->pop;

	x = eventtime_single(proposal,world, actualpop, tentry->timeslice, tentry->lineages, /*tentry->*/ age, &event, &to, &from);
	//if (x>=HUGE)
	//  {
	//    fprintf(stdout,"%i> X Proposal failed: time=%f origin time=%f nexttime=%f, deltatime=%f\n", myID, proposal->time, proposal->origin->tyme, nexttime, age);
	//    loopcleanup(assign, world, proposal, oldpop, timevector);
	//    return 0;
        //}
	  
	x += EPSILON*UNIF_RANDUM();
	proposal->time = age + x;
	//fprintf(world->options->logfile,"#proposal %li %li %c %f\n",from, to, event, proposal->time);
	if(proposal->time < 0.0 || isnan(proposal->time))
	  {
	    error("abort");
	  }
        if(proposal->time < proposal->origin->tyme)
        {
	  fprintf(stdout,"%i> Proposal failed because of unordered entry in time list, abort this sample\nProposed time=%f origin time=%f nexttime=%f, deltatime=%f\n", myID, proposal->time, proposal->origin->tyme, nexttime, age);
            // we end up here when the migration events exceed the upper limit
	  loopcleanup(assign, world, proposal, oldpop, newpop, timevector);
	  //warning("Too many migrations?\n");
	  return 0;
        }
        if (proposal->time < nexttime)
        {
            if (event == 'm' || event == 'd')
            {
	      if (!migrate (proposal, proposal->origin, event, from, to))
                {
		  // we end up here when the migration events exceed the upper limit
		  loopcleanup(assign, world, proposal, oldpop, newpop, timevector);
		  return 0;
                }
	      age = proposal->time;
	      //printf("%c: to=%li from=%li time=%f\n", event,
	      //       proposal->migr_table[proposal->migr_table_counter - 1].to,
	      //       proposal->migr_table[proposal->migr_table_counter - 1].from,
	      //       proposal->migr_table[proposal->migr_table_counter - 1].time);
	      continue;
            }
            else
            {   /*coalesce */

	      //printf("%i> %c: pop=%li time=%f\n", myID, event, actualpop, proposal->time); 
	      //if (event != 'c')
	      //	error("Abort in new_treeupdate(), line 1024");
	      //printf("**\n");
                chooseTarget (proposal, timevector, proposal->bordernodes,
                              &bordernum);
                if(bordernum == 0)
                {
		  if (!migrate (proposal, proposal->origin, event, from, to))
                    {
		      // we end up here when the migration events exceed the upper limit
		      loopcleanup(assign, world, proposal, oldpop, newpop, timevector);
		      //warning("too many migrations\n");
		      return 0;
                    }
                    age = proposal->time;
                    continue;
                }
		//if(proposal->divlist->div_elem>0)
		//  {
		//    memcpy(proposal->world->param0, proposal->param0save,  sizeof(MYREAL)* (size_t) np);
		//  }
                pretendcoalesce1p (proposal);
                coalesced = TRUE;
                break;
            }
        }   /*end if proposal->time < nextime */
        age = nexttime;
        zz++;
        if(zz >= timevector->T)
        {
            break;
        }
        tentry = &(*timevector).tl[zz]; /*next entry in timelist */
	// skyparam
	if (tentry->age > proposal->world->times[skyz])
	  {
	    zz--;
	    tentry = &(*timevector).tl[zz]; /*next entry in timelist */
	    tentry->timeslice = skyz;
	    nexttime =  proposal->world->times[skyz];
	    skyz++;
	  }
	else
	  {
	    nexttime = tentry->age;
	    tentry->timeslice = skyz-1;
	  }
    }
    if (!coalesced)
    {
        if (!beyond_last_node(proposal, (*timevector).tl, (*timevector).T - 1, &slider))
        {
	  loopcleanup(assign, world, proposal, oldpop, newpop, timevector);
	  //warning("problem with beyond last node\n");
	  return 0;
        }
	//if(proposal->divlist->div_elem>0)
	//{
	//  memcpy(proposal->world->param0, proposal->param0save,  sizeof(MYREAL)* (size_t) np);
	//}
        pretendcoalesce1p (proposal);
    }
    if (acceptlike (world, proposal, g, timevector, assign))
      {
        if (proposal->time > world->root->tyme)
        {   /*saveguard */
            world->root->tyme += proposal->time;
        }
#ifdef DEBUGXX
	int i;
	for(i=0;i<proposal->migr_table_counter; i++)
	  {
	    printf("%i> PMT:  %li %li, %f %c \n",myID,proposal->migr_table[i].from,proposal->migr_table[i].to,proposal->migr_table[i].time,proposal->migr_table[i].event);
	  }
	for(i=0;i<proposal->migr_table_counter2; i++)
	  {
	    printf("%i> PMT2: %li %li, %f %c\n",myID,proposal->migr_table2[i].from,proposal->migr_table2[i].to,proposal->migr_table2[i].time,proposal->migr_table2[i].event);
	  }
	printf("%i> PM: %li %li, %f %c ",myID,proposal->target->pop,proposal->target->actualpop,proposal->time,'c');
	printf("@\n");
#endif  /*end DEBUGXX */
        coalesce1p (proposal);
	if(world->cold && assign)
	  {
	    proposal->world->unassigned[assignindex]->accept[world->locus] += 1;
	    record_assignment(world->locus, world);
	  }
        if(world->options->bayes_infer)
        {
            world->bayes->starttime = -world->bayes->starttime;
            world->bayes->stoptime = proposal->time;
            world->treelen = 0.0;
            calc_treelength (world->root->next->back, &world->treelen);
        }
        //
        world->likelihood[g] = treelikelihood (world);
        /* create a new timelist */
        construct_tymelist (world, &world->treetimes[0]);
        world->migration_counts = 0;
        /* report the number of migration on the tree */
        count_migrations (world->root->next->back, &world->migration_counts);
        free_timevector (timevector);
#ifndef TESTING2
	free_masterproposal (proposal);
#endif
        return 1;   /* new tree accepted */
      }
    else
      {
        // record the time interval that was used for the lineage.
        if(world->options->bayes_infer)
	  {
            world->bayes->stoptime = -proposal->time;
	  }
	loopcleanup(assign, world, proposal, oldpop, newpop, timevector);
	//warning("acceptlike rejected\n");
	return 0;
      }
    error("do not go here in new_tree_update");
    //return -1;
}

void set_things(MYREAL t1, MYREAL t2, char e1, char e2, long to1, long from1, long to2, long from2, 
		MYREAL *time, char *event, long *to, long *from)
{
  if (t1 < t2)
    {
      *time = t1 + EPSILON * UNIF_RANDUM();
      *event = e1;
      *to = to1;
      *from = from1;
    }
  else
    {
      *time = t2 + EPSILON * UNIF_RANDUM();
      *event = e2;
      *to = to2;
      *from = from2;
    }
}

// calculates the time=age and event for each lineage (of the final 2) 
// in beyond_....()
// augmented for mittag-leffler alpha
void calc_time_per_line(proposal_fmt * proposal, char type, long popto, boolean same,MYREAL *time, char *event, long *to, long *from)
{
  (void) type;
  (void) popto;
  char events;
  long i,j;
  long   pop = *to;
  MYREAL age = *time;
  //long f; //from
  long t; //to
  long frompop=-1;
  long froms = -1;
  long tos = -1;
  double temptime;
  double time_coal = (double) HUGE;
  double time_mig  = (double) HUGE;
  double time_mig_old;
  double time_spec = (double) HUGE;
  world_fmt * world = proposal->world;
  double  r0 = LOG(UNIF_RANDUM());
  double  r1;
  long msta = world->mstart[pop];
  long msto = world->mend[pop];
  //species_fmt *s = NULL; 
	  
  // coalescence time if same is true
  if(same)
    time_coal = -r0 *  (proposal->param0[pop]) * 0.5;
  else
    time_coal = (double) HUGE;
  // migration
  time_mig_old = (double) HUGE;
  time_mig = (double) HUGE;
  for (i=msta;i<msto;i++)
    {
      if(shortcut(i, world, &j))
	continue;
      else
	{
	  r1  = LOG(UNIF_RANDUM());
	  temptime = -r1 /(world->data->geo[j] * world->param0[j]);
	  if (temptime < time_mig_old)
	    {
	      time_mig_old = temptime;
	      time_mig = temptime;
	      m2mm(j,world->numpop,&frompop,&t);
	      if (t!=pop)
		error("calc_time_per_line problem with migration to");	      
	    }
	}
    }
  // speciation
  if (world->has_speciation)
    {
      time_spec = time_to_speciation(world, pop, age, &events, &tos, &froms);
    }
  else
    time_spec = (double) HUGE;
  // adding fuzz so that we do not end up with identical times: 
  age += EPSILON * UNIF_RANDUM();
  if (time_coal < time_mig && time_coal < time_spec)
    {
      *time = age + time_coal;
      *from = pop;
      *event='c';
    }
  else if (time_mig < time_coal && time_mig < time_spec)
    {
      *time = age + time_mig;
      *from = frompop;
      *event = 'm';
    }
  else if (time_spec < time_coal && time_spec < time_mig)
    {
      *time = age + time_spec;
      *from = froms;
      *event = 'd';
    }
  else
    {
      //warning("no event was suitable for calc_per_line");
      *event = ' ';
      *from = -1;
      *time = (double) HUGE;
    }
}



// simulating beyond the last node in the tree
// replaces pre_population in mcmc1.c
//
int beyond_last_node(proposal_fmt* proposal, vtlist *tentry, long gte, long *slider)
{
  (void) slider;
  long pop1, pop2;
  char type1;
  char type2;
  long popto1, popto2;
  MYREAL age1, horizon;
  char event;
  MYREAL time;
  long to, from;
  long pmc1, pmc2;
  boolean dangling;
  MYREAL time1, time2;
  char event1, event2;
  long to1, to2, from1, from2;
  //some standard stuff copied from pre_population
  if (gte > 0)
    proposal->realtarget = tentry[gte - 1].eventnode;
  else
    proposal->realtarget = tentry[0].eventnode->next->back; //?????gte
  if (proposal->realtarget == proposal->oback)
    {
      proposal->realtarget = crawlback (proposal->osister)->back;
    }
  if (proposal->realtarget->type == 'm' || proposal->realtarget->type == 'd')
    {
      proposal->target = crawlback (proposal->realtarget->next);
      if (proposal->target == proposal->oback)
        {
	  proposal->target = proposal->osister;
        }
    }
  else
    {
      proposal->target = proposal->realtarget;
    }
  proposal->tsister = NULL;
  pop2 = proposal->realtarget->pop;
  pop1 = proposal->migr_table_counter > 0 ? proposal->migr_table[proposal->migr_table_counter -
								 1].from : proposal->origin->pop;
  //type1 = proposal->migr_table_counter > 0 ? proposal->migr_table[proposal->migr_table_counter -
  //								  1].event : proposal->origin->type; 
  age1 = MAX (proposal->realtarget->tyme,
	      proposal->migr_table_counter >
	      0 ? proposal->migr_table[proposal->migr_table_counter -
				       1].time : proposal->origin->tyme);
  horizon = MAX (proposal->oback->tyme, age1);

  //the residual tree is only one lineages after the last node, and the simulated line
  //can be in the same population, or different populations, it also could be that the
  // lineages have not crossed into the ancestor yet.
  // one population
  event = ' ';
  // pop1 is the actual population at proposal->time in the dangling lineage
  // pop2 is the last lineage in the residual tree 
  time = age1;
  proposal->time = time;
  // after the last node in the residual tree but before the last node in the 
  // original tree
  while (time < horizon)
    {
      time = proposal->time;
      pop2 = proposal->realtarget->pop;
      pop1 = proposal->migr_table_counter > 0 ? proposal->migr_table[proposal->migr_table_counter -
								     1].from : proposal->origin->pop;
      type1 = proposal->migr_table_counter > 0 ? proposal->migr_table[proposal->migr_table_counter -
								      1].event : proposal->origin->type;
      popto1 = proposal->migr_table_counter > 0 ? proposal->migr_table[proposal->migr_table_counter -
								       1].to : proposal->origin->actualpop;
      // calculate only for the dangling lineage, gets the timeinterval (not the full age)
      to = pop1;
      calc_time_per_line(proposal, type1, popto1, pop1 == pop2, &time, &event, &to, &from);
      // start proposal->time pop1, pop2
      // proposal time        to/from pop2
      if (time >= horizon)
        break;

      proposal->time = time;
      if (event == 'c')
	{
	  //printf("@%c: pop=%li time=%f\n", event, to , proposal->time); 
	  return 1;
	}
      if (!migrate (proposal, proposal->origin, event, from, to))
	{
	  return 0; // too many migrations! abort this trial
	}
    }
  proposal->time =  horizon;
  to1 = pop1;
  to2 = pop2;
  from1 = -1;
  from2 = -1;
  event = '@';
  // went past last node information in the old tree
  while (event != 'c')
    {
      time1 = time2 = proposal->time;
      pmc1 = proposal->migr_table_counter;
      pop1 = pmc1 > 0 ? proposal->migr_table[pmc1 - 1].from : proposal->origin->pop;
      type1= pmc1 > 0 ? proposal->migr_table[pmc1 - 1].event : proposal->origin->type;
      popto1= pmc1 > 0 ? proposal->migr_table[pmc1 - 1].to : proposal->origin->actualpop;
      pmc2 = proposal->migr_table_counter2;
      pop2 = pmc2 > 0 ? proposal->migr_table2[pmc2 - 1].from : proposal->realtarget->pop;
      type2= pmc2 > 0 ? proposal->migr_table2[pmc2 - 1].event : proposal->realtarget->type;
      popto2= pmc2 > 0 ? proposal->migr_table2[pmc2 - 1].to : proposal->realtarget->actualpop;
      to1 = pop1;
      to2 = pop2;
      calc_time_per_line(proposal, type1, popto1, pop1 == pop2, &time1, &event1, &to1, &from1);
      calc_time_per_line(proposal, type2, popto2, pop1 == pop2, &time2, &event2, &to2, &from2);
      //printf("%i> _both() interval line1: (%f,%c) line 2: (%f, %c) \n",myID,time1 - proposal->time, event1,time2 - proposal->time, event2);
      if (time1 < time2)
	{
	  dangling = TRUE;
	  time = time1;
	  event = event1;
	  to = to1;
	  from = from1;
	}
      else
	{
	  dangling = FALSE;
	  time = time2;
	  event = event2;
	  to = to2;
	  from = from2;
	}
      if (time >= 1e10)// HUGE, but HUge seems not to catch all DEBUG TODO PROBLEM
	return 0;
      proposal->time = time;
      if (event == 'c')
	{
	  //printf("@@%c: to=%li from=%li time=%f\n", event, to, from, proposal->time); 
	  return 1;
	}
      if(dangling)
	{
	  if (!migrate (proposal, proposal->origin, event, from, to))
	    {
	      return 0;  /* migration limit reached */
	    }
	  //pmc1 = proposal->migr_table_counter;
	  //pop1 = pmc1 > 0 ? proposal->migr_table[pmc1 - 1].from : proposal->origin->pop;
	}
      else
        {
	  if (!migrateb (proposal, proposal->realtarget, event, from, to))
	    {
	      return 0;  /* migration limit reached */
	    }
	      //pmc2 = proposal->migr_table_counter2;
	      //pop2 = pmc2 > 0 ? proposal->migr_table2[pmc2 - 1].from : proposal->realtarget->pop;
	}
    }
  return 0;
}

long get_species_record(world_fmt * world,long which)
{
  //long f,t;
  //m2mm(which,world->numpop,&f,&t);
  species_fmt * s = get_which_species_model(which, world->species_model, world->species_model_size);
  return s->id;
}


long propose_new_spec_mu(world_fmt * world, long which, boolean *is_mu, MYREAL *newparam)
{
  species_fmt *s;
  long t;
  MYREAL r;
  //long x = which-world->numpop2-world->bayes->mu;
  //long remainder = x % 2; // even are the mu and odd are the sigma [0,1,2,3,....]
  //x = (long) x/2; //[0,1,2,3,4,5,6,7,8,9]==>[0,0,1,1,2,2,3,3,4,4]

  s = get_which_species_model(which, world->species_model, world->species_model_size);
  r = RANDUM();
  //if (remainder==0) // pick mu or sigma
  if (which == s->paramindex_mu)
    {
      *is_mu = TRUE;
      do{
	*newparam = (*propose_new[which])(world->param0[which],which,world,&r);
	check_min_max_param(newparam,(double) s->min, (double) s->max);
      } while (*newparam <= (double) s->min && *newparam >= (double) s->max);
      s->mu = (double) (*newparam);
    }
  else
    {
      if (which != s->paramindex_sigma)
	error("wrong index for change if species std variable");
      *is_mu = FALSE;
      do{
      	*newparam = (*propose_new[which])(world->param0[which],which,world,&r);
      	check_min_max_param(newparam, (double) s->sigmamin,(double) s->sigmamax);
      } while (*newparam <= (double) s->sigmamin && *newparam >= (double) s->sigmamax);
      s->sigma = (double) (*newparam);
    }
  if (s->type == 't')
    {
      for (t=0;t<world->species_model_size;t++)
	{
	  species_fmt *other = &world->species_model[t];
	  if(s!=other && s->from == other->from && other->type == 't')
	    {
	      other->mu = s->mu;
	      other->sigma = s->sigma;
	      world->param0[other->paramindex_mu] = (double) s->mu;
	      world->param0[other->paramindex_sigma] = (double) s->sigma;
	    }
	}
    }

#ifdef DEBUG
  //printf("%i> propose_newmu_sigma: %li <- %li: specmu=%f (param[%li]=%f) sigma=%f mumin=%f mumax=%f\n",myID, s->to, s->from, s->mu,which,world->param0[which], s->sigma, s->min, s->max);
#endif
  return s->id;
}
void adjust_averagediv(world_fmt * world, long which, MYREAL *newparam)
{
  species_fmt *s;
  long t;
  s = get_which_species_model(which, world->species_model, world->species_model_size);
  if (which==s->paramindex_mu)
    {
      check_min_max_param(newparam,(double) s->min, (double) s->max);
      s->mu = *newparam;
    }
  else
    {
      check_min_max_param(newparam, (double) s->sigmamin,(double) s->sigmamax);
      s->sigma = *newparam;
    }
  if (s->type == 't')
    {
      for (t=0;t<world->species_model_size;t++)
	{
	  species_fmt *other = &world->species_model[t];
	  if(s->from == other->from && other->type == 't')
	    {
	      other->mu = s->mu;
	      other->sigma = s->sigma;
	      world->param0[other->paramindex_mu] = (double) s->mu;
	      world->param0[other->paramindex_sigma] = (double) s->sigma;
	    }
	}
    }
}

// used in 
MYREAL wait_event_species(world_fmt *world, vtlist*tli, MYREAL t0, MYREAL t1, boolean waitonly, MYREAL *eventprob)
{
  const long numpop = world->numpop;
  species_fmt * s;
  long ssize = world->species_model_size;
  long i,j;
  //long npp = world->numpop2 + world->bayes->mu;
  MYREAL specw = 0.0;
  long *lineages = tli->lineages;
  node *eventnode = tli->eventnode;
  long tlit = tli->to;
  MYREAL mu;
  MYREAL sigma;
  //*eventprob = 0.0;
  if (!world->has_speciation)
    return 0.0;
  for (i=0;i<numpop;i++)
    {
      for (j=0;j<numpop;j++)
	{
	  if(i==j)
	    continue;
	  s = get_fixed_species_model(j, i, world->species_model, ssize);
	  if (s==NULL)
	    continue;
	  else
	    {
	      mu = world->param0[s->paramindex_mu];
	      sigma = world->param0[s->paramindex_sigma];
	      //if (mu<EPSILON)
	      //  warning("%i> heat=%f speciation mu is zero\n",myID, world->heat);
	      specw += lineages[s->to] * (*log_prob_wait_speciate)(t0,t1,mu,sigma,s);
	      if (!waitonly && eventnode->type == 'd' && tlit==i)
		*eventprob = (*log_point_prob_speciate)(t1,mu,sigma,s)-log(lineages[s->to]);
	    }
	}
    }
  return specw;
}

// for skyline plots calculate waiting for speciation event
MYREAL wait_D(long pop, MYREAL t0, MYREAL t1, long *lineages, world_fmt *world)
{
  //warning("wait_D: needs fix to allow two or more ancestors");
  //const long numpop = world->numpop;
  species_fmt * s;
  long ssize = world->species_model_size;
  //long npp = world->numpop2 + world->bayes->mu;
  MYREAL specw = 0.0;
  MYREAL mu;
  MYREAL sigma;
  if (!world->has_speciation)
    return 0.0;
  s = get_species_model(pop, world->species_model, ssize);
  if (s==NULL)
    return 0.0;
  else
    {
      mu = world->param0[s->paramindex_mu];
      sigma = world->param0[s->paramindex_sigma];
      specw = lineages[s->to] * (*log_prob_wait_speciate)(t0,t1,mu,sigma,s);
    }
  return specw;
}


void set_first_speciestree(node *mrca, world_fmt *world)
{
    long ssize = world->species_model_size;
    long i, j;
    long from;
    long the_from;
    long to;
    boolean found=FALSE;
    char type;
    node *p,*q;
    node *tmp;
    the_from = world->species_model[0].from;
    for (i=0;i<ssize;i++)
    {
        from = world->species_model[i].from;
        for (j=0;j<ssize;j++)
        {
            to = world->species_model[i].to;
            if (from == to)
            {
                found = TRUE;
                break;
            }
            found = FALSE;
        }
        if (found)
            the_from = from;
    }
    mrca->pop = the_from;
    mrca->actualpop = the_from;
    p = mrca->next->back;
    q = mrca->next->next->back;
    type = (char) set_type(world, p->actualpop,mrca->actualpop, world->options->custm2, world->numpop);
    tmp = add_migration (world, p, type, mrca->actualpop, p->actualpop,
                         (MYREAL) RANDDOUBLE(0.0, mrca->tyme - p->tyme));
    mrca->next->back = tmp;
    tmp->back = mrca->next;
    type = set_type(world, q->actualpop,mrca->actualpop, world->options->custm2, world->numpop);
    tmp = add_migration (world, q, type, mrca->actualpop, q->actualpop,
                         (MYREAL) RANDDOUBLE(0.0, mrca->tyme - q->tyme));
    mrca->next->next->back = tmp;
    tmp->back = mrca->next->next;
}


// assembles the histogram from the sampled parameter values in db->data
// there is a similar function that reads from the bayesallfile mdimfile 
void construct_species_hist(species_fmt *db, world_fmt *world, long locus, 
			    long npa, long pai, long offset, long numbin, 
			    MYREAL *mini, MYREAL *maxi, float **results,
			    float *total, float *themean, float *thestd)
{
  (void) numbin;
  (void) npa;
  bayes_fmt *bayes = world->bayes;
  //char *custm2 = world->options->custm2;
  long      j, j0;
  long      i;
  long      bin;
  long      nb=0;

  MYREAL    delta;
  double    value;
  double    *values = db->data+offset;
  //long    floorindex;
  long size = db->size;
  long halfsize = size / 2;
  float *p = (float *) mycalloc(halfsize,sizeof(float));
  long z = 0;
  for(j0=0;j0<pai;j0++)
    {
      if(shortcut(j0,world,&j))//1229 || custm2[j0]=='d')
	{
	  continue;
	}
      //      if (j0>=world->numpop2)
      //	j = j0;
      nb += bayes->histogram[locus].bins[j];
    }
  //delta = bayes->deltahist[pai];
  for(i=0;i < size; i+=2)
    {
      value = values[i];
      p[z] =  (float) value;
      z++;
    }
  qsort(p, (size_t) halfsize, sizeof(float), floatcmp);
  for(i=0;i < halfsize; i++)
    {
      delta = bayes->deltahist[pai];
      value =  (double) p[i];
      *themean += (float) value;
      *thestd += (float) (value * value);
      if(value < mini[pai])
	{
	  warning("%i> TOO SMALL value=%f, mini[%li]=%f\n",myID,value,pai,mini[pai]);
	}
      if(value > maxi[pai])
	{
	  warning("%i> TOO LARGE value=%f, maxi[%li]=%f\n",myID,value,pai,maxi[pai]);
	  bin = bayes->histogram[locus].bins[pai] - 1;
	}
      else
	{
	  bin = (long) ((value - mini[pai])/delta);
	  if(bin<0)
	    bin=0;
	}
      if((bin) > bayes->histogram[locus].bins[pai])
	{
	  warning("%i> value not counted for histogram: bin=%li > histbins=%li\n", myID,bin, bayes->histogram[locus].bins[pai]);
	  *total +=1;
	} 
      (*results)[nb+bin] += (float) 1.;
      *total += 1;
    }
  //
  // adjusting the integral =1
  const float ttt = (float)(1.0/((double) *total));
  for(bin=0;bin < bayes->histogram[locus].bins[pai]; bin++)
    {
      (*results)[nb + bin] *= ttt;
    }
  myfree(p);
}


void construct_locusspecies_histogram(world_fmt *world, long locus, MYREAL *mini, MYREAL *maxi, float **results)
{
  bayes_fmt *bayes = world->bayes;
  species_fmt * adb = world->species_model;
  species_fmt * db;
  //long s;
  long i;
  long j0,j;
  long np = world->numpop2 + bayes->mu;
  long npa = np + 2 * world->species_model_size;
  long pa = np;
  long pai;
  float themean;
  float thestd;
  float total;
  long pai2;
  float themean2;
  float thestd2;
  float total2;
  long numbin = 0;
  for(j0=0; j0 < np; j0++)
    {
      if(shortcut(j0,world,&j))
	{
	  continue;
	}
      //if(j0 == world->numpop2)
      //j=world->numpop2;
      //1229else
      //1229{
	  //1229if(world->options->custm2[j]=='d')
	  //1229  continue;
      //1229}
      numbin += bayes->histogram[locus].bins[j];
    }    
  pai = pa;
  for (i=0;i<world->species_model_size;i++)
    {
      db = &(adb[i]);
      //s = db->size;
      themean = 0.0f;
      thestd = 0.0f;
      total = 0.0f;
      construct_species_hist(db,world,locus, npa, pai, 0,numbin, mini, maxi, 
			     results, &total,&themean,&thestd);
      world->bayes->histogram[locus].means[pai] = (double) (themean/total);
      world->bayes->histogram[locus].stds[pai]  = (double) (thestd / total);
      world->bayes->histtotal[locus*npa+pai] = (MYREAL) total;
      pai2= pai+1;
      themean2 = 0.0;
      thestd2 = 0.0;
      total2 = 0;
      numbin += bayes->histogram[locus].bins[pai];
      pai = pai2+1;

      construct_species_hist(db,world,locus, npa, pai2, 1, numbin, mini, maxi, 
			     results, &total2,&themean2,&thestd2);
      world->bayes->histogram[locus].means[pai2] = (double) (themean2/total2);
      world->bayes->histogram[locus].stds[pai2]  = (double) (thestd2 / total2);
      world->bayes->histtotal[locus*npa+pai2] = (MYREAL) total2;
      numbin += bayes->histogram[locus].bins[pai2];
    }
}

void get_div_from_tree(node *root, divtime_fmt **div, long *divtime_alloc, long *divtime_num)
{
  if (*divtime_num >= *divtime_alloc - 1)
    {
      *divtime_alloc += HUNDRED;
      *div = (divtime_fmt *) realloc(*div, sizeof(divtime_fmt) * (size_t) (*divtime_alloc));
    }
  if (root == NULL)
    return;
  switch(root->type)
    {
    case 't':
      return;
    case 'r':
      get_div_from_tree(root->next->back, div, divtime_alloc, divtime_num);
      break;
    case 'i':
      get_div_from_tree(root->next->back, div, divtime_alloc, divtime_num);
      get_div_from_tree(root->next->next->back, div, divtime_alloc, divtime_num);
      break;
    case 'm':
      get_div_from_tree(root->next->back, div, divtime_alloc, divtime_num);
      break;
    case 'd':
      (*div)[*divtime_num].age = root->tyme;
      (*div)[*divtime_num].from = root->pop;
      (*div)[*divtime_num].to = root->actualpop;
      (*divtime_num) += 1;
      get_div_from_tree(root->next->back, div, divtime_alloc, divtime_num);
      break;
    default:
      error("missing type information in get_div_from_tree");
      //break;
    }
}

void record_parameters(world_fmt *world)
{
  //extract ages of coalescences, migrations, divergences
  //
  if (world->cold)
    {
      world->divtime_num = 0;
      get_div_from_tree(world->root, &world->divtime, &world->divtime_alloc, &world->divtime_num);
      if (world->divtime_num>0)
	{
#ifdef MPI
	  send_divtime("%li\n",world->divtime_num);
#else
	  fprintf(world->divtimefile,"%li\n",world->divtime_num);
#endif
	  int i;
	  species_fmt *s = NULL; 	 
	  for (i=0; i<world->divtime_num; i++)
	    {
	      s = get_fixed_species_model(world->divtime[i].from, world->divtime[i].to, world->species_model, world->species_model_size); 
#ifdef MPI
	      send_divtime("%f %li %li %f %f\n",world->divtime[i].age,world->divtime[i].from,world->divtime[i].to, s->mu, s->sigma);
#else
	      fprintf(world->divtimefile,"%f %li %li %f %f\n",world->divtime[i].age,world->divtime[i].from,world->divtime[i].to, (double) s->mu, (double) s->sigma);
#endif

	    }
	}
    }
  //save to file or send to master
  // tree# age thetas1,2,3 Ms21 12 ...  Ds21 12.. Stds21 12 ...
  // - better: tree# age paramID value
  
}

/*void adjust_priorsfor species
@@@@@@
      if(world->options->custm2[i]=='d')
	{
	  long f;
	  long t;
	  m2mm(i,world->numpop,&f,&t);
	  s = get_species_model(t, world->species_model, world->species_model_size);
	  bayes->minparam[i] = s->min;
	  bayes->maxparam[i] = s->max;
	  bayes->meanparam[i] = s->mu;
	  if(options->prior->alpha[i]
	  bayes->deltahist[z++] = (s->max - s->min)/ options->bayespriorm[world->numpop].bins;//use M number bins TODO FIX DEBUG
	  bayes->deltahist[z++] = (s->sigmamax - s->sigmamin)/ options->bayespriorm[world->numpop].bins;//use M number bins TODO FIX DEBUG
	}
*/
#endif /* #ifndef SPECIATETESTER */
// testing particular parts of the speciate.c files

#ifdef SPECIATETESTER
//clang -g speciate.c -DSPECIATETESTER -DNOJPEG -DNOPNG  -I../SFMT-1.4.1 -I../haru -DPRETTY -DLETTERPAPER -o tester
// call like this tester mu sigma | mean

double inverse_cumstd_normal(double p)
{
  // Coefficients in rational approximations
  const double a[] = {-3.969683028665376e+01,  2.209460984245205e+02,
		      -2.759285104469687e+02,  1.383577518672690e+02,
		      -3.066479806614716e+01,  2.506628277459239e+00};

  const double b[] = {-5.447609879822406e+01,  1.615858368580409e+02,
		      -1.556989798598866e+02,  6.680131188771972e+01,
		      -1.328068155288572e+01 };

  const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01,
		      -2.400758277161838e+00, -2.549732539343734e+00,
		      4.374664141464968e+00,  2.938163982698783e+00};

  const double d[] = {7.784695709041462e-03, 3.224671290700398e-01,
		     2.445134137142996e+00,  3.754408661907416e+00};

  // Define break-points.
  const double plow  = 0.02425;
  const double phigh = 1 - plow;

  // Rational approximation for lower region:
  if ( p < plow ) {
    double q  = sqrt(-2.*log(p));
    return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
      ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }

  // Rational approximation for upper region:
  if ( phigh < p ) {
    double q  = sqrt(-2.*log(1-p));
    return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
      ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }

  // Rational approximation for central region:
  double q = p - 0.5;
  double r = q*q;
  return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
    (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
}

// z = erfinv(x) where x is between -1 and 1 returns a value z so that x=erf(z)
double erfinv(double x)
{
  if (x < -1.0 || x > 1.0)
    {
      //warning("erfinv received illegal value");
      return (double) -HUGE;
    }
      return inverse_cumstd_normal((x+1.0)/2.0)/SQRT2;
    }

double randum()
{
  return ((double) random())/RAND_MAX;
}

long random_integer(long start, long stop)
{
  return start + (long) (random() * (stop-start+1));
}

// weibull
MYREAL time_to_speciate_weibull(world_fmt *world, long pop, MYREAL t0, char *event, long *to, long *from)
{
  //  f[a_, b_, t0_] := 
  //        Module[{value = 0, count = 0}, 
  //           While[value < t0 && count++ < 100, 
  //                 value = b (-Log[RandomReal[]])^(1/a)]; Return[value]]
  long ssize = world->species_model_size;
  MYREAL mu;
  MYREAL sigma;
  MYREAL interval;
  double tmax;
  double r;
  double xx, yy;
  double t1 = (double) HUGE;
  double tnew;
  boolean changed=FALSE;
  species_fmt *s = NULL;
  
  for (from=0; from < world->numpop; from++)
    {
      if (from == pop)
	continue;
      s = get_fixed_species_model(from, pop, world->species_model, ssize);
      if(s==NULL)
	continue;
      lambda = world->param0[s->paramindex_mu];
      k = world->param0[s->paramindex_sigma];
      tmax = (double) s->max;
      r = log(UNIF_RANDUM());
      xx = log(t0/lambda) * k;
      yy = exp(xx) - r;
      tnew = lambda * pow(yy,1.0/k);
      if (tnew < t1)
	{
	  t1 = tnew;
	  changed=TRUE;
	}
    }
  if(changed==FALSE)
    return (double) HUGE;
  if (t1<t0)
    t1 = t0 + EPSILON;
  *event = 'd';
  *from = s->from;
  *to = s->to;//should be the same as pop
  if (t0 <= t1 && t1 < tmax)
    {
      interval =  t1-t0;  
      return interval;
    }
  else 
    return (double) HUGE;
}




int main(int argc, char **argv)
{
  char event = 'd';
  long to = 1;
  long from = 0;
  double age = 0.0;
  long pop = 1;
  world_fmt *world = calloc(1,sizeof(world_fmt));
  long ssize = 1;
  world->species_model_size = ssize;
  long numpop=2;
  world->param0 = calloc(6,sizeof(double));
  world->param0[0] = 0.01;
  world->param0[1] = 0.01;
  world->param0[2] = 100;
  world->param0[3] = 100;
  world->param0[4] = atof(argv[1]);
  world->param0[5] = atof(argv[2]);
  world->species_model = (species_fmt *) mycalloc((ssize+1),sizeof(species_fmt));
  world->species_model[0].type='t';
  world->species_model[0].min = 0.0;
  world->species_model[0].max=10.0;
  world->species_model[0].sigmamin = 0.0;
  world->species_model[0].sigmamax=10.0;
  world->species_model[0].from = from;
  world->species_model[0].to = to;
  world->species_model[0].paramindex_mu=4;
  world->species_model[0].paramindex_sigma=5;
  world->species_model[0].size = 0;
  world->species_model[0].allocsize = 10;
  world->species_model[0].data = calloc(10,sizeof(double));
  for (long i=0; i<1000; i++)
    fprintf(stdout,"%f\n",time_to_speciate_normal(world, pop, age, &event, &to, &from));
  return 0;
}

#endif



