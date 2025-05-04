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
/*! \file priors.c 

priors.c supplies prior distributions
*/
#include "migration.h"
#include "random.h"
#include "definitions.h"
#include "tools.h" 
#include "sighandler.h"

#ifdef HAVE_LGAMMA
#define LGAMMA lgamma
#else
#define LGAMMA mylgamma
#endif

extern int myID;
void is_priorkind(prior_fmt *p, char *priorkind);
MYREAL cdf_gamma(MYREAL a, MYREAL b, MYREAL x);
MYREAL trunc_gamma_rand(MYREAL alpha, MYREAL beta, MYREAL lower, MYREAL upper);
MYREAL normal_rand(MYREAL mean, MYREAL std);
double trunc_random_normal2(double mi, double ma, double mean, double std);

MYREAL propose_uni_newparam (MYREAL param, long which, world_fmt *world, MYREAL *r);
MYREAL propose_exp_newparam (MYREAL param, long which, world_fmt *world, MYREAL *r);
MYREAL propose_expb_newparam (MYREAL param,long which, world_fmt *world, MYREAL *r);
MYREAL propose_mult_newparam (MYREAL param,long which, world_fmt *world, MYREAL *r);
MYREAL propose_normal_newparam (MYREAL param,long which, world_fmt *world, MYREAL *r);
MYREAL propose_gamma_newparam (MYREAL param, long which, world_fmt *world, MYREAL *r);

MYREAL log_prior_ratio_uni  (MYREAL newparam, MYREAL oldparam, bayes_fmt *bayes, long which);
MYREAL log_prior_ratio_exp  (MYREAL newparam, MYREAL oldparam, bayes_fmt *bayes, long which);
MYREAL log_prior_ratio_wexp (MYREAL newparam, MYREAL oldparam, bayes_fmt *bayes, long which);
MYREAL log_prior_ratio_mult (MYREAL newparam, MYREAL oldparam, bayes_fmt *bayes, long which);
MYREAL log_prior_ratio_normal (MYREAL newparam, MYREAL oldparam, bayes_fmt *bayes, long which);
MYREAL log_prior_ratio_gamma(MYREAL newparam, MYREAL oldparam, bayes_fmt * bayes, long which);


MYREAL log_prior_uni  (world_fmt * world, long numparam);//for heating
MYREAL log_prior_exp  (world_fmt * world, long numparam);//for heating
MYREAL log_prior_wexp (world_fmt * world, long numparam);//for heating
MYREAL log_prior_mult (world_fmt * world, long numparam);//for heating
MYREAL log_prior_normal (world_fmt * world, long numparam);//for heating
MYREAL log_prior_gamma(world_fmt *world, long numparam);

MYREAL log_prior_uni1  (world_fmt * world, long numparam, MYREAL);
MYREAL log_prior_exp1  (world_fmt * world, long numparam, MYREAL);
MYREAL log_prior_wexp1 (world_fmt * world, long numparam, MYREAL);
MYREAL log_prior_mult1 (world_fmt * world, long numparam, MYREAL);
MYREAL log_prior_normal1 (world_fmt * world, long numparam, MYREAL);
MYREAL log_prior_gamma1(world_fmt *world, long numparam, MYREAL val);


MYREAL hastings_ratio_uni(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam);
MYREAL hastings_ratio_exp(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam);
MYREAL hastings_ratio_expb(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam);
MYREAL hastings_ratio_mult(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam);
MYREAL hastings_ratio_normal(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam);
MYREAL hastings_ratio_gamma(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam);

void init_hyperpriorrecord(hyper_fmt ** hyperp,long numparam);
void hyper_gamma_record(MYREAL alpha, MYREAL beta, hyper_fmt * hg);
void hyper_report(world_fmt *world, bayes_fmt *bayes);
MYREAL propose_uniform(MYREAL param, MYREAL minparam, MYREAL maxparam, MYREAL *r, MYREAL delta);
void set_option_prior(prior_fmt **p, int type, MYREAL mini, MYREAL maxi, MYREAL mean, long bins);
prior_fmt * copy_option_prior(prior_fmt *pmodel, option_fmt *options, MYREAL ratemin);
void copy_prior(prior_fmt *target, prior_fmt *source);
prior_fmt * find_prior_menu(long from, long to, long priortype, option_fmt * options);
void find_prior(long from, long to, long priortype, option_fmt * options, prior_fmt *result);
prior_fmt * set_default_prior(option_fmt *options, long from, long to, long priortype);
long set_all_default_priors(option_fmt *options);
prior_fmt * insert_prior(prior_fmt *p, prior_fmt **plist,long *z,option_fmt * options);
prior_fmt * insert_prior_old(option_fmt *options, long index, prior_fmt *p, long type);
void check_bayes_priors(option_fmt *options, data_fmt *data, world_fmt *world);
prior_fmt * get_prior_list(prior_fmt **list, int type);
MYREAL find_beta_truncgamma(MYREAL mean, MYREAL alpha, MYREAL lower, MYREAL upper);
long * which_type(long which, world_fmt * world, long *mainindex);
/*
void setup_type_matrix(world_fmt *world)
{
  long typecounts = 2 + world->species_model_size * 2 + world->bayes->mu;
  world->bayes->rules->indicesnum = (long *) mycalloc(typecounts,sizeof(long));
  world->bayes->rules->indices = (long **) mycalloc(typecounts,sizeof(long * ));
  long i;
  world->bayes->rules->indices[0] = (long *) mycalloc(world->numpop,sizeof(long));
  world->bayes->rules->indicesnum[0] = world->numpop;
  for (i=0;i<world->numpop;i++)
    {
      if (world->options->custm2[i]=='*')
	world->bayes->rules->indices[0][i] = i;
      else
	world->bayes->rules->indices[0][i] = -1;
    }
  if (world->numpop>1)
    {
      world->bayes->rules->indicesnum[1] = world->numpop2 - world->numpop;
      for (i=world->numpop; i < world->numpop2;i++)
	{
	  if (world->options->custm2[i]=='*')
	    world->bayes->rules->indices[1][i] = i;
	  else
	    world->bayes->rules->indices[1][i] = -1;
	}
      long start = world->numpop2+world->bayes->mu;
      world->bayes->rules->indicesnum[3] = world->species_model_size;
      world->bayes->rules->indicesnum[4] = world->species_model_size;
      for (i=start;i < start + world->species_model_size*2; i+=2)
	{
	  world->bayes->rules->indices[3][i] = i;
	  world->bayes->rules->indices[4][i+1] = i+1;
	}
    }
  if (world->bayes->mu)
    {
      world->bayes->rules->indicesnum[2] = 1;
      world->bayes->rules->indices[2][i] = world->numpop2;
    }
  else
    {
      world->bayes->rules->indicesnum[2] = 1;
      world->bayes->rules->indices[2][i] = -1;
    }
}
*/

long * which_type(long which, world_fmt * world, long *mainindex)
{
  long ** indices = world->bayes->rules->indices;
  if (which < world->numpop)
    {
      *mainindex = 0;
      return indices[0]; //world->bayes->rules->indices[0];
    }
  else if (which < world->numpop2)
    {
      *mainindex = 1;
      return indices[1];
    }
  else if (which == world->numpop2 && world->bayes->mu)
    {
      *mainindex = 2;
      return indices[2];
    }
  else if (which >= world->numpop2 && which % 2 == 0)
    {
      *mainindex = 3;
      return indices[3];
    }
  else
    {
      *mainindex = 4;
      return indices[4];
    }
}

// check whether the value drawn from the prior fits into the rule set
// specified in bayes->rules->custm2rules[i] with i=0 .. bayes->rules->rulesnum
// currently accepts only a single rule
/* not used currently
static boolean rule_check(MYREAL value, long which, world_fmt *world)
{
  char alfabet[] = "abcdefghijklmnopqrstuvwz";
  if (world->bayes->rules == NULL)
    return TRUE;
  char ** rules = world->bayes->rules->custm2rules;
  //long numrules = world->bayes->rules->rulesnum;
  char myrule = rules[0][which];
  long i=0;
  if(strchr("SsMm0c", myrule))
    return TRUE;
  
  long *indices;
  long mainindex;
  indices = which_type(which, world,&mainindex);
  for (i=0; i< world->bayes->rules->indicesnum[mainindex]; i++)
    {
      long ind = indices[i];
      if (myrule < alfabet[ind] && value >= world->param0[ind])
	{
	  myfree(indices);
	  return FALSE;
	}
    }
  return TRUE;
}
*/


///
/// Normal distribution
/// this generates only one random number, but could be easily changed to supply two numbers
MYREAL normal_rand(MYREAL mean, MYREAL std)
{
  MYREAL u = UNIF_RANDUM();
  MYREAL v = UNIF_RANDUM();
  return (sqrt(-2 * log(u)) * cos(TWOPI*v)) * std + mean;
  //return sqrt(-2 * log(u)) * sin(TWOPI*v);
}

// for prior_fmt
static float trunc_random_normal(float *p)
{
  float r = (float) normal_rand((double) p[2],(double) p[3]);
  while (r < p[0] || r > p[1])
    r = (float) normal_rand((double) p[2], (double) p[3]);
  return r;
}

double trunc_random_normal2(double mi, double ma, double mean, double std)
{
  (void) ma;
  long count=0;
  double r = normal_rand(mean,std);
  while (count++ < 100 && (r < mi)) // was || r > ma
    r = normal_rand(mean, std);
  if(count>=100)
    {
      //warning("random normal failed: %f, %f, %f\n",mi,r, ma);
      return (double) -1.0;
    }
  else
    return r;
}

// for prior_fmt
static float cdf_normal(float x, float *p)
{
  return (float) (0.5 * (1.0 + erf(((double) (x-p[2]))/(sqrt(2.0) * (double) p[3]))));
}

static float trunc_cdf_normal(float x, float *p)
{
  float aa =  cdf_normal((p[0]-p[2])/p[3],p);
  float Z = cdf_normal((p[1]-p[2])/p[3],p) - aa;
  return (cdf_normal(x,p) - aa)/ Z;
}

// uniform for prior_fmt
static float trunc_random_uni(float *p)
{
  MYREAL u = UNIF_RANDUM();
  return (float) ((double) p[0] + (double) (p[1]-p[0]) * u);
}

// uniform for prior_fmt
static float trunc_cdf_uni(float x, float *p)
{
  //x==(y-p[0])/(p[1]-p[0])
  return x * (p[1] - p[0]) + p[0];
}

// uniform for prior_fmt
static float trunc_random_exp(float *p)
{
  MYREAL u = UNIF_RANDUM();
  float r = (float) (-log(u) * (double) p[2]);
    while(r < p[0] || r > p[1])
      {
	u = UNIF_RANDUM();
	r = (float) (-log(u) * (double) p[2]);
      }
    return r;
}

// uniform for prior_fmt
static float trunc_cdf_exp(float x, float *p)
{
  double lambda = 1.0/ (double) p[2];
  //float cdfx = 1.0 - exp(-lambda * x);
  double cdfa =  (1.0 - exp(-lambda * (double) p[0]));
  double cdfb =  (1.0 - exp(-lambda * (double) p[1]));
  //return cdfx/(cdfb-cdfa);
  return (float) (-log((double) -x * (cdfb-cdfa) + 1.0)/lambda);
}

static float trunc_random_gamma(float *p)
{
  MYREAL alpha = (double) p[2];
  MYREAL beta = (double) p[3];
  return  (float) trunc_gamma_rand(alpha,beta,(double) p[0],(double) p[1]);
}

static double bisection(int n, double a, double b, double tol, float p[], float x) {
  MYREAL alpha = (double) p[2];
  MYREAL beta = (double) p[3];
  if (beta < SMALLEPSILON)
    {
      beta = SMALLEPSILON;
    }
  MYREAL Z = cdf_gamma(alpha, beta, (double) p[1]) -cdf_gamma(alpha, beta, (double) p[0]);
  float fa = (float) (cdf_gamma(alpha, beta, a)/Z - (double) x);
  //float fb = (float) cdf_gamma(alpha, beta, b)/Z - x;
  double     c = (a + b) / 2.0;
  float fc;
  int count = 0;
  
  while (((b - a) > tol) && (count++ < n)) 
    {
      c = (a + b) / 2.0;
      fc = (float) (cdf_gamma(alpha, beta, c)/Z - (double) x);
      if (fc < 0.0f) 
	{
	  if (fa < 0.0f) 
	    {
	      a = c;
	      fa = fc;
	    } 
	  else 
	    {
	      b = c;
	      //fb = fc;
	    }
	} 
      else 
	{
	  if (fa >= 0.0f) 
	    {
	      a = c;
	      fa = fc;
	    } 
	  else 
	    {
	      b = c;
	      //fb = fc;
	    }
	}
    }
  return c;
}


static float trunc_cdf_gamma(float x, float *p)
{
  return  (float) (bisection(100, (double) p[0], (double) p[1], 0.0001, p, x)); 
}


static MYREAL logpdf_gamma(MYREAL a, MYREAL b, MYREAL x)
{
  if (a>0.0 && b> 0.0)
    return -(x/b) - a*log(b) + (-1 + a)*log(x) - LGAMMA(a);
  else 
    return (double) -HUGE;
}

MYREAL cdf_gamma(MYREAL a, MYREAL b, MYREAL x)
{
  //Gamma[a, 0, xmax/b]/Gamma[a]
  return incompletegamma(x/b,a);
}

static MYREAL logpdf_truncgamma(MYREAL a, MYREAL b, MYREAL  xmin, MYREAL xmax, MYREAL x)
{
  if((x > xmin) && (x <= xmax))
    return logpdf_gamma(a,b,x)-log(cdf_gamma(a,b,xmax)-cdf_gamma(a,b,xmin));
  else 
    return (double) -HUGE;
}


static MYREAL mean_truncgamma(MYREAL alpha, MYREAL beta, MYREAL lower, MYREAL upper)
{
  MYREAL n1, n2, d1, d2;
  MYREAL m=0;
  MYREAL nom;
  MYREAL denom;
  //1>  beta*(Gamma(1+alpha,lower/beta,upper/beta)/Gamma(alpha,lower/beta,upper/beta)
  //2> (beta*(Gamma(1 + alpha,upper/beta) - Gamma(1 + alpha,lower/beta)))/
  //   (Gamma(alpha,upper/beta) - Gamma(alpha,lower/beta));
  // 1 and 2 should be the same
  // 
  n1 = logincompletegamma(upper/beta,1+alpha);
  n2 =  logincompletegamma(lower/beta,1+alpha);
  d1 = logincompletegamma(upper/beta,alpha);
  d2 =  logincompletegamma(lower/beta,alpha);
  nom = (1.0 - exp(n2-n1));
  denom = exp(d1-n1) - exp(d2-n1);
  if ((denom == 0.0) && (nom == 0.0))
    return beta * alpha;
  else
    m = beta * alpha * nom / denom;
  return m;
}


MYREAL find_beta_truncgamma(MYREAL mean, MYREAL alpha, MYREAL lower, MYREAL upper)
{
  /* from wikipedia:
     INPUT: Function f, endpoint values a, b, tolerance TOL, maximum iterations NMAX
     CONDITIONS: a < b, either f(a) < 0 and f(b) > 0 or f(a) > 0 and f(b) < 0
     OUTPUT: value which differs from a root of f(x)=0 by less than TOL
     N ← 1
     While N ≤ NMAX { limit iterations to prevent infinite loop
     c ← (a + b)/2 new midpoint
     If (f(c) = 0 or (b – a)/2 < TOL then { solution found
     Output(c)
     Stop
     }
     N ← N + 1 increment step counter
     If sign(f(c)) = sign(f(a)) then a ← c else b ← c new interval
     }
     Output("Method failed.") max number of steps exceeded  
  */
  MYREAL tolerance = EPSILON;
  MYREAL nmax = 1000;
  long n=1;
  MYREAL beta = (upper+lower)* 0.5; /*=c in readme*/
  MYREAL fbeta;
  MYREAL fa;
  MYREAL a = 0.0;  //was lower
  MYREAL b = upper*2.0; // was upper
  //return mean/alpha;
  if (a==0.0)
    a += SMALLEPSILON;
  while (n<= nmax)
    {
      beta = (a+b)/2.0;
      fbeta = mean_truncgamma(alpha,beta,0.0, (double) HUGE) - mean; // see fa   
      if(fbeta  == 0.0 || (b-a)/2. < tolerance)
	{
	  //printf("beta=%f\n",beta);
	  return beta;
	}
      n += 1;
      fa = mean_truncgamma(alpha,a,0.0,(double) HUGE) - mean; //lower,upper replaced 
      // with 0 and large value to allow better mean recovery with tight lower,upper bounds    
      if(((fbeta > 0.0) && (fa < 0.0)) || ((fbeta<0.0) && (fa>0.0)))
	b = beta;
      else
	a = beta;
    }
  warning("limit reached without achieving necessary accuracy");
  //printf("beta=%f --------------\n",beta);
  return beta;
}


/*void use_hyperprior_mean(world_fmt * world)
{
  long numpop2 = world->numpop2;
  long numparam = numpop2 + bayes->mu+ world->species_model_size * 2;
  bayes_fmt * bayes = world->bayes;
  double hyperp = bayes->hyperp;
  for (i=0; i < numparam; i++)
    {
      
    }
    }*/



void init_hyperpriorrecord(hyper_fmt ** hyperp,long numparam)
{
  long i;
  (*hyperp) = (hyper_fmt *) mycalloc(numparam, sizeof(hyper_fmt));
  for(i=0;i<numparam;i++)
    {
      onepass_mean_std_start(&(*hyperp)[i].alpha, &(*hyperp)[i].alphastd, &(*hyperp)[i].alphan);
      onepass_mean_std_start(&(*hyperp)[i].mean, &(*hyperp)[i].meanstd, &(*hyperp)[i].meann);
    }
}

static void hyper_normal(MYREAL *mean, MYREAL *std, MYREAL origmean, MYREAL origstd, MYREAL lower, MYREAL upper, MYREAL factor)
{
  // use a normal to draw new means, alpha, and then recalculate beta
  // using the original prior mean as the mean and std of the normal
  MYREAL nmean = normal_rand(origmean,factor*origmean);
  MYREAL nstd = normal_rand(origstd,factor*origstd);
  while (nstd <= EPSILON)
    {
      nstd = normal_rand(origstd,origstd);
    }
  while ((nmean <= lower) || (nmean > upper))
    {
      nmean = normal_rand(origmean,origmean);
    }
  *mean = nmean;
  *std = nstd;
}

static void hyper_gamma(MYREAL *alpha, MYREAL *beta, MYREAL origmean, MYREAL origalpha, MYREAL lower, MYREAL upper, MYREAL factormean, MYREAL factoralpha)
{
  // use a normal to draw new means, alpha, and then recalculate beta
  // using the original prior mean as the mean and std of the normal
  MYREAL nmean = normal_rand(origmean,factormean*origmean);
  MYREAL nalpha = origalpha;
  if(factoralpha>0.0)
    {
      nalpha = normal_rand(origalpha,factoralpha*origalpha);
      while (nalpha < lower || nalpha > upper)
	{
	  nalpha = normal_rand(origalpha, factoralpha*origalpha);
	    }
    }
  while (nmean < lower || nmean > upper)
    {
      nmean = normal_rand(origmean,factormean*origmean);
    }
  MYREAL nbeta = find_beta_truncgamma(nmean,nalpha,lower,upper);
  //printf("%i> %c nmean=%f origmean=%f alpha=%f, beta=%f\n", myID, nmean < origmean ? '-' : '+', nmean, origmean, nalpha, nbeta);
  *alpha = nalpha;
  *beta = nbeta;

}

void hyper_gamma_record(MYREAL alpha, MYREAL beta, hyper_fmt * hg)
{
  onepass_mean_std_calc(&hg->alpha,&hg->alphastd,&hg->alphan,alpha);
  onepass_mean_std_calc(&hg->mean,&hg->meanstd,&hg->meann,alpha*beta);
}

void hyper_report(world_fmt *world, bayes_fmt *bayes)
{
  hyper_fmt *hyperp = bayes->hyperp;
  long i;
  fprintf(world->outfile, "Hyperpriors are active\n");
  fprintf(world->outfile, "----------------------\n");
  fprintf(world->outfile, "Hyper prior is distributed: Normal(priormean,%f*std)\n",bayes->hyperfactormean);  
  fprintf(world->outfile, "                            Normal(prioralpha,%f*std)\n",bayes->hyperfactoralpha);
  fprintf(world->outfile, "Parameter     PriorMean Std  PriorAlpha   Std     N\n");
  for (i=0;i<bayes->numparams;i++)
    {
      fprintf(world->outfile,"%li        %f  %f  %f %f %li\n", 
	      i+1, hyperp->mean, hyperp->meanstd, 
	      hyperp->alpha, hyperp->alphastd, 
	      hyperp->alphan);
    }
}

///
/// Gamma prior retrieve a new value from a truncated gamma between lower and upper
/// currently does not use the old parameter, but this 
MYREAL propose_gamma_newparam (MYREAL param, long which, world_fmt *world, MYREAL *r)
{
  (void) param;
  (void) r;
  MYREAL rr;
  bayes_fmt *bayes = world->bayes;
  // test with delta
  // const MYREAL delta = bayes->delta[which];
  // const MYREAL l = param - delta;
  // const MYREAL u = param + delta;
  // const MYREAL lower = l < bayes->minparam[which] ? bayes->minparam[which] : l;
  // const MYREAL upper = u > bayes->maxparam[which] ? bayes->maxparam[which] : u;
  const MYREAL lower = bayes->minparam[which];
  const MYREAL upper = bayes->maxparam[which]; 
  MYREAL alpha = bayes->alphaparam[which];
  MYREAL beta = bayes->betaparam[which];
  if (bayes->hyperprior && bayes->hyperp[which].count++ % bayes->hyperinterval == 0)
    {
      const MYREAL mean = bayes->meanparam[which];
      hyper_gamma(&alpha, &beta, mean, bayes->alphaorigparam[which], lower, upper,bayes->hyperfactormean,bayes->hyperfactoralpha);
      bayes->alphaparam[which] = alpha;
      bayes->betaparam[which] = beta;
      hyper_gamma_record(alpha,beta, &bayes->hyperp[which]);
    }
  //rr = trunc_gamma_rand(alpha, beta, param-delta, param+delta);
  rr = trunc_gamma_rand(alpha, beta, lower, upper);
  //while (!rule_check(rr,which, world))
  //    rr = trunc_gamma_rand(alpha, beta, lower, upper);
  return rr;
}


///
/// Hastings ratio calculator for gamma distribution
/// P(new -> old)    P(old)
/// ------------- = -------
/// P(old -> new)    P(new)
/// cancels with log_prior_gamma -> 0.0
MYREAL hastings_ratio_gamma(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam)
{
  (void) oldparam;
  (void) delta;
  (void) r;
  if((newparam > bayes->maxparam[whichparam]) || (newparam < bayes->minparam[whichparam]))
    return (double) -HUGE;
  else
      return 0.;
    }


///
/// Log Prior gamma distribution ratios between old and new parameter:
/// cancels with hastings ratio
MYREAL log_prior_ratio_gamma(MYREAL newparam, MYREAL oldparam, bayes_fmt * bayes, long which)
{
  (void) oldparam;
  if((newparam > bayes->maxparam[which]) || (newparam < bayes->minparam[which]))
    return (double) -HUGE;
  else
    return 0.;
}

///
/// Gamma prior distribution for theta or migration rate used in heating acceptance
/// uses pdf_truncgamma(a,b,min,max,x)
MYREAL log_prior_gamma(world_fmt *world, long numparam)
{
  //  long frompop;
  //long topop;
  MYREAL p0;
  long numpop = world->numpop;
  long start = ((numparam <= numpop || numpop==1) ? 0 : numpop);
  long stop = ((start == numpop) ? world->numpop2 : numpop);
  long i;
  MYREAL * param0 = world->param0;
  bayes_fmt * bayes = world->bayes;
  MYREAL a;
  MYREAL b;
  MYREAL val=0.0;
  for(i = start; i < stop; i++)
    {
      if(!strchr("0c", world->options->custm2[i]))
	{
	  //@@if(i>=numpop)//@@ && !world->options->usem)
	  //@@  {
#ifndef PRIORTEST 
	      //@@     m2mm(i,numpop,&frompop,&topop);
#endif
	      //@@p0 = param0[i] * param0[topop];
	  //@@ }
	  //@@else
	  //@@{
	      p0 = param0[i];
	      //@@  }

	  if((p0 > bayes->maxparam[i]) || (p0 < bayes->minparam[i]))
	    return (double) -HUGE;
	  else
	    {
	      a = bayes->alphaparam[i];
	      b = bayes->betaparam[i];
	      val += logpdf_truncgamma(a,b,bayes->minparam[i],bayes->maxparam[i],p0);
	      //old!! val += -p0 * ib + a * log(ib) + (a - 1.) * log(p0) - LGAMMA(a);
	    }
	}
    }
  return val;
}

///
/// 
MYREAL log_prior_gamma1(world_fmt *world, long numparam, MYREAL val)
{
  bayes_fmt * bayes = world->bayes;
  MYREAL retval;
  MYREAL a = bayes->alphaparam[numparam];
  MYREAL b = bayes->betaparam[numparam]; 
  retval =  logpdf_truncgamma(a,b,bayes->minparam[numparam],bayes->maxparam[numparam],val);
  return retval;
}

///////////////////////////////////////////////////////////////////////
// uniform prior
///
/// Uniform flat prior, coding based on Rasmus Nielsen, veryfied with mathematica
/// the correlation among values is dependent on delta
MYREAL propose_uniform(MYREAL param, MYREAL minparam, MYREAL maxparam, MYREAL *r, MYREAL delta)
{
  MYREAL thisdelta;
  MYREAL rr;
  MYREAL np=0.0;
  *r = UNIF_RANDUM();
  rr = (*r) - 0.5;
  thisdelta = rr * delta;
  np = param + thisdelta;
  if (np > maxparam)
    {
      np =  2. * maxparam - np;
    }
  else
    {
      if (np < minparam)
	{
	  np = 2. * minparam - np;
	}
    }
  return np;
}

MYREAL
propose_uni_newparam (MYREAL param, long which, world_fmt *world, MYREAL *r)
{
  bayes_fmt *bayes = world->bayes;
  MYREAL minparam = bayes->minparam[which];
  MYREAL delta = bayes->delta[which];
  MYREAL maxparam = bayes->maxparam[which];
  MYREAL np = propose_uniform(param, minparam, maxparam, r, delta);
  //while (!rule_check(np,which, world))
  //  {
  //    np = propose_uniform(param, minparam, maxparam, r, delta);
  //  }
  return np;
}



///
/// Hastings ratio calculator for uniform distribution
MYREAL hastings_ratio_uni(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam)
{
  (void) newparam;
  (void) oldparam;
  (void) delta;
  (void) r;
  (void) bayes;
  (void) whichparam;

    //  return  (newparam - oldparam)/bayes->priormean[whichparam];
    return 0.;
}


///////////////////////////////////////////////////////////////////////
// exponential prior
/// \brief Exponential prior distribution
/// exponential prior distribution using the gamma prior with alpha=1
///
MYREAL
propose_exp_newparam (MYREAL param, long which, world_fmt *world,  MYREAL * r)
{
  MYREAL rr = propose_gamma_newparam (param, which, world, r);
  //while (!rule_check(rr,which, world))
  //  {
  //    *r = UNIF_RANDUM();
  //    rr = propose_gamma_newparam (param, which, world, r);
  //  }
  return rr;
}

///
/// Hastings ratio calculator for exponential distribution
MYREAL hastings_ratio_exp(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam)
{
  (void) newparam;
  (void) oldparam;
  (void) delta;
  (void) r;
  (void) bayes;
  (void) whichparam;
    return 0.;
}


///
/// uses an exponential prior with boundaries minparm and maxparm and uses also a window around the old parameter
/// the window is of arbitrary size
/// (see propose_exp_newparam() for details)
MYREAL
propose_expb_newparam (MYREAL param, long which, world_fmt *world, MYREAL *r)
{
  MYREAL rr = propose_gamma_newparam (param, which, world, r);
  //while (!rule_check(rr,which, world))
  //  {
  //    *r = UNIF_RANDUM();
  //    rr = propose_gamma_newparam (param, which, world, r);
  //  }
  return rr;
}

///
/// Hastings ratio calculator for exponential distribution
MYREAL hastings_ratio_expb(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam)
{
  (void) newparam;
  (void) oldparam;
  (void) delta;
  (void) r;
  (void) bayes;
  (void) whichparam;

    return 0.;
}

///////////////////////////////////////////////////////////////////////
// Multiplier prior (not working?)
///
/// Multiplier move, the parameter is multiplied with the rate 1/(lambda*m)
///
/// Requirements:
///     a < m < b
///     a = 1/b
///
///
MYREAL
propose_mult_newparam (MYREAL param, long which, world_fmt *world, MYREAL *r)
{
    const MYREAL delta = world->bayes->delta[which];
    const MYREAL maxparam = world->bayes->maxparam[which];
    MYREAL multiplier;  // minmult < multiplier < maxmult
    MYREAL maxmult = delta;
    //MYREAL minmult = 1/maxmult;
    MYREAL lambda = 2. * log(maxmult); // tuning parameter \lambda = 2 ln(b)
    MYREAL np;
    multiplier = EXP(lambda * ((*r) - 0.5));	
    np = multiplier * param;
    (*r) = multiplier;
    
    while(np > maxparam)
      {
	np = maxparam * maxparam / np;
      }
    //while (!rule_check(np,which, world))
    //{
    //	*r = UNIF_RANDUM();
    //	multiplier = EXP(lambda * ((*r) - 0.5));	
    //	np = multiplier * param;
    //	(*r) = multiplier;
    //	
    //	while(np > maxparam)
    //	  {
    //	    np = maxparam * maxparam / np;
    //	  }
    //}
    return np;
}

///
/// Hastings ratio calculator for exponential distribution with multiplicator
/// using the old parameter as mean value
/// simply the jacobian from x->param, matrix if derivatives [should be seocnd right?]
/// [check this with the green paper]
/// first derivatives
MYREAL hastings_ratio_mult(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam)
{
    // q(x|x')     1/(lambda(1/m))    | dx'/dx       dx'/du  |
    // ------------      -----------------------------    |                      |   =
    // q(x'|x)     1/(lambda m)       | du'/dx       du'/du  |
    //
    //          |  dmx/dx       dmx/dm     |
    // = m^2    |                          | = m^2  |(m (-1/m^2) - x zero| = m^2 m/m^2 = m
    //          |  d(1/m)/dx    d(1/m)/dm  |
  (void) newparam;
  (void) oldparam;
  (void) delta;
  (void) bayes;
  (void) whichparam;
  return log(r); // rate multiplier
}


///////////////////////////////////////////////////////////////////////
// ratio normal prior

MYREAL log_prior_ratio_normal(MYREAL newparam,
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

MYREAL log_prior_normal(world_fmt *world, long numparam)
{
    MYREAL val = world->param0[numparam];
    return log_prior_normal1(world,numparam,val);
}
///
/// normal prior distribution, returns 1/(sqrt(Pi*1) exp(-(x-mu)^2/2)
MYREAL log_prior_normal1(world_fmt *world, long numparam, MYREAL val)
{
    bayes_fmt * bayes = world->bayes;
    long i = numparam;
    MYREAL retval = (double) -HUGE;
    MYREAL p0 = val;
    MYREAL mean = bayes->meanparam[i];
    MYREAL x;
    if((p0 <= bayes->maxparam[i]) && (p0 >= bayes->minparam[i]))
    {
        x = p0 - mean;
        retval = LOG2PIHALF - x*x/2.0;
        // bayes->maxparam[i] - bayes->minparam[i];
        if (retval > 0.0)
            return retval;
    }
    return retval;
}

MYREAL hastings_ratio_normal(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam)
{
  (void) newparam;
  (void) oldparam;
  (void) delta;
  (void) r;
  (void) bayes;
  (void) whichparam;
    return 0.;
}

MYREAL
propose_normal_newparam (MYREAL param, long which, world_fmt *world, MYREAL *r)
{
  (void) param;
  bayes_fmt *bayes = world->bayes;
    MYREAL rr = (*r) - 0.5;
    MYREAL lower = bayes->minparam[which];
    //MYREAL meanparam = bayes->meanparam[which];
    MYREAL upper = bayes->maxparam[which];
    MYREAL std = 1.0;
  MYREAL mean = bayes->meanparam[which];
  if (bayes->hyperprior && bayes->hypercount++ % bayes->hyperinterval == 0)
    {
      hyper_normal(&mean, &std, mean, 1.0, lower, upper, 1.0);
    }
    rr = normal_rand(mean,std);
  return rr;
}


///////////////////////////////////////////////////////////////////////

void set_option_prior(prior_fmt **p, int type, MYREAL mini, MYREAL maxi, MYREAL mean, long bins)
{
  MYREAL b = 10.;
  //(*p) = calloc(1,sizeof(prior_fmt));
  (*p)->next = NULL;
  (*p)->number = -1;
  (*p)->type = type;
  switch(type)
    {
    case THETAPRIOR:
      strcpy((*p)->ptypename,"THETA"); break;
    case MIGPRIOR:
      strcpy((*p)->ptypename,"MIG"); break;
    case SPECIESTIMEPRIOR:
      strcpy((*p)->ptypename,"SPLIT"); break;
    case SPECIESSTDPRIOR:
      strcpy((*p)->ptypename,"SPLITSTD"); break;
    case RATEPRIOR:
      strcpy((*p)->ptypename,"RATE"); break;
    case GROWTHPRIOR:
      strcpy((*p)->ptypename,"GROWTH"); break;
    default:
      error("unknown prior in check_bayes_priors()");
      //break;
    }
  (*p)->delta = (maxi-mini)/b;
  (*p)->min = mini;
  (*p)->mean = mean;
  (*p)->max = maxi;
  (*p)->bins = bins;

  //  p->delta = (p->max - p->min)/b;
  (*p)->alpha = -1.0;
}

void is_priorkind(prior_fmt *p, char *priorkind)
{
  const char text[7][20] = {TUNIFORMPRIOR, TEXPPRIOR, TWEXPPRIOR, TMULTPRIOR, TGAMMAPRIOR, TNORMALPRIOR, TOTHER};
  const int numkind = NUMPRIORKIND;
  if (p->kind < NUMPRIORKIND)
    strncpy(priorkind,text[p->kind],PRIORKINDLENGTH);
  else
    strncpy(priorkind,text[NUMPRIORKIND-1],PRIORKINDLENGTH);
  return;
}
  

prior_fmt * copy_option_prior(prior_fmt *pmodel, option_fmt *options, MYREAL ratemin)
{
  (void) options;
  prior_fmt * p = (prior_fmt *) calloc(1,sizeof(prior_fmt));
  p->next = NULL;
  p->kind = pmodel->kind;
  p->type = pmodel->type;
  strcpy(p->ptypename, pmodel->ptypename);
  p->number = pmodel->number;
  p->min = pmodel->min+ratemin;
  p->mean = pmodel->mean;
  if (p->mean < p->min)
    p->mean = p->min;
  p->std = pmodel->std;
  p->max = pmodel->max;
  p->delta = pmodel->delta;
  p->alpha = pmodel->alpha;
  p->beta = pmodel->beta;
  p->bins = pmodel->bins;
  return p;
}


prior_fmt * insert_prior(prior_fmt *p, prior_fmt **plist,long *z,option_fmt * options)
{
  prior_fmt *q;
  if (plist[*z]->type == p->type)
    {
      q = copy_option_prior(plist[*z],options,0.0);
      p->next = q;
      p = q;
      (*z)++;
    }
  else
    {
      q = copy_option_prior(p,options,0.0);
      p->next = q;
      p = q;
    }
  return p;
}

prior_fmt * insert_prior_old(option_fmt *options, long index, prior_fmt *p, long type)
{
  (void) index;
  prior_fmt *q, *tmp;
  if (p==NULL)
    return NULL;
  if ((p->next==NULL) || (p->next != NULL &&  p->next->type != type))
    {
      q = copy_option_prior(p, options,0.0);
      tmp = p->next;
      p->next = q;
      q->next = tmp;
      p = q;
    }
  else
    return p->next;
  return p;
}
/*
long count_d(char *s)
{
  long count=0;
  char *t = s;
  while (*t != '\0')
    {
      if (*t == 'd' || *t == 'D')
	count++;
      t++;
    }
  return count;
  }*/


void copy_prior(prior_fmt *target, prior_fmt *source)
{
      target->next = source->next;
      target->kind = source->kind;
      target->type = source->type;
      strcpy(target->ptypename,source->ptypename);
      target->number = source->number;
      target->min = source->min;
      target->mean = source->mean;
      target->std = source->std;
      target->max = source->max;
      target->delta = source->delta;
      target->alpha = source->alpha;
      target->beta = source->beta;
      target->bins = source->bins;
      target->from = source->from;
      target->to = source->to;
      memcpy(target->v,source->v,sizeof(float)*4);
      target->random = source->random;
      target->cdf = source->cdf;
    }

 // find_prior returns an available prior with type this version returns a pointer to bayes_priors
prior_fmt * find_prior_menu(long from, long to, long priortype, option_fmt * options)
{
  prior_fmt * r = NULL;
  prior_fmt * p = options->bayes_priors;
  long i=0;
  for (i=0;i<options->bayes_priors_num;i++)
    {
      if (priortype == p[i].type)
	{
	  if(from == p[i].from)
	    {
	      if (to == p[i].to)
		{
		  r = &p[i];
		  return r;
		}
	    }
	  else
	    {
	      if(-1 == p[i].from)
		{
		  if(-1 == p[i].to)
		    {
		      r = &p[i];
		      return r;
		    }
		}
	    }
	}
    }
  // default if no prior found add default:
  //
  options->bayes_priors_num +=1;
  options->bayes_priors = (prior_fmt *) myrealloc(options->bayes_priors, options->bayes_priors_num * sizeof(prior_fmt));
  r = &options->bayes_priors[options->bayes_priors_num-1];
  r->from = from;
  r->to = to;
  switch(priortype)
    {
    case THETAPRIOR:
	set_option_prior(&r, THETAPRIOR, SMALLEST_THETA, DNA_GUESS_THETA*10.0, DNA_GUESS_THETA, BAYESNUMBIN);
      break;
    case MIGPRIOR:
	set_option_prior(&r,MIGPRIOR, SMALLEST_MIGRATION, DNA_GUESS_MIG*10.0, DNA_GUESS_MIG, BAYESNUMBIN);
      break;
    case RATEPRIOR:
      set_option_prior(&r,RATEPRIOR, SMALLEST_RATE, BIGGEST_RATE, 1.0, BAYESNUMBIN);	
      break;
    case SPECIESTIMEPRIOR:
      set_option_prior(&r,SPECIESTIMEPRIOR, SMALLEST_DNASPECIES, BIGGEST_DNASPECIES, DNA_GUESS_THETA , BAYESNUMBIN);
      break;
    case SPECIESSTDPRIOR:
      set_option_prior(&r,SPECIESSTDPRIOR, SMALLEST_DNASPECIES, BIGGEST_DNASPECIES, DNA_GUESS_THETA , BAYESNUMBIN);
      break;
    case GROWTHPRIOR:
      set_option_prior(&r, GROWTHPRIOR, LOWERGROWTH, 0.0, UPPERGROWTH, BAYESNUMBIN);	   
      break;
    default:
      error("no default prior found");
    }
  return r;
}

long set_all_default_priors(option_fmt *options)
{
  if (options->bayes_priors_num == 0)
    {
      set_default_prior(options, -1, -1, THETAPRIOR);
      if(options->numpop>1 || strlen(options->custm)>1)
	{
	  set_default_prior(options, -1, -1, MIGPRIOR);
	  //check options->custm for off-diagonal dD0
	  if((strchr(options->custm, 'd')!=NULL) || (strchr(options->custm, 'D')!=NULL))
	    {
	      set_default_prior(options, -1, -1, SPECIESTIMEPRIOR);
	      set_default_prior(options, -1, -1, SPECIESSTDPRIOR);
	    }
	}
      if(options->rate)
	set_default_prior(options, -1, -1, RATEPRIOR);
    }
  return options->bayes_priors_num;
}

prior_fmt * set_default_prior(option_fmt *options, long from, long to, long priortype)
{
  options->bayes_priors_num +=1;
  options->bayes_priors = (prior_fmt *) myrealloc(options->bayes_priors, options->bayes_priors_num * sizeof(prior_fmt));
  prior_fmt * r = &options->bayes_priors[options->bayes_priors_num-1];
  r->from = from;
  r->to = to;
  int datatype = (int) options->datatype;
  r->kind = UNIFORMPRIOR;
  switch(priortype)
    {
    case THETAPRIOR:
      if(strchr(SEQUENCETYPES, datatype))
	set_option_prior(&r, THETAPRIOR, SMALLEST_THETA, DNA_GUESS_THETA*10.0, DNA_GUESS_THETA, BAYESNUMBIN);
      else
	set_option_prior(&r, THETAPRIOR, SMALLEST_THETA, ALLELE_GUESS_THETA*10.0, ALLELE_GUESS_THETA, BAYESNUMBIN);	   
      break;
    case MIGPRIOR:
      if(strchr(SEQUENCETYPES, datatype))
	set_option_prior(&r,MIGPRIOR, SMALLEST_MIGRATION, DNA_GUESS_MIG*10.0, DNA_GUESS_MIG, BAYESNUMBIN);
      else
	set_option_prior(&r,MIGPRIOR, SMALLEST_MIGRATION, ALLELE_GUESS_MIG*10.0, DNA_GUESS_MIG, BAYESNUMBIN);
      break;
    case RATEPRIOR:
      set_option_prior(&r,RATEPRIOR, SMALLEST_RATE, BIGGEST_RATE, 1.0, BAYESNUMBIN);	
      break;
    case SPECIESTIMEPRIOR:
      if(strchr(SEQUENCETYPES, datatype))
	set_option_prior(&r,SPECIESTIMEPRIOR, SMALLEST_SPLIT, DNA_GUESS_SPLIT * 50.0, DNA_GUESS_SPLIT , BAYESNUMBIN);
      else
	set_option_prior(&r,SPECIESTIMEPRIOR, SMALLEST_SPLIT, ALLELE_GUESS_SPLIT*50.0, ALLELE_GUESS_SPLIT , BAYESNUMBIN);
      break;
    case SPECIESSTDPRIOR:
      if(strchr(SEQUENCETYPES, datatype))
	set_option_prior(&r,SPECIESSTDPRIOR, SMALLEST_SPLITSTD, DNA_GUESS_SPLITSTD * 50.0, DNA_GUESS_SPLITSTD , BAYESNUMBIN);
      else
	set_option_prior(&r,SPECIESSTDPRIOR, SMALLEST_SPLITSTD, ALLELE_GUESS_SPLITSTD * 50.0, ALLELE_GUESS_SPLITSTD , BAYESNUMBIN);
      break;
    default:
      error("no default prior found");
    }
  return r;
}


 // find_prior returns an available prior with type 
void find_prior(long from, long to, long priortype, option_fmt * options, prior_fmt *result)
{
  prior_fmt * r = NULL;
  prior_fmt * p = options->bayes_priors;
  long i=0;
  for (i=0;i<options->bayes_priors_num;i++)
    {
      if (priortype == p[i].type)
	{
	  if(from == p[i].from)
	    {
	      if (to == p[i].to)
		{
		  r = &p[i];
		  copy_prior(result,r);
		  return;
		}
	    }
	  else
	    {
	      if(-1 == p[i].from)
		{
		  if(-1 == p[i].to)
		    {
		      r = &p[i];
		      copy_prior(result,r);
		      return;
		    }
		}
	    }
	  //	  copy_prior(result,r);
	  //return; 
	}
    }
  // default if no prior found add default:
  //
  result->from = from;
  result->to = to;
  switch(priortype)
    {
    case THETAPRIOR:
      set_option_prior(&result, THETAPRIOR, SMALLEST_THETA, DNA_GUESS_THETA*10.0, DNA_GUESS_THETA, BAYESNUMBIN);	   
      break;
    case MIGPRIOR:
      set_option_prior(&result,MIGPRIOR, SMALLEST_MIGRATION, DNA_GUESS_MIG*10.0, DNA_GUESS_MIG, BAYESNUMBIN);
      break;
    case RATEPRIOR:
      set_option_prior(&result,RATEPRIOR, SMALLEST_RATE, BIGGEST_RATE, 1.0, BAYESNUMBIN);	
      break;
    case SPECIESTIMEPRIOR:
      set_option_prior(&result,SPECIESTIMEPRIOR, SMALLEST_DNASPECIES, BIGGEST_DNASPECIES, DNA_GUESS_THETA , BAYESNUMBIN);
      break;
    case SPECIESSTDPRIOR:
      set_option_prior(&result,SPECIESSTDPRIOR, SMALLEST_DNASPECIES, BIGGEST_DNASPECIES, DNA_GUESS_THETA , BAYESNUMBIN);
      break;
    case GROWTHPRIOR:
      set_option_prior(&result, GROWTHPRIOR, LOWERGROWTH, 0.0, UPPERGROWTH, BAYESNUMBIN);	   
      break;
    default:
      error("no default prior found");
    }
}

/// checks the settings of the number of long an short chain for bayes options and resets useless settings
void check_bayes_priors(option_fmt *options, data_fmt *data, world_fmt *world)
{
  (void) data;
  const long numpop = world->numpop;
  const long numpop2 = numpop * numpop;
  const int  has_mu = (int) options->bayesmurates;
  const long a = world->species_model_size * 2;
  const long np = numpop2 + has_mu + a + world->options->growpops_numalloc;
  //MYREAL ratemin = 0.0;
  long       i;
  long       j=0;
  //long       pos;
  double      beta=0.0;
  //longpair   *typelist; 
  //long oldnp;
  long z;
  long w=0;
  long from;
  long to;
  prior_fmt  *plist = (prior_fmt *) mycalloc(np, sizeof(prior_fmt));
  for (i=0; i<numpop;i++)
    {
      find_prior(i, i, THETAPRIOR, options, &plist[w]);//uses options->bayes_priors [uses a return ptr!]
      plist[w++].bins = options->bayes_posterior_bins[THETAPRIOR];
    }
  for (z=0; z<numpop;z++)
    for (i=z+1; i<numpop;i++)
    {
      find_prior(i, z, MIGPRIOR, options, &plist[w]);
      plist[w++].bins = options->bayes_posterior_bins[MIGPRIOR];
      find_prior(z, i, MIGPRIOR, options, &plist[w]);
      plist[w++].bins = options->bayes_posterior_bins[MIGPRIOR];
    }
  if(has_mu)
    {
      find_prior(numpop2, numpop2, THETAPRIOR, options, &plist[w]);
      plist[w++].bins = options->bayes_posterior_bins[RATEPRIOR];
    }
  z = numpop2 + has_mu;
  if (options->has_speciation)
    {
      for (i=0; i < numpop2; i++)
	{
	  if (strchr("tdDT", options->custm2[i]))
	    {
	      m2mm(i,numpop,&from,&to);
	      find_prior(from, to, SPECIESTIMEPRIOR, options, &plist[w]);
	      plist[w++].bins = options->bayes_posterior_bins[SPECIESTIMEPRIOR];
	      find_prior(from, to, SPECIESSTDPRIOR, options, &plist[w]);
	      plist[w++].bins = options->bayes_posterior_bins[SPECIESSTDPRIOR];
	    }
	}
    }
  if (world->has_growth)
    {
      for (i=0; i<numpop;i++)
	{
	  //printf("@w %li %li\n",i, w);
	  find_prior(i, i, GROWTHPRIOR, options, &plist[w]);//uses options->bayes_priors [uses a return ptr!]
	  plist[w++].bins = options->bayes_posterior_bins[GROWTHPRIOR];
	}
    }
  myfree(options->bayes_priors);
  options->bayes_priors = plist;
  options->bayes_priors_num = np;
  for(j=0;j<np;j++)
    {
      options->bayes_priors[j].v[0]= (float) options->bayes_priors[j].min;
      options->bayes_priors[j].v[1]= (float) options->bayes_priors[j].max;
      //options->bayes_priors[j].v[2]= (float) options->bayes_priors[j].mean;
      //options->bayes_priors[j].v[3]= (float) options->bayes_priors[j].std;
      switch(options->bayes_priors[j].kind)
	{
	case EXPPRIOR:
	case WEXPPRIOR:
	  options->bayes_priors[j].v[2]= (float) options->bayes_priors[j].mean;
	  options->bayes_priors[j].v[3]= (float) options->bayes_priors[j].std;
	  options->bayes_priors[j].random = trunc_random_exp; 
	  options->bayes_priors[j].cdf = trunc_cdf_exp;
	  break;
	case GAMMAPRIOR:
	  beta = (find_beta_truncgamma(options->bayes_priors[j].mean, options->bayes_priors[j].alpha, 
				       options->bayes_priors[j].min, options->bayes_priors[j].max));
	  options->bayes_priors[j].beta = beta;
	  options->bayes_priors[j].v[2]= (float) options->bayes_priors[j].alpha;
	  options->bayes_priors[j].v[3]= (float) options->bayes_priors[j].beta;
	  options->bayes_priors[j].random = trunc_random_gamma; 
	  options->bayes_priors[j].cdf = trunc_cdf_gamma;
	  break;
	case NORMALPRIOR:
	  options->bayes_priors[j].v[2]= (float) options->bayes_priors[j].mean;
	  options->bayes_priors[j].v[3]= (float) options->bayes_priors[j].std;
	  options->bayes_priors[j].random = trunc_random_normal; 
	  options->bayes_priors[j].cdf = trunc_cdf_normal;
	  break;
	case UNIFORMPRIOR:
	default:
	  options->bayes_priors[j].v[2]= (float) options->bayes_priors[j].mean;
	  options->bayes_priors[j].v[3]= (float) options->bayes_priors[j].std;
	  options->bayes_priors[j].random = trunc_random_uni; 
	  options->bayes_priors[j].cdf = trunc_cdf_uni;
	  break;
	}
    }
}


prior_fmt * get_prior_list(prior_fmt **list, int type)
{
  prior_fmt * p = list[0];
  prior_fmt *oldp = p;
  
  while(p->type != type)
    {
      p = p->next;
      if (p==NULL)
	return NULL;
    }
  while(p->type == type)
    {
      oldp = p;
      p = p->next;
      if (p==NULL)
	break;
    }
  return oldp;
}

// test section
// compile using this:
// gcc -o priortest -g priors.c random.c sighandler.c tools.c -DPRIORTEST -DMERSENNE_TWISTER -DNOJPEGLIB -DMEXP=19937
// then call by  
// priortest alpha beta
// it will print 10000 numbers
// the mean and standard deviation
//
#ifdef PRIORTEST
#include <stdio.h>
#include "sighandler.h"

char * generator;
int myID;
long *seed;

 int main(long argc, char **argv)
 {
   long i;
   MYREAL xx;
   MYREAL pxx;
   MYREAL a;
   MYREAL b;
   MYREAL mean=0.0;
   MYREAL var=0.0;
   
   world_fmt *world;
   world = (world_fmt*) calloc(1,sizeof(world_fmt));
   world->bayes = (bayes_fmt*) calloc(1,sizeof(bayes_fmt));
   world->bayes->maxparam = calloc(1,sizeof(MYREAL));
   world->bayes->minparam = calloc(1,sizeof(MYREAL));
   world->bayes->alphaparam = calloc(1,sizeof(MYREAL));
   world->bayes->betaparam = calloc(1,sizeof(MYREAL));
   world->bayes->meanparam = calloc(1,sizeof(MYREAL));
   world->numpop=1;

   generator = (char *) mycalloc (1,sizeof(char) * 80);

   init_gen_rand(123789);

   a = atof(argv[1]);
   b = atof(argv[2]);
   world->bayes->alphaparam[0] = a;
   world->bayes->meanparam[0] = a*b; 
   world->bayes->minparam[0] = 1.0; 
   world->bayes->maxparam[0] = 50.0; 
   world->bayes->betaparam[0] = find_beta_truncgamma(a*b, a, world->bayes->minparam[0],world->bayes->maxparam[0]);
   printf("Truncated gamma distribution with alpha=%f, beta=%f, lower=%f, upper=%f\n",
	  a,world->bayes->betaparam[0],world->bayes->minparam[0],world->bayes->maxparam[0]);
   b = world->bayes->betaparam[0];
   for (i=0;i<10000;i++)
     {
       xx=trunc_gamma_rand(a,b,world->bayes->minparam[0], world->bayes->maxparam[0]);
       pxx = log_prior_gamma1(world, 0, xx);
       mean += (xx - mean)/(i+1);
       var  += (xx - mean) * (xx-mean);
       printf("%f %f\n",xx, pxx);
     }
   printf("results of random gamma truncated 0 and 500, using a=%f, b=%f\n",a,b);
   printf("Mean=%f Standard deviation = %f\n",mean, sqrt(var/10000));
   printf("Expected=%f expected standard deviation = %f\n",a*b, sqrt(a)*b);
   MYREAL v = (world->bayes->maxparam[0]-world->bayes->minparam[0]) * 0.2;
   MYREAL fx = cdf_gamma(a,b,0.5);
   MYREAL fx2 = logpdf_gamma(a,b,v);
   MYREAL fx3 = logpdf_truncgamma(a,b,world->bayes->minparam[0],world->bayes->maxparam[0],v);
   printf("v=%f CDF(%f,%f,%f)=%f PDF(%f,%f,%f)=%f TPDF(%f,%f,%f,%f,%f)=%f\n",v,a,b,0.5,fx,a,b,v,fx2,a,b,
	  world->bayes->minparam[0],world->bayes->maxparam[0],v,fx3);
 }


float CDF(prior_fmt *prior, float p)
{
  float x = prior->cdf(prior->values);
  return x;
}

float random_prior(prior_fmt *prior)
{
  float x = prior->random(prior->values);
  return x;
}

#endif
