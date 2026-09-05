/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 R e p o r t e r   R O U T I N E  reports things progress=True or verbose
 
 Peter Beerli
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2012 Peter Beerli, Tallahassee FL
 
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
 
 
$Id: reporter.c 2158 2013-04-29 01:56:20Z beerli $
 
-------------------------------------------------------*/
/* \file reporter.c 
Routines that report progress and also calculates Gelman-Rubin convergence statistic
*/

#include "migration.h"
#include "mcmc.h"


#include "random.h"
#include "tools.h"
#include "options.h"
#include "sighandler.h"
#include "migrate_mpi.h"
#include "speciate.h"

#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif


void both_chain_means (MYREAL *mc, MYREAL *lc, MYREAL *tc, long len,
                       long lastn, long n);
void calc_gelmanb (MYREAL *gelmanb, MYREAL *mc, MYREAL *tc, MYREAL *lc,
                   long len, long lastn, long n);
void calc_gelmanw (MYREAL *gelmanw, world_fmt * world, MYREAL *mc, MYREAL *tc,
                   long len, long lastn, long n);
void calc_gelmanr (MYREAL *gelmanr, MYREAL *gelmanw, MYREAL *gelmanb,
                   long len, long lastn, long n);
void calc_average_biggest_gelmanr (MYREAL *gelmanr, long len, MYREAL *meanR,
                                   MYREAL *bigR);
void print_gelmanr (MYREAL average, MYREAL biggest);
MYREAL calc_s (long tthis, MYREAL *tc, world_fmt * world);
MYREAL calc_s_bayes (long tthis, MYREAL *tc, world_fmt * world);
void chain_means (MYREAL *thischainmeansm, world_fmt * world);

void calc_gelmanw2 (MYREAL *gelmanw, MYREAL *s1, MYREAL *s2, long len);
void all_chain_means (MYREAL *mc, MYREAL *chainmeans, long *nmeans, long len, long maxreplicate);
void calc_allgelmanb (MYREAL *gelmanb, MYREAL *mc, MYREAL *chainmeans, long *nmeans, long len, long maxreplicate);
void calc_allgelmanw2 (MYREAL *gelmanw, MYREAL *chain_s, long *nmeans, long len, long maxreplicate);
void calc_allgelmanr (MYREAL *gelmanr, MYREAL *gelmanw, MYREAL *gelmanb, long *nmeans, long len, long maxreplicate);

void collect_acceptance(world_fmt *world);
void collect_ess_values(world_fmt *world);
void print_bayes_ess(FILE *file, world_fmt *world, MYREAL *autocorr, MYREAL *effsample);
MYREAL single_chain_var(world_fmt *world, unsigned long T, MYREAL *variance, MYREAL *autoc, MYREAL *effsample);
void calculate_ess_frombayes(world_fmt *world, long T, MYREAL *params, long locus, MYREAL *autoc, MYREAL *effsample);
boolean max_ess(const MYREAL * ess, long n, const MYREAL minimum, MYREAL *miness);
void chain_means_ml (MYREAL *thischainmeans, world_fmt * world);
void chain_means_bayes (MYREAL *thischainmeans, world_fmt * world);
void  calc_chain_s(MYREAL *cs, MYREAL *cm, world_fmt *world, long replicate);
void convergence_progress(FILE *file, world_fmt *world);
void methods(world_fmt *world);
void citations(world_fmt *world);
void convergence_check (world_fmt * world, boolean progress);
void convergence_check_bayes (world_fmt *world,  long maxreplicate);
//##



//public functions
///
/// convergence indicator for ML runs

void convergence_check (world_fmt * world, boolean progress)
{
    static MYREAL *lastchainmeans, *chainmeans, *thischainmeans;
    static MYREAL *gelmanw, *gelmanb, *gelmanr;

    static boolean done = FALSE;
    static long len = 0;
    static long lastn = 0;
    static boolean first = TRUE;

    long maxreplicate = (world->options->replicate
                         && world->options->replicatenum >
                         0) ? world->options->replicatenum : 1;

    long n = 0;
    if (world->chains == 1 && maxreplicate <= 1)
        return;

    if (world->start)
        first = TRUE;
    if (progress || world->options->gelman)
    {
        if (!done)
        {
            done = TRUE;
            // len defines the length of arrays that
            // have to hold all km, kt, p, and mindex means (ml) or parametes (bayes)
	    len = world->numparam;

            lastchainmeans = (MYREAL *) mycalloc (len, sizeof (MYREAL));
            thischainmeans = (MYREAL *) mycalloc (len, sizeof (MYREAL));
            chainmeans = (MYREAL *) mycalloc (len, sizeof (MYREAL));
            gelmanw = (MYREAL *) mycalloc (len, sizeof (MYREAL));
            gelmanb = (MYREAL *) mycalloc (len, sizeof (MYREAL));
            gelmanr = (MYREAL *) mycalloc (len, sizeof (MYREAL));
        }
        n = -1; // a dummy value because the gelman_ functions expect a dummy last value
        memset (thischainmeans, 0, sizeof (MYREAL) * (size_t) len);
        if (first)
        {
            first = FALSE;
            chain_means (lastchainmeans, world);
            return;
        }
        else
        {
            chain_means (thischainmeans, world);
            both_chain_means (chainmeans, lastchainmeans, thischainmeans, len,
                              lastn, n);
            calc_gelmanb (gelmanb, chainmeans, thischainmeans, lastchainmeans,
                          len, lastn, n);
            calc_gelmanw (gelmanw, world, thischainmeans, lastchainmeans, len,
                          lastn, n);
            calc_gelmanr (gelmanr, gelmanw, gelmanb, len, lastn, n);
            calc_average_biggest_gelmanr (gelmanr, len, &world->convergence->gelmanmeanRall,
                                          &world->convergence->gelmanmaxRall);
            memcpy (lastchainmeans, thischainmeans, sizeof (MYREAL) * (size_t) len);
            lastn = n;
        }
    }
}

///
/// convergence indicator for Bayesian runs
void convergence_check_bayes (world_fmt *world,  long maxreplicate)
{
    MYREAL *gelmanw, *gelmanb, *gelmanr;
    long len = 0;
    long n_i = 0;
    long n_j = 0;
    long i;
    long j;

    MYREAL *chain_i;
    MYREAL *chain_j;
    MYREAL *s_i;
    MYREAL *s_j;
    MYREAL *chain_averages;
    MYREAL *chain_means = world->convergence->chain_means;
    long   *chain_counts = world->convergence->chain_counts;
    MYREAL *chain_s = world->convergence->chain_s;
    long   *nmeans;

    if (world->chains == 1 && maxreplicate <= 1)
        return;
    // len defines the length of arrays that
    /* BUG FIX: this used to be numparam+1, disagreeing with
       chain_means_bayes()'s own numparam-wide slots (world.c's/main.c's
       allocation used yet a third formula, numpop2+1 -- fixed there too)
       -- the stray "+1" never corresponded to a value chain_means_bayes
       ever wrote, so gelmanw/gelmanb/gelmanr's last slot was always a
       divide-by-zero-derived NaN that poisoned the averaged R statistic. */
    len = world->numparam;
    nmeans  = (long *) mycalloc (maxreplicate, sizeof (long));
    gelmanw = (MYREAL *) mycalloc (len, sizeof (MYREAL));
    gelmanb = (MYREAL *) mycalloc (len, sizeof (MYREAL));
    gelmanr = (MYREAL *) mycalloc (len, sizeof (MYREAL));
    chain_averages = (MYREAL *) mycalloc (len, sizeof (MYREAL));

    // can be done independently before this function
    // evaluate the within_variances for each monitored parameter: 
    // wp[i] = 1/(n-1) Sum[p[i]- mean[p[i]]]
    // evaluate the total_within_variance w = 1/m Sum[wp[i]]
    // needs to be done here
    // evaluate between_variance for each monitored parameter: 
    // bp[i] =  Sum[p[j]-mean[p[j]]]
    // where m is the number of indpendent chains, m = 2 ... many, shall we evaluate all pairs?
    // evaluate between_variance for all parameters 
    // b = n/(m-1) .... [check program?
    // calculate expected variance V = (1-1/n) W + 1/n B
    // Gelman-Rubin statistic sqrt(R) = Sqrt[V/W]

    for(i = 0; i < maxreplicate; i++)
      {
	chain_i = &chain_means[i*len];
	s_i = &chain_s[i*len];
	n_i = chain_counts[i];
	nmeans[i] = n_i;
	//report_values(chain_i, len, "chain_i");
	// setting whether this replicate is already recorded or not
	// would allow to use asynchronous MPI stuff
	if(chain_i[0] <= 0.)
	  continue;
	for(j = 0; j < i; j++)
	  {
	    chain_j = &chain_means[j*len];
	    s_j = &chain_s[j*len];
	    n_j = chain_counts[j];
	    //report_values(chain_j, len, "chain_j");
	    if(chain_j[0] <= 0.)
	      continue;
	    both_chain_means (chain_averages, chain_i, chain_j, len, n_i, n_j);
	    //report_values(chain_averages, len, "averages");
	    calc_gelmanb (gelmanb, chain_averages, chain_i, chain_j, len, n_i, n_j);
	    //report_values(gelmanb, len, "gelmanb");
	    calc_gelmanw2 (gelmanw, s_i, s_j, len);
	    //report_values(gelmanw, len, "gelmanw");
	    calc_gelmanr (gelmanr, gelmanw, gelmanb, len, n_i, n_j);
	    //report_values(gelmanr, len, "gelmanr");
	    calc_average_biggest_gelmanr (gelmanr, len, 
					  &world->convergence->gelmanmeanmaxR[j * maxreplicate + i],
					  &world->convergence->gelmanmeanmaxR[i * maxreplicate + j]);
	  }
      }
    all_chain_means(chain_averages,chain_means, nmeans, len, maxreplicate);
    calc_allgelmanb(gelmanb,chain_averages, chain_means, nmeans, len, maxreplicate);
    calc_allgelmanw2 (gelmanw, chain_s, nmeans, len, maxreplicate);
    calc_allgelmanr (gelmanr, gelmanw, gelmanb, nmeans, len, maxreplicate);
    calc_average_biggest_gelmanr (gelmanr, len, 
				  &world->convergence->gelmanmeanRall,
				  &world->convergence->gelmanmaxRall);
    
    myfree(gelmanw);
    myfree(gelmanb);
    myfree(gelmanr);
    myfree(chain_averages);
    myfree(nmeans);
}



void
both_chain_means (MYREAL *mc, MYREAL *lc, MYREAL *tc, long len, long lastn,
                  long n)
{
    long i;

    for (i = 0; i < len; i++)
    {
        mc[i] = (lc[i] * lastn + tc[i] * n) / (n + lastn);
    }
}

///
/// calculates the overall means of all replicate chains
void
all_chain_means (MYREAL *mc, MYREAL *chainmeans, long *nmeans, long len, long maxreplicate)
{
    long i;
    long j;
    long nsum = 0;
    MYREAL sum = 0.;
    for (i = 0; i < len; i++)
    {
      sum = 0.;
      nsum = 0;
      for(j=0; j < maxreplicate; j++)
	{
	  //if(nmeans[j]>1)
	  //  {
	      sum += chainmeans[j * len + i] * nmeans[j];
	      nsum += nmeans[j];
	      //  }
	}
      mc[i] = sum / nsum;
    }
}


void
calc_gelmanb (MYREAL *gelmanb, MYREAL *mc, MYREAL *tc, MYREAL *lc, long len,
              long lastn, long n)
{
  (void) lastn;
  (void) n;
    long i;

    for (i = 0; i < len; i++)
    {
      gelmanb[i] = /* 1/(2-1)* */  (pow ((lc[i] - mc[i]), 2.) + (pow ((tc[i] - mc[i]), 2.)));
    }

}
void
calc_allgelmanb (MYREAL *gelmanb, MYREAL *mc, MYREAL *chainmeans, long *nmeans, long len, long maxreplicate)
{
    long i;
    long j;
    //long nsum = 0;;
    MYREAL sum = 0.;
    MYREAL val;
    for (i = 0; i < len; i++)
      {
	sum = 0.;
	//nsum = 0;
	for(j=0; j < maxreplicate; j++)
	  {
	    if(nmeans[j]>1)
	      {
		val = chainmeans[j * len + i] -  mc[i];
		sum += val * val;
		//nsum += nmeans[j];
	      }
	}
      gelmanb[i] = sum / (maxreplicate -1);

      //gelmanb[i] = nn * (pow ((lc[i] - mc[i]), 2.) + (pow ((tc[i] - mc[i]), 2.)));
    }

}


///
/// collect the autocorrelation values and the ESS values
void collect_acceptance(world_fmt *world)
{
  const long npp = world->numparam;
  long i;
  // parameter
  for(i=0;i<npp; i++)
    {
       if(world->bayes->map[i][1] != INVALID)
	 {
	   world->accept_archive[i] += world->bayes->accept[i];
	   world->trials_archive[i] += world->bayes->trials[i];
	 }
    }
  // genealogy
  world->accept_archive[npp] += world->bayes->accept[npp];
  world->trials_archive[npp] += world->bayes->trials[npp];
}

///
/// collect the autocorrelation values and the ESS values
void collect_ess_values(world_fmt *world)
{
  const long nnn = world->numparam;
  long i;
  //parameter
  for(i=0;i<nnn; i++)
    {
      if(world->bayes->map[i][1] != INVALID)
	{
	  // onepass mean of autocorrelation
	  world->auto_archive[i] += (world->autocorrelation[i] - world->auto_archive[i])/world->archive_n;
	  // summing ess values
	  world->ess_archive[i] += world->effective_sample[i];
	}
    }
  //genealogy
  // onepass mean of autocorrelation for genealogy
  world->auto_archive[nnn] += (world->autocorrelation[nnn] - world->auto_archive[nnn])/world->archive_n;
  // summing ess values for genealogy
  world->ess_archive[nnn] += world->effective_sample[nnn];  
  //n++;
  world->archive_n += 1;
}


//
// print out the acceptance ratios for all the different Bayesian updates
void
print_bayes_ess(FILE * file,  world_fmt *world, MYREAL *autocorr, MYREAL *effsample)
{
  long j0, j;          
  long topop    =0;  
  long frompop  =0;  
  char *stemp;       
  long trials   =0;
  long numpop = world->numparamcumvec[THETAPRIOR];
  long numpop2 = world->numparamcumvec[MIGPRIOR];
  long npr = world->numparamcumvec[RATEPRIOR];
  long nps = world->numparamcumvec[SPLITSTDPRIOR];
  long npg = world->numparamcumvec[GROWTHPRIOR];
  long npa = world->numparam;
  bayes_fmt *bayes = world->bayes;
  stemp = (char *) mycalloc(LINESIZE,sizeof(char));
  
  //species_fmt *s = NULL;
  long z=0;
  FPRINTF(file,"\n\n\nAutocorrelation for all parameters and the genealogies\n");
  FPRINTF(file,"-------------------------------------------------------------------\n\n");
  FPRINTF(file,"Parameter           Autocorrelation           Effective Sample size\n");
    // population sizes
    for(j0=0; j0 < npa; j0++)
    {
      //        if(!strchr("0c", bayes->custm2[j]))
      if(shortcut(j0,world,&j))
        {
	  continue;
	}
      else
	{
	  if (j<numpop) //theta
	    {
	      FPRINTF(file,"Theta_%-3li              %8.3f         %17.3f\n", j+1, autocorr[j], effsample[j]);
	      if(effsample[j]<ESSMINIMUM && file==world->outfile)
		record_warnings(world,"Param %li: Effective sample size of run seems too short! ",j+1);	  
	    }
	  else if (j < numpop2) // migration rates
	    {	      
	      m2mm (j, world->numpop, &frompop, &topop);
	      if(world->options->usem)
		{
		  mysnprintf(stemp,LINESIZE, "M_%li->%li", frompop+1, topop+1);
		}
	      else
		{
		  mysnprintf(stemp,LINESIZE, "xN_%lim_%li->%li", topop+1, frompop+1, topop+1);
		}
	      FPRINTF(file, "%-12.12s           %8.3f         %17.3f\n", stemp, autocorr[j],effsample[j]);
	      if(effsample[j]<ESSMINIMUM && file==world->outfile)
		record_warnings(world,"Param %li: Effective sample size of run seems too short! ",j+1);		
	    }
	  else if(bayes->mu && j==world->numpop2)
	    {	      
	      FPRINTF(file, "Rate of mutation rate (%li) %8.3f         %17.3f\n", j+1,autocorr[j],effsample[j]);
	      if(effsample[j]<ESSMINIMUM && file==world->outfile)
		record_warnings(world,"Param %li: Effective sample size of run seems too short! ",j+1);
	    }
	  else if (world->has_speciation && j < nps)
	    {
	      z=0;
	      species_fmt * s = get_which_species_model(j, world->species_model, world->species_model_size);
	      long from = s->from;
	      long to = s->to;
	      trials=world->trials_archive[j];
	      if (j == s->paramindex_mu)
		{
		  mysnprintf(stemp,LINESIZE,"_%li->%li",1+from,1+to);
		  FPRINTF(file, "D%-12.12s          %8.3f         %17.3f\n", stemp, autocorr[j],effsample[j]);
		}
	      else
		{
		  FPRINTF(file, "S%-12.12s          %8.3f         %17.3f\n", stemp, autocorr[j], effsample[j]);
		}
	      if(effsample[j]<ESSMINIMUM && file==world->outfile)
		record_warnings(world,"Param %li: Effective sample size of run seems too short! ",j+1);
	    }
	  else if(world->has_growth && j < npg)
	    {
	      trials=world->trials_archive[j];
	      long d = j - nps;
	      if (d < world->options->growpops_numalloc && world->options->growpops[d]==0)
		d++;
	      mysnprintf(stemp,LINESIZE,"Growth_%li",d+1);
	      FPRINTF(file, "%-12.12s           %8.3f         %17.3f\n", stemp, autocorr[j],effsample[j]);
	      if(effsample[j]<ESSMINIMUM && file==world->outfile)
		record_warnings(world,"Param %li: Growth: Effective sample size of run seems too short! ",j+1);
	    }
	  else if(world->has_mlalpha && world->tri_mlalpha != FIXED)
	    {
	      trials=world->trials_archive[j];
	      long d = j - npg;
	      if (d < world->options->mlalphapops_numalloc && world->options->mlalphapops[d]==0)
		d++;
	      mysnprintf(stemp,LINESIZE,"ML-alpha_%li",d+1);
	      FPRINTF(file, "%-12.12s           %8.3f         %17.3f\n", stemp, autocorr[j],effsample[j]);
	      if(effsample[j]<ESSMINIMUM && file==world->outfile)
		record_warnings(world,"Param %li: ML-alpha: Effective sample size of run seems too short! ",j+1);
	    }
	}
    }
    // accepted trees
    FPRINTF(file, "Genealogies            %8.3f         %17.3f\n", autocorr[npa],effsample[npa]);
    if(effsample[npa]<ESSMINIMUM && file==world->outfile)
      record_warnings(world,"Genealogies %li: Effective sample size of run seems too short! ",npa);		      
    if(world->loci>1 && file!=stdout)
      {
	FPRINTF(file,"(*) averaged over loci.\n");
      }
    myfree(stemp);
}




///
/// returns variance and other measurement of a single chain used in Bayesian burnin_chain() 
/// the standard deviation and the mean are calculated using a result of
/// B. P. Welford 1962 (from D. E. Knuth Art of computer programming Vol 2 page 232
/// the recursion form of the autocorrelation r was developed with Koffi Sampson
/// April 2007, see mathematica code below
/// newff[xx_] := Module[{},
///		     nt = Length[xx];
///		     r = 0;
///		     x1 = xx[[1]];
///		     mo = x1;
///		     mn = x1;
///		     n = 2;
///		     s = 0;
///		     While[n <= nt, 
///			   xo = xx[[n - 1]];
///			   xn = xx[[n]];
///			   mn = mn  + (xn - mn)/n; 
///			   s = s + (xn - mo) (xn - mn);
///KOFFI       r = r  + xo xn + mo (n mo  - x1 - xo) - mn ((n + 1) mn - x1 - xn);
///PETER       r = r + (xo - mo)(xn-mn)
///			   mo = mn;
///			   n++;
///		     ];
///  r/s
///    ]
/// The effective sample is calculated as the length the sample n * (1-r)/(1+r)
/// 
MYREAL single_chain_var(world_fmt *world, unsigned long T, MYREAL *variance, MYREAL *autoc, MYREAL *effsample)
{
  static double *mean;
  static double *S;
  static double *r;
  static double *xold;
  static double *xstart;
  static boolean done=FALSE;
  static long n = 0;
  static long nn = 0;
  long nn1 = 0;
  static long nps = 0;
  static long npg = 0;
  long i;
  long j;
  long start;

  double x  = 0.0;
  double xo = 0.0;
  double x1 = 0.0;
  double m  = 0.0;
  double mo = 0.0;
  double v  = 0.0;
  double *delta;
  double temp;
  if(!done)
    {
      nn = world->numparam;
      nn1 = nn + 1;
      nps = world->numparamcumvec[SPLITSTDPRIOR];
      npg = world->numparamcumvec[GROWTHPRIOR];
      mean = (double *) mycalloc(nn1,sizeof(double));
      S    = (double *) mycalloc(nn1,sizeof(double));
      r    = (double *) mycalloc(nn1,sizeof(double));
      xold = (double *) mycalloc(nn1,sizeof(double));
      xstart = (double *) mycalloc(nn1,sizeof(double));
      done=TRUE;
    }
  // reset static variable for each chain
  if(world==NULL)
    {
      n=0;
      memset(mean,0,(size_t) nn * sizeof(double));
      memset(S,0,(size_t) nn * sizeof(double));
      memset(r,0,(size_t) nn * sizeof(double));
      memset(xold,0,(size_t) nn * sizeof(double));
      memset(xstart,0,(size_t) nn * sizeof(double));
      return 0;
    }
  // syntax so that an init really works
  if(T==0)
    return 0;

  delta = (double *) mycalloc(nn+1,sizeof(double));
  //handles BAYES or ML
  start = 0; //(world->options->bayes_infer ? 0 : nn-1);
  if(n==0) //initialization of mean calculator
    {
      n += 1;
      for(i=start; i < nn; i++)
	{
	  if(shortcut(i,world,&j))
	    continue;
	  else
	    {
	      if (i<nps)
		{
		  j  = world->bayes->map[i][1];
		  x         = world->param0[j];
		}
	      else if (i < npg)
		{
		  j = i;
		  x = world->growth[world->options->growpops[j-nps]-1];
		}
	      else if (i < nn)
		{
		  j = i;
		  x = world->mlalpha[world->options->mlalphapops[j-npg]-1];
		}
	    }
	  mean[j]   = x;
	  S[j]      = 0.0;
	  r[j]      = 0.0;
	  xstart[j] = x;
	  xold[j] = x;
	}
      //genealogy
      x = world->likelihood[world->G];
      mean[nn]   = x;
      S[nn]      = 0.0;
      r[nn]      = 0.0;
      xstart[nn] = x;
      xold[nn] = x;
      v = 0.0;
      myfree(delta);
    }
  else
    {
      // n is at least 1
      n += 1;
      for(i=start; i < nn; i++)
	{
	  if(world->bayes->map[i][1] == INVALID)
	    continue;
	  else
	    {
	      j  = world->bayes->map[i][1];
	      x  = (double) world->param0[j];
	    }
	  xo        = xold[j];
	  x1        = xstart[j];
	  delta[j]  = x - mean[j];
	  mo        = mean[j];
	  mean[j]  += delta[j]/n;
	  m         = mean[j];
	  S[j]     += delta[j] * (x - m);
	  temp = xo * x + mo * (n * mo - x1 - xo) - m * ((n+1) * m - x1 - x);
	  r[j]     += temp;
	  autoc[j]  = S[j] > 0.0 ? (MYREAL) (r[j]/S[j]): (MYREAL) 1.0;
	  effsample[j] = (MYREAL) (n * (1. - autoc[j])/(1. + autoc[j]));
	  v        += S[j];
	  //printf("[%li] n=%li effsample=%g atoc=%g r=%g (%f) s=%g mean=%g\n", j, n,effsample[j],autoc[j],r[j], temp, S[j],mean[j]);
	  xold[j] = x;
	}
      x = world->likelihood[world->G];
      xo        = xold[nn];
      x1        = xstart[nn];
      delta[nn]  = x - mean[nn];
      mo        = mean[nn];
      mean[nn]  += delta[nn]/n;
      m         = mean[nn];
      S[nn]     += delta[nn] * (x - m);
      temp = xo * x + mo * (n * mo - x1 - xo) - m * ((n+1) * m - x1 - x);
      r[nn]     += temp;
      autoc[nn]  = S[nn] > 0.0 ? (MYREAL) (r[nn]/S[nn]): (MYREAL) 1.0;
      effsample[nn] = (MYREAL) (n * (1. - autoc[nn])/(1. + autoc[nn]));
      v        += S[nn];
      //printf("%i> [%li] n=%li effsample=%g atoc=%g r=%g (%f) s=%g mean=%g\n", myID, nn, n,effsample[nn],autoc[nn],r[nn], temp, S[nn],mean[nn]);
      xold[nn] = x;
      
      myfree(delta);
    }
  //  printf("@single_chain_var() [%li]",(long) world->cold);
  //for(i=0; i < nn; i++)
  //printf("%f ",autoc[i]);
  //printf("\n@ ");
  //for(i=0; i < nn; i++)
  //  printf("%li ",world->bayes->accept[i]);
  //printf("\n");
  *variance =  (MYREAL) v / n;
  return *variance;
}
///
/// this function mimics single_chain_var() but is used only in the reading bayesallfile function read_bayes_fromfile()
/// taking into account that there may be a mixture of different loci and some variables are not where one expect them
/// during the standard run
void calculate_ess_frombayes(world_fmt *world, long T, MYREAL *params, long locus, MYREAL *autoc, MYREAL *effsample)
{
  static double *mean;
  static double *S;
  static double *r;
  static double *xold;
  static double *xstart;
  static boolean done=FALSE;
  static long *n;
  static long nn = 0;
  long nnbase = 0;
  long nnbase1 = 0 ;
  long i;
  long j;
  long z;
  long start;

  double x  = 0.0;
  double xo = 0.0;
  double x1 = 0.0;
  double m  = 0.0;
  double mo = 0.0;
  //static double *v;
  double *delta;
  double temp;
  long nn1=0;
  if(!done)
    {
      nnbase = world->numparam;
      nnbase1 = nnbase + 1;
      nn = world->loci * (nnbase+1);
      n = (long *) mycalloc(world->loci,sizeof(long));
      mean = (double *) mycalloc(nn,sizeof(double));
      S    = (double *) mycalloc(nn,sizeof(double));
      r    = (double *) mycalloc(nn,sizeof(double));
      xold = (double *) mycalloc(nn,sizeof(double));
      xstart = (double *) mycalloc(nn,sizeof(double));
      done=TRUE;
    }
  // syntax so that an init really works
  if(T==0)
    return;

  delta = (double *) mycalloc(nn,sizeof(double));
  //handles BAYES or ML
  start = 0;//(world->options->bayes_infer ? 0 : nnbase - 1 );
  if(n[locus]==0) //initialization of mean calculator
    {
      n[locus] += 1;
      for(i=start; i < nnbase1; i++)
	{
	  if(i<nnbase)
	    {
	      if(world->bayes->map[i][1] == INVALID)
		continue;
	      else
		{
		  j  = world->bayes->map[i][1];
		  z  = locus * nnbase1 + j;
		  x  = params[j+2];
		}
	    }
	  else
	    {
	      z = locus * nnbase1 + nnbase;
	      x = params[1];
	    }
	  mean[z]   = x;
	  S[z]      = 0.0;
	  r[z]      = 0.0;
	  xstart[z] = x;
	  xold[z] = x;
	}
    }
  else
    {
      // n is at least 1
      n[locus] += 1;
      for(i=start; i < nnbase1; i++)
	{
	  if(i<nnbase)
	    {
	      if(world->bayes->map[i][1] == INVALID)
		continue;
	      else
		{
		  j  = world->bayes->map[i][1];
		  z  = locus * nnbase + j;
		  x  = params[j+2];
		}
	    }
	  else
	    {
	      z = locus * nnbase1 + nnbase;
	      x = params[1];
	    }	  
	  xo        = xold[z];
	  x1        = xstart[z];
	  delta[z]  = x - mean[z];
	  mo        = mean[z];
	  mean[z]  += delta[z]/n[locus];
	  m         = mean[z];
	  S[z]     += delta[z] * (x - m);
	  temp = xo * x + mo * (n[locus] * mo - x1 - xo) - m * ((n[locus]+1) * m - x1 - x);
	  r[z]     += temp;
	  autoc[z]  = S[z] > 0.0 ? (MYREAL) (r[z]/S[z]): (MYREAL) 1.0;
	  effsample[z] = (MYREAL) (n[locus] * (1. - autoc[z])/(1. + autoc[z]));
	  //v        += S[z];
	  //printf("[%li] n=%li effsample=%g atoc=%g r=%g (%f) s=%g mean=%g\n", j, n[locus],effsample[j],autoc[j],r[j], temp, S[j],mean[j]);
	  xold[z] = x;
	}
    }
  myfree(delta);    
  //  *variance =  (MYREAL) v[locus] / n[locus];
  //return *variance;
}

boolean max_ess(const MYREAL * ess, long n, const MYREAL minimum, MYREAL *miness)
{
  long i;
  long z=0;
  *miness = (double) HUGE;
  for(i=0; i < n; i++)
    {
      if (ess[i]< *miness)
	*miness = ess[i];
      z += (long) (ess[i] >= minimum);
    }
  if(z==n)
    {
      return TRUE;
    }
  return FALSE;
}

/// calculate effective sample size
/// which is the N' = N (1-r)/(1+r)
/// where r is the autocorrelation coefficient for lag 1.
/// with r(1) = sum[(x[i]-xmean)(x[i+1]-xmean)/standarddeviation(x),i,0,n-1]

void
calc_gelmanw (MYREAL *gelmanw, world_fmt * world, MYREAL *mc, MYREAL *tc,
              long len, long lastn, long n)
{
  (void) lastn;
  (void) n;
    long i;
    MYREAL s1, s2;
    for (i = 0; i < len; i++)
      {
	s1 = calc_s (i, tc, world);
	s2 = calc_s (i, mc, world);
	gelmanw[i] = 0.5 * (s1 + s2);
      }
}

void
calc_gelmanw2 (MYREAL *gelmanw, MYREAL *s1, MYREAL *s2, long len)
{
    long i;

    for (i = 0; i < len; i++)
    {
        gelmanw[i] = 0.5 * (s1[i] + s2[i]);
    }
}

///
/// calculates the overall s^2 of all replicate chains
void
calc_allgelmanw2 (MYREAL *gelmanw, MYREAL *chain_s, long *nmeans, long len, long maxreplicate)
{
    long i;
    long j;
    MYREAL sum = 0.;
    for (i = 0; i < len; i++)
    {
      sum = 0.;
      for(j=0; j < maxreplicate; j++)
	{
	  if(nmeans[j]>1)
	    {
	      sum += chain_s[j * len + i];
	    }
	}
      gelmanw[i] = sum / maxreplicate;
    }
}



void
calc_gelmanr (MYREAL *gelmanr, MYREAL *gelmanw, MYREAL *gelmanb, long len,
              long lastn, long n)
{
    long i;
    MYREAL nn = (n + lastn) / 2.;
    MYREAL sqplus = 0.;
    //    MYREAL v;
    for (i = 0; i < len; i++)
    {
      sqplus = ((nn - 1.) / nn * gelmanw[i] + gelmanb[i]);
      //v = 1.5 * sqplus;
      gelmanr[i] = sqrt (sqplus / gelmanw[i]);// - (nn-1)/(2. * nn));
    }
}

void
calc_allgelmanr (MYREAL *gelmanr, MYREAL *gelmanw, MYREAL *gelmanb, long *nmeans, long len, long maxreplicate)
{
    long i;
    long j;
    long nn=0;
    MYREAL sqplus;
    //    MYREAL v;

    for(j=0; j < maxreplicate; j++)
      {
	nn += nmeans[j];
      }
    nn = (long) (ceil(nn/maxreplicate));
    for (i = 0; i < len; i++)
    {
      sqplus = (nn - 1.) / nn * gelmanw[i] + gelmanb[i];
      //      v = (maxreplicate+1)/maxreplicate * sqplus;
      gelmanr[i] = sqrt (sqplus / gelmanw[i]);// - (nn-1)/(maxreplicate * nn));
    }
}

void
calc_average_biggest_gelmanr (MYREAL *gelmanr, long len,
                              MYREAL *meanR, MYREAL *bigR)
{
    long i;
    MYREAL average = 0;
    MYREAL biggest = 0.;
    for (i = 0; i < len; i++)
    {
        if (biggest < gelmanr[i])
            biggest = gelmanr[i];
        average += gelmanr[i];
    }
    if (len > 0)
        *meanR = average / len;
    else
        *meanR = average;
    *bigR = biggest;
}

///
/// calculate the variance of a chain tc is the parameter average 
MYREAL
calc_s (long tthis, MYREAL *tc, world_fmt * world)
{
  timearchive_fmt * atl     = &world->atl[world->rep][world->locus];
  long                    T = atl->T;
  long              numpop  = world->numpop;
  //  long              numpop2 = numpop * numpop;
  long              startp;
  long              startl;
  long              startkm;
  long              i;
  long              j;
  MYREAL            xx;
  MYREAL            s = 0.0;

  startkm = numpop;
  startp  = startkm + numpop;
  startl  = startp + numpop;

  if (tthis < startkm)
    {
      i = tthis;
      for (j = 0; j < T; j++)
	{
	  xx = atl->tl[j].kt[i] - tc[i];
	  s += xx * xx;
	}
      if(s>0.)
	s /= T - 1.;
      return s;
    }
    else
    {
        if (tthis < startp)
        {
            i = tthis - startkm;
	    for (j = 0; j < T; j++)
	      {
		xx = atl->tl[j].km[i] - tc[i];
		s += xx * xx;
	      }
	    s /= T - 1.;
	    return s;
        }
        else
        {
            if (tthis < startl)
	      {
                i = tthis - startp;
		for (j = 0; j < T; j++)
		  {
		    xx = atl->tl[j].p[i] - tc[i];
		    s += xx * xx;
		  }
		s /= T - 1.;
		return s;
	      }
            else
	      {
                i = tthis - startl;
		for (j = 0; j < T; j++)
		  {
		    xx = atl->tl[j].mindex[i] - tc[i];
		    s += xx * xx;
		  }
		s /= T - 1.;
		return s;
	      }
        }
    }
}

///
/// calculate the variance of a chain tc is the parameter average 
/// by supplying the length of the chain
/*
static MYREAL calc_s2 (long tthis, MYREAL *tc, world_fmt * world, long T)
{
  timearchive_fmt * atl     = &world->atl[world->rep][world->locus];
  long              numpop  = world->numpop;
  //  long              numpop2 = numpop * numpop;
  long              startp;
  long              startl;
  long              startkm;
  long              i;
  long              j;
  MYREAL            xx;
  MYREAL            s = 0.0;

  startkm = numpop;
  startp  = startkm + numpop;
  startl  = startp + numpop;

  if (tthis < startkm)
    {
      i = tthis;
      for (j = 0; j < T; j++)
	{
	  xx = atl->tl[j].kt[i] - tc[i];
	  s += xx * xx;
	}
      if(s>0.)
	s /= T - 1.;
      return s;
    }
    else
    {
        if (tthis < startp)
        {
            i = tthis - startkm;
	    for (j = 0; j < T; j++)
	      {
		xx = atl->tl[j].km[i] - tc[i];
		s += xx * xx;
	      }
	    s /= T - 1.;
	    return s;
        }
        else
        {
            if (tthis < startl)
	      {
                i = tthis - startp;
		for (j = 0; j < T; j++)
		  {
		    xx = atl->tl[j].p[i] - tc[i];
		    s += xx * xx;
		  }
		s /= T - 1.;
		return s;
	      }
            else
	      {
                i = tthis - startl;
		for (j = 0; j < T; j++)
		  {
		    xx = atl->tl[j].mindex[i] - tc[i];
		    s += xx * xx;
		  }
		s /= T - 1.;
		return s;
	      }
        }
    }
}
*/

///
/// calculate the variance of a chain tc is the parameter average 
MYREAL
calc_s_bayes (long tthis, MYREAL *tc, world_fmt * world)
{
  //long T            = world->convergence->chain_counts[world->rep];
  long nn           = 2+world->numparam;
  long pnum         = world->bayes->numparams;
  
  MYREAL  * params     = world->bayes->params;
  long              i;
  long              j;
  MYREAL            xx;
  MYREAL            s = 0.0;

  /* BUG FIX: bayes_save_parameter() (bayes.c) lays each recorded step out
     as [oldval, likelihood, param_0, param_1, ..., param_{n-1}] (an
     nn = 2+numparam -wide row) -- so the real parameters live at raw
     column tthis+2, not tthis. tc[] (the chain's own mean, from
     chain_means_bayes()) is already 0-indexed per real parameter,
     matching tc[tthis] to real parameter tthis. Reading
     params[j*nn+tthis] instead of params[j*nn+tthis+2] compared every
     parameter's raw samples against a mean two columns over (or, for
     tthis 0/1, against oldval/likelihood entirely), corrupting the
     within-chain variance for essentially every parameter. */
  i = tthis + 2;
  for (j = 0; j < pnum /*T*/; j++)
    {
      xx = params[j * nn + i] - tc[tthis];
      s += xx * xx;
    }
  if(s>0.)
    s /= pnum - 1.;
  return s;
}


///
/// calculates means for the different events on a tree for the Gelman-Rubin statistic
/// this version of the statistics operates on the compressed treedata and not the parameters
/// 
void chain_means_ml (MYREAL *thischainmeans, world_fmt * world)
{
  timearchive_fmt * atl     = &world->atl[world->rep][world->locus];
  long              T       = atl->T;
  long              numpop  = world->numpop;
  long              numpop2 = numpop * numpop;
  long              startp;
  long              startl;
  long              startkm;
  long              i;
  long              j;

  startkm = numpop;
  startp  = startkm + numpop;
  startl  = startp + numpop;
  
  for (j = 0; j < T; j++)
    {
      for (i = 0; i < numpop; i++)
        {
	  thischainmeans[i]           += atl->tl[j].kt[i];
	  thischainmeans[i + startkm] += atl->tl[j].km[i];
	  thischainmeans[i + startp]  += atl->tl[j].p[i];
        }
      for (i = 0; i < numpop2; i++)
	{
	  thischainmeans[i + startl]  += atl->tl[j].mindex[i];
	}
    }
  for (i = 0; i < numpop; i++)
    {
      thischainmeans[i]           /= T;
      thischainmeans[i + startkm] /= T;
      thischainmeans[i + startp]  /= T;
    }
  for (i = startl; i < numpop2 + startl; i++)
    {
      thischainmeans[i] /= T;
    }
}
///
/// calculates means for the different parameters for the Gelman-Rubin statistic
/// for bayesian inference
void chain_means_bayes (MYREAL *thischainmeans, world_fmt * world)
{
  //  long              T       = world->convergence->chain_counts[world->rep];
  long              T       = world->bayes->numparams;
  long              i;
  long              j;
  MYREAL           *params  = world->bayes->params;
  long              nn      = 2+world->numparam;

  for (j = 0; j < T; j++)
    {
      for (i = 2; i < nn; i++)
        {
	  thischainmeans[i-2]           += params[j * nn + i];
        }
    }
  for (i = 2; i < nn; i++)
    {
      thischainmeans[i-2]           /= T;
    }
}

///
/// chain_means
void
chain_means (MYREAL *thischainmeans, world_fmt * world)
{
  /* BUG FIX: this used to set chain_counts[world->rep] to its own
     current value if already nonzero, else 1 -- so once the very first
     call left it at 1, every later call was a no-op, and it could never
     become anything but 1 for the life of the run, regardless of how
     many posterior samples were actually recorded. chain_counts[] is
     read as a genuine per-replicate sample count/weight everywhere
     downstream (all_chain_means()'s weighted average, calc_allgelmanb()'s/
     calc_allgelmanw2()'s "nmeans[j] > 1" sanity guard) -- with the old
     logic that guard could never pass (nmeans[j] was always exactly 1),
     so convergence_check_bayes()'s combined-across-all-replicates
     within-/between-chain variances (gelmanw/gelmanb) were silently
     always 0, making the reported "Mean sqrt(R)"/"Maximum sqrt(R)"
     NaN/0 on every run that used replicates, independent of the
     numparam/offset fixes above. world->bayes->numparams is the same
     real per-locus/replicate sample count chain_means_bayes()/
     calc_s_bayes() already average over (their own T/pnum) -- using it
     here too makes chain_counts[] agree with what was actually
     recorded. */
  world->convergence->chain_counts[world->rep] = world->bayes->numparams;
  chain_means_bayes (thischainmeans, world);
}

///
/// calculates the average over 
void  calc_chain_s(MYREAL *cs, MYREAL *cm, world_fmt *world, long replicate)
{
  long len;
  long start;
  long stop;
  long i;
  /* BUG FIX: was numpop2+1, disagreeing with chain_means_bayes()'s own
     numparam-wide per-replicate slots in cm/cs (and with world.c's/
     main.c's matching allocation fix). */
  len = world->numparam;
  start = replicate * len;
  stop = start + len;
  for(i = start; i < stop; i++)
    {
      cs[i] = calc_s_bayes (i-start, &cm[start], world);
    }
}

///
/// report convergence statistic to screen or file
void convergence_progress(FILE *file, world_fmt *world)
{
  long i;
  long j;
  long maxreplicate = (world->options->replicate
		       && world->options->replicatenum >
		       0) ? world->options->replicatenum : 1;
  
  MYREAL *gelmanmeanR = world->convergence->gelmanmeanmaxR; //upper half of matrix
  MYREAL *gelmanmaxR = world->convergence->gelmanmeanmaxR;  //lower half of matrix

  fprintf(file,"\n\nGelman-Rubin convergence statistic\n");
  fprintf(file,"----------------------------------\n\n");
  fprintf(file,"Values close to 1.0, especially values < 1.2 are a sign of convergence of the\n");
  fprintf(file,"chains. On very short chains this statistic does not work well\n\n");
  if(world->options->gelmanpairs)
    {
      fprintf(file,"Pairwise evaluation of replicates\n[above diagonal: mean; below diagonal: maximum value]\n");
      fprintf(file,"%4li ", 1L);
      for(i=1; i < maxreplicate; i++)
	{
	  fprintf(file,"%6li ", i+1);
	}
      fprintf(file,"\n");
      for(i=0; i < maxreplicate; i++)
	{
	  for(j=0; j < maxreplicate; j++)
	    {
	      if(i==j)
		fprintf(file,"------ ");
	      else
		{
		  if(i < j) 
		    {
		      if(gelmanmeanR[i * maxreplicate + j] <=0.)
			{
			  fprintf(file,"  -N-  ");
			}
		      else
			{
			  fprintf(file," %1.3f ",gelmanmeanR[i * maxreplicate + j]);
			}
		    }
		  else
		    {
		      if(gelmanmaxR[i * maxreplicate + j] <=0.)
			{
			  fprintf(file," -N-   ");
			}
		      else
			{
			  fprintf(file,"*%1.3f ",gelmanmaxR[i * maxreplicate + j]); 
			}
		    }
		}
	    }
	  fprintf(file,"\n");
	}
      fprintf(file,"\n\n");
    }
  fprintf(file,"Overall convergence:\n");
  fprintf(file,"Mean sqrt(R)    = %f\n", world->convergence->gelmanmeanRall);
  fprintf(file,"Maximum sqrt(R) = %f\n\n\n", world->convergence->gelmanmaxRall);
}

void methods(world_fmt *world)
{
    char title[] = "Method description suggestion:";
    fprintf(world->outfile,"%s\n",title);
    fprintf(world->outfile,"Migrate was run in %s with the following run time settings:",
	    world->options->bayes_infer ? "Bayesian mode (Beerli and Felsenstein 1999, 2001; Beerli 1998,2006)" : 
	    "Maximum Likelihood mode (Beerli 1998, Beerli and Felsenstein 1999, 2001)");
    if (world->options->bayes_infer)
      {
	fprintf(world->outfile,"A total of x replicates were run. Each replicate represents a full\n");
	fprintf(world->outfile,"MCMC run of one cold chain with x-1 heated connected chains. Each replicate\n");
	fprintf(world->outfile,"visited a total of x states (Genealogy were visited with frequency of x and\n");
	fprintf(world->outfile,"each parameter was visited with this average frequency x. The run parameters\n");
	fprintf(world->outfile,"were set so that each locus sampled x observations from the visited states using\n");
	fprintf(world->outfile,"an increment of y. We used uniform prior distributions with ranges of a to b\n" );
	fprintf(world->outfile,"and c to d for Theta and M, respectively");
      }
}

void citations(world_fmt *world)
{
  char title[] = "Citation suggestions:";
  char texttitle[5][50]  = {"General method:", "Maximum likelihood:", "Bayesian inference:", "Likelihood ratio test:", "Bayes factors:"};
  char citations[9][320]  = {"[1] P. Beerli, 1998. Estimation of migration rates and population sizes in geographically structured populations., in Advances in Molecular Ecology, G. R. Carvalho, ed., vol. 306 of NATO sciences series, Series A: Life sciences, ISO Press, Amsterdam, pp. 39–53.",\
    "[2] P. Beerli and J. Felsenstein, 1999. Maximum-likelihood estimation of migration rates and effective popu- lation numbers in two populations using a coalescent approach, Genetics, 152:763–773.",\
    "[3] P. Beerli and J. Felsenstein, 2001. Maximum likelihood estimation of a migration matrix and effective population sizes in n subpopulations by using a coalescent approach, Proceedings of the National Academy of Sciences of the United States of America, 98: p. 4563–4568.",\
    "[4] P. Beerli, 2004. Effect of unsampled populations on the estimation of population sizes and migration rates between sampled populations, Molecular Ecology, 13: 827–836.",\
      "[5] P. Beerli, 2006. Comparison of Bayesian and maximum-likelihood inference of population genetic parameters. Bioinformatics 22:341-345",\
    "[6] P. Beerli, 2007. Estimation of the population scaled mutation rate from microsatellite data, Genetics, 177:1967–1968.",\
    "[7] P. Beerli, 2009. How to use migrate or why are Markov chain Monte Carlo programs difficult to use?, in Population Genetics for Animal Conservation, G. Bertorelle, M. W. Bruford, H. C. Hauffe, A. Rizzoli, and C. Vernesi, eds., vol. 17 of Conservation Biology, Cambridge University Press, Cambridge UK, pp. 42–79.",\
    "[8] P. Beerli and M. Palczewski, 2010. Unified framework to evaluate panmixia and migration direction among multiple sampling locations, Genetics, 185: 313–326."};

    fprintf(world->outfile,"%s\n",title);
//general
    fprintf(world->outfile,"%s\n",texttitle[0]);
    fprintf(world->outfile,"%s\n%s\n%s\n",citations[0],citations[6],citations[5]);
//likelihood
    fprintf(world->outfile,"%s\n",texttitle[1]);
    fprintf(world->outfile,"%s\n%s\n",citations[1],citations[2]);
//Bayesian inference
    fprintf(world->outfile,"%s\n",texttitle[2]);
    fprintf(world->outfile,"%s\n%s\n",citations[4],citations[8]);
//likelihood
    fprintf(world->outfile,"%s\n",texttitle[4]);
    fprintf(world->outfile,"%s\n",citations[7]);

}
