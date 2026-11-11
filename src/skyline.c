/*------------------------------------------------------
Maximum likelihood estimation
of migration rate  and effectice population size
using a Metropolis-Hastings Monte Carlo algorithm
-------------------------------------------------------
    variation over time routines   R O U T I N E S

    Peter Beerli 2006, Tallahassee
    beerli@fsu.edu

    Copyright 2006 Peter Beerli, Tallahassee

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

 $Id$
    */
/*! \file skyline.c 

this file contains functions that calculate the expected parameters based on 
individual time intervals. Results will be printed as a histogram over time (similar to the 
skyline plots of Rambaut, Strimmer etc) and also will print a table/file with values for each parameter

*/
#include <stdlib.h>
#include "skyline.h"
#include "random.h"
#include "tools.h"
#include "sighandler.h"
#include "world.h"
#include "speciate.h"
#ifdef PRETTY
#include "pretty.h"
#endif
#ifdef MPI
#include "migrate_mpi.h"
#else
extern int myID;
#endif
MYREAL waiting_theta(long pop, 
		     MYREAL *param, 
		     MYREAL * mig0list, 
		     long *lineages, 
		     MYREAL t0, MYREAL t1, 
		     MYREAL inv_mu_rate,
		     long numpop, world_fmt *world);
MYREAL waiting_M(long pop,long to, 
		 MYREAL *param, 
		 MYREAL * mig0list, 
		 long *lineages, 
		 MYREAL t0, MYREAL t1, 
		 MYREAL inv_mu_rate,
		 long numpop, world_fmt *world, boolean usem);

MYREAL waiting_D(long pop,long to, 
		 MYREAL *param, 
		 MYREAL * mig0list, 
		 long *lineages, 
		 MYREAL t0, MYREAL t1, 
		 MYREAL inv_mu_rate,
		 long numpop, world_fmt *world);
void print_expected_values_list(FILE *file, long locus, tetra **eventbins, MYREAL eventbinsize, long *eventbinnum, world_fmt *world);
void print_expected_values_tofile(FILE *file,  world_fmt *world);
void prepare_expected_values(world_fmt *world);
void print_expected_values_title(FILE *file, boolean progress);
//##


///
/// calculate the expected Theta for a coalescence event and returns the expected value
MYREAL waiting_theta(long pop, 
		     MYREAL *param, 
		     MYREAL * mig0list, 
		     long *lineages, 
		     MYREAL t0, MYREAL t1, 
		     MYREAL inv_mu_rate,
		     long numpop, world_fmt *world)
{
  long i;
  long kpop;
  long line;
  MYREAL interval = t1-t0;
  MYREAL expected = 0.0;
  //MYREAL specw=0.0;
  for(i=0;i<numpop; i++)
    {
      line = lineages[i];
      if(line>0)
	{
	  //specw += wait_D(i, t0, t1, lineages,world);
	  expected += wait_D(i, t0, t1, lineages,world)/interval;
	  expected += line * (line-1)/param[i] +  mig0list[i] * line;
	}
    }
  line = lineages[pop];
  kpop = line * ( line -1);
  expected = 1./interval - (expected - kpop / param[pop])*inv_mu_rate;
  //expected = (specw - (expected - kpop / param[pop]))*inv_mu_rate;
  // 
  if(expected > EPSILON)
    //if(interval > EPSILON)
    //return 0.5*exp(-specw + interval * expected) * kpop * interval * kpop * interval; 
    return kpop/expected;
  else
    return -90.;
}

///
/// calculate the expected M for a coalescence event and returns the expected value
MYREAL waiting_M(long pop,long to, 
		 MYREAL *param, 
		 MYREAL * mig0list, 
		 long *lineages, 
		 MYREAL t0, MYREAL t1, 
		 MYREAL inv_mu_rate,
		 long numpop, world_fmt *world, boolean usem)
{
  long i;
  long kpop;
  MYREAL expected = 0.0;
  MYREAL line;
  MYREAL interval = t1 - t0;
  MYREAL pr = usem ? 1.0 : param[to];
  for(i=0;i<numpop; i++)
    {
      line = lineages[i];
      if(line>0)
	{
	  expected += wait_D(i, t0, t1, lineages,world)/interval;
	  expected += line * (line-1) /param[i] +  mig0list[i] * line;
	}
    }
  kpop = lineages[to];
  expected = 1./interval - (expected - kpop * param[pop]/pr) * inv_mu_rate;
  if(kpop > EPSILON)
    return pr * expected / kpop;
  else
    return -9.;
}
///
/// calculate the expected D for an event and returns the expected value
MYREAL waiting_D(long pop,long to, 
		 MYREAL *param, 
		 MYREAL * mig0list, 
		 long *lineages, 
		 MYREAL t0, MYREAL t1, 
		 MYREAL inv_mu_rate,
		 long numpop, world_fmt *world)
{
  (void) pop;
  (void) to;
  (void) inv_mu_rate;;
  long i;
  //long kpop;
  MYREAL expected = 0.0;
  MYREAL line;
  MYREAL specw;//, specwpop;
  MYREAL interval = t1 - t0;
  for(i=0;i<numpop; i++)
    {
      line = lineages[i];
      if(line>0)
	{
	  specw = wait_D(i, t0, t1, lineages,world)/interval;
	  //if (i==to)
	  //  specwpop = specw;
	  expected += line * (line-1)/param[i] +  mig0list[i] * line + specw * line;
	}
    }
  //kpop = lineages[to];
  //expected = 1./interval - (expected - specwpop) * inv_mu_rate;
  return t1;
  //if(kpop > EPSILON)
  //  return expected / kpop;
  //else
  //  return -900.;
}

/* ---- empirical skyline (skyline=YES): generalized skyline estimator ----
   For every time bin and parameter the sampled genealogies contribute
     exposure  E = sum over intervals of rate-per-unit-parameter x overlap
     events    n = number of events of that kind inside the bin
   Theta_p: rate k_p(k_p-1)/Theta_p, so E = sum k_p(k_p-1) dt and Theta = E/n;
   M_qp:    rate k_p M (or k_p xNm/Theta_p), so E = sum k_p dt (/Theta_p) and
            M = n/E.
   This is the maximum-likelihood value of a parameter that is constant within
   the bin, given the genealogies (the classic skyline of Pybus et al. 2000
   pooled over a bin and over genealogies). The earlier version averaged the
   per-interval values k(k-1)dt over the bin weighted by the time the interval
   spent in the bin, which favours long intervals (about 2 Theta for constant
   size), and dropped values above 10x the running mean.
   Raw slots while sampling: [0] exposure (compensated sum, correction in [3]),
   [1] events, [4] genealogies that reached the bin, [6] id of the last such
   genealogy; [0]-[5] are summed when MPI workers send their bins.
   prepare_expected_values() turns them into the reported values. */
static long skyline_genealogy_id = 0;

/// called once for each sampled genealogy before its intervals are recorded
void skyline_next_genealogy(void)
{
  skyline_genealogy_id++;
}

static void skyline_ensure_bin(tetra **eventbins, long *eventbinnum, long pop, long bin)
{
  long i, j, allocsize;
  if (bin < eventbinnum[pop])
    return;
  allocsize = bin + 10;
  eventbins[pop] = (tetra *) myrealloc(eventbins[pop], (size_t) allocsize * sizeof(tetra));
  for (i = eventbinnum[pop]; i < allocsize; i++)
    for (j = 0; j < 9; j++)
      eventbins[pop][i][j] = 0.0f;
  eventbinnum[pop] = allocsize;
}

static void skyline_add_exposure(tetra **eventbins, long *eventbinnum, long pop,
                                 double t0, double t1, double binsize, double rate)
{
  long bin, first, last;
  double lo, hi, x, y, t;
  const float gid = (float) skyline_genealogy_id;
  if (rate <= 0.0 || t1 <= t0)
    return;
  first = (long) (t0 / binsize);
  last = (long) (t1 / binsize);
  skyline_ensure_bin(eventbins, eventbinnum, pop, last);
  for (bin = first; bin <= last; bin++)
    {
      lo = MAX(t0, binsize * (double) bin);
      hi = MIN(t1, binsize * (double) (bin + 1));
      if (hi <= lo)
        continue;
      x = rate * (hi - lo);
      /* compensated (Kahan) summation: the bins are floats */
      y = x - (double) eventbins[pop][bin][3];
      t = (double) eventbins[pop][bin][0] + y;
      eventbins[pop][bin][3] = (float) ((t - (double) eventbins[pop][bin][0]) - y);
      eventbins[pop][bin][0] = (float) t;
      if (eventbins[pop][bin][6] != gid)
        {
          eventbins[pop][bin][6] = gid;
          eventbins[pop][bin][4] += 1.0f;
        }
    }
}

static void skyline_add_event(tetra **eventbins, long *eventbinnum, long pop,
                              double age, double binsize)
{
  long bin = (long) (age / binsize);
  skyline_ensure_bin(eventbins, eventbinnum, pop, bin);
  eventbins[pop][bin][1] += 1.0f;
}

/// records one interval of a sampled genealogy (from age-interval to age,
/// lineages during the interval) and the event that ends it
void calculate_expected_values(tetra **eventbins, 
			       long *eventbinnum, 
			       MYREAL eventinterval, 
			       MYREAL interval, 
			       MYREAL age,
			       char type,
			       long from, 
			       long to, 
			       long * lineages, 
			       long numpop, 
			       world_fmt *world)
{
  const boolean usem = world->options->usem;
  const double inv_mu_rate = 1.0 / (double) world->options->mu_rates[world->locus];
  const double binsize = (double) eventinterval;
  double t0 = (double) (age - interval);
  const double t1 = (double) age;
  long p, q, m;
  double k, pr;
  if (age > 1000.0)
    return;
  if (t0 < 0.0)
    t0 = 0.0;
  for (p = 0; p < numpop; p++)
    {
      k = (double) lineages[p];
      if (k <= 0.0)
        continue;
      if (k > 1.0 && world->bayes->map[p][1] != INVALID)
        skyline_add_exposure(eventbins, eventbinnum, p, t0, t1, binsize, k * (k - 1.0) * inv_mu_rate);
      for (q = 0; q < numpop; q++)
        {
          if (q == p)
            continue;
          m = mm2m(q, p, numpop);   /* lineages in p move to q (backward) */
          if (world->bayes->map[m][1] == INVALID)
            continue;
          pr = usem ? 1.0 : (double) world->param0[p];
          if (pr > 0.0)
            skyline_add_exposure(eventbins, eventbinnum, m, t0, t1, binsize, k * inv_mu_rate / pr);
        }
    }
  if (from == to)
    {
      if (type != 't' && type != 'T' && type != 'd' && type != 'D' && type != 'r'
          && world->bayes->map[to][1] != INVALID)
        skyline_add_event(eventbins, eventbinnum, to, t1, binsize);
    }
  else if (type != 'd' && type != 'D' && type != 't' && type != 'T')
    {
      m = mm2m(from, to, numpop);
      if (world->bayes->map[m][1] != INVALID)
        skyline_add_event(eventbins, eventbinnum, m, t1, binsize);
    }
}


///
/// set up the skyline plot histogram containers, only when also the migration histograms are recorded
/// this function needs to be called AFTER the setup_mighist() function
/// the main vehicle is the evenbins that contain 9 basic bins for each time bin, each parameter
/// 0) Time
/// 1) Parameter value (average)
/// 2) standard deviation  (average standard deviation)
/// 3) standard deviation of weights (average weight)
/// 4) Number of values put into bin (average count)
/// 5) autocorrelation (average autocorrelation coeff)
/// 6)
/// 7)
/// 8)
void setup_expected_events (world_fmt * world, option_fmt * options)
{
  (void) options;
  long locus, i, j,z;
    long allocsize = 1;
    long npall = world->numparam; // MPI pack_skyline_buffer() sends numparam slots
    MYREAL  binsize = (double) world->options->eventbinsize; // in mutational units (=time scale)
    if (world->options->mighist && world->options->skyline)
    {
      for (locus = 0; locus < world->loci; locus++)
        {
	  world->mighistloci[locus].eventbinsize = binsize;
	  world->mighistloci[locus].eventbins = (tetra **) mycalloc (npall, sizeof (tetra *));
	  world->mighistloci[locus].eventbinnum = (long *) mycalloc (npall, sizeof (long));
	  for (i = 0; i < npall; i++)
            {
	      world->mighistloci[locus].eventbinnum[i] = allocsize;
	      world->mighistloci[locus].eventbins[i] = (tetra *) mycalloc (allocsize, sizeof (tetra));
	      for(j=0;j<allocsize;j++)
		{
		  for(z=0;z<9;z++)
		    world->mighistloci[locus].eventbins[i][j][z] = 0.;
		}
            }
        }
    }
}

///
/// Destroy the skyline plot histogram container
void
destroy_expected_events (world_fmt * world)
{
    long locus, i;
    long npall = world->numparam; // MPI pack_skyline_buffer() sends numparam slots
    if (world->options->mighist && world->options->skyline)
    {
      for (locus = 0; locus < world->loci; locus++)
        {	  
	  for (i = 0; i < npall; i++)
            {
	      myfree(world->mighistloci[locus].eventbins[i]);
            }
	  myfree(world->mighistloci[locus].eventbins);
	  myfree(world->mighistloci[locus].eventbinnum);
        }
    }
}


void print_expected_values_list(FILE *file, long locus, tetra **eventbins, MYREAL eventbinsize, long *eventbinnum, world_fmt *world)
{
  long i;
  long pop;
  long frompop;
  long topop;
  long numpop  = world->numpop;
  long numpop2 = world->numpop2;
  long npall   = numpop2 + world->species_model_size * 2 + world->bayes->mu;
  MYREAL age;
  
  for(pop = 0; pop < npall; pop++)
    {  
      age = 0.;
      if(world->bayes->map[pop][1] == INVALID)
	continue;	
      if(pop < numpop)
	{
	  fprintf(file,"\nLocus: %li   Parameter: %s_%li\n", locus+1, "Theta",pop+1);  
	}
      else
	{
	    if(pop<world->numpop2)
	      {
		m2mm(pop,numpop,&frompop,&topop);
		fprintf(file,"\nLocus: %li   Parameter: %s_(%li,%li)\n", locus+1, "M", frompop+1, topop+1);  
	      }
	    else
              {
                if(pop==world->numpop2 && world->bayes->mu)
                  continue;
                d2mm(pop,world,&frompop,&topop);
                fprintf(file,"\nLocus: %li   Parameter: %s_(%li,%li)\n", locus+1, "D", frompop+1, topop+1);
              }
	}
      fprintf(file,"Time        Parameter       Frequency of visit\n");
      fprintf(file,"----------------------------------------------\n");
      for(i = 0; i < eventbinnum[pop]; i++)
	{
	  age += eventbinsize;
	  if(eventbins[pop][i][1] < (float) SMALL_VALUE)
	    continue;
	  if(eventbins[pop][i][0] < 0)
	    error("nono");
	  fprintf(file,"%10.10f  %10.10f     %10.10f\n", age, (double) eventbins[pop][i][0], (double) eventbins[pop][i][1]);
	}
    } 
}

#define SOME_ELEMENTS 100

void print_expected_values_tofile(FILE *file,  world_fmt *world)
{
  long i;
  long pop;
  long frompop;
  long topop;
  long locus;
  long sumloc;
  long numpop = world->numpop;
  //long numpop2 = world->numpop2;
  long npall = world->numpop2+world->species_model_size * 2 + world->bayes->mu;
  MYREAL age;
  tetra **eventbins;
  long *eventbinnum;
  MYREAL eventbinsize;

  fprintf(file,"# Raw record of the skyline histogram for all parameters and all loci\n");  
  fprintf(file,"# The time interval is set to %f\n", (double) world->options->eventbinsize);  
  fprintf(file,"# produced by the program %s (http://popgen.csit.fsu.edu/migrate.hml)\n",
	  MIGRATEVERSION);  
  fprintf(file,"# written by Peter Beerli 2006, Tallahassee,\n");
  fprintf(file,"# if you have problems with this file please email to beerli@fsu.edu\n");  
  fprintf(file,"#\n");
  fprintf(file,"# Order of the parameters:\n");
  fprintf(file,"# Parameter-number Parameter\n");
  for(pop=0;pop<npall;pop++)
    {
      if(world->bayes->map[pop][1] == INVALID)
	continue;	
      if(pop < numpop)
	{
	  fprintf(file,"# %6li    %s_%li\n", pop+1, "Theta",pop+1);  
	}
      else
	{
	  if(pop<world->numpop2)
	    {
	      m2mm(pop,numpop,&frompop,&topop);
	      fprintf(file,"# %6li    %s_(%li,%li)\n", pop+1, (world->options->usem ? "M" : "xNm"), frompop+1, topop+1);  
	    }
	  else
	    {
	      if(pop==world->numpop2 && world->bayes->mu)
		continue;
	      d2mm(pop,world,&frompop,&topop);
	      fprintf(file,"# %6li    %s_(%li,%li)\n", pop+1, "D", frompop+1, topop+1);
	      pop++; //this advances over the standard deviation
	    }
	}
    }
  fprintf(file,"#\n#----------------------------------------------------------------------------\n");
  fprintf(file,"# Locus Parameter-number Bin Age Parameter-value Events-per-genealogy\n");
  fprintf(file,"#        Standard-error Exposure-per-genealogy Events Autocorrelation(unused, 0)\n");
  fprintf(file,"#----------------------------------------------------------------------------\n");
  fprintf(file,"# Generalized skyline: in every time bin (Age = bin middle) the value is the\n");
  fprintf(file,"# maximum-likelihood value of a parameter that is constant in the bin, given\n");
  fprintf(file,"# the sampled genealogies: Theta = sum k(k-1) dt / number of coalescences,\n");
  fprintf(file,"# M = number of migrations / sum k dt (k lineages, dt time in the bin).\n");
  fprintf(file,"# Standard error = value / sqrt(events per genealogy). A bin without\n");
  fprintf(file,"# coalescences (Theta) or without lineages has all values 0; the list ends\n");
  fprintf(file,"# with the last bin that has an estimate. No skyline for divergence.\n");
  if(world->loci>1)
    {
      sumloc =1;
      fprintf(file,"# Locus %li is sum over all loci, when there are more than 1 locus\n", world->loci+1);
    }
  else
    {
      sumloc = 0;
    }

  for(locus=0; locus < world->loci + sumloc; locus++)
    {
      if(!world->data->skiploci[locus])
	{
	  eventbins =  world->mighistloci[locus].eventbins;
	  eventbinnum = world->mighistloci[locus].eventbinnum;
	  eventbinsize = world->mighistloci[locus].eventbinsize;
	  for(pop = 0; pop < npall; pop++)
	    {  
	      if(world->bayes->map[pop][1] == INVALID)
		continue;
	      if (pop>world->numpop2)
		{
		  if ((pop-world->numpop2-world->bayes->mu)%2 == 1)
		    continue;
		}
	      age = eventbinsize / 2.;
	      for(i = 0; i < eventbinnum[pop]; i++)
		{
		  fprintf(file,"%li\t%li\t%li\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%li\t%10.10f\n", 
			  locus+1, pop+1, i+1, age, (double) eventbins[pop][i][0], (double) eventbins[pop][i][1], (double) eventbins[pop][i][2],(long) eventbins[pop][i][4], (double) eventbins[pop][i][5]);
		  age += eventbinsize;
		}
	    } 
	}
    }
}

/// genealogies sampled for a locus: the largest genealogy count of any bin
/// (every sampled genealogy covers the first bin of every Theta)
static float skyline_genealogies(tetra **eventbins, long *eventbinnum, long npall)
{
  long pop, i;
  float g = 0.0f;
  for (pop = 0; pop < npall; pop++)
    for (i = 0; i < eventbinnum[pop]; i++)
      if (eventbins[pop][i][4] > g)
        g = eventbins[pop][i][4];
  return g;
}

/// raw sums of one bin -> reported values:
/// [0] value, [1] events per genealogy, [2] standard error, [3] exposure per
/// genealogy, [4] number of events, [5] 0 (no autocorrelation)
static void skyline_finish_bin(tetra *bin, double exposure, double events,
                               double events_per_genealogy, double exposure_per_genealogy,
                               boolean is_theta)
{
  long j;
  double value, se;
  for (j = 0; j < 9; j++)
    (*bin)[j] = 0.0f;
  if (exposure <= 0.0 || exposure_per_genealogy <= 0.0)
    return;
  if (is_theta)
    {
      if (events <= 0.0)
        return;   /* no coalescence: no estimate */
      value = exposure / events;
      se = value / sqrt(events_per_genealogy);
    }
  else
    {
      value = events / exposure;
      /* Poisson error of the events of one genealogy; at least one event */
      se = sqrt(MAX(events_per_genealogy, 1.0)) / exposure_per_genealogy;
    }
  (*bin)[0] = (float) value;
  (*bin)[1] = (float) MAX(events_per_genealogy, (double) FLT_MIN);
  (*bin)[2] = (float) se;
  (*bin)[3] = (float) exposure_per_genealogy;
  (*bin)[4] = (float) events;
}

static boolean skyline_reported(world_fmt *world, long pop)
{
  if (world->bayes->map[pop][1] == INVALID)
    return FALSE;
  return pop < world->numpop2;   /* Theta and M; no empirical skyline for D */
}

void prepare_expected_values(world_fmt *world)
{
  long locus, pop, i, n;
  long npall = world->numpop2 + world->bayes->mu + world->species_model_size * 2;
  long allmax = 0, used = 0;
  tetra **eventbins;
  tetra **all = NULL;
  long *eventbinnum;
  float g;
  double **sumE = NULL, **sumN = NULL, **sumNpg = NULL, **sumEpg = NULL;

  /* longest skyline over loci */
  for (locus = 0; locus < world->loci; locus++)
    {
      if (world->data->skiploci[locus])
        continue;
      for (pop = 0; pop < npall; pop++)
        if (world->mighistloci[locus].eventbinnum[pop] > allmax)
          allmax = world->mighistloci[locus].eventbinnum[pop];
    }
  if (world->loci > 1 && allmax > 0)
    {
      sumE = (double **) mycalloc(npall, sizeof(double *));
      sumN = (double **) mycalloc(npall, sizeof(double *));
      sumNpg = (double **) mycalloc(npall, sizeof(double *));
      sumEpg = (double **) mycalloc(npall, sizeof(double *));
      for (pop = 0; pop < npall; pop++)
        {
          sumE[pop] = (double *) mycalloc(allmax, sizeof(double));
          sumN[pop] = (double *) mycalloc(allmax, sizeof(double));
          sumNpg[pop] = (double *) mycalloc(allmax, sizeof(double));
          sumEpg[pop] = (double *) mycalloc(allmax, sizeof(double));
        }
    }
  /* per locus: pool the raw sums for the combination, then finish the locus */
  for (locus = 0; locus < world->loci; locus++)
    {
      if (world->data->skiploci[locus])
        continue;
      eventbins = world->mighistloci[locus].eventbins;
      eventbinnum = world->mighistloci[locus].eventbinnum;
      g = skyline_genealogies(eventbins, eventbinnum, npall);
      if (g <= 0.0f)
        continue;
      used++;
      for (pop = 0; pop < npall; pop++)
        {
          if (!skyline_reported(world, pop))
            {
              eventbinnum[pop] = 0;
              continue;
            }
          n = 0;
          for (i = 0; i < eventbinnum[pop]; i++)
            {
              const double e = (double) eventbins[pop][i][0];
              const double k = (double) eventbins[pop][i][1];
              if (sumE != NULL)
                {
                  sumE[pop][i] += e;
                  sumN[pop][i] += k;
                  sumNpg[pop][i] += k / (double) g;
                  sumEpg[pop][i] += e / (double) g;
                }
              skyline_finish_bin(&eventbins[pop][i], e, k, k / (double) g, e / (double) g,
                                 pop < world->numpop);
              if (eventbins[pop][i][1] > 0.0f)
                n = i + 1;
            }
          eventbinnum[pop] = n;   /* drop the empty tail */
        }
    }
  /* all loci: the pooled estimate (loci share the parameter); the events per
     genealogy add over loci, so the error is that of one multi-locus sample */
  if (sumE != NULL)
    {
      if (world->mighistloci[world->loci].eventbins == NULL)
        world->mighistloci[world->loci].eventbins = (tetra **) mycalloc(npall, sizeof(tetra *));
      if (world->mighistloci[world->loci].eventbinnum == NULL)
        world->mighistloci[world->loci].eventbinnum = (long *) mycalloc(npall, sizeof(long));
      world->mighistloci[world->loci].eventbinsize = world->mighistloci[0].eventbinsize;
      all = world->mighistloci[world->loci].eventbins;
      for (pop = 0; pop < npall; pop++)
        {
          all[pop] = (tetra *) mycalloc(allmax, sizeof(tetra));
          n = 0;
          if (skyline_reported(world, pop) && used > 0)
            for (i = 0; i < allmax; i++)
              {
                skyline_finish_bin(&all[pop][i], sumE[pop][i], sumN[pop][i], sumNpg[pop][i],
                                   sumEpg[pop][i], pop < world->numpop);
                if (all[pop][i][1] > 0.0f)
                  n = i + 1;
              }
          world->mighistloci[world->loci].eventbinnum[pop] = n;
          myfree(sumE[pop]);
          myfree(sumN[pop]);
          myfree(sumNpg[pop]);
          myfree(sumEpg[pop]);
        }
      myfree(sumE);
      myfree(sumN);
      myfree(sumNpg);
      myfree(sumEpg);
    }
}
  
void print_expected_values_title(FILE *file, boolean progress)
{
  if(progress)
    {      
      fprintf(file,"\n\nParameter changes over time\n");
      fprintf(file,"---------------------------\nSEE IN PDF FILE AND SKYLINE FILE\n");
    }
}

void print_expected_values(world_fmt * world, option_fmt *options)
{
  (void) options;
  long locus;
  long sumloc = (world->loci > 1) ? 1 : 0;
  long npall = world->numpop2 + world->bayes->mu + world->species_model_size * 2;
  if(world->options->skyline)
    {
	  // prepare skyline histogram for printing
	  prepare_expected_values(world);
	  // print skyline to file
	  print_expected_values_tofile(world->skylinefile, world);
      // print title to screen
      print_expected_values_title(stdout, world->options->progress);
      // print title to ascii-outfile
      print_expected_values_title(world->outfile, TRUE);
      // print content to stdout and to ascii - outfile
      for(locus=0; locus < world->loci+sumloc; locus++)
      	{
	  if(world->options->verbose)
	    print_expected_values_list(stdout, locus, world->mighistloci[locus].eventbins, 
	    			       world->mighistloci[locus].eventbinsize, 
	    			       world->mighistloci[locus].eventbinnum, world);
	  print_expected_values_list(world->outfile, locus, world->mighistloci[locus].eventbins, 
	    			     world->mighistloci[locus].eventbinsize, 
	    			     world->mighistloci[locus].eventbinnum, world);
	}
#ifdef PRETTY
	  pdf_skyline_histogram(world->loci, npall,  world, FALSE);
	  if(world->options->bayes_infer)
	    pdf_skyline_histogram(world->loci, npall,  world, TRUE);
#endif
    }
}


void debug_skyline(world_fmt *world, char text[])
{
  long i;
  long pop;
  long locus;
  tetra ** eventbins;
  long *eventbinnum;
  MYREAL eventbinsize;
  MYREAL age;
  FILE *file = stdout;
  long numpop2 = world->numpop2;
  fprintf(file,"#%i -------- %s ----------\n",myID,text);
  fprintf(file,"#\n#----------------------------------------------------------------------------\n");
  fprintf(file,"# Locus Parameter-number Bin Age Parameter-value(*) Parameter-Frequency(*)\n");
  fprintf(file,"#----------------------------------------------------------------------------\n");
  fprintf(file,"# (*) values with -1 were NEVER visited\n");
  for(locus=0; locus < world->loci; locus++)
    {
      if(!world->data->skiploci[locus])
	{
	  eventbins =  world->mighistloci[locus].eventbins;
	  eventbinnum = world->mighistloci[locus].eventbinnum;
	  eventbinsize = world->mighistloci[locus].eventbinsize;
	  for(pop = 0; pop < numpop2; pop++)
	    {  
	      age = eventbinsize / 2.;
	      for(i = 0; i < eventbinnum[pop]; i++)
		{
		  fprintf(file,"%li %li %li %10.10f %10.10f %10.10f %10.10f\n", 
			  locus+1, pop+1, i+1, age, (double) eventbins[pop][i][0], (double) eventbins[pop][i][1], (double) (eventbins[pop][i][0]/ eventbins[pop][i][1]));
		  age += eventbinsize;
		}
	    } 
	}
    }
  fprintf(file,"#%i >>>>>>>>>>>>>>> %s <<<<<<<<<<<<<<<<<END\n",myID,text);
}

