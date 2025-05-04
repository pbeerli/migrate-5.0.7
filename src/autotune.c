/*
 Bayesian inference or maximum likelihood inference of structured population models
 
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
 
 
 */

#include "migration.h"
#include "tools.h"
#include "sighandler.h"
#include "migrate_mpi.h"
#include "random.h"
#include "bayes.h"
#include "mcmc.h"
#include "reporter.h"
#include "seqerror.h"
#include "haplotype.h"
#include "skyparam.h"

extern void print_menu_equilib (world_fmt * world);
extern int myID;
// functions
void autotune_proposal(world_fmt *world, long which);
void present_burnin_info(world_fmt *world, MYREAL ess, MYREAL acceptance, MYREAL var, MYREAL oldvar,  long step);
MYREAL mean_acceptance_rate(world_fmt * world);
MYINLINE boolean auto_stop_burnin(world_fmt *world,  long step,  long stop, MYREAL * var, MYREAL *autocorrelation, MYREAL * effective_sample);
long  expected_end_burnin(world_fmt *world, MYREAL percent, long starttime, char *text);
void burnin_progress(MYREAL percent);
void burnin_bayes(world_fmt * world);
void burnin_chain (world_fmt * world);


// functions implementation
void autotune_proposal(world_fmt *world, long which)
{
    MYREAL ratio;
    worldoption_fmt *wopt = world->options;
    bayes_fmt *bayes = world->bayes;
    MYREAL delta = bayes->delta[which];
    long space = MAX(10, (long) world->numpop2+1);
    MYREAL ma = bayes->maxparam[which];
    MYREAL mi = bayes->minparam[which];
    MYREAL mindelta = (ma-mi)/200.;

    if(world->in_burnin && world->cold && !world->options->prioralone)
    {
        // formula: delta = (pR - 1) * (min-max)
        // for pR=0.44 and min=0 and max=0.1 --> delta= -0.66 * -0.1 = 0.066
        // pR = delta/(min-max) + 1
        // Binom(n,k) pR^k (1-pR)^(n-k)
        // n=1:k=1: delta/(min-max)+1
        //      delta/(-10*delta)+1 ==> 1-1/10 = 0.9
        // n=1:k=0: (1-delta/(min-max)-1)) ==> -1/10 = -0.1
        if(wopt->has_autotune)
        {
	  if(bayes->trials[which] % space == 0)
	    {
	      ratio = bayes->accept[which]/(1.0+bayes->trials[which]);
	      if(wopt->autotune < ratio)
		{
		  delta *= 1.0101; /*0.99^(r/(r-1))*/
		}
	      else
		{
		  delta *= 0.990;
		}
	      if(delta > ma)
                delta = ma;
	      if(delta < mindelta)
		delta = mindelta;
	      bayes->delta[which] = delta;
	    }
	}
    }
}

void present_burnin_info(world_fmt *world, MYREAL ess, MYREAL acceptance, MYREAL var, MYREAL oldvar,  long step)
{
   long z=0;
#ifdef MPI
    char p[LINESIZE];
    char p1[LINESIZE];
    long bufsize;
    if(myID!=MASTER)
    {
      bufsize = sprintf(p,"%li %f %f %f %f %li\n",world->locus, ess, acceptance, var,oldvar,step);
        sprintf(p1,"B%li",bufsize);
        MYMPISEND (p1, SMALLBUFSIZE, MPI_CHAR, (MYINT) MASTER, (MYINT) myID+BURNTAG, comm_world);
        MYMPISEND (p, bufsize, MPI_CHAR, (MYINT) MASTER, (MYINT) myID+BURNTAG, comm_world);
    }
#endif
    z = world->burnin_z;
    if(world->burnin_stops_alloc <= z)
      {
	world->burnin_stops_alloc += z;
	world->burnin_stops = (burnin_record_fmt *) myrealloc(world->burnin_stops, world->burnin_stops_alloc * sizeof(burnin_record_fmt));  
      }
    world->burnin_stops[z].locus     =  world->locus;
    world->burnin_stops[z].replicate =  world->rep;
    world->burnin_stops[z].stopstep  = step;
    world->burnin_stops[z].ess       = ess;
    world->burnin_stops[z].accept  = acceptance;
    world->burnin_stops[z].variance  = var;
    world->burnin_stops[z].oldvariance = oldvar;
    world->burnin_stops[z].worker    = myID;
    world->burnin_z++;
}

MYREAL mean_acceptance_rate(world_fmt * world)
{
   long which;
   long i;
  bayes_fmt * bayes = world->bayes;
  long nng = world->numpop2 + ( long) bayes->mu + 1 + world->species_model_size * 2 + world->grownum; 
  MYREAL sum = 0.0;
  long count = 0;
  if (!world->options->bayes_infer)
    return -1.0;

  for (i = 0; i < nng; i++)
    {
      if(shortcut(i,world,&which))
	continue;
      else
	{
	  sum += bayes->accept[which]/(1.0+bayes->trials[which]);
	  count += 1;
	}
    }
  return sum/count;
}

MYINLINE boolean auto_stop_burnin(world_fmt *world,  long step,  long stop, MYREAL * var, MYREAL *autocorrelation, MYREAL * effective_sample)
{
    char autostop = world->options->burnin_autostop;
    const  long delta = ((stop >= 10000) ? 1000 : (stop / 10));
    const  long nn = world->numpop2 + ( long) world->bayes->mu + 1 + world->species_model_size * 2 + world->grownum;
    MYREAL oldvar= *var;
    MYREAL acceptance;
    MYREAL ess;
    boolean done = FALSE;
    boolean acceptanceOK = FALSE;
    boolean essOK = FALSE;
    boolean vardiffOK = FALSE;
    single_chain_var (world, (unsigned long) step, var, autocorrelation, effective_sample);
    if (step > delta && oldvar>0.0)
      vardiffOK = (fabs(*var/oldvar - 1.0) < world->varheat);
    else
      vardiffOK = FALSE;
    essOK = max_ess(effective_sample,nn,world->essminimum, &ess);
    acceptance = mean_acceptance_rate(world);
    acceptanceOK = acceptance > world->options->autotune - 0.05 && acceptance < world->options->autotune + 0.05;  
    switch(autostop)
    {
    case 'a':
      done = vardiffOK;
      break;
    case 't':
      done = acceptanceOK;
      break;
    case 'e':
      done = essOK;
      break;
    case ' ':
    default:
      break;
    }
    if(done || step==stop)
    {
      if(world->cold && world->options->burnin_autostop !=' ')
	present_burnin_info(world,ess, acceptance, *var,oldvar,step);
      return done;
    }
    return FALSE;
}

long  expected_end_burnin(world_fmt *world, MYREAL percent, long starttime, char *text)
{
    char nowstr[STRSIZE]; // will hold seconds since epoch
    get_time (nowstr, "%s");
    long mytime = atol(nowstr);
    if (starttime==0)
      return mytime;
    else
      {
	mytime -= starttime;
	//printf("burningtime: %li %li\n",starttime,mytime);
	mytime = (long) ((100 * mytime / (double) percent) * world->options->heated_chains);
	get_time (nowstr, "%H:%M:%S");
#ifdef MPI
	FPRINTF(stdout, "[%3i] %s   %s in %li seconds\n",myID,nowstr, text, mytime); 
#else
	FPRINTF(stdout, "%s   %s in %li seconds\n",nowstr, text, mytime); 
#endif
      }
    return mytime;
}

void burnin_progress(MYREAL percent)
{
    char nowstr[STRSIZE]; // will hold time of day
    get_time (nowstr, "%H:%M:%S");
#ifdef MPI
    FPRINTF(stdout, "[%3i] %s   Burn-in %3.1f%% complete\n",myID, nowstr, percent); 
#else
    FPRINTF(stdout, "%s   Burn-in %3.1f%% complete\n",nowstr, percent); 
#endif
}


void burnin_bayes(world_fmt * world)
{
#ifndef MPI
  const boolean progress = world->options->progress && world->cold;
  long starttime;
  boolean reportdone=FALSE;
#endif
  const  long stop = world->options->burn_in * world->increment;
  const  long nn = world->numpop2 + ( long) world->bayes->mu+ world->species_model_size * 2 + world->grownum;
   long delta = ((stop > 100000) ? (stop / 10) : (stop/3));
  MYREAL * autocorrelation;
  MYREAL * effective_sample;
  //MYREAL * acceptances;
  MYREAL var= (MYREAL) HUGE;
   long step;
  int choice;
  boolean success=FALSE;
  boolean done=FALSE;
  if(stop==0)
    return;
  if(delta==0)
    delta=1;

  autocorrelation  = (MYREAL *) mycalloc((size_t) (3 + 3*nn),sizeof(MYREAL));
  effective_sample = autocorrelation + nn + 1;
  //acceptances      = effective_sample + nn + 1;

  single_chain_var (NULL, 0, &var, NULL, NULL);
  double r = -1.0; 
  double *choices = world->options->choices; //cumulative probability of choice
  // between treeupdate, haplotypeing, and parameter updates
  //update ratio for trees versus parameters
#ifndef MPI
  starttime = expected_end_burnin(world,0.0, 0, " " );
#endif
  for (step = 1; step <= stop; step++)
    {
#ifndef MPI
      // too much writing for MPI therefore we do not report progress on burnin
      if (!reportdone && progress && step % delta == 0)
	{
	  reportdone=TRUE;
	  expected_end_burnin(world,100*(step-1)/stop * world->options->heated_chains, starttime, "Burn-in complete ");
	}
#endif
      r = RANDUM();
      choice = 0;
      while(r>choices[choice])
	{
	  choice++;
	}
      switch(choice)
	{
	case 0:
	  success = (boolean) tree_update (world, 0, NOASSIGNING);
	  break;
	case 1:
	  world->bayesaccept = bayes_update (world);
	  world->param_like = world->bayes->oldval;
	  break;
	case 2:
	  success = (boolean) swap_haplotypes(world);
	  break;
	case 3:
	  // time update
	  world->param_like = bayes_update_timeparam(world, &success);
	  break;
	case 4:
	  // assignment of individuals
	  success = (boolean) tree_update (world, 0, ASSIGNING);
	  break;
	case 5:
	  change_freq(world);
	  break;
	default:
	  error("failure in choosing among updates -- Bayes");
	}
	        if((step % delta) == 0)
        {
	  done = auto_stop_burnin(world, step, stop, &var, autocorrelation, effective_sample);
          }
	if(done)
	  break;
        world->bayes->count = 0;
	world->bayes->hypercount=0;
    }
  world->treesdone += stop;
  myfree(autocorrelation);
}

void burnin_chain (world_fmt * world)
{
  char autostop = world->options->burnin_autostop;
  const  long nng = world->numpop2 + ( long) world->bayes->mu + 1 + world->species_model_size * 2 + world->grownum; //to set zero including the genealogy
  const boolean treeprint = (boolean) world->options->treeprint;
     long z=0;
    world->burnin_z=0;        
    world->options->treeprint = myNONE;
    
    world->in_burnin = TRUE;
    if (world->cold)
    {
        print_menu_equilib (world);
    }
    burnin_bayes(world);
    memset(world->bayes->accept,0,sizeof(long) * (size_t) nng);
    memset(world->bayes->trials,0,sizeof(long) * (size_t) nng);
    if (world->cold)
      {
	z = world->burnin_z-1;
	switch(autostop)
	  {
	  case 'a':
	    FPRINTF(stdout,"[%3i]            stopped at step %li with var-ratio=%.2f/%.2f=%.3f (min(ESS)=%.2f, avg(acceptance)=%.2f)\n",
		    myID,world->burnin_stops[z].stopstep,world->burnin_stops[z].variance,world->burnin_stops[z].oldvariance,
		    world->burnin_stops[z].variance/world->burnin_stops[z].oldvariance,
		    world->burnin_stops[z].ess, world->burnin_stops[z].accept);
	    break;
	  case 't':
	    FPRINTF(stdout,"[%3i]            stopped at step %li with avg(acceptance)=%.2f (var-ratio=%.2f/%.2f=%.3f, min(ESS)=%.2f)\n",
		    myID,world->burnin_stops[z].stopstep,
		    world->burnin_stops[z].accept,
		    world->burnin_stops[z].variance,world->burnin_stops[z].oldvariance,
		    world->burnin_stops[z].variance/world->burnin_stops[z].oldvariance,
		    world->burnin_stops[z].ess);
	    break;
	  case 'e':
	    FPRINTF(stdout,"[%3i]            stopped at step %li with min(ess)=%.2f (var-ratio=%.2f/%.2f=%.3f, avg(acceptance)=%.2f)\n",
		    myID,world->burnin_stops[z].stopstep,
		    world->burnin_stops[z].ess,
		    world->burnin_stops[z].variance,world->burnin_stops[z].oldvariance,
		    world->burnin_stops[z].variance/world->burnin_stops[z].oldvariance,
		    world->burnin_stops[z].accept);
	    break;
	  default:
	    break;
	  }
      }
    world->options->treeprint = treeprint;
    world->in_burnin = FALSE;
}
