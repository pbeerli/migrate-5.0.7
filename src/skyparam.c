// skyline plot parametrized
#include "migration.h"
#include "tree.h"
#include "slice.h"
#include "tools.h"
#include "skyparam.h"
#include "random.h"
#include "sighandler.h"
#include "bayes.h"
#include "priors.h"
#include "pretty.h"
#include <math.h>

extern MYREAL log_prior_ratio_uni(MYREAL newparam, 
			   MYREAL oldparam, 
			   bayes_fmt * bayes, 
				  long which);

long slice_bayes_update_timeparam(world_fmt * world);
#ifdef MPI
void      skyline_param_record(float *temp, long *z, world_fmt * world);
#else
void      skyline_param_record(char *temp, long *c, world_fmt * world);
#endif
void  create_recording_times(float ***times, long numtimesalloc, long loci, long elements);
void  recreate_recording_times(float ***times, long numtimesalloc, long loci, long elements);

//##
///
/// standard metropolis proposal
MYREAL bayes_update_timeparam(world_fmt * world, boolean *success)
{
  //uniform_proposal(long which, world_fmt * world, MYREAL *oldparam, boolean *success)
  long npx = world->numpop2 + (world->bayes->mu ? 1 : 0);
  long blocks = world->timeelements;
  long addition=0;
  long end = world->timeelements +  world->timeelements * (world->numpop2 + addition);
  long whichblock = RANDINT (0, blocks-1);
  long which = RANDINT(0,npx-1);
  long w,wb;
  MYREAL newval;
  MYREAL newparam;
  MYREAL *oldtimeparam;
  MYREAL *skyparam;
  //const MYREAL      murate = world->options->mu_rates[world->locus];  
  bayes_fmt * bayes = world->bayes;
  /* todo set a map so that we can not only map out invalids but allow to have different models for different times*/
    while(bayes->map[which][1] == INVALID)
    {
      which  = RANDINT(0,npx-1);
    }
    doublevec1d(&oldtimeparam,end);
    memcpy(oldtimeparam, world->times,sizeof(MYREAL) * (size_t) end);
    w = bayes->map[which][1];
    wb = whichblock*npx + w;  
    skyparam = world->times + blocks + wb - w;
    //MYREAL mean = bayes->priormean[w];
    MYREAL oldval = probg_treetimes(world);
    MYREAL r = UNIF_RANDUM();
    // draw a new parameter from the prior distribution
    // for migration parameters we need to distinguish whether the 
    // prior is in terms of M or xNm
    newparam = (MYREAL) propose_uni_newparam(skyparam[w], w, world, &r);
    MYREAL hastingsratio = 1.0;
    //printf("param[%li/%li]=%f/%f\n",wb, which, newparam,ne);
    skyparam[w] = newparam;
    newval = probg_treetimes(world);
    //Acceptance or rejection of the new value
    *success = bayes_accept(newval, oldval,world->heat, hastingsratio);
    if(*success)
      {
	return newval;
      }
    else
      {
	if(world->options->prioralone)
	  {
	    *success = TRUE;
	    return newval;
	  }
	memcpy(world->times, oldtimeparam,sizeof(MYREAL) * (size_t) end);
	return oldval;
      }
}


long slice_bayes_update_timeparam(world_fmt * world)
{
  long npx = world->numpop2 + (world->bayes->mu ? 1 : 0);
  long blocks = world->timeelements;
  long addition=0;
  long end = world->timeelements + world->timeelements * (world->numpop2 + addition);
  long whichblock = RANDINT (0, blocks-1);
  long which = RANDINT(0,npx-1);
  long ba,w,wb;
  MYREAL newval;
  MYREAL *oldtimeparam;
  /* todo set a map so that we can not only map out invalids but allow to have different models for different times*/
    while(world->bayes->map[which][1] == INVALID)
    {
      which  = RANDINT(0,npx-1);
    }
    doublevec1d(&oldtimeparam,end);
    memcpy(oldtimeparam, world->times,sizeof(MYREAL) * (size_t) end);
    w = world->bayes->map[which][1];
    wb = whichblock*npx + w;  
    newval = sliceRatio(&world->timek[wb], wb, world, log_prior_ratio_uni);
    world->logprior = calculate_prior(world);
    world->bayes->oldval = newval;
    world->param_like = newval;
    ba = 1;
    world->bayes->accept[which] += ba;
    world->bayes->trials[which] += 1;
    // put in here the prior distribution recorder
    // record_prior(world->bayes->priorhist, world->locus, which, 
    // 
    myfree(oldtimeparam);
    return ba;
}

void insert_time_boundaries(timelist_fmt *timevector, world_fmt *world)
{
  MYREAL *times = world->times;
  long timeelements = world->timeelements;
  long start = timevector->T;
  long z;
  long i;
  if(timeelements>1)
    {
      if (timevector->allocT <= timevector->T + timeelements-1)
	{
	  increase_timelist2 (&timevector, timeelements); 
	}
      timevector->T += timeelements;
      z = 1;
      for(i=start; i < timevector->T; i++)
	{
	  timevector->tl[i].age = times[z];
	  timevector->tl[i].eventnode = NULL;
	  timevector->tl[i].slice = i;
	  timevector->tl[i].from = -1;
	  timevector->tl[i].to = -1;
	} 
    }
}

#ifdef MPI
void      skyline_param_record(float *temp, long *z, world_fmt * world)
{
  long i;
  long addition=0;
  long end = world->timeelements +  world->timeelements * (world->numpop2 + addition);
  for (i=0; i < end; i++)
    {
      temp[*z] = (float) world->times[i]; //times and timek 
      *z += 1;
    }
}
#else
void      skyline_param_record(char *temp, long *c, world_fmt * world)
{
  long i;
  long addition=0;
  long end = world->timeelements +  world->timeelements * (world->numpop2 + addition);
  for (i=0; i < end; i++)
    *c += sprintf(temp+ *c,"\t%f", world->times[i]);//times and timek 
}
#endif

void  create_recording_times(float ***times, long numtimesalloc, long loci, long elements)
{
  long i;
  (*times) = (float **) mycalloc(loci, sizeof(float *));
  for (i=0;i<loci;i++)
    {
      (*times)[i] = (float *) mycalloc(elements* (size_t) numtimesalloc,sizeof(float));
    }
}
void  recreate_recording_times(float ***times, long numtimesalloc, long loci, long elements)
{
  long i;
  for (i=0;i<loci;i++)
    {
      (*times)[i] = (float *) myrealloc((*times)[i], (size_t) (elements*numtimesalloc) * sizeof(float));
    }
}

void
skyline_param_reader(world_fmt *world, long step, long locus, char **input)
{
  (void) step;
  // ignore step, we simply append all parameters chronologically per locus, since we combine replicates
  // and the order for the summary does not matter this should be OK, for now I leave the parameter step in 
  long addition=0;
  long i;
  long end = world->timeelements +  world->timeelements * (world->numpop2 + addition);
  if(world->recording_times==NULL)
    {
      world->numtimes = 0;
      world->numtimesalloc = 1000; //DEBUG
      create_recording_times(&world->recording_times, world->numtimesalloc, world->loci,end);
    }
  if(world->numtimes >= world->numtimesalloc)
    {
      world->numtimesalloc += 1000; //DEBUG
      recreate_recording_times(&world->recording_times, world->numtimesalloc, world->loci,end);
    }
  float * recording_times = world->recording_times[locus] + world->numtimes * end; 
  for (i=0; i < end; i++)
    {
      recording_times[i] = (float) atof(strsep(input,"\t")); 
    }
}
