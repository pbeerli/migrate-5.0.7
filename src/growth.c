 /*------------------------------------------------------
 inference of population parameters
 using a Metropolis-Hastings Monte Carlo algorithm
 -------------------------------------------------------
 growth routines
 
 Peter Beerli 2013-2020, Tallahassee
 beerli@fsu.edu
 
 Copyright 2017 Peter Beerli
 
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
#include "growth.h"
#include "sighandler.h"
#include "bayes.h"
boolean is_notin(long x, long *order, long orderlen);
void reset_growth(world_fmt * world);
#if defined(MPI) && !defined(PARALIO) /* */
void print_growth_record(float *temp, long *z, world_fmt *world);
#else
void  print_growth_record(char *temp, long *c, world_fmt * world);
#endif
boolean init_growpop(worldoption_fmt * wopt, option_fmt *options, long numpop);


void reset_growth(world_fmt * world)
{
  long i;
  if (world->has_growth)
    {
      for(i=0;i<world->options->growpops_numalloc;i++)
	{
	  if (world->options->growpops[i]!=0)
	    world->growth[world->options->growpops[i]-1]=0.0;
	}
    }
}


boolean init_growpop(worldoption_fmt * wopt, option_fmt *options, long numpop)
{
  long i;
  boolean use_growth = FALSE;
  for (i=0;i<options->growpops_numalloc;i++)
    {
      if (options->growpops[i]>0)
	{
	  use_growth = TRUE;
	  break;
	}
    }
  if (use_growth)
    {
      wopt->growpops = (long*) mycalloc(numpop,sizeof(long));
      memcpy(wopt->growpops,options->growpops, sizeof(double) * (size_t) options->growpops_numalloc);
      if (options->growpops_numalloc < numpop)
	{
	  long last = options->growpops_numalloc - 1;
	  if (options->growpops[last]!=0)
	    {
	      for (i=last+1;i<numpop;i++)
		wopt->growpops[i]=i;
	    }
	}
      wopt->growpops_numalloc = numpop; //options->growpops_numalloc;
    }
  else
    {
      wopt->growpops = NULL;
      wopt->growpops_numalloc = 0;
    }
  return use_growth;
}

boolean is_notin(long x, long *order, long orderlen)
{
  long i;
  for (i=0;i<orderlen;i++)
    {
      if(x==order[i])
	return FALSE;
    }
  return TRUE;
}

void init_growth(world_fmt * world, long numpop)
{
  //boolean is_first= FALSE;
  long i;
  long j;
  long *tmp=NULL;
  //long z=0;
  world->grownum = 0;
  long order[1000] = {0}; // potentially 1000 populations
  long last=0;
  if(world->options->growpops != NULL)
    {
      world->has_growth = TRUE;
      tmp = (long*) mycalloc(numpop,sizeof(long));

      for (i=0; i < numpop; i++)
	tmp[i]=i+1;//1,2,3,4,5,6,7,8,9

      last = 0;
      order[last]=0;
      //evaluate the order of combined estimator
      for (i=0; i < world->options->growpops_numalloc; i++)
	{
	  long x = world->options->growpops[i];	  // for example: 1,0,1,1,2,0,0,2
	  if (x != 0)
	    {
	      if (is_notin(x,order,numpop))
		{
		  last++;
		  order[last]=x;//0,1,2,
		}
	    }
	  else
	    {
	      tmp[i]=0; // for example: 1,0,3,4,5,0,0,8
	    }
	}
      // no find the first occurence of the elements in order
      for (j=1;j<last+1; j++) //the first element in order is zero, we skip that
	{
	  //  is_first=FALSE;
	  for (i=0; i < world->options->growpops_numalloc; i++)
	    {
	      if (world->options->growpops[i]==order[j])	       
		{
		  //if (!is_first)
		  //{
		  //  first = i;
		  //  is_first = TRUE;
		  //  tmp[i]=i; //is this needed? I guess not!
		  //}
		  //else
		  //{
		  //    tmp[i] = first; //1: 1,0,1,1,5,0,0,8,9
		  tmp[i] = order[j]; //1: 1,0,1,1,5,0,0,8,9
		  //2: 1,0,1,1,2,0,0,2,9
		  //3: 1,0,1,1,2,0,0,2,3
		  //		}
		}
	    }
	}
      j = 0;
      for(i= world->options->growpops_numalloc;i<numpop;i++)
	{
	  j++;
	  tmp[i]=order[last]+j;
	}
      last += j;
      myfree(world->options->growpops);
      world->options->growpops=tmp;
      world->growth = (double*) mycalloc(last,sizeof(double));
      world->savegrowth = (double*) mycalloc(last,sizeof(double));
      world->grownum = last; //we do not calculate fo zero
    }
  else
    {
      world->has_growth = FALSE;
    }
  //myfree(world->options->growpops);
  //world->options->growpops = tmp;
}


#if defined(MPI) && !defined(PARALIO) /* */
void print_growth_record(float *temp, long *z, world_fmt *world)
{
  long i;
  if (world->has_growth)
    {
      for (i=0; i<world->numpop;i++)
	{
	  temp[(*z)++] = (float) world->growth[i];
	}
    }
}
#else /*not MPI or MPI & PARALIO*/
void  print_growth_record(char *temp, long *c, world_fmt * world)
{
  long i;
  if (world->has_growth)
    {
      for (i=0; i<world->numpop;i++)
	{
	  *c += sprintf(temp+ *c,"\t%f", world->growth[i]); 
	} 
    }
}
#endif

void construct_locusgrowth_histogram(world_fmt *world, long locus, MYREAL *mini, MYREAL *maxi, double **results)
{
  bayes_fmt *bayes = world->bayes;
  long i;
  long j0;
  double themean;
  double thestd;
 
  long numbin = 0;
  long pa;
 
  long rpa;
  long total=0;
  long *bins = bayes->histogram[locus].bins;
  boolean *visited;
  
  long grownum = world->options->growpops_numalloc;
  long np = world->numpop2 + bayes->mu + 2 * world->species_model_size;
  
  visited = (boolean *) mycalloc(world->options->growpops_numalloc, sizeof(boolean));
  for(i=0;i<np;i++)
    numbin += bins[i];
  
  for(j0=0; j0 < grownum; j0++)
    {
      long pick = world->options->growpops[j0];
      if (pick == 0)
	continue;
      else
	rpa = pick-1;
      if(!visited[rpa])
        {
	  themean = 0.0;
	  thestd = 0.0;
	  total = 0;
	  pa = rpa + np; //this points into the histograms (!! all other parameters, and then growth!!)	  
	  //np=#all param before growth, rpa= #growhindex, numbin=#bins before #grwothindex
	  construct_param_hist(world,locus,np, rpa,numbin, mini, maxi, results, &total,&themean,&thestd);
	  world->bayes->histogram[locus].means[pa] = themean/total;
	  world->bayes->histogram[locus].stds[pa]  = thestd / total;
	  world->bayes->histtotal[locus*np+pa] = (MYREAL) total;
	  visited[rpa] = TRUE;
	  numbin += bins[pa];
        }
    }
}

  
