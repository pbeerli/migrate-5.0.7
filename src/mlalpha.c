 /*------------------------------------------------------
 inference of population parameters
 using a Metropolis-Hastings Monte Carlo algorithm
 -------------------------------------------------------
 mlalpha routines
 
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
#include "mlalpha.h"
#include "sighandler.h"
#include "options.h"
#include "bayes.h"
#include "tools.h"
void reset_mlalpha(world_fmt * world);
#if defined(MPI) && !defined(PARALIO) /* */
void print_mlalpha_record(float *temp, long *z, world_fmt *world);
#else
void  print_mlalpha_record(char *temp, long *c, world_fmt * world);
#endif
boolean init_mlalphapop(worldoption_fmt * wopt, option_fmt *options, long numpop);

void print_parm_mlalphapops(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options, data_fmt *data);
void set_mlalpha(char **value, char **tmp, option_fmt *options);

void reset_mlalpha(world_fmt * world)
{
  long i;
  if (world->has_mlalpha)
    {
      for(i=0;i<world->options->mlalphapops_numalloc;i++)
	{
	  if (world->options->mlalphapops[i]!=0)
	    world->mlalpha[world->options->mlalphapops[i]-1]=1.0;
	}
    }
}


boolean init_mlalphapop(worldoption_fmt * wopt, option_fmt *options, long numpop)
{
  long i;
  boolean use_mlalpha = FALSE;
  for (i=0;i<options->mlalphapops_numalloc;i++)
    {
      if (options->mlalphapops[i]>0)
	{
	  use_mlalpha = TRUE;
	  break;
	}
    }
  if (use_mlalpha)
    {
      //world->has_mlalpha = TRUE;
      wopt->mlalphapops = (long*) mycalloc(numpop,sizeof(long));
      memcpy(wopt->mlalphapops,options->mlalphapops, sizeof(double) * (size_t) options->mlalphapops_numalloc);
      if (options->mlalphapops_numalloc < numpop)
	{
	  long last = options->mlalphapops_numalloc - 1;
	  if (options->mlalphapops[last]!=0)
	    {
	      for (i=last+1;i<numpop;i++)
		wopt->mlalphapops[i]=i;
	    }
	}
      wopt->mlalphapops_numalloc = numpop; //options->mlalphapops_numalloc;
    }
  else
    {
      //world->has_mlalpha = FALSE;
      wopt->mlalphapops=NULL;
      wopt->mlalphapops_numalloc = 0;
    }
  return use_mlalpha;
}


void init_mlalpha(world_fmt * world, long numpop)
{
  //boolean is_first= FALSE;
  long i;
  long j;
  long *tmp=NULL;
  //long z=0;
  world->mlalphanum = 0;
  long order[1000] = {0}; // potentially 1000 populations
  long last=0;
  if(world->options->mlalphapops != NULL)
    {
      world->has_mlalpha = TRUE;
      tmp = (long*) mycalloc(numpop,sizeof(long));

      for (i=0; i < numpop; i++)
	tmp[i]=i+1;//1,2,3,4,5,6,7,8,9

      last = 0;
      order[last]=0;
      //evaluate the order of combined estimator
      for (i=0; i < world->options->mlalphapops_numalloc; i++)
	{
	  long x = world->options->mlalphapops[i];	  // for example: 1,0,1,1,2,0,0,2
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
	  for (i=0; i < world->options->mlalphapops_numalloc; i++)
	    {
	      if (world->options->mlalphapops[i]==order[j])	       
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
      for(i= world->options->mlalphapops_numalloc;i<numpop;i++)
	{
	  j++;
	  tmp[i]=order[last]+j;
	}
      last += j;
      myfree(world->options->mlalphapops);
      world->options->mlalphapops=tmp;
      world->mlalpha = (double*) mycalloc(last,sizeof(double));
      world->savemlalpha = (double*) mycalloc(last,sizeof(double));
      world->mlalphanum = last; //we do not calculate for zero
    }
  else
    {
      world->has_mlalpha = FALSE;
    }
  //myfree(world->options->mlalphapops);
  //world->options->mlalphapops = tmp;
}


#if defined(MPI) && !defined(PARALIO) /* */
void print_mlalpha_record(float *temp, long *z, world_fmt *world)
{
  long i;
  if (world->has_mlalpha)
    {
      for (i=0; i<world->numpop;i++)
	{
	  temp[(*z)++] = (float) world->mlalpha[i];
	}
    }
}
#else /*not MPI or MPI & PARALIO*/
void  print_mlalpha_record(char *temp, long *c, world_fmt * world)
{
  long i;
  if (world->has_mlalpha)
    {
      for (i=0; i<world->numpop;i++)
	{
	  *c += mysnprintf(temp+ *c,LINESIZE,"\t%f", world->mlalpha[i]); 
	} 
    }
}
#endif

void construct_locusmlalpha_histogram(world_fmt *world, long locus, MYREAL *mini, MYREAL *maxi, double **results)
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
  
  long mlalphanum = world->options->mlalphapops_numalloc;
  long np = world->numparamcumvec[GROWTHPRIOR];
  //was with growth before world->numpop2 + bayes->mu + 2 * world->species_model_size;
  
  visited = (boolean *) mycalloc(world->options->mlalphapops_numalloc, sizeof(boolean));
  for(i=0;i<np;i++)
    numbin += bins[i];
  
  for(j0=0; j0 < mlalphanum; j0++)
    {
      long pick = world->options->mlalphapops[j0];
      if (pick == 0)
	continue;
      else
	rpa = pick-1;
      if(!visited[rpa])
        {
	  themean = 0.0;
	  thestd = 0.0;
	  total = 0;
	  pa = rpa + np; //this points into the histograms (!! all other parameters, and then mlalpha!!)	  
	  //np=#all param before mlalpha, rpa= #mlalpha hindex, numbin=#bins before #grwothindex
	  construct_param_hist(world,locus,np, rpa,numbin, mini, maxi, results, &total,&themean,&thestd);
	  world->bayes->histogram[locus].means[pa] = themean/total;
	  world->bayes->histogram[locus].stds[pa]  = thestd / total;
	  world->bayes->histtotal[locus*np+pa] = (MYREAL) total;
	  visited[rpa] = TRUE;
	  numbin += bins[pa];
        }
    }
}

  
///
/// print the parmfile entry for the population mlalpha parameter setting and labeling
/// 0=no mlalpha, other labels are marking mlalpha either individual populations or combinations
/// for example:
/// population-mlalpha={0} # all populations use mlalpha=1
/// population-mlalpha={1} # all populations have the same mlalpha rate
/// population-mlalpha={1 2 3} # all 3 populations have individual rates, if there are more than
///                           # than 3 populations, pop 4 etc will be in lockstep with 3
/// population-mlalpha={1 0 1} # population 1 and 3 are in lockstep, population 2 mlalpha=1
void print_parm_mlalphapops(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options, data_fmt *data)
{
  (void) data;
  long pos;
  long i;
  char *input;
  print_parm_br(bufsize, buffer, allocbufsize);
  print_parm_comment(bufsize, buffer, allocbufsize, "      print the parmfile entry for the population mlalpha parameter setting and labeling");
  print_parm_comment(bufsize, buffer, allocbufsize, "      0=no mlalpha, other labels are marking mlalpha either individual populations or combinations");
  print_parm_comment(bufsize, buffer, allocbufsize, "      for example:");
  print_parm_comment(bufsize, buffer, allocbufsize, "      population-mlalpha={0} # all populations are constant");
  print_parm_comment(bufsize, buffer, allocbufsize, "      population-mlalpha={1} # all populations are exponentially mlalphaing with the same mlalpha rate");
  print_parm_comment(bufsize, buffer, allocbufsize, "      population-mlalpha={1 2 3} # all 3 populations are mlalphaing with individual rates, if there are more than");
  print_parm_comment(bufsize, buffer, allocbufsize, "                                # than 3 populations, pop 4 etc will mlalpha in lockstep with 3");
  print_parm_comment(bufsize, buffer, allocbufsize, "      population-mlalpha={1 0 1} # population 1 and 3 mlalpha in lockstep, population 2 is constant in size.");
  input = (char *) mycalloc(options->mlalphapops_numalloc * 50, sizeof(char));
  pos = mysnprintf(input,LINESIZE, "population-mlalpha={%li", options->mlalphapops[0]);
  for(i=1; i < options->mlalphapops_numalloc; i++)
    {
      pos += mysnprintf(input + pos,LINESIZE,", %li", options->mlalphapops[i]);
    }
  mysnprintf(input + pos,LINESIZE,"}");
  print_parm_mutable(bufsize, buffer, allocbufsize, "%s", input);
  print_parm_br(bufsize, buffer, allocbufsize);
  myfree(input); 
}

void set_mlalpha(char **value, char **tmp, option_fmt *options)
{
      long z=0;
      options->mlalphapops[0] = 0;
      get_next_word(value,"{}:,; ",tmp);
      while(*tmp != NULL)
	{
	  if(z >= options->mlalphapops_numalloc)
	    {
	      options->mlalphapops_numalloc = z+1;
	      options->mlalphapops = (long*) myrealloc(options->mlalphapops, options->mlalphapops_numalloc * sizeof(long));
	    }
	  options->mlalphapops[z] = atol(*tmp);
	  //printf("%i> population-mlalphath relabel [%li]=%li (input: |%s|)\n",myID,z, options->mlalphapops[z],*tmp);
	  z++;
	  get_next_word(value,"{}:,; ",tmp);
	}
}
