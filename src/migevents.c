/*------------------------------------------------------
Maximum likelihood estimation
of migration rate  and effectice population size
using a Metropolis-Hastings Monte Carlo algorithm
-------------------------------------------------------
    mighistogram   R O U T I N E S that use less memory

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
/*! \file migevents.c 

this file contains functions that repport the numbers and times of coalescence and migration events
Results will be printed as a histogram over time and also will print the histogram to file and
print a summary table
*/
#include <stdlib.h>

#include "migevents.h"
#include "skyline.h"
#include "random.h"
#include "tools.h"
#include "sighandler.h"
#include "world.h"
#include "pretty.h"
#include <string.h>
#ifdef MPI
#include "migrate_mpi.h"
#else
extern int myID;
#endif

void setup_mig_coal_events (world_fmt * world, option_fmt * options);
void minmax (histogram_fmt * hist, float *tempmin, float *tempmax);
void destroy_mig_coal_events (world_fmt * world);
void print_event_values_list(FILE *file, long locus, duo **eventbins, MYREAL eventbinsize, long *eventbinnum, world_fmt* world);
void read_event_values_fromfile(FILE *file, world_fmt *world);
void print_event_values_tofile(FILE *file,  world_fmt *world);
void prepare_event_values(world_fmt *world, MYREAL *sums0,MYREAL *sums1);
void print_event_values_title(FILE *file, boolean progress);
void store_events (world_fmt * world, timelist_fmt * ltl, long np, long rep);
//##
///
/// allocates the necessary memory for recording of the events through time
void
setup_mighist (world_fmt * world, option_fmt * options)
{
    if (world->options->mighist)
      {
        world->mighistloci = (mighistloci_fmt *)
        mycalloc (world->loci+1, sizeof (mighistloci_fmt));
        world->mighistlocinum = 0;
        setup_mig_coal_events(world,options);
      }
    else
        world->mighistloci = NULL;    
}

///
/// destructor of the container that records events through time
void
destroy_mighist (world_fmt * world)
{
    long locus;
    long pop;
    if (world->options->mighist)
      {
        for(locus=0; locus < world->loci; locus++)
          {
	    for(pop=0;pop < world->numpop2; pop++)
	      {
		myfree(world->mighistloci[locus].migeventbins[pop]);
	      }
            myfree(world->mighistloci[locus].migeventbins);
            //  if(world->options->skyline)
            //  myfree(world->mighistloci[locus].eventbins);
          }
        myfree(world->mighistloci);
      }
}

void minmax (histogram_fmt * hist, float *tempmin, float *tempmax)
{
    long i;
    MYREAL tmp1, tmp2;
    MYREAL tmpmin = MYREAL_MAX;
    MYREAL tmpmax = -MYREAL_MAX;
    
    for (i = 0; i < hist->count; i++)
      {
        if ((tmp1 = hist->time[i]) < tmpmin)
            tmpmin = tmp1;
        if ((tmp2 = hist->time[i]) > tmpmax)
            tmpmax = tmp2;
      }
    *tempmax = (float) tmpmax;
    *tempmin = (float) tmpmin;
}

///
/// print and calculate means, medians, of migration and coalescence event times
/// over all times and of the last (the MRCA) event.
void
print_mighist_output (FILE * out, world_fmt * world, MYREAL *sums, boolean mrca) 
{
    const long column   = (mrca ? 1 : 0 );
    const long numpop2  = world->numpop2;
    const long npp      =  numpop2 + world->bayes->mu + world->species_model_size * 2 + world->grownum;
    const long startpop = (world->options->mighist_all ? 0 : world->numpop);
    const long endpop   = (mrca ? world->numpop : npp);
    const long loci1    = ((world->loci == 1) ? 1 : world->loci + 1);
    const long loci1pop = loci1 * npp;
    const long sumloc   = ((world->loci == 1) ? 0 : world->loci);
    
    long   locus; 
    long   pop;
    long   i;
    long   lp;
    long   frompop;
    long   topop;
    long * eventbinnum = NULL;
    long * eventbinmax = NULL;
    
    duo ** eventbins;
    
    double *total;
    double *meantime;
    double *mediantime;
    double *stdtime;
    double *freq;
    double *values;
    double **allvalues=NULL;
    double ntotal;
    double eventbinsize;
    double age;
    double binfreq;

    eventbinmax = (long *) mycalloc(npp, sizeof(long));

    
    if(loci1>1)
      {
        for (pop = startpop; pop < endpop; pop++)
          {
            for (locus = 0; locus < world->loci; locus++)
              {
                eventbinnum = world->mighistloci[locus].migeventbinnum;
                if(eventbinnum[pop] > eventbinmax[pop])
                    eventbinmax[pop] = eventbinnum[pop];
              }  
          }
        allvalues = (double **) mycalloc(npp,sizeof(double*));
        for (pop = startpop; pop < endpop; pop++)
          { 
            allvalues[pop] = (double *) mycalloc(eventbinmax[pop]+1,sizeof(double));
          }
      }
    // initialize the vectors used for the table
    values = (double *) mycalloc(1, sizeof(double));
    total = (double *) mycalloc(loci1, sizeof(double));
    meantime = (double *) mycalloc(loci1pop * 4, sizeof(double));
    mediantime = meantime + (loci1pop);
    stdtime = mediantime + (loci1pop);
    freq = stdtime + (loci1pop);
    
    for (locus = 0; locus < world->loci; locus++)
      {
        eventbins = world->mighistloci[locus].migeventbins;
        eventbinnum = world->mighistloci[locus].migeventbinnum;
        eventbinsize = (double) world->mighistloci[locus].eventbinsize;
        
        for (pop = startpop; pop < endpop; pop++)
          {
            age = eventbinsize / 2.0;
            values = (double *) myrealloc(values, eventbinnum[pop] * sizeof(double));
            
            lp = locus * npp + pop;
            
            meantime[lp] = 0.;
            stdtime[lp] = 0.;
            freq[lp] = 0.;
            
            for (i = 0; i < eventbinnum[pop]; i++)
              {
                binfreq = (double) (eventbins[pop][i][column]);
                freq[lp] += binfreq; // total frequency of events per parameter 
                total[locus] += binfreq; // total of frequency this should add up to world->loci
                values[i] = binfreq * age; // frequency times age [how many times was this age recorded.]
                meantime[lp] += values[i];
                stdtime[lp] += age * age * binfreq;
                age += eventbinsize;
                if(loci1 > 1)
                  {
                    allvalues[pop][i] += values[i];
                  }
                if(freq[lp] < 0.499999)
                  {
                    mediantime[lp] = age;
                  }
              }
            // freq[lp] should be one, but the divisions are done to fix accumulated rounding issues.
	    if(freq[lp]>0.0)
	      {
		meantime[lp] /= freq[lp];
		stdtime[lp] = ((stdtime[lp]/freq[lp]) - (meantime[lp]*meantime[lp]));
	      }
            if(loci1>1)
              {
                meantime[sumloc * npp + pop] += meantime[lp]/world->loci;
		mediantime[sumloc * npp + pop] += mediantime[lp]/world->loci;
                stdtime[sumloc * npp + pop] += stdtime[lp]/world->loci;
                mediantime[sumloc * npp + pop] += mediantime[lp]/world->loci;
              }
          }
        // calculating the frequences over all parameters using the sums
        ntotal = 0.;
        for (pop = startpop; pop < endpop; pop++)
          {
            lp = locus * npp + pop;
            ntotal += (double) sums[lp];
            freq[lp] =  (double) sums[lp];  
            stdtime[lp] = (double) sqrt((double) stdtime[lp]);
          }
        for (pop = startpop; pop < endpop; pop++)
          {
            lp = locus * npp + pop;
            freq[lp] /= ntotal;
            if(loci1>1)
	      freq[sumloc * npp + pop] += freq[lp]/world->loci;
          }
	if(loci1>1)
	  total[sumloc] += total[locus];
      }
    
    myfree(values);
    
    if(loci1 > 1)
      {
        for (pop = startpop; pop < endpop; pop++)
          {
            lp = sumloc * npp + pop;
            //mediantime[lp] = -1.; //DEBUG: allvalues[pop][(long) (eventbinnum[pop]/2)]/freq[lp];
            stdtime[lp] = (double) sqrt((double) stdtime[lp]);
          }
      }
    
    if(mrca)
      {
        FPRINTF (out, "\n\nTime and probability of most recent common ancestor\n");
        FPRINTF (out,     "===================================================\n\n");
      }
    else
      {
        FPRINTF (out, "\n\nApproximate Summary of %s events\n", 
                 world->options->mighist_all ? "coalescence and migration" : "migration");
        FPRINTF (out, "===============%s================\n\n", 
                 world->options->mighist_all ? "=========================" : "=========");
      }
    for (locus = 0; locus < loci1; locus++)
      {    /* Each locus + Summary */
        if (locus != world->loci)
            FPRINTF (out, "Locus %li\n", locus + 1);
        else
            FPRINTF (out, "Over all loci\n");
        FPRINTF (out,
                 "---------------------------------------------------------\n");
        FPRINTF (out,
                 "Population   Time                             Frequency\n");
        FPRINTF (out, "             -----------------------------\n");
        FPRINTF (out, "From    To   Average    Median     Std\n");
        FPRINTF (out,
                 "---------------------------------------------------------\n");
	//        for (pop = (world->options->mighist_all ? 0 : world->numpop); pop <  (mrca ? world->numpop : world->numpop2); pop++)
        for (pop = 0; pop < npp; pop++)
          {
	    if(world->bayes->map[pop][1] == INVALID)
	      continue;
	    if (pop>world->numpop2)
	      {
		if ((pop-world->numpop2-world->bayes->mu)%2 == 1)
		  continue;
		d2mm(pop,world,&frompop,&topop);
	      }
	    else
	      {
		m2mm(pop,world->numpop,&frompop,&topop);
	      }
            lp = locus * npp + pop;
            if(freq[lp]<=0.0)
              {
                FPRINTF (out,
                         "%4li %4li    -------    -------    -------    %3.5f\n",
                         frompop + 1,
                         topop + 1,
                         0.0);
              }
            else
              {
                FPRINTF (out,
                         "%4li %4li    %3.5f    %3.5f    %3.5f    %3.5f\n",
                         frompop + 1,
                         topop + 1,
                         (double) meantime[lp], (double) mediantime[lp],
                         (double) stdtime[lp],
                         (double) freq[lp]);
              }
          }
        FPRINTF (out,
                 "---------------------------------------------------------\n");
        FPRINTF (out, "\n");
      }
#ifdef PRETTY
pdf_print_time_table (world, meantime, mediantime, stdtime, freq, mrca);
#endif

 myfree(total);
 myfree(meantime);
 myfree(allvalues);
 myfree(eventbinmax);
}



///
/// gather the events from the genealogy for each timeinterval.
/// This function constructs a count-histogram for all event types
/// and a specific histogram for the last coalescent (the MRCA)
void calculate_event_values(duo **eventbins, long *eventbinnum, MYREAL eventinterval, 
                            MYREAL interval, MYREAL age, char type, long from, long to, 
                            long * lineages, long numpop, boolean is_last, world_fmt *world)
{
  (void) lineages;
    MYREAL inv_eventinterval = 1. / eventinterval;
    long bin = (long) (age * inv_eventinterval); // this bin
                                                 //  long k;
    long i;
    long pop;
    MYREAL val;
    long lastbin = (long) ((age - interval) * inv_eventinterval);
    long diff;
    duo *bins=NULL;
    long allocsize=1;
    if(lastbin < 0)
      { 
        lastbin=0;
      }
    //    if(bin > 10000 || bin < 0) // tooo many bins and it overflows probably 
    //                           //I should stop with some bins like 10000
    //  {
    //    warning("Bin (%li) is smaller than zero or bigger than 10000", bin);
    //    return;
    //  }
    diff = bin - lastbin;
    val = 1.0;
    if(from == to)
      {
        pop = to;
        //fprintf(stdout,"%i> c bin=%li ebin[%li]=%li\n",myID, bin, pop, eventbinnum[pop]);
      } 
    else
      {
	if(!strchr("dDtT",type))
	  pop = mm2m(from, to, numpop);
	else
	  pop = mm2d(from,to,world);
        //fprintf(stdout,"%i> m bin=%li ebin[%li]=%li\n",myID, bin, pop, eventbinnum[pop]);
      }
    //fflush(stdout);
    if(bin>=eventbinnum[pop])
      {
        allocsize = bin+1;
        eventbins[pop] = (duo *) myrealloc(eventbins[pop], allocsize * sizeof(duo));
        //	  fprintf(stdout,"%i> ral: as=%li, old=%li, last val=(%f,%f)\n",myID,allocsize, eventbinnum[pop],
        //  eventbins[pop][eventbinnum[pop]-1][0],eventbins[pop][eventbinnum[pop]-1][1]);
        for(i=eventbinnum[pop];i < allocsize; i++)
          {
            eventbins[pop][i][0] = 0.;
            eventbins[pop][i][1] = 0.;
          }
        eventbinnum[pop] = allocsize;
      }
    bins = eventbins[pop];
    //  fprintf(stdout,"m diff=%li ebin[pop]=%li, from=%i to=%li val=%f interval=%f bin=%li\n",diff , eventbinnum[pop], from, pop, val, interval, bin);
    
    if(diff < 0)
        error("time difference is negative -- not possible");
    bins[bin][0] += (float) val;
    if(is_last)
      bins[bin][1] += (float) val;
}

///
/// set up the event plot histogram containers, only when also the migration histograms are recorded
/// this function needs to be called by setup_mighist() function
void
setup_mig_coal_events (world_fmt * world, option_fmt * options)
{
  (void) options;
    long locus, i, j;
    long allocsize = 1;
    long npall = world->numpop2+ world->bayes->mu + world->species_model_size * 2 + world->grownum;
    MYREAL  binsize = (double) world->options->eventbinsize; // in mutational units (=time scale)
    if (world->options->mighist)
      {
        for (locus = 0; locus < world->loci; locus++)
          {
            world->mighistloci[locus].eventbinsize = binsize;
	    if(world->mighistloci[locus].migeventbins == NULL)
	      world->mighistloci[locus].migeventbins = (duo **) mycalloc (npall, sizeof (duo *));
	    if(world->mighistloci[locus].migeventbinnum == NULL)
	      world->mighistloci[locus].migeventbinnum = (long *) mycalloc (npall, sizeof (long));
            
            for (i = 0; i < npall; i++)
              {
                world->mighistloci[locus].migeventbinnum[i] = allocsize;
		if(world->mighistloci[locus].migeventbins[i] == NULL)
		  world->mighistloci[locus].migeventbins[i] = (duo *) mycalloc (allocsize, sizeof (duo));
                for(j=0;j<allocsize;j++)
                  {
                    world->mighistloci[locus].migeventbins[i][j][0] = 0.;
                    world->mighistloci[locus].migeventbins[i][j][1] = 0.;
                  }
              }
          }
      }
}

///
/// Destroy the event plot histogram container
void destroy_mig_coal_events (world_fmt * world)
{
    long locus, i;
    long npall = world->numpop2  + world->bayes->mu + world->species_model_size * 2 + world->grownum;
    if (world->options->mighist)
      {
        for (locus = 0; locus < world->loci; locus++)
          {	  
            for (i = 0; i < npall; i++)
              {
                myfree(world->mighistloci[locus].migeventbins[i]);
              }
            myfree(world->mighistloci[locus].migeventbins);
            myfree(world->mighistloci[locus].migeventbinnum);
          }
      }
}


void print_event_values_list(FILE *file, long locus, duo **eventbins, MYREAL eventbinsize, long *eventbinnum, world_fmt* world)
{
    long i;
    long pop;
    long frompop;
    long topop;
    long numpop = world->numpop;
    MYREAL age;
    long npall = world->numpop2+world->species_model_size * 2 + world->bayes->mu;
    for(pop = 0; pop < npall; pop++)
      {  	
	if(world->bayes->map[pop][1] == INVALID)
	  continue;	
	if (pop>world->numpop2)
	  {
	    if ((pop-world->numpop2-world->bayes->mu)%2 == 1)
	      continue;
	  }
	
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
        fprintf(file,"Time        Frequency of visit   Frequency of MRCA\n");
        fprintf(file,"--------------------------------------------------\n");
        age = eventbinsize/2.;
        for(i = 0; i < eventbinnum[pop]; i++)
          {
            fprintf(file,"%10.10f  %10.10f     %10.10f\n", age, (double) eventbins[pop][i][0], (double) eventbins[pop][i][1]);
            age += eventbinsize;
          }
      } 
}

///
/// read events from mighistfile and push them into the mighistlocus structure
#define SOME_ELEMENTS 100
void read_event_values_fromfile(FILE *file, world_fmt *world)
{
  long locus;
  long pop;
  long i;
//  MYREAL age;
  char *input;
  char *inptr;
  duo **eventbins;
  long *eventbinnum;
  //MYREAL eventbinsize;
  long spacer=0;
  long *allocsize;
  long npall = world->numpop2+world->species_model_size * 2 + world->bayes->mu;

  input = (char *) mycalloc(SUPERLINESIZE, sizeof(char));
  allocsize = (long *) mycalloc((npall*world->loci),sizeof(long));
  for(i=0;i<npall*world->loci;i++)
    {
      allocsize[i]=1;
    }
  while(FGETS(input,SUPERLINESIZE,file) != EOF)
    {
      // grab the commentlines
      while(input[0]=='#')
	{
	  FGETS(input,LINESIZE,file);
	}
      // read the mighistfile
      if(input !=NULL)
	{
	  inptr = input;
	  locus      = atol(strsep(&inptr,"\t"))-1;
	  spacer = locus * npall;
	  if(locus == -1)
	    error("help");
	  pop       = atol(strsep(&inptr,"\t"))-1;
	  i          = atol(strsep(&inptr,"\t"))-1;
	  //age        = atof(strsep(&inptr,"\t"));
      (void) strsep(&inptr,"\t");
	  eventbins =  world->mighistloci[locus].migeventbins;
	  eventbinnum = world->mighistloci[locus].migeventbinnum;
	  eventbinnum[pop] = i+1;
	  if(allocsize[spacer + pop] < eventbinnum[pop])
	    {
	      allocsize[spacer + pop] += SOME_ELEMENTS;
	      eventbins[pop] = (duo *) myrealloc(eventbins[pop],allocsize[spacer + pop] * sizeof(duo));
	    }
	  eventbins[pop][i][0] = (float) atof(strsep(&inptr,"\t"));
	  eventbins[pop][i][1] = (float) atof(strsep(&inptr,"\t"));
	}
    }
  myfree(allocsize);
}


void print_event_values_tofile(FILE *file,  world_fmt *world)
{
    long i;
    long pop;
    long frompop;
    long topop;
    long locus;
    long sumloc;
    long numpop = world->numpop;

    long npall = world->numpop2+world->species_model_size * 2 + world->bayes->mu;
    MYREAL age;
    duo **eventbins;
    long *eventbinnum;
    MYREAL eventbinsize;
    
    fprintf(file,"# Raw record of the events histogram for all parameters and all loci\n");  
    fprintf(file,"# The time interval is set to %f\n", (double) world->options->eventbinsize);  
    fprintf(file,"# produced by the program %s (http://popgen.scs.fsu.edu/migrate.hml)\n",
            MIGRATEVERSION);  
    fprintf(file,"# written by Peter Beerli 2006-2014, Tallahassee,\n");
    fprintf(file,"# if you have problems with this file please email to beerli@fsu.edu\n");  
    fprintf(file,"#\n");
    fprintf(file,"# Order of the parameters:\n");
    fprintf(file,"# Parameter-number Parameter\n");
    for(pop=0;pop<npall;pop++)
      {
	if(world->bayes->map[pop][1] == INVALID)
	  continue;
	if (pop>world->numpop2)
	  {
	    if ((pop-world->numpop2-world->bayes->mu)%2 == 1)
	      continue;
	  }
	
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
	      }
	  }
      }
    fprintf(file,"#\n#----------------------------------------------------------------------------\n");
    fprintf(file,"# Locus Parameter-number Bin Age Events-frequency(*) MRCA-frequency(*)\n");
    fprintf(file,"#----------------------------------------------------------------------------\n");
    fprintf(file,"# (*) values with -1 were NEVER visited\n");
    fprintf(file,"# @@@@@Locus param bin age freq mrcafreq\n");
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
	    eventbins =  world->mighistloci[locus].migeventbins;
	    eventbinnum = world->mighistloci[locus].migeventbinnum;
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
		    fprintf(file,"%li\t%li\t%li\t%10.10f\t%10.10f\t%10.10f\n", 
			    locus+1, pop+1, i+1, age, (double) eventbins[pop][i][0], (double) eventbins[pop][i][1]);
		    age += eventbinsize;
		  }
	      } 
	  }
      }
}

///
/// calculates frequencies for all event bins (over all and MRCA)
void prepare_event_values(world_fmt *world, MYREAL *sums0,MYREAL *sums1)
{
    MYREAL sum0;
    MYREAL sum1;
    long   z = 0;
    long   locus;
    long   pop;
    long   i;
    //long   numpop2 = world->numpop2;
    long npall = world->numpop2+world->species_model_size * 2 + world->bayes->mu;
    duo    **eventbins;
    duo    **eventbins_all;
    long   *eventbinnum;
    long   eventbinnum_allmax=0;
    
    for(locus=0; locus < world->loci; locus++)
      {
	if(!world->data->skiploci[locus])
	  {
	    eventbins =  world->mighistloci[locus].migeventbins;
	    eventbinnum = world->mighistloci[locus].migeventbinnum;
	    for(pop = 0; pop < npall; pop++)
	      {
		sum0 = (MYREAL) 0.0;
		sum1 = (MYREAL) 0.0;
		
		if(eventbinnum[pop] > eventbinnum_allmax)
		  eventbinnum_allmax = eventbinnum[pop];
		
		for(i = 0; i < eventbinnum[pop]; i++)
		  {
		    if(eventbins[pop][i][0] <= 0.0f)
		      {    
			continue;
		      }
		    else
		      {
			sum0 += (double) eventbins[pop][i][0];                    
		      }
		    
		    if(eventbins[pop][i][1] <= 0.0f)
		      { 
			continue;
		      }
		    else
		      {
			sum1 += (double) eventbins[pop][i][1];
		      }
		  }
		// calculate frequency from bincount data
		for(i = 0; i < eventbinnum[pop]; i++)
		  {
		    if(eventbins[pop][i][0] <= 0.0f)
		      continue;
		    eventbins[pop][i][0] /= (float) sum0;
		    if(eventbins[pop][i][1] <= 0.0f)
		      continue;
		    eventbins[pop][i][1] /= (float) sum1;
		  }
		sums0[z]=sum0;
		sums1[z++]=sum1;
	      } 
	  }
      }
    // for all loci, even when there is only one
    if(world->mighistloci[world->loci].migeventbins == NULL)
    world->mighistloci[world->loci].migeventbins = (duo **) mycalloc(npall,sizeof(duo *));
    if(world->mighistloci[world->loci].migeventbinnum == NULL)
    world->mighistloci[world->loci].migeventbinnum = (long *) mycalloc(npall,sizeof(long));
    world->mighistloci[world->loci].eventbinsize = world->mighistloci[0].eventbinsize;
    for(pop=0; pop< npall; pop++)    
      {
	if(world->mighistloci[world->loci].migeventbins[pop] == NULL)
	  world->mighistloci[world->loci].migeventbins[pop] = (duo *) mycalloc(eventbinnum_allmax,sizeof(duo));
        world->mighistloci[world->loci].migeventbinnum[pop] = eventbinnum_allmax;  
      }
    eventbins_all = world->mighistloci[world->loci].migeventbins;
    for (locus = 0; locus < world->loci; locus++)
      {
	if(!world->data->skiploci[locus])
	  {
	    eventbinnum = world->mighistloci[locus].migeventbinnum;
	    eventbins = world->mighistloci[locus].migeventbins;
	    for(pop=0; pop< npall; pop++)
	      {
		for(i=0; i < eventbinnum[pop]; i++)
		  {
		    if(eventbins[pop][i][0] > 0.0f)
		      {
			eventbins_all[pop][i][0] += eventbins[pop][i][0]/world->loci;
		      }
		    if(eventbins[pop][i][1] > 0.0f)
		      {
			eventbins_all[pop][i][1] += eventbins[pop][i][1]/world->loci;
		      }
		  }
	      }
	  }
      }
}

///
/// print titles for events over time
void print_event_values_title(FILE *file, boolean progress)
{
    if(progress)
      {      
        fprintf(file,"\n\nEvents over time\n");
        fprintf(file,"--------------------\n\n");
      }
}

///
/// print event statistics to outfile, mighistfile and also PDF
void print_event_values(world_fmt * world)
{
    long locus;
    long sumloc = (world->loci > 1) ? 1 : 0;
    MYREAL *sums0;
    MYREAL *sums1;
    long npall = world->numpop2+world->species_model_size * 2 + world->bayes->mu;
    if(world->options->mighist)
      {
    
	sums0 = (MYREAL *) mycalloc((world->loci * npall * 2),sizeof(MYREAL));
	sums1 = sums0 + (world->loci * npall);

	if(world->options->datatype != 'g')
	  {
	    // prepare event histogram for printing and get total of observations
	    // for all (sums0) and for the mrcas only (sums1)
	    prepare_event_values(world,sums0, sums1);
	    
	    // print event to file
	    print_event_values_tofile(world->mighistfile, world);
	  }
	else
	  {
	    read_event_values_fromfile(world->mighistfile,world);
	    prepare_event_values(world,sums0, sums1);
	  }
        // print title to screen
        print_event_values_title(stdout, world->options->progress);
        
        // print title to ascii-outfile
        print_event_values_title(world->outfile, TRUE);
#ifdef PRETTY
        if(world->options->verbose)
            pdf_print_eventtime_table(world);
#endif
        // print content to stdout and to ascii - outfile [it seems not useful in the output file
        if(world->options->verbose)
          {
            for(locus=0; locus < world->loci+sumloc; locus++)
              {
		if(!world->data->skiploci[locus])
		  {
		    print_event_values_list(stdout, locus, world->mighistloci[locus].migeventbins, 
					    world->mighistloci[locus].eventbinsize, 
					    world->mighistloci[locus].migeventbinnum, world);
		    print_event_values_list(world->outfile, locus, world->mighistloci[locus].migeventbins, 
					    world->mighistloci[locus].eventbinsize, 
					    world->mighistloci[locus].migeventbinnum, world);
		  }
	      }
          }
#ifdef PRETTY
        pdf_event_histogram(world->loci, npall,  world);
#endif
        // this include also printing to the pretty interface
        print_mighist_output (world->outfile, world, sums0, FALSE);
        print_mighist_output (world->outfile, world, sums1, TRUE);
	myfree(sums0);
      }
}

/// store the event statistics
void store_events (world_fmt * world, timelist_fmt * ltl, long np, long rep)
{
  (void) rep;
  (void) np;
    long j;
    long T;
    MYREAL t;
    mighistloci_fmt *bb = NULL;
    if (world->in_last_chain && world->options->mighist)
      {
	world->options->mighist_counter = 0;

	if(!world->data->skiploci[world->locus])
	  {
	    bb = &(world->mighistloci[world->locus]);
	    
	    if((*ltl).tl[0].eventnode->type != 't')
	      {
		calculate_event_values(bb->migeventbins, bb->migeventbinnum, bb->eventbinsize, 
				       (*ltl).tl[0].age, 
				       (*ltl).tl[0].age, 
				       (*ltl).tl[0].eventnode->type, 
				       (*ltl).tl[0].from, (*ltl).tl[0].to,
				       (*ltl).tl[0].lineages, world->numpop, FALSE,
				       world);
		// calculate skyline values
		if(world->options->skyline)
		  {
		    if((*ltl).tl[0].eventnode->visited)
		      {
			calculate_expected_values(bb->eventbins, bb->eventbinnum, bb->eventbinsize, 
						  (*ltl).tl[0].age, 
						  (*ltl).tl[0].age, 
						  (*ltl).tl[0].eventnode->type,
						  (*ltl).tl[0].from, (*ltl).tl[0].to,
						  (*ltl).tl[0].lineages, world->numpop, world);
		      }
		  }
	      }
	    T = (*ltl).T - 1;
	    for (j = 1; j < T; j++)
	      {
		t = (*ltl).tl[j].age - (*ltl).tl[j - 1].age;
		if((*ltl).tl[j].eventnode->type != 't')
		  {
		    // calculate event values
		    calculate_event_values(bb->migeventbins, bb->migeventbinnum, bb->eventbinsize, 
					   t, 
					   (*ltl).tl[j].age, 
					   (*ltl).tl[j].eventnode->type, 
					   (*ltl).tl[j].from, (*ltl).tl[j].to,
					   (*ltl).tl[j].lineages, world->numpop, is_same(j,T-1),
					   world);
		    
		    // calculate skyline values
		    if(world->options->skyline)
		      {
			if((*ltl).tl[j].eventnode->visited)
			  {
			    calculate_expected_values(bb->eventbins, bb->eventbinnum, bb->eventbinsize, t, 
						      (*ltl).tl[j].age,
						      (*ltl).tl[j].eventnode->type, 
						      (*ltl).tl[j].from, (*ltl).tl[j].to,
						      (*ltl).tl[j].lineages, world->numpop, world);
			  }	  
		      }
		  }
	      }
	  }
      }
}
