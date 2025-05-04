/* histogrammer for bayes histogram data
   takes the file bayesallfile and reads it into the histogram structure
   to use the calculate_hpd etc and also call the pretty printer functions
  
   (c) Peter Beerli 2007-2012

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
#include "definitions.h"
#include "migration.h"
#include "migrate_mpi.h"
#include "bayes.h"
#include "tools.h"
#include "sighandler.h"
#include "reporter.h"
#include "correlation.h"
#include "pretty.h"
#include "skyparam.h"
#include "speciate.h"
#include <errno.h>
#include <string.h>
#include <assert.h>

extern int myID;
extern int numcpu;
extern MYREAL scaling_prior(world_fmt *world, long numparam, MYREAL val);

// functions
#ifdef ZNZ
void read_bayes_fromfile(znzFile fmdimfile, world_fmt *world,option_fmt *options, char **files, long fnum);
void read_from_bayesmdim_minimal_info(znzFile mdimfile, world_fmt *world,option_fmt *options, data_fmt *data);
#else
void read_bayes_fromfile(FILE *fmdimfile, world_fmt *world,option_fmt *options, char **files, long fnum);
void read_from_bayesmdim_minimal_info(FILE *mdimfile, world_fmt *world,option_fmt *options, data_fmt *data);
#endif
long get_fullbinsum(MYREAL *lowerbound, MYREAL *upperbound, world_fmt *world, option_fmt *options, long locus);
boolean checking_bayesallfile(world_fmt *world, option_fmt *options, data_fmt *data, long ***unfinished);

// functions implementation
#ifdef ZNZ
void read_from_bayesmdim_minimal_info(znzFile mdimfile, world_fmt *world,option_fmt *options, data_fmt *data)
#else
void read_from_bayesmdim_minimal_info(FILE *mdimfile, world_fmt *world,option_fmt *options, data_fmt *data)
#endif
{
  char *input;
  //long nrep;
  long pop;
  long tmp;
  boolean done=FALSE;
  boolean recordedusem=TRUE;

  input = (char *) mycalloc(LINESIZE , sizeof(char));
  while(done==FALSE && FGETS3(input,LINESIZE,mdimfile) != EOF)
    {
      if(input[0] == '#' && strstr(input,"begin"))
	{
	  FGETS3(input,LINESIZE,mdimfile);
	  printf("%i>>>>>>> read from bayesallfile <<<<<<<<<<<<<<<<<<<<<<\n",myID);
	  options->custm = (char *) myrealloc(options->custm, sizeof(char) * (strlen(input)+1));
	  strcpy(options->custm,input+9);
	  printf("%i> custom       = %s\n",myID, options->custm);
	  FGETS3(input,LINESIZE,mdimfile);
	  options->custm2 = (char *) myrealloc(options->custm2, sizeof(char) * (strlen(input)+1));
	  strcpy(options->custm2,input+9);
	  printf("%i> custom2      = %s\n",myID, options->custm2);
	  FGETS3(input,LINESIZE,mdimfile);
	  sscanf (input+3, "%li %li %li %li %li %i", &world->loci, &world->numpop,
		  &world->numpop2, &tmp, &options->replicatenum,&recordedusem);
	  printf("%i> loci         = %li\n",myID, world->loci);
	  printf("%i> numpop       = %li\n",myID, world->numpop);
	  printf("%i> numpop^2     = %li\n",myID, world->numpop2);
	  printf("%i> replicate    = %li\n",myID, tmp);
	  printf("%i> replicatenum = %li\n",myID, options->replicatenum);
	  printf("%i> use_M        = %li\n",myID, (long) recordedusem);
	  // fill some more...
	  data->numpop = world->numpop;
	  options->newpops_numpop = world->numpop;
	  //xcode   nrep = options->replicatenum;
	  //xcode   if (nrep == 0)
	  //xcode     nrep = 1;
	  options->replicate = (boolean) tmp;
	  if(options->usem != recordedusem)
	    options->recordedusem = recordedusem;
	  data->popnames = (char **) mymalloc (sizeof (char *) * (size_t) world->numpop);
        for (pop = 0; pop < world->numpop; pop++)
	  {
            data->popnames[pop] = (char *) mycalloc (1, sizeof (char) * LINESIZE);
	    FGETS3(input,LINESIZE,mdimfile);
	    sscanf (input+3, "%s", data->popnames[pop]);
	    //printf("%i> population = %s\n",myID, data->popnames[pop]);
	  }
	done=TRUE;
	}
    }
  options->muloci = data->loci = world->loci;
  data->skiploci =
    (boolean *) myrealloc (data->skiploci,
			   sizeof (boolean) * (size_t) (data->loci + 1));
  memset (data->skiploci, 0, sizeof (boolean) * (size_t) (data->loci + 1));
  data->numpop = world->numpop;
  printf("%i>>>>>>> end read from bayesallfile <<<<<<<<<<<<<<<<<<\n",myID);
  myfree(input);
}
			    

long get_fullbinsum(MYREAL *lowerbound, MYREAL *upperbound, world_fmt *world, option_fmt *options, long locus)
{
  (void) locus;
  long temp=0;
  long i;
  long n=world->numpop2;
  for(i = 0; i < world->numpop; i++)
    {
      temp += options->bayes_priors[i].bins;
      lowerbound[i] = options->bayes_priors[i].min;
      upperbound[i] = options->bayes_priors[i].max;
    }
  for(i = world->numpop; i < world->numpop2; i++)
    {
      temp += options->bayes_priors[i].bins;
      lowerbound[i] = options->bayes_priors[i].min;
      upperbound[i] = options->bayes_priors[i].max;
    }
  if (world->bayes->mu)
    {
      i = world->numpop2;
      n += 1;
      temp += options->bayes_priors[i].bins;
      lowerbound[i] = options->bayes_priors[i].min;
      upperbound[i] = options->bayes_priors[i].max;
    }
  if (world->has_speciation)
    {
      n +=  2* world->species_model_size;
      for(i = world->numpop2+world->bayes->mu; i < world->numpop2+world->bayes->mu+ 2* world->species_model_size; i++)
	{
	  temp += options->bayes_priors[i].bins;
	  lowerbound[i] = options->bayes_priors[i].min;
	  upperbound[i] = options->bayes_priors[i].max;
	}
    }
  if(world->has_growth)
    {
      for(i = n; i < n+world->grownum; i++)
	{
	  temp += options->bayes_priors[i].bins;
	  lowerbound[i] = options->bayes_priors[i].min;
	  upperbound[i] = options->bayes_priors[i].max;
	}

    }
  return temp;
}

///
/// check_bayesallfile 
/// used for checkpointing, to figure out what was done and what failed
///
boolean checking_bayesallfile(world_fmt *world, option_fmt *options, data_fmt *data, long ***unfinished)
{
  boolean done=FALSE;
  boolean recover_needed = FALSE;
  long locus, replicate;
  char * input = (char *) mycalloc(SUPERLINESIZE, sizeof(char));
  char *inptr = input;
  //FILE * mdimfile = fopen(options->bayesmdimfilename,"r");
  read_from_bayesmdim_minimal_info(world->bayesmdimfile, world, options, data);
  long repmax = number_replicates2(options);
  intvec2d(unfinished, world->loci,repmax);
  //long mysize = world->loci * repmax;

#ifdef ZNZ
  znzFile mdimfile = world->bayesmdimfile;
#else
  FILE *mdimfile = world->bayesmdimfile;
#endif
  while(done==FALSE && FGETS3(input,SUPERLINESIZE,mdimfile) != EOF)
      {
	// grab the commentlines
	if (input[0] == '#' || input[0]=='S' || input[0]=='\0')
	  continue;
	inptr = input;	      
	/*step       =*/ atol(strsep(&inptr,"\t"));
	locus      = atol(strsep(&inptr,"\t"))-1;
	if(locus < 0)
	  error("help");
	replicate  = atol(strsep(&inptr,"\t"))-1;
	(*unfinished)[locus][replicate] += 1;
      }
  // report 
  printf("Checking bayesallfile, to find last loci and replicates worked on\n");
  long maxsample1 = options->lsteps - 1;
  for (locus=0; locus < world->loci; locus++)
    {
      for (replicate=0; replicate < repmax; replicate++)
	{
	  long xx = (*unfinished)[locus][replicate];
	  printf("%5li ", xx);
	  if (xx < maxsample1)
	    recover_needed = TRUE;
	} 
      printf("\n");
    }
  printf("\n");
  //
  myfree(input);
  if(recover_needed)
    {
#ifdef ZNZ
      znzclose(mdimfile);
      world->bayesmdimfile = znzopen(options->bayesmdimfilename, "a", options->use_compressed);      
      znzwrite("\n", (size_t) 1, (size_t) sizeof(char), world->bayesmdimfile);
#else
      fclose(mdimfile);
      world->bayesmdimfile = fopen(options->bayesmdimfilename, "a");
      FPRINTF(world->bayesmdimfile,"\n");
#endif
    }
  else
    {
#ifdef ZNZ
      znzclose(mdimfile);
      world->bayesmdimfile = znzopen(options->bayesmdimfilename, "r", options->use_compressed);      
#else
      fclose(mdimfile);
      world->bayesmdimfile = fopen(options->bayesmdimfilename, "r");
#endif
    }
  return recover_needed;
}  


///
/// read bayesallfile from disk and creates all the needed parts to 
/// recreate the output and pdf output
#ifdef ZNZ
void read_bayes_fromfile(znzFile fmdimfile, world_fmt *world,option_fmt *options, char **files, long fnum)
#else
  void read_bayes_fromfile(FILE *fmdimfile, world_fmt *world,option_fmt *options, char **files, long fnum)
#endif
{
  //  const long nn = world->numpop2 + world->bayes->mu * world->loci + 1;// One is for Log(Prob(Data|Model)
  long *n = NULL;
  long j0, j, z0, z;
  long step=0;
  long locus=0;
  //  long frompop;
  //long topop;
  long nnn=0;
  long t;
  long bin;
  //const long numpop = world->numpop;
  const long numpop2 = world->numpop2;
  const long np = world->numpop2 + world->bayes->mu + world->species_model_size * 2;
  const long npg = np + world->grownum;
  const long hc = world->options->heated_chains; 
  long numbins = 0;
  long numbinsall = 0;
  char *input;
  char *inptr;
  bayes_fmt * bayes = world->bayes;
  bayeshistogram_fmt *hist;
  MYREAL post;
  MYREAL like;
  MYREAL *params;
  MYREAL *delta = bayes->deltahist;
  MYREAL *lowerbound;
  MYREAL *upperbound;
  MYREAL *autocorrelation;
  MYREAL *ess;
  boolean *done;
  MYREAL *oldmeans;
  //MYREAL theta;
  //boolean notusem = !world->options->usem;
#ifndef ZNZ
  FILE *mdimfile = fmdimfile;
#else
  znzFile mdimfile = fmdimfile;
#endif
  long f;
  done = (boolean *) mycalloc(world->loci, sizeof(boolean));
  params = (MYREAL *) mycalloc((2 + npg), sizeof(MYREAL));
  oldmeans = (MYREAL *) mycalloc((2+ npg), sizeof(MYREAL));
  autocorrelation = (MYREAL *) mycalloc((2 * world->loci * npg), sizeof(MYREAL));
  ess = autocorrelation + world->loci * npg;
  lowerbound = (MYREAL *) mycalloc(npg, sizeof(MYREAL));
  upperbound = (MYREAL *) mycalloc(npg, sizeof(MYREAL));
  n = (long *) mycalloc((npg * world->loci), sizeof(long));
  input = (char *) mycalloc(SUPERLINESIZE, sizeof(char));
#ifdef DEBUG  
  printf("Begin reading the bayesallfile back into the system\n");
#endif
  // files has always at least one entry
  if (files!=NULL)
    {
      mdimfile = NULL;
    }
  else
    {
      error("Filelist should be filled\n");
    }
#ifdef ZNZ 
  unsigned long bytes = SUPERLINESIZE > ONEMEGABYTE ? SUPERLINESIZE : ONEMEGABYTE;
#endif
  for(f=0;f<fnum;f++)
    {
      input[0]='\0';
#ifdef ZNZ
      assert(files!=NULL);
      mdimfile = znzopen(files[f], "r", options->use_compressed);      
#else
	  mdimfile = fopen(files[f],"r");      
#endif
#ifdef DEBUG
	  printf("%s opened\n",files[f]);
#endif
	  if (mdimfile == NULL)
	    {
	      printf("errno = %s (%d).\n", strerror(errno), errno);
	      exit(1);
	    }
#ifdef ZNZ
      znzbuffer(mdimfile,bytes);
#endif
      while(FGETS3(input,SUPERLINESIZE,mdimfile) != EOF)
	{
	  // grab the commentlines
	  if (input[0] == '#' || input[0]=='S' || input[0]=='\0')
	    continue;
	  inptr = input;	      
	  step       = atol(strsep(&inptr,"\t"));
	  //printf("%li\n",step);
	  locus      = atol(strsep(&inptr,"\t"))-1;
	  if(locus == -1)
		error("help");
	      //  replicate  = atol(strsep(&inptr,"\t"))-1;
	      if (inptr!=NULL)
		(void) strsep(&inptr,"\t");
	      else
		continue;
	      if (inptr!=NULL)
		{
		  post       = atof(strsep(&inptr,"\t"));
		  params[0] = post;
		}
	      else
		continue;
	      if (inptr!=NULL)
		{
		  like       = atof(strsep(&inptr,"\t"));
		  params[1] = like;
		}
	      else
		continue;
	      if (inptr!=NULL)
		//  probg      = atof( strsep(&inptr,"\t"));
		(void) strsep(&inptr,"\t");
	      else
		continue;
	      if (inptr!=NULL)
		//  prior      = atof( strsep(&inptr,"\t"));
		(void) strsep(&inptr,"\t");
	      else
		continue;
	      if (inptr!=NULL)
		//  T = atol(strsep(&inptr,"\t"))+1;
		(void) strsep(&inptr,"\t");
	      else
		continue;
	      if (inptr!=NULL)
		//  treelength = atof(strsep(&inptr,"\t"));
		(void) strsep(&inptr,"\t");
	      else
		continue;

	      if(!done[locus])
		{
		  
		  done[locus] = TRUE;
		  // allocate the number of bins for the histogram
		  bayes->histogram[locus].binsum = get_fullbinsum(lowerbound, upperbound, world, options, locus);
		  bayes->histogram[locus].results = (double *) mycalloc(bayes->histogram[locus].binsum + 1, sizeof(double));
		  bayes->histogram[locus].set95 = (char *) mycalloc(bayes->histogram[locus].binsum* 2 + 2, sizeof(char));
		  bayes->histogram[locus].set50 = world->bayes->histogram[locus].set95 + bayes->histogram[locus].binsum + 1;
		  if(bayes->histogram[locus].covariance==NULL)
		    doublevec2d(&bayes->histogram[locus].covariance,npg,npg);
		}
	      hist = &bayes->histogram[locus];
	      numbinsall = 0;
	      n[locus] += 1; // we use the same n for all variables
	      for(j0=0;j0 < numpop2; j0++)
		{
		  if(shortcut(j0,world,&j))
		    {
		      continue;
		    }
		  if (inptr!=NULL)
		    params[j+2] =  atof(strsep(&inptr,"\t"));
		  else
		    continue;
		  // n[j] += 1;
		  oldmeans[j] = hist->means[j];
		  hist->means[j] += (params[j+2] - hist->means[j]) / n[locus];
		  numbinsall += hist->bins[j];
		  numbins = numbinsall - hist->bins[j];
		  
		  if ((upperbound[j] - params[j+2]) < -EPSILON)
		    {
		      warning("above upper bound: %f\n",params[j+2]);
		      continue;
		    }
		  bin = (long) ((params[j+2]-lowerbound[j]) / delta[j]);
		  hist->minima[j0] = lowerbound[j];
		  hist->maxima[j0] = upperbound[j];
		  hist->results[numbins + bin] += 1.;
		  bayes->histtotal[locus * npg + j] += 1;
		}
	      if(bayes->mu && j0==numpop2)
		{
		  numbins = numbinsall;
		  if (inptr!=NULL)
		    params[j0+2] = atof(strsep(&inptr,"\t"));
		  else
		    continue;
		  //n[j0+locus] += 1;
		  hist->means[j0] += (params[j0+2] - hist->means[j0]) / n[locus];//n[j0+locus];
		  bin = (long) ((params[j0+2]-lowerbound[j0]) / delta[j0]); 
		  hist->minima[j0] = lowerbound[j0];
		  hist->maxima[j0] = upperbound[j0];
		  hist->results[numbins + bin] += 1.;
		  bayes->histtotal[locus * npg + j0] += 1;
		}
	      if(world->has_speciation)
		{
		  for(j0=numpop2+bayes->mu;j0 < numpop2+bayes->mu+world->species_model_size*2; j0++)
		    {
		      if(shortcut(j0,world,&j))
			{
			  continue;
			}
		      if(inptr!=NULL)
			params[j+2] =  atof(strsep(&inptr,"\t"));
		      else
			continue;
		      //n[j] += 1;
		      oldmeans[j] = hist->means[j];
		      hist->means[j] += (params[j+2] - hist->means[j]) / n[locus];//n[j];
		      numbinsall += hist->bins[j];
		      numbins = numbinsall - hist->bins[j];
		      
		      if (params[j+2]>upperbound[j])
			{
			  warning("above upper bound: %f\n",params[j+2]);
			  continue;
			}
		      bin = (long) ((params[j+2]-lowerbound[j]) / delta[j]);
		      hist->minima[j0] = lowerbound[j];
		      hist->maxima[j0] = upperbound[j];
		      hist->results[numbins + bin] += 1.;
		      bayes->histtotal[locus * npg + j] += 1;
		    }
		}
	      if (world->has_growth)
		{
		  long grownum = world->options->growpops_numalloc;
		    for(j0=0; j0 < grownum; j0++)
		      {
			long pick = world->options->growpops[j0];
			if (pick == 0)
			  continue;
			else
			  j = pick + np - 1;
			if(inptr!=NULL)
			  params[j+2] =  atof(strsep(&inptr,"\t"));
			else
			  continue;
			//n[j] += 1;
			oldmeans[j] = hist->means[j];
			hist->means[j] += (params[j+2] - hist->means[j]) / n[locus];//n[j];
			numbinsall += hist->bins[j];
			numbins = numbinsall - hist->bins[j];
			
			if (params[j+2]>upperbound[j])
			  {
			    warning("above upper bound: %f\n",params[j+2]);
			    continue;
			  }
		      bin = (long) ((params[j+2]-lowerbound[j]) / delta[j]);
		      hist->minima[j] = lowerbound[j];
		      hist->maxima[j] = upperbound[j];
		      hist->results[numbins + bin] += 1.;
		      bayes->histtotal[locus * npg + j] += 1;
		    }
		}
	      hist->n = n[locus]; //assumes that all are the same (should be!)
	      for(j0=0;j0 < numpop2; j0++)
		{
		  if(shortcut(j0,world,&j))
		    continue;
		  else
		    {
		      for(z0=0;z0 < numpop2; z0++)
			{
			  if(shortcut(z0,world,&z))
			    continue;
			  else
			    {
			      hist->covariance[j][z] += (params[z+2] - hist->means[z]) * (params[j+2]-oldmeans[j]); 
			    }
			}
		    }
		}	  
	      //	      if(world->options->datatype == 'g')
	      if(options->checkpointing)
		{
		  for(t=0;t<hc;t++)
		    {
		      if(inptr!=NULL)
			world->bf[locus * hc + t] = atof(strsep(&inptr,"\t"));
		      else
			continue;
		    }
		  // dummy read of thermo sum up to this point
		  //lsum = atof(strsep(&inptr,"\t"));
		  if(inptr!=NULL)
		    (void ) strsep(&inptr,"\t");
		  else
		    continue;
		  // harmonic mean: scaler contains the log value, hm contains 1.
		  if(inptr!=NULL)
		    {
		      world->hmscale[locus] = atof(strsep(&inptr,"\t"));
		      world->hm[locus] = 1.; 
		    }
		  else
		    continue;
		  //calculate_ess_frombayes (world, step, params, locus, autocorrelation, ess);
		  //covariance_bayes(world,locus);
		}
	      else
		{
		  if(inptr!=NULL)
		    (void) strsep(&inptr,"\n");
		  else
		    continue;		  
		}
	  }
#ifdef ZNZ
      znzclose(mdimfile);
#else
      fclose(mdimfile);
#endif
      /*long i;
      for(i=0;i<world->loci;i++)
	{
	  for(j0=0;j0<world->numpop2+world->bayes->mu+2*world->species_model_size;j0++)
	    {
	      if(shortcut(j0,bayes,&j))
		{
		  continue;
		}	      
	      printf("@@2@@@ locus=%li parameter=%li mean=%f \n",i, j0, bayes->histogram[locus].means[j]);
	    }
	}
      */
      if(options->mdimdelete)
	{
	  remove(options->bayesmdimfilename);
	}
      // skyline parameters if present
      if (world->options->skyline_param)
	skyline_param_reader(world, step, locus, &inptr);
      // end skyline parameters
#ifdef DEBUG
      fprintf(stdout,"%s closed\n",files[f]);
#endif
      mdimfile=NULL;
    }
#ifdef DEBUG
  printf("End reading the bayesallfile back into the system\n");
#endif
  if(world->options->datatype == 'g')
    {
      // reset the archiving machinery
      memset(world->auto_archive,0, sizeof(MYREAL) * (size_t) (2 * (npg + 1)));
      nnn = 1;
	     //for(j=0;j<world->loci;j++)
	     //{
	  
	  for(t=0;t<npg; t++)
	    {
	      // onepass mean of autocorrelation
	      world->auto_archive[t] += (autocorrelation[t] - world->auto_archive[t])/nnn;
	      // summing ess values
	      world->ess_archive[t] += ess[t];
	      //printf("j=%li t=%li %f\n", j, t, world->ess_archive[t]);
	    }
	  nnn++;
	     //}
    }
  myfree(params);
  myfree(oldmeans);
  myfree(autocorrelation);
  myfree(lowerbound);
  myfree(upperbound);
  myfree(n);
  myfree(done);
  myfree(input);
}
