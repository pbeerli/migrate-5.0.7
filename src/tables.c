// Tables printing in migrate
// MIT opensource license
// consolidation from other files to simplify future changes
// (c) Peter Beerli 2021
//

#include "migration.h"
#include "bayes.h"
#include "marginallike.h"
#include "pretty.h"
#include "reporter.h"
#include "migrate_mpi.h"

void print_bayesfactor(world_fmt **universe, option_fmt * options);
void print_burnin_autostop(world_fmt * world);
void print_heatingreport(world_fmt **universe, option_fmt* options);
void print_mcmc_run_character(world_fmt * world);
void print_finaltree(world_fmt *world, option_fmt *options);

void print_bayesfactor(world_fmt **universe, option_fmt * options)
{
  world_fmt * world = universe[0];//EARTH
  FILE * file = world->outfile;
  long t;
  const long hc = world->options->heated_chains;
  long locus;
  MYREAL heat0 = 1.;
  MYREAL heat1 = 1.;
  MYREAL heat2 = 1.0;
  MYREAL val0  = 0.;
  MYREAL val1  = 0.;
  MYREAL sval0  = 0.;
  MYREAL sval1  = 0.;
  MYREAL val2  = 0.;
  MYREAL bfsum = 0.;
  MYREAL bfsum2 = 0.;
  MYREAL approxlsum = 0.;
  MYREAL hsum = 0.;
  MYREAL lsum;
  MYREAL lsum0;
  MYREAL ratio = 0.0;
  MYREAL sratio = 0.0;
  MYREAL allratio = 0.0;
  MYREAL sallratio = 0.0;
  MYREAL scaling_factor = 0.0;
  MYREAL *locusweight = world->data->locusweight;//invariant loci treatment
  double min_bf = -HUGE;
  double min_bfb = -HUGE;
  double min_ratio = -HUGE;
  double max_bf = -HUGE;
  double max_bfb = -HUGE;
  double max_ratio = -HUGE;
    
  // calculate the harmonic mean score
  for(locus=0;locus < world->loci; locus++)
    {
      if(world->data->skiploci[locus])
	continue;
      hsum +=  locusweight[locus] * (world->hmscale[locus] - log (world->hm[locus]));//invariant loci treatment
    }

  fprintf(file,"\n\n\nLog-Probability of the data given the model (marginal likelihood = log(P(D|thisModel))\n");
  fprintf(file,"--------------------------------------------------------------------\n[Use this value for Bayes factor calculations:\nBF = Exp[log(P(D|thisModel) - log(P(D|otherModel)]\nshows the support for thisModel]\n\n");
  if(options->heating)
    {
      pdf_bayes_factor_header(world,options);
      //fprintf(file,"\n\nLocus          TI(1a)       BTI(1b)         SS(2)         HS(3)\n");
      //fprintf(file,"---------------------------------------------------------------\n");
      fprintf(file,"\n\nLocus          TI(1a)       BTI(1b)         HS(2)\n");
      fprintf(file,    "-------------------------------------------------\n");
      pdf_bayes_factor_rawscores_header(world,options);
      allratio = 0.0;
      for(locus=0;locus<world->loci;locus++)
	{
	  if(world->data->skiploci[locus])
	    continue;
	  ratio = 0.0;
	  sratio = 0.0;
	  lsum = 0.;
	  lsum0 = 0.;
	  for(t=1; t < hc; t++)
	    {
	      heat2 = heat0;
	      val2 = val0;
#ifdef MPI
	      heat0 = 1./options->heat[t-1];
	      heat1 = 1./options->heat[t];
#else
	      if(options->checkpointing)
		{
		  if(options->adaptiveheat!=NOTADAPTIVE)
		    {
		      heat0 = 1./ universe[t-1]->averageheat;
		      heat1 = 1./ universe[t]->averageheat;
		    }
		  else
		    {
		      heat0 = 1./options->heat[t-1];
		      heat1 = 1./options->heat[t];
		    }
		}
	      else
		{
		  heat0 = 1./options->heat[t-1];
		  heat1 = 1./options->heat[t];
		}
#endif
	      //Simpson's rule and trapezoidal are the same when I only have function values at a and b
	      val0 = locusweight[locus] * world->bf[locus * hc + t-1];
	      val1 = locusweight[locus] * world->bf[locus * hc + t];
	      //printf("\"log mL:\", %i, %f, %f, %f, %f\n", myID, heat0, heat1, val0, val1); 
	      ratio += val0 - val1;

	      // stepping stones
	      sval0 = locusweight[locus] * log(world->steppingstones[locus * hc + t-1]) + world->steppingstone_scalars[locus * hc + t-1];
	      sval1 = locusweight[locus] * log(world->steppingstones[locus * hc + t]) + world->steppingstone_scalars[locus * hc + t];
	      //printf("\"log mL:\", %i, %f, %f, %f, %f\n", myID, heat0, heat1, sval0, sval1); 
	      sratio += sval0 - sval1;
	      //if(isnan(sratio))
	      //{
	      //  fprintf(stderr,"Steppingstone calculation failed:\n");
	      //  fprintf(stderr,"[%li,%li] sratio= (log(%f)+%f)-(log(%f)+%f=NaN\n",locus,t,
	      //	  world->steppingstones[locus * hc + t-1],
	      //	  world->steppingstone_scalars[locus * hc + t-1],
	      //	  world->steppingstones[locus * hc + t],
	      //	  world->steppingstone_scalars[locus * hc + t]);
	      //}
	      //we keep last element to adjust for Bezier approximation
	      lsum0 = (heat0 - heat1) * (val0 + val1) * 0.5;
	      lsum += lsum0; 
	    }
	  //	  (x2 y1 - x1 y2)/(x1 - x2)
	  // this last addition to the lsum calculates the chunk between the last temperature and 
	  // the infinitely hot temperature as a linear approximation of the the second hottest temperature
	  // this is certainly rough, but in simulations with 3 populations one can see that with large number
	  // of temperatures this looks reasonable, and one can approximate the integral more accurately with 
	  // with only a few columns.
	  // using Bezier to approximate nice curve between the last two points to mimick the curve that
	  // can be found with 16 or 32 heated chains, handle points are calculated using adhoc decisions
	  // (comparison with 32 heated chains) using 0.8 of the interval for handle_y1 and the intercept
	  // between the first and the second last point to calculate the handle_y2 see sumbezier()
	  // for implementation. Currently this is not tunable.
	  //	  approxlsum = sumbezier(100, heat1, world->bf[locus * hc + t-1], 
	  //			 heat0, world->bf[locus * hc + t-2], 
	  //			 1.0, world->bf[locus * hc]);
	  MYREAL ratio2=0.0;
	  approxlsum = sumbezier(100L, heat1, val1, 
				 heat0, val0, 
				 heat2, val2, &ratio2);
	  ratio += val1;

	  //fprintf(file,"  %5li  %12.2f  %12.2f  %12.2f  %12.2f\n", locus + 1, lsum, lsum-lsum0+approxlsum,
	  //	  sratio, ratio);
	  if (options->tersepdf)
	    {
	      if (min_bf > (lsum - lsum0+approxlsum))
		{
		  min_bf  = lsum;
		  min_bfb = lsum - lsum0+approxlsum;
		  min_ratio = ratio;
		}
	      if (max_bf < (lsum - lsum0+approxlsum))
		{
		  max_bf  = lsum;
		  max_bfb = lsum - lsum0+approxlsum;
		  max_ratio = ratio;
		}
	      //lsum, lsum-lsum0+approxlsum, ratio 
	    }
	  else
	    {
	      fprintf(file,"  %5li  %12.2f  %12.2f  %12.2f\n", locus + 1, lsum, lsum-lsum0+approxlsum,
		      ratio); 
	      pdf_bayes_factor_rawscores(locus, lsum, lsum-lsum0+approxlsum, sratio, world->hmscale[locus] - log (world->hm[locus]));
	    }
	  bfsum2 += approxlsum + lsum - lsum0;  	  
	  bfsum += lsum; //+ world->bfscale[locus];
	  allratio += ratio;
	  sallratio += sratio;
	}

      if (options->tersepdf)
	{
	  fprintf(file,"  %s %12.2f  %12.2f  %12.2f\n", "Lowest", min_bf, min_bfb, min_ratio);
	  fprintf(file,"  %s %12.2f  %12.2f  %12.2f\n", "Highest", max_bf, max_bfb, max_ratio);
	  pdf_bayes_factor_rawscores_minmax(BFMIN, min_bf, min_bfb, min_ratio);
	  pdf_bayes_factor_rawscores_minmax(BFMAX, max_bf, max_bfb, max_ratio);
	}
      // print out "ALL" row for both multiloci and single loci run (the single locus run has the same
      // form as the multilocus run so that I can grep the results more easily for model comparison
      scaling_factor = combine_scaling_factor(world);
      bfsum  += scaling_factor;
      bfsum2 += scaling_factor;
      hsum   += scaling_factor;
      allratio += scaling_factor;
      sallratio += scaling_factor;
      if(world->loci>1)
	{
	  fprintf(file,"---------------------------------------------------------------\n");
	  //	  fprintf(file,"  All    %12.2f  %12.2f  %12.2f  %12.2f\n[Scaling factor = %f]\n",
	  //	  bfsum, bfsum2, sallratio, hsum, scaling_factor);
	  fprintf(file,"  All    %12.2f  %12.2f  %12.2f\n[Scaling factor = %f]\n",
		  bfsum, bfsum2,  hsum, scaling_factor);
	}
      if(world->loci>1)
	{
	  pdf_bayes_factor_rawscores(-1L, bfsum, bfsum2, sallratio,hsum);
	}
    }
  else   // -----------------------------------------not heating
    {
      fprintf(file,"No model selection evaluation is possible without heating\n");
    }
  fprintf(file,"\n\n(1a) TI: Thermodynamic integration: log(Prob(D|Model)): Good approximation with many temperatures\n");
  fprintf(file,"(1b) BTI: Bezier-approximated Thermodynamic integration: when using few temperatures USE THIS!\n");
  //future fprintf(file,"(2)  SS: Steppingstone Sampling (Xie et al 2011)\n");
  fprintf(file,"(2)  HS: Harmonic mean approximation: Overestimates the marginal likelihood, poor variance\n\n");
  
  pdf_bayes_factor_comment(world, scaling_factor);
}

void print_burnin_autostop(world_fmt * world)
{
  long z;
  FILE * file = world->outfile;
  worldoption_fmt * options = world->options;
  long maxreplicate = (options->replicate
		       && options->replicatenum >
		       0) ? options->replicatenum : 1;
  fprintf(file,"\n\n\nStop of burn-in phase due to convergence\n");
  fprintf(file,"----------------------------------------\n");
  switch(options->burnin_autostop)
    {
    case 'a':
      fprintf(file,"[Stopping criteria was: Variance ratio between the last two groups of 1000 steps < %f]\n\n",world->varheat); 
      break;
    case 't':
      fprintf(file,"[Stopping criteria was: reached prescribed acceptance ratio of %f]\n\n",world->options->autotune); 
      break;
    case 'e':
      fprintf(file,"[Stopping criteria was: All effective MCMC sample sizes > %f]\n\n",world->essminimum); 
      break;
    case ' ':
      fprintf(file,"\n\n");
      break;
    }
  
  fprintf(file,"Locus  Replicate  Steps  ESS*    Accept* Variance ratio (new/old variance)\n");
  fprintf(file,"-----  ---------  ------ ------- ------- ---------------------------------\n");
  for(z=0; z < world->loci * maxreplicate; z++)
    {
      burnin_record_fmt * bz = &world->burnin_stops[z];
      if(bz->oldvariance > 0.0)
	{
	  fprintf(file,"%5li  %5li  %10li   %6.1f %6.4f %10.4f (%f/%f)\n",1 + bz->locus,
		  bz->replicate,
		  bz->stopstep,
		  bz->ess,
		  bz->accept,
		  bz->variance/world->burnin_stops[z].oldvariance,
		  bz->variance,
		  bz->oldvariance);
	}
      //	  bz->worker = myID;
    }
  fprintf(file,"(*=worst)\n\n");
  pdf_burnin_stops(world, maxreplicate);
}

void print_heatingreport(world_fmt **universe, option_fmt * options)
{
  long locus;
  // if adaptive heating print a table with the average temperatures
  if(options->heating)
    {
      long t;
      world_fmt * world = universe[0];
      const long hc = world->options->heated_chains;
      if(options->adaptiveheat!=NOTADAPTIVE)
	{
	  fprintf(world->outfile,"\n\n\nAverage temperatures during the run using %s\n",
		  (options->adaptiveheat==STANDARD) ? "standard adaptive heating scheme" : "bounded adaptive heating scheme" );
	  fprintf(world->outfile,"===========================================================================\n\n");
	  fprintf(world->outfile,"Chain Temperature\n");
	  for(t = 0; t < options->heated_chains; t++)
	    {
	      fprintf(world->outfile,"%5li %10.5f\n",t+1,universe[t]->averageheat);
	    }
	  fprintf(world->outfile,"Adaptive heating often fails, if the average temperatures are very close together\n");
	  fprintf(world->outfile,"try to rerun using static heating! If you want to compare models using marginal\n");
	  fprintf(world->outfile,"likelihoods then you MUST use static heating\n");
	  pdf_print_averageheat(universe,options);
	}
      else
	{
	  // print heating table for static heating
	  fprintf(world->outfile,"\n\n\nTemperatures during the run using the standard heating scheme\n" );
	  fprintf(world->outfile,"===========================================================================\n\n");
	  fprintf(world->outfile,"Chain Temperature               log(marginal likelihood)  log(mL_steppingstone)\n");
	  // locus means indicator for chain
	  for(t = 0; t < options->heated_chains; t++)
	    {
	      double nloc=0.0;
	      double bfsum = 0.0;
	      double ssum = 0.0;
	      for(locus = 0; locus < world->loci; locus++)
		{
		  if(world->data->skiploci[locus])
		    {
		      continue;
		    }
		  else
		    {
		      nloc += world->data->locusweight[locus];
		    }
		  bfsum += world->data->locusweight[locus] * world->bf[locus * hc + t];
		  ssum += log(world->steppingstones[locus * hc + t]) + world->steppingstone_scalars[locus * hc + t];
		}
	      fprintf(world->outfile,"%5li %10.5f          %10.5f  %10.5f\n",t+1,universe[t]->heat, bfsum/nloc, ssum/nloc);
	    }
	  pdf_print_averageheat(universe,options);
	}
    }
}


void print_mcmc_run_character(world_fmt * world)
{
  // printing of MCMC run characteristics
  fprintf(world->outfile,"\n\nMCMC run characteristics\n");
  fprintf(world->outfile,"========================\n\n");
  bayes_print_hyperprior(world->outfile,world);
  bayes_print_accept(world->outfile,world);
  pdf_bayes_print_hyperpriors(world);
  pdf_bayes_print_accept(world);
  print_bayes_ess(world->outfile, world, world->auto_archive, world->ess_archive);
  if(world->options->progress)
    {
      print_bayes_ess(stdout,world, world->auto_archive, world->ess_archive);
    }
  pdf_bayes_print_ess(world);
  if(strchr("aet",world->options->burnin_autostop))
    {
      print_burnin_autostop(world);
    }
}


void print_finaltree(world_fmt *world, option_fmt *options)
{
  long locus;
  long i;
  if(options->treeprint)// && options->checkpointing)
    {
      if(world->options->treeinmemory)
	{
	  for(locus=0;locus<world->loci; locus++)
	    {
	      fprintf(world->treefile,"%s", world->treespace[locus]);
	    }
	}
#ifdef NEXUSTREE
	    FPRINTF(world->treefile,"\nend;\n");
	    //FPRINTF(world->treefile,"begin paup;\n");
	    //FPRINTF(world->treefile,"gsi /taxsets=(");
	    //for (i=0;i<world->numpop-1;i++)
	    //  FPRINTF(world->treefile,"deme%ld ",i);
	    //FPRINTF(world->treefile,"deme%ld",i);
	    //FPRINTF(world->treefile,") nperms=100000\n");
	    //FPRINTF(world->treefile,"end;\n\n");
#endif
    }
}

