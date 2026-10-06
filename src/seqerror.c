// sequencing error estimation for each nucleotide independently
// Fall 2012
// PB
#include "migration.h"
#include "mcmc.h"
#include "random.h"
#include "tools.h"
#include "tree.h"
#include "sighandler.h"
#include "sequence.h"
#include "migrate_mpi.h"
#include "seqerror.h"
#include "bayes.h"
extern int myID;

//#define INDIX(a,b,c) ((a)*(b)+(c))

// functions
//void fill_world_seqerror(world_fmt *world, option_fmt *options);
void destroy_seqerror(world_fmt* world);
void change_freq_tip(world_fmt *world, node *tip, MYREAL *errorrates);
//void change_freq(world_fmt *world);
//##


// function implementations
void fill_world_seqerror(world_fmt *world, option_fmt *options)
{
  long locus;
  long mult = options->seqerrorcombined ? 1 : 4;
  world->seqerrorcombined = options->seqerrorcombined;
  world->seqerrorallocnum = (long *) mycalloc(world->loci,sizeof(long));
  world->seqerrorratesnum = (long *) mycalloc(world->loci,sizeof(long));
  world->seqerrorcount = (long *) mycalloc(world->loci,sizeof(long));
  world->seqerrorsteps = (long *)  mycalloc(world->loci, sizeof(long));
  world->seqerrorrates = (MYREAL **)  mycalloc(world->loci, sizeof(MYREAL *));
  for (locus=0;locus<world->loci; locus++)
    {
      world->seqerrorcount[locus] = 0;
      world->seqerrorratesnum[locus] = 1;
      world->seqerrorallocnum[locus] = 10;
      world->seqerrorrates[locus] = (MYREAL *) mycalloc((world->seqerrorallocnum[locus] * mult), sizeof(MYREAL));
      world->seqerrorrates[locus][0] = options->seqerror[0];
      /* four rates only without combining (the condition was inverted: the
         separate rates started at 0, the combined one wrote unused slots) */
      if(!world->seqerrorcombined)
	{ 
	  world->seqerrorrates[locus][1] = options->seqerror[1]; 
	  world->seqerrorrates[locus][2] = options->seqerror[2];
	  world->seqerrorrates[locus][3] = options->seqerror[3];
	}
    }
}

void destroy_seqerror(world_fmt* world)
{
  long locus;
  myfree(world->seqerrorallocnum);
  myfree(world->seqerrorratesnum);
  myfree(world->seqerrorcount);
  myfree(world->seqerrorsteps);
  for (locus=0;locus<world->loci; locus++)
    {
      myfree(world->seqerrorrates[locus]);
    }
  myfree(world->seqerrorrates);
}


void change_freq_tip(world_fmt *world, node *tip, MYREAL *errorrates)
{
  long sublocus;
  long k, l;
  //long xs;
  const long locus = world->locus;
  mutationmodel_fmt * s;
  
  const long sublocistart = world->sublocistarts[locus];
  const long sublociend   = world->sublocistarts[locus+1];
  for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
    {
      s = &world->mutationmodels[sublocus];
      const long xs = sublocus - sublocistart;
      
      long numpatterns = s->numpatterns;
      for (k = 0; k < numpatterns; k++)
	{
	  //j = s->alias[k] - 1;
	  for (l = 0; l < s->numsiterates; l++)
	    {
	      set_nucleotide(tip->x[xs].s[k][l], tip->sequence[k],errorrates);
	    }
	}
    }
}




// log density of the Beta(1,10) prior of a sequencing error rate, up to a
// constant: mean 1/11, highest at 0
static MYREAL seqerror_logprior(MYREAL e)
{
  if (e <= 0.0 || e >= 1.0)
    return (MYREAL) -HUGE;
  return 9.0 * log(1.0 - e);
}

/// Metropolis-Hastings update of a sequencing error rate at a fixed
/// genealogy: one rate (the combined one, or one of the four nucleotide
/// rates) is moved by a reflected window on [0,1], the tips are rewritten,
/// and the move is accepted with the heated data-likelihood ratio times the
/// Beta(1,10) prior ratio. It used to draw from Beta(10,1) (mean 0.91), had
/// no prior or proposal term, rode on a tree move and restored only one of
/// the four rates after a rejection.
void change_freq(world_fmt *world)
{
  const boolean is_combined = world->seqerrorcombined;
  const long mult = is_combined ? 1 : 4;
  const long seqpn = mult;
  const long locus = world->locus;
  const long sumtips = world->sumtips;
  const MYREAL width = 0.05;
  node ** tips = world->nodep;
  MYREAL cur[4], prop[4];
  MYREAL oldlike, newlike, h, e, enew;
  long i, k, end, before;
  boolean success;

  end    = world->seqerrorratesnum[locus]*seqpn;
  before = end - seqpn;
  for (i = 0; i < 4; i++)
    cur[i] = world->seqerrorrates[locus][before + (is_combined ? 0 : i)];
  memcpy(prop, cur, sizeof(cur));
  k = is_combined ? 0 : random_integer(0,3);
  e = cur[k];
  enew = e + (UNIF_RANDUM() - 0.5) * width;
  if (enew < 0.0)
    enew = -enew;
  if (enew > 1.0)
    enew = 2.0 - enew;
  if (is_combined)
    prop[0] = prop[1] = prop[2] = prop[3] = enew;
  else
    prop[k] = enew;
  h = seqerror_logprior(enew) - seqerror_logprior(e);

  oldlike = world->likelihood[world->G];
  for (i = 0; i < sumtips; i++)
    change_freq_tip(world, tips[i], prop);
  set_all_dirty(world->root->next, crawlback (world->root->next), world, world->locus);
  first_smooth(world,world->locus);
  newlike = treelikelihood(world);
  if (world->options->prioralone)
    success = bayes_accept(0.0, 0.0, world->heat, h);
  else
    success = bayes_accept(newlike, oldlike, world->heat, h);
  if (success)
    {
      world->likelihood[world->G] = newlike;
      memcpy(cur, prop, sizeof(cur));
      world->seqerrorcount[locus] += 1;
    }
  else
    {
      for (i = 0; i < sumtips; i++)
        change_freq_tip(world, tips[i], cur);
      set_all_dirty(world->root->next, crawlback (world->root->next), world, world->locus);
      first_smooth(world,world->locus);
      world->likelihood[world->G] = oldlike;
    }

  // record the current rates as before: a new row every increment of the
  // cold chain while sampling, otherwise the last row holds the state
  if(!world->in_burnin)
    world->seqerrorsteps[locus] += 1;
  if (!world->in_burnin && world->cold && world->seqerrorsteps[locus] % world->increment == 0)
    {
      if (world->seqerrorratesnum[locus] + 100 > world->seqerrorallocnum[locus])
        {
          world->seqerrorallocnum[locus] += 100;
          world->seqerrorrates[locus]= (MYREAL *) myrealloc(world->seqerrorrates[locus],sizeof(MYREAL) * (size_t) (mult*world->seqerrorallocnum[locus]));
        }
      end = world->seqerrorratesnum[locus]*seqpn;
      for (i = 0; i < seqpn; i++)
        world->seqerrorrates[locus][end + i] = cur[i];
      world->seqerrorratesnum[locus] += 1;
    }
  else
    {
      end    = world->seqerrorratesnum[locus]*seqpn;
      before = end - seqpn;
      for (i = 0; i < seqpn; i++)
        world->seqerrorrates[locus][before + i] = cur[i];
    }
}


void seqerror_report(world_fmt *world, char *seqerrorfile)
{
  boolean is_combined = world->seqerrorcombined;
  long i;
  long locus;
  // overall loci
  long na=0, nc=0, ng=0, nt=0;
  MYREAL fa=0.0, fc=0.0, fg=0.0, ft=0.0;
  MYREAL sfa=0.0, sfc=0.0, sfg=0.0, sft=0.0;
  long na1=0, nc1=0, ng1=0, nt1=0;
  MYREAL fa1=0.0, fc1=0.0, fg1=0.0, ft1=0.0;
  MYREAL sfa1=0.0, sfc1=0.0, sfg1=0.0, sft1=0.0;
  FILE *file = NULL;
  long ii;
  if (seqerrorfile!=NULL)
    file = fopen(seqerrorfile,"w");
  
  fprintf(world->outfile,"\n\nEstimation of sequencing error\n");
  fprintf(world->outfile,"------------------------------------------------\n\n");
  fprintf(world->outfile,"Locus   Nucleotide   Mean error    Standard deviation    n\n");
  fprintf(world->outfile,"-----------------------------------------------------------\n");
  
  if(is_combined)
    {
      onepass_mean_std_start(&fa, &sfa, &na);
      for (locus=0;locus<world->loci;locus++)
	{
	  onepass_mean_std_start(&fa1, &sfa1, &na1);
	  for (i=0;i<world->seqerrorratesnum[locus];i++)
	    {
	      onepass_mean_std_calc(&fa, &sfa, &na, world->seqerrorrates[locus][i]);
	      onepass_mean_std_calc(&fa1, &sfa1, &na1, world->seqerrorrates[locus][i]);
	      if (file)
		{
		  fprintf(file,"%li %f \n", locus+1, world->seqerrorrates[locus][i]);
		}
	    }
	  onepass_mean_std_end(&fa1, &sfa1, &na1);
	  fprintf(world->outfile,"% 5li     N        %4.6f            %4.6f    %8li\n", locus+1,fa1, sfa1, na1);
	}
      onepass_mean_std_end(&fa, &sfa, &na);
      fprintf(world->outfile,"%5.5s     N        %4.6f            %4.6f    %8li\n", "All" ,fa, sfa, na);
    }
  else
    {
      onepass_mean_std_start(&fa, &sfa, &na);
      onepass_mean_std_start(&fc, &sfc, &nc);
      onepass_mean_std_start(&fg, &sfg, &ng);
      onepass_mean_std_start(&ft, &sft, &nt);
      for (locus=0;locus<world->loci;locus++)
	{
	  onepass_mean_std_start(&fa1, &sfa1, &na1);
	  onepass_mean_std_start(&fc1, &sfc1, &nc1);
	  onepass_mean_std_start(&fg1, &sfg1, &ng1);
	  onepass_mean_std_start(&ft1, &sft1, &nt1);
	  for (i=0;i<world->seqerrorratesnum[locus];i++)
	    {
	      ii = i * 4;
	      onepass_mean_std_calc(&fa, &sfa, &na, world->seqerrorrates[locus][ii]);
	      onepass_mean_std_calc(&fc, &sfc, &nc, world->seqerrorrates[locus][ii+1]);
	      onepass_mean_std_calc(&fg, &sfg, &ng, world->seqerrorrates[locus][ii+2]);
	      onepass_mean_std_calc(&ft, &sft, &nt, world->seqerrorrates[locus][ii+3]);		
	      onepass_mean_std_calc(&fa1, &sfa1, &na1, world->seqerrorrates[locus][ii]);
	      onepass_mean_std_calc(&fc1, &sfc1, &nc1, world->seqerrorrates[locus][ii+1]);
	      onepass_mean_std_calc(&fg1, &sfg1, &ng1, world->seqerrorrates[locus][ii+2]);
	      onepass_mean_std_calc(&ft1, &sft1, &nt1, world->seqerrorrates[locus][ii+3]);		
	      if (file)
		{
		  fprintf(file,"%li %f %f %f %f\n", locus, world->seqerrorrates[locus][ii],
			  world->seqerrorrates[locus][ii+1],world->seqerrorrates[locus][ii+2],
			  world->seqerrorrates[locus][ii+4]);
		}
	    } 
	  onepass_mean_std_end(&fa1, &sfa1, &na1);
	  onepass_mean_std_end(&fc1, &sfc1, &nc1);
	  onepass_mean_std_end(&fg1, &sfg1, &ng1);
	  onepass_mean_std_end(&ft1, &sft1, &nt1);
	  fprintf(world->outfile,"% 5li         A        %4.6f        %4.6f    %8li\n", locus+1,fa1, sfa1, na1);
	  fprintf(world->outfile,"% 5li         C        %4.6f        %4.6f    %8li\n", locus+1,fc1, sfc1, nc1);
	  fprintf(world->outfile,"% 5li         G        %4.6f        %4.6f    %8li\n", locus+1,fg1, sfg1, ng1);
	  fprintf(world->outfile,"% 5li         T        %4.6f        %4.6f    %8li\n", locus+1,ft1, sft1, nt1);

	}
      onepass_mean_std_end(&fa, &sfa, &na);
      onepass_mean_std_end(&fc, &sfc, &nc);
      onepass_mean_std_end(&fg, &sfg, &ng);
      onepass_mean_std_end(&ft, &sft, &nt);
      fprintf(world->outfile,"  All         A        %4.6f        %4.6f    %8li\n", fa, sfa, na);
      fprintf(world->outfile,"  All         C        %4.6f        %4.6f    %8li\n", fc, sfc, nc);
      fprintf(world->outfile,"  All         G        %4.6f        %4.6f    %8li\n", fg, sfg, ng);
      fprintf(world->outfile,"  All         T        %4.6f        %4.6f    %8li\n", ft, sft, nt);
    }
  FClose(file);
}


///
/// save all population assignments from the worker nodes
void
get_seqerror (world_fmt * world, option_fmt * options)
{
#ifdef MPI
    long maxreplicate = (options->replicate
                         && options->replicatenum >
                         0) ? options->replicatenum : 1;
    
    if (myID == MASTER && world->has_estimateseqerror)
    {
        mpi_results_master (MIGMPI_SEQERROR, world, maxreplicate,
                            unpack_seqerror_buffer);
    }
#else
    (void) world;
    (void) options;
#endif
}
