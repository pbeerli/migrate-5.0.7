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
      if(world->seqerrorcombined)
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




// changing a frequency
void change_freq(world_fmt *world)
{
  boolean is_combined = world->seqerrorcombined;
  long mult = is_combined ? 1 : 4;
  long i;
  long count;
  long sumtips = world->sumtips;
  node ** tips = world->nodep;
  long locus = world->locus;
  // use a beta distribution with heavy weight towards zero
  // but that code fails currently for testing a uniform
  MYREAL old;
  long newpos;
  long seqpn;
  long end;
  long before;
  MYREAL errorrates[4];
  MYREAL mynew = random_beta(10.0,1.0);
  if(is_combined)
    {
      seqpn = 1;
      newpos = 0;
    }
  else
    {
      seqpn = 4;
      newpos = random_integer(0,3);
    }
  end    = world->seqerrorratesnum[locus]*seqpn;
  before = end - seqpn; 
  old = world->seqerrorrates[locus][before+newpos];
  if(is_combined)
    { 
      errorrates[0] = mynew;
      errorrates[1] = errorrates[0];
      errorrates[2] = errorrates[0];
      errorrates[3] = errorrates[0];
    }
  else
    { 
      errorrates[0] = world->seqerrorrates[locus][before];
      errorrates[1] = world->seqerrorrates[locus][before+1]; 
      errorrates[2] = world->seqerrorrates[locus][before+2];
      errorrates[3] = world->seqerrorrates[locus][before+3];
      errorrates[newpos] = mynew;
    }


  for (i=0;i<sumtips;i++)
    {
      change_freq_tip(world, tips[i], errorrates);
    }

  set_all_dirty(world->root->next, crawlback (world->root->next), world, world->locus);
  first_smooth(world,world->locus);
  count = tree_update (world, world->G, NOASSIGNING);
  world->seqerrorcount[world->locus] += count;
  if(count == 0)
    {
      // reset because changes were not accepted
      errorrates[newpos] = old;
      for (i=0;i<sumtips;i++)
	{
	  change_freq_tip(world, tips[i], errorrates);
	}
      set_all_dirty(world->root->next, crawlback (world->root->next), world, world->locus);
      first_smooth(world,world->locus);
      tree_update (world, world->G, NOASSIGNING);
    }    
  if(!world->in_burnin)
    world->seqerrorsteps[locus] += 1;
  if (!world->in_burnin && world->seqerrorsteps[locus] % world->increment == 0)
    {
      if (world->cold)
	{
#ifdef DEBUG
	  if(is_combined)
	    printf("%i> change_freq() %li %f %li %f\n",myID, world->locus, world->likelihood[world->G], count, errorrates[0]);
	  else
	    printf("%i> change_freq() %li %f %li %f %f %f %f\n",myID, world->locus, world->likelihood[world->G], count, errorrates[0], errorrates[1], errorrates[2], errorrates[3]);
#endif       
	  if (world->seqerrorratesnum[locus] + 100 > world->seqerrorallocnum[locus])
	    {
	      world->seqerrorallocnum[locus] += 100;
	      world->seqerrorrates[locus]= (MYREAL *) myrealloc(world->seqerrorrates[locus],sizeof(MYREAL) * (size_t) (mult*world->seqerrorallocnum[locus]));
	    }
	  end    = world->seqerrorratesnum[locus]*seqpn;
	  world->seqerrorrates[locus][end]   = errorrates[0]; 
	  if(!is_combined)
	    {
	      world->seqerrorrates[locus][end+1] = errorrates[1]; 
	      world->seqerrorrates[locus][end+2] = errorrates[2];
	      world->seqerrorrates[locus][end+3] = errorrates[3];
	    }
	  world->seqerrorratesnum[locus] += 1;
	}
      else
	{
	  world->seqerrorrates[locus][0] = errorrates[0]; 
	  if(!is_combined)
	    {
	      world->seqerrorrates[locus][1] = errorrates[1]; 
	      world->seqerrorrates[locus][2] = errorrates[2];
	      world->seqerrorrates[locus][3] = errorrates[3];
	    }
	}
    }
  else
    {
      // overwrite old record (aka do not record change)
      end    = world->seqerrorratesnum[locus]*seqpn;
      before = end - seqpn; 
      world->seqerrorrates[locus][before]   =  errorrates[0]; 
      if(!is_combined)
	{
	  world->seqerrorrates[locus][before+1] = errorrates[1]; 
	  world->seqerrorrates[locus][before+2] = errorrates[2];
	  world->seqerrorrates[locus][before+3] = errorrates[3];
	}
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
