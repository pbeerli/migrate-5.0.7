/* \file world.c */
/*------------------------------------------------------
Maximum likelihood estimation
of migration rate  and effectice population size
using a Metropolis-Hastings Monte Carlo algorithm
-------------------------------------------------------
    W O R L D   R O U T I N E S

    creates tree structures

    prints results,
    and finally helps to destroy itself.

    Peter Beerli, started 1996, Seattle
    beerli@fsu.edu

    Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
    Copyright 2003-2007 Peter Beerli, Tallahassee FL

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


$Id: world.c 2169 2013-08-24 19:02:04Z beerli $
-------------------------------------------------------*/
#include "migration.h"
#include "sighandler.h"
#include "assignment.h"
#include "menu.h"
#include "mcmc.h"
#include "random.h"
#include "reporter.h"
#include "tools.h"
#include "growth.h"
#include "bayes.h"
#include "laguerre.h"
#include "options.h"
#include "priors.h"
#include "tree.h"
#include "data.h"
#include "migrate_mpi.h"
#include "mittag_leffler.h"
#include "sequence.h"
#include "seqerror.h"
#include "migevents.h"
#include "speciate.h"
#include "skyline.h"
#include "skyparam.h"
#include "haplotype.h"

//#ifdef PRETTY
#include "pretty.h"
//#endif
#ifdef UEP
#include "uep.h"
#endif

#include "mutationmodel.h"
#ifdef BEAGLE
#include "calculator.h"
#endif

#ifdef ALTIVEC
#include "altivec.h"
#endif

#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif

#ifndef EARTH
#define EARTH universe[0]
#endif

//#define PLANESIZE 36
#define PLANEBIGTICKS 6
//#define PLANEBIGTICKVALUES {1,2,3,4,5,6}
//#define CONTOURLEVELS 8
//#define CONTOURS_LOCUS {0.0,-3.35670/2., -9.48773/2., -13.2767/2.,	\
//    0.0,-3.35670/2., -9.48773/2., -13.2767/2.}
//#define CONTOURS_LOCI  {0.0,-4.35146/2., -11.0705/2., -15.0863/2.,	\
//    0.0,-4.35146/2., -11.0705/2., -15.0863/2.}

extern int myID;
static long * mdimfilecount;

/* prototypes ------------------------------------------- */
void create_world (world_fmt ** world, long loci);
void init_world (world_fmt * world, data_fmt * data, option_fmt * options);
//void calc_simple_param (world_fmt * world, data_fmt * data);
void print_menu_locus (FILE *file, world_fmt * world, long locus);
void print_menu_chain (char type, long chain, long steps, world_fmt * world,
                       option_fmt * options, long rep);
void set_bounds (long *increment, long *steps, long *chain,
                 const option_fmt * options, const char type);
//void burnin_chain (world_fmt * world);
void print_list (world_fmt ** universe, option_fmt * options, data_fmt * data);
//void plot_surface (world_fmt * world, option_fmt * options,
//                   data_fmt * data, char ****plane, long x);
//void plot_surface_header2 (FILE * outfile, long locus, char plotvar[]);
//void plot_surface2 (FILE * outfile, long x, long locus, long numpop,
//                    char ****plane, plotmax_fmt ** plotmax,
//                    char **popnames, char plotvar[], long migrmodel);
void free_universe (world_fmt ** worlds, long numworlds, option_fmt *options);
void free_world(world_fmt *world, option_fmt * options);
void precalc_world (world_fmt * world);
void reprecalc_world (world_fmt * world, long that);
extern MYREAL prob_tree (world_fmt * world, timelist_fmt * tyme);
/* private functions */
void create_timearchive (timearchive_fmt *** atl, long loci,
                         long samples, long numpop, long replicates);
void init_plotplane (world_fmt * world);
void create_cov (world_fmt * world);
void print_menu_equilib (world_fmt * world);
void print_finish (world_fmt * world, long filepos);
void copy_time (world_fmt * world, timelist_fmt * ltl, long from,
                long to, long np, long rep);
void archive_timelist (tarchive_fmt * atl, vtlist * tl, long T,
                       long np, world_fmt * world);
void copy_timelist (tarchive_fmt * from, tarchive_fmt * to, long np);
long whichp (long from, long to, long pop);
/* long whichl (long from, long to, long pop); */
/* long calc_T(timearchive_fmt * old, long steps); */
void print_results (world_fmt ** universe, option_fmt * options, data_fmt * data);
void print_fst (world_fmt * world, option_fmt * options,
                data_fmt * data, MYREAL **fstparam);
void prepare_print_nu (MYREAL nu, char *str);
void prepare_print_nm (MYREAL nm, MYREAL nmu, char *strllike);
void print_menu_coalnodes (FILE * file, world_fmt * world, long G, long rep);
void print_menu_createplot (boolean progress);
void print_mathematica (world_fmt * world, MYREAL ***plane, long x, long y);
void print_cov (world_fmt * world, long numpop, long loci, MYREAL ***cov);
void print_cov_table (FILE * outfile, long locus, world_fmt * world,
                      MYREAL ***cov, MYREAL *corr, long addvar);
void
print_cov_table2 (FILE * outfile, long locus, world_fmt * world,
		  MYREAL ***cov, MYREAL *corr, long addvar);
void free_seqx (node * p, world_fmt * world);
void free_x (node * p, world_fmt * world);
void increase_timearchive (world_fmt * world, long locus,
                           long sample, long numpop, long rep);
void print_CV (world_fmt * world);
void print_gelmanr (MYREAL average, MYREAL biggest);
void print_progress(worldoption_fmt * options, world_fmt * world,
                    long rep, long visited, long accepted);
void print_param (char ** file, long *bufsize, long *allocbufsize, boolean usem, world_fmt *world, long nn,
                  char spacer[]);
void set_contours (MYREAL **contours, long df, long numpop);
void print_result_header (char *titletext, world_fmt * world);
void print_result_population (long pop, world_fmt * world,
                              option_fmt * options, data_fmt * data);
void print_result_param (FILE * file, MYREAL *param, long numpop, long pop,
                         boolean usem);
void print_result_fst (long pop, world_fmt * world, data_fmt * data);
void prognose_time (char *nowstr, world_fmt * world,
                    long options_increment, long steps, char *spacer, boolean tobuffer);
void klone (world_fmt * original, world_fmt * kopie,
            option_fmt * options, data_fmt * data, MYREAL temperature);
void klone_part (world_fmt * original, world_fmt * kopie,
                 option_fmt * options, data_fmt * data, MYREAL temperature);
void clone_polish (world_fmt * original, world_fmt * kopie);
long chance_swap_tree (world_fmt * tthis, world_fmt * that);
void advance_clone_like (world_fmt * world, long accepted, long *j);
void polish_world (world_fmt * world);
void fill_worldoptions (worldoption_fmt * wopt, option_fmt * options, long numpop);
void fill_worlddata (worlddata_fmt * wdata, world_fmt * world,
                     data_fmt * data, long numpop, boolean readsum);
void print_replicate(world_fmt *world, long maxrep, long rep, long locus);
#ifdef LONGSUM
void print_fluctuate_header(world_fmt *world);
void setup_fluctuate(world_fmt *world, option_fmt * options);
void print_popstring(long pop, world_fmt *world, option_fmt *options, data_fmt *data);
void print_fluctuate(world_fmt **universe, option_fmt *options, data_fmt *data);
#endif /*LONGSUM*/


boolean auto_updating_bayes(world_fmt *world, boolean autostop, boolean essstop,
		     long *z, long delta,MYREAL *var,MYREAL *oldvar,MYREAL *oldvar2,
		     MYREAL *autocorrelation,MYREAL *effective_sample);

boolean auto_updating_ML(world_fmt *world, boolean autostop, boolean essstop,
		     long *z, long delta,MYREAL *var,MYREAL *oldvar,MYREAL *oldvar2,
		     MYREAL *autocorrelation,MYREAL *effective_sample);

boolean updating(world_fmt *world);

boolean autostopping(world_fmt *world, boolean autostop, boolean essstop,
		     long *z,long step, long delta,MYREAL *var,MYREAL *oldvar,MYREAL *oldvar2,
		     MYREAL *autocorrelation,MYREAL *effective_sample);
void set_ticks (MYREAL **ticks, MYREAL *plotrange, short type);
boolean is_same(long i, long j);
void print_popstring(long pop, world_fmt *world, option_fmt *options, data_fmt *data);
void klone_tree_setup(world_fmt *kopie, option_fmt *options);
void set_replicates (world_fmt * world, long repkind, long rep, long *repstart,
		     long *repstop);
long get_numparam(world_fmt* world);

//##
extern void   unset_penalizer_function(boolean inprofiles);

extern MYREAL (*log_prior_theta) (world_fmt *, long);//for heating
extern MYREAL (*log_prior_mig) (world_fmt *, long);//for heating
extern MYREAL (*log_prior_rate) (world_fmt *, long);//for heating

extern double page_height;

/* Functions ++++++++++++++++++++++++++++++++++++++++++++++++ */


/// \brief allocation of memory for the world structure
///
/// Allocation of the structure world, this structure is permanent throught
/// the program and is used almost in all functions of migrate, it allocates 
/// itself, the times of events on the working genealogy.
void
create_world (world_fmt ** world, long loci)
{
  (void) loci;
    (*world) = (world_fmt *) mycalloc (1, sizeof (world_fmt));
    (*world)->treetimes = (timelist_fmt *) mycalloc (1, sizeof (timelist_fmt));
    (*world)->param_like = -MYREAL_MAX;
}

/// \brief assign parameters from options to world->options
/// 
/// copy the options into the world structure, this function also allocates
/// also the memory for the arrays copied
/// [removed the custm copy operation because init_world is doing this already -- 
/// speciation related changes in wolrd that are not shown in options, thus options and world->options
/// are different in custm and custm2
void
fill_worldoptions (worldoption_fmt * wopt, option_fmt * options, long numpop)
{
    long i;
    long optnumpop2;
    long numpop2 = numpop * numpop;
   // MYREAL total;
    wopt->smoothing[0]=options->smoothing[0];
    wopt->smoothing[1]=options->smoothing[1];
    wopt->smoothing_window[0]=options->smoothing_window[0];
    wopt->smoothing_window[1]=options->smoothing_window[1];
    wopt->allposteriors = options->allposteriors;
    wopt->bayes_infer = options->bayes_infer;
    wopt->gamma = options->gamma;
    wopt->alphavalue = options->alphavalue;
    wopt->murates = options->murates;
    wopt->murates_fromdata = options->murates_fromdata;
    // variable mutation rates are not necessarily filled in 
    // if they are empty than a default value (1.0) is filled in
    if(options->mu_rates == NULL)
      options->mu_rates = (MYREAL *) mycalloc (options->muloci, sizeof (MYREAL));
    // the inheritance scalars are typically not filled in
    // and are set to 1.0 in this function
    if(options->inheritance_scalars == NULL)
      {
	options->inheritance_scalars = (MYREAL *) mycalloc (options->muloci, sizeof (MYREAL));
      }
    else
      {
	options->inheritance_scalars = (MYREAL *) myrealloc (options->inheritance_scalars, options->muloci * sizeof (MYREAL));
      }
    if(options->inheritance_scalars_numalloc < options->muloci)
      {
	double last  = options->inheritance_scalars[options->inheritance_scalars_numalloc-1];
	for (i = options->inheritance_scalars_numalloc; i < options->muloci; i++)
	  {
	    options->inheritance_scalars[i] = last;
	  }
      }
    // allocation of inheritance scalars for world->options.
    if(wopt->inheritance_scalars == NULL)
      wopt->inheritance_scalars = (MYREAL *) mycalloc (options->muloci, sizeof (MYREAL));
    //
    // fill in growth using growpops in options
    init_growpop(wopt,options, numpop);
    //
    // allow to run the prior only
    wopt->prioralone = options->prioralone;
    //
    // world->option  uses als mu_rates and need to allocate them
    if(wopt->mu_rates == NULL)
      {
	wopt->mu_rates = (MYREAL *) mycalloc ((2 * options->muloci), sizeof (MYREAL));
	//	printf("%i> murate size %li\n",myID, 2 * options->muloci * sizeof (MYREAL));
	wopt->lmu_rates = wopt->mu_rates + options->muloci;
      }

    //    if (options->murates)
    // {
    //    memcpy (wopt->mu_rates, options->mu_rates,
    //            sizeof (MYREAL) * options->muloci);
    //    for (i = 0; i < options->muloci; i++)
    //    {
    //        wopt->lmu_rates[i] = LOG (options->mu_rates[i]);
    //     }
    //}
    //else
    //{
    for (i = 0; i < options->muloci; i++)
      {
	if(!options->bayesmurates)
	  {
	    if(options->mu_rates[i] <= 0.0)
	      {
		wopt->mu_rates[i] = 1.0;
		wopt->lmu_rates[i] = 0.0;
	      }
	    else
	      {
		wopt->mu_rates[i] = options->mu_rates[i];
		wopt->lmu_rates[i] = log(options->mu_rates[i]);
	      }
	  }
	wopt->inheritance_scalars[i] = options->inheritance_scalars[i];
      }
    wopt->migration_model = options->migration_model;
    if(wopt->custm==NULL)
      {
	wopt->custm = (char *) mycalloc ((2 * (numpop2 + 2)), sizeof (char));
	wopt->custm2 = wopt->custm + (numpop2 + 2);
	//1229 }
	optnumpop2 = (long) strlen(options->custm);
	if(optnumpop2>numpop2)
	  optnumpop2=numpop2;
	sprintf (wopt->custm, "%-*.*s", (int) optnumpop2, (int) optnumpop2, options->custm);
	sprintf (wopt->custm2, "%-*.*s", (int) optnumpop2, (int) optnumpop2, options->custm2);
      }//1229
    if(wopt->thetag==NULL)
      wopt->thetag = (MYREAL *) mycalloc (numpop2 + 2, sizeof (MYREAL));
    memcpy (wopt->thetag, options->thetag,
            sizeof (MYREAL) * (size_t) options->numthetag);
    if (wopt->mg==NULL)
      wopt->mg = (MYREAL *) mycalloc (numpop2 + 2, sizeof (MYREAL));
    memcpy (wopt->mg, options->mg, sizeof (MYREAL) * (size_t) options->nummg);
    wopt->zeron = 0; //=options->zeron;
    if(wopt->zeroparam==NULL)
      wopt->zeroparam = (long *) mycalloc ((1 + wopt->zeron), sizeof (long));
    //memcpy (wopt->zeroparam, options->zeroparam, sizeof (long) * wopt->zeron);
    wopt->constn = 0; //=options->constn;
    if(wopt->constparam==NULL)
      wopt->constparam = (long *) mycalloc ((1 + wopt->constn), sizeof (long));
    //    memcpy (wopt->constparam, options->constparam,
    //        sizeof (long) * wopt->constn);
    wopt->symn = 0; //=options->symn;
    if (wopt->symparam==NULL)
      wopt->symparam = (twin_fmt *) mycalloc ((1 + wopt->symn), sizeof (twin_fmt));
    // memcpy (wopt->symparam, options->symparam, sizeof (twin_fmt) * wopt->symn);
    wopt->sym2n = 0; //=options->sym2n;
    if (wopt->sym2param==NULL)
      wopt->sym2param = (quad_fmt *) mycalloc ((1 + wopt->sym2n), sizeof (quad_fmt));
    //memcpy (wopt->sym2param, options->sym2param,
    //        sizeof (quad_fmt) * wopt->sym2n);
    wopt->mmn = 0; //=options->mmn;
    wopt->tmn = 0; //options->tmn;
    if (wopt->mmparam==NULL)
      wopt->mmparam = (long *) mycalloc ((1 + wopt->mmn), sizeof (long));
    //memcpy (wopt->mmparam, options->mmparam, sizeof (long) * wopt->mmn);
    wopt->tersepdf = options->tersepdf;
    wopt->mixplot = options->mixplot;
    wopt->mixfile = options->mixfile;
    wopt->progress = options->progress;
    wopt->writelog = options->writelog;
    wopt->logfile = options->logfile;
    wopt->plotnow = options->plotnow;
    wopt->verbose = options->verbose;
    wopt->replicate = options->replicate;
    wopt->gelman = options->gelman;
    wopt->gelmanpairs = options->gelmanpairs;
    wopt->lcepsilon = options->lcepsilon;
    wopt->simulation = options->simulation;
    wopt->datatype = options->datatype;
    if(wopt->bayes_infer)
      wopt->lsteps = options->lsteps;
    else
      wopt->lsteps = MAX(options->ssteps,options->lsteps);
    wopt->lincr = options->lincrement;
    wopt->loglsteps = LOG ((double) wopt->lsteps);
    wopt->treeprint = options->treeprint;
    wopt->movingsteps = options->movingsteps;
    wopt->acceptfreq = options->acceptfreq;
    wopt->rcategs = options->rcategs;
    wopt->categs = options->categs;
    wopt->heating = options->heating;
    wopt->heated_chains = options->heated_chains;
    wopt->heating_interval = options->heating_interval;
    wopt->heating_count = 0;
    wopt->adaptiveheat = options->adaptiveheat;
    if (wopt->heat==NULL)
      wopt->heat = (MYREAL *) mycalloc (wopt->heated_chains, sizeof (MYREAL));
    if (wopt->averageheat==NULL)
      wopt->averageheat = (MYREAL *) mycalloc (wopt->heated_chains, sizeof (MYREAL));
    memcpy (wopt->heat, options->heat, sizeof (MYREAL) * (size_t) wopt->heated_chains);
    wopt->profilemethod = options->profilemethod;
    wopt->printprofile = options->printprofile;
    wopt->printprofsummary = options->printprofsummary;
    wopt->profileparamtype = options->profileparamtype;
    wopt->df = options->df;
    wopt->lchains = options->lchains;
    wopt->replicatenum = options->replicatenum;
    for(i=0;i<options->muloci;i++)
      wopt->micro_threshold[i] = options->micro_threshold;
    wopt->micro_stepnum = options->micro_stepnum;
    wopt->msat_tuning[0] = options->msat_tuning[0];
    wopt->msat_tuning[1] = options->msat_tuning[1];
    if (wopt->rrate==NULL)
      wopt->rrate = (MYREAL *) mycalloc (wopt->rcategs, sizeof (MYREAL));
    memcpy (wopt->rrate, options->rrate, sizeof (MYREAL) * (size_t) wopt->rcategs);
    if (wopt->rate==NULL)
      wopt->rate = (MYREAL *) mycalloc (wopt->categs, sizeof (MYREAL));
    memcpy (wopt->rate, options->rate, sizeof (MYREAL) * (size_t) wopt->categs);
    if (wopt->probcat==NULL)
      wopt->probcat = (MYREAL *) mycalloc (wopt->rcategs, sizeof (MYREAL));
    memcpy (wopt->probcat, options->probcat, sizeof (MYREAL) * (size_t) wopt->rcategs);
    wopt->fastlike = options->fastlike;
    wopt->pluschain = options->pluschain;
    wopt->mighist = options->mighist;
    wopt->mighist_all = options->mighist_all;
    wopt->mighist_counter = options->mighist_counter;
    wopt->mighist_increment = options->mighist_increment;
    wopt->skyline = options->skyline;
    wopt->skyline_param = options->skyline_param;
    wopt->timeelements = options->timeelements;
    wopt->eventbinsize = options->eventbinsize;
    wopt->burn_in = options->burn_in;
    wopt->burnin_autostop = options->burnin_autostop;
    wopt->usem = options->usem;
    wopt->minmigsumstat = options->minmigsumstat;
    wopt->lambda = options->lambda;
    set_updating_choices(wopt->choices, options, STANDARD);
    memcpy(wopt->slice_sampling,options->slice_sampling,sizeof(boolean) * PRIOR_SIZE);
    wopt->has_bayesfile = options->has_bayesfile;
    wopt->has_bayesmdimfile = options->has_bayesmdimfile;
    wopt->bayesmdiminterval = options->bayesmdiminterval;
    wopt->has_datefile = options->has_datefile;
    // we use the inverse of the number of generations per year
    wopt->generation_year = options->generation_year;//inverse of generations per year
    // this will be used to standardize the interal timescale (expected mutations) with
    // an external timescale of "years" 
    wopt->mutationrate_year_numalloc = options->mutationrate_year_numalloc;
    if(wopt->mutationrate_year==NULL)
      wopt->mutationrate_year = (MYREAL *) mycalloc(options->mutationrate_year_numalloc,sizeof(MYREAL));
    //total = 0;
    for(i = 0; i < options->mutationrate_year_numalloc; i++)
      {
	wopt->mutationrate_year[i] = options->mutationrate_year[i];
      }
    set_meanmu(wopt,options,options->muloci);
#ifdef UEP
    wopt->uep = options->uep;
    wopt->ueprate = options->ueprate;
    wopt->uepmu = options->uepmu;
    wopt->uepnu = options->uepnu;
    wopt->uepfreq0 = options->uepfreq0;
    wopt->uepfreq1 = options->uepfreq1;
#endif
    wopt->heatedswap_off = options->heatedswap_off;
    if (options->bayes_infer)
      {
	wopt->has_autotune = options->has_autotune;
	wopt->autotune = options->autotune;
      }
    else
      {
	wopt->has_autotune = FALSE;
	wopt->autotune = ' ';
      }
  wopt->onlyvariable = options->onlyvariable; 	    
  wopt->has_variableandone = options->has_variableandone; 
  wopt->firstinvariant = options->firstinvariant;     
}


/// \brief copies data related variables into world
///
/// copies data related variables into a structure wolrd->data
void
fill_worlddata (worlddata_fmt * wdata, world_fmt * world,
                data_fmt * data, long numpop, boolean readsum)
{
    long numpop2 = numpop * numpop;
    wdata->skiploci =
        (boolean *) mycalloc (1, sizeof (boolean) * (size_t) (data->loci + 1 + data->loci));
    wdata->geo = (MYREAL *) mycalloc (1, sizeof (MYREAL) * (size_t) numpop2);
    memcpy (wdata->geo, data->geo, sizeof (MYREAL) * (size_t) numpop2);
    wdata->lgeo = (MYREAL *) mycalloc (1, sizeof (MYREAL) * (size_t) numpop2);
    memcpy (wdata->lgeo, data->lgeo, sizeof (MYREAL) * (size_t) numpop2);
    if (!readsum)
    {
      wdata->maxalleles = (long *) mycalloc ((2 * data->loci), sizeof (long));
        memcpy (wdata->maxalleles, data->maxalleles,
                sizeof (long) * (size_t) data->loci);
        wdata->seq[0]->sites = wdata->maxalleles + data->loci;
	if(data->seq[0]==NULL)
	  {
	    long sublocus;
	    for(sublocus=0; sublocus < data->allsubloci; sublocus++)
	      {
		wdata->seq[0]->sites[sublocus] = world->mutationmodels[sublocus].numsites;
	      }
	    // perhaps need to deal with mutiple addon? wdata->seq[0]->addon = data->seq[0]->addon;
	  }
	else
	  {
	    memcpy (wdata->seq[0]->sites, data->seq[0]->sites,
		   sizeof (long) * (size_t) data->loci);
	    wdata->seq[0]->addon = data->seq[0]->addon;
	  }
        wdata->sumfile = data->sumfile;
        wdata->sampledates = data->sampledates;
	wdata->maxsampledate = data->maxsampledate;
	wdata->numind = data->numind;
#ifdef UEP        
        wdata->uepsites = data->uepsites;
#endif
    }
    wdata->locusweight = (MYREAL *) mycalloc(data->loci+1, sizeof(MYREAL));
    memcpy(wdata->locusweight,data->locusweight,sizeof(MYREAL)*(size_t) data->loci);//invariant loci treatment
#ifdef DEBUG
    //    long locus;
    //for (locus=0; locus<data->loci;locus++)
    //  fprintf(stdout,"%i> invariant loci treatment: locus=%li weight=%f\n",myID, locus, wdata->locusweight[locus]);
#endif
}

/// initialization of world, Step II allocation, filling values
void
init_world (world_fmt * world, data_fmt * data, option_fmt * options)
{
#ifdef MPI
  long maxnum = 0;
#endif
  long numpop;
  //long atemp;
  long locus, i, rep, step;
  long size;
  long custmlen = 0;
  long addition;
  long sumloc = data->loci > 1 ? 1 : 0;
  long convergence_len;
  long convergence_len2;
  long maxreplicate;
  long numchains=1;
  long npall;
  //mittag-leffler
  world->mlalpha = options->mlalpha;
  world->mlinheritance = options->mlinheritance;
  world->options = (worldoption_fmt *) mycalloc (1, sizeof (worldoption_fmt));
  world->buffer = (char *) mycalloc (LINESIZE, sizeof (char));
  if(options->murates_fromdata)
    options->muloci = data->loci;
  world->root = NULL;
  world->contribution = NULL;
  world->options->treeinmemory=FALSE; // if set TRUE (see later) this will record best tree to a string
#ifdef MPI
  if(myID==MASTER && world->cold)
    {
      maxnum = (options->replicate == TRUE) ? ((options->replicatenum > 0) ? options->replicatenum  : options->lchains ) : 1;
      maxnum = data->loci * maxnum;
      //      printf("%i> MASTER: mpi request_stack allocation is %li deep (numcpu=%i)\n",myID, maxnum, numcpu);
      if(maxnum < numcpu)
        maxnum = numcpu;
      world->mpistacknum = 0;
      world->mpistack_requestnum = 0;
      world->mpistack_numalloc = maxnum;
      world->mpistack = (int *) mycalloc(maxnum,sizeof(int));
      world->mpistack_request_numalloc = maxnum;
      world->mpistack_request = (mpirequest_fmt *) mycalloc(maxnum,sizeof(mpirequest_fmt));
      for(i=0;i<maxnum;i++)
	world->mpistack_request[i].tempstr = (char *) mycalloc(SMALLBUFSIZE,sizeof(char));
      //printf("%i=================> stack allocation is maxnum=%li\n",myID,maxnum);
    }

  world->who = (int *) mycalloc (data->loci, sizeof (int));
#ifdef SLOWNET
  world->profilewho = (int *) mycalloc (options->newpops_numpop * options->newpops_numpop + 2, sizeof (int));
  //OLD   world->profilewho = (int *) mycalloc (data->numpop * data->numpop + 2, sizeof (int));
#endif
#endif
      
      world->data = (worlddata_fmt *) mycalloc (1, sizeof (worlddata_fmt));
      world->data->seq = (seqmodel_fmt **) mycalloc (1, sizeof (seqmodel_fmt *));
      world->data->seq[0] = (seqmodel_fmt *) mycalloc (1, sizeof (seqmodel_fmt));
      world->data->seq[0]->done = FALSE;
      options->muloci = data->loci;
      world->options->micro_threshold = (long *) mycalloc(data->loci,sizeof(long));
      numpop = options->newpops_numpop;
      if (numpop==0)
	{
	  numpop = data->numpop;
	  options->newpops_numpop = numpop;
	}
      world->data->allsubloci = data->allsubloci;
      fill_worldoptions (world->options, options, numpop);
      world->has_unassigned = options->has_unassigned;
      world->has_unassignedfreq = options->has_unassignedfreq;
      world->has_unassignedpoplist = options->has_unassignedpoplist;
      if (options->has_unassignedpoplist)
	{
	  world->unassignedpoplistnum = options->unassignedpoplistnum;
	  world->unassignedpoplist = (long *) mycalloc(options->unassignedpoplistnum,sizeof(long));
	  memcpy(world->unassignedpoplist,options->unassignedpoplist,sizeof(long)*world->unassignedpoplistnum);
	}
      if (world->has_unassigned)
	fill_world_unassigned(world);
      fill_worlddata (world->data, world, data, numpop, options->readsum);
      if (!options->readsum)
	{
	  world->loci = data->loci;
	  world->skipped = 0;
	  world->numpop = numpop;
	  world->numpop2 = numpop * numpop;
	  if (world->numpop == 1)
	    options->plot = FALSE;
	  if(world->cold)
	    {
	      world->options->treeinc = options->treeinc;
	      switch(world->options->treeprint)
		{
		case ALL:
		case LASTCHAIN:
		case BEST:
		  world->options->treeinmemory=TRUE;
		  world->treespace = (char**) mycalloc(world->loci, sizeof(char *));
		  world->treespacenum = (long *) mycalloc(world->loci,sizeof(long));
		  world->treespacealloc = (long *) mycalloc(world->loci,sizeof(long));
		  world->besttreelike = (MYREAL *) mycalloc(world->loci,sizeof(MYREAL));
		  for(i=0;i<world->loci;i++)
		    {
		      world->treespace[i] = (char *) mycalloc(LONGLINESIZE, sizeof(char));
		      world->treespacealloc[i] = LONGLINESIZE;
		      world->besttreelike[i] = (double) -HUGE;
		    }
		  break;
		default:
		  world->options->treeinmemory = FALSE;
		  world->options->treeprint = myNONE;
		  break;
		}
	    }
	  else
	    {
	      world->options->treeinmemory = FALSE;
	      world->options->treeprint = myNONE;
	    }
	}
      custmlen = (long) strlen (options->custm);
      fillup_custm (custmlen, world, options);

      world->has_estimateseqerror = options->has_estimateseqerror;
      if (world->has_estimateseqerror)
	{
	  fill_world_seqerror(world,options);
	}
      world->bayesaccept = -1L; // default setting for all modes to do only tree changes and tree heat swaps.

      // bayes init is deferred after relabeling of populations

      if(options->bayes_infer)
	{
	  world->bayes = (bayes_fmt *) mycalloc(1, sizeof(bayes_fmt));
	  world->bayes->count=0;
	  world->has_speciation = options->has_speciation;
	  world->has_migration = options->has_migration;
	  init_speciesvector(world, options);
	  // growth
	  init_growth(world,numpop);
	  bayes_init(world->bayes,world, options);
	  if(world->cold)
	    {
	      check_bayes_priors(options, data, world);
	      bayes_init_histogram(world, options);
	      // speciation
	      if (options->recorddivtime)
		{
		  world->divtime_alloc = HUNDRED;
		  world->divtime_num = 0;
		  world->divtime = (divtime_fmt *) mycalloc(world->divtime_alloc,sizeof(divtime_fmt));
		}
	    }
	  fill_speciesvector(world, options);
	  bayes_fill(world, options);
	}
      //else
      //	{
      //  // we use the world->bayes->map for the skylineplots
	//  world->bayes = (bayes_fmt *) mycalloc(1, sizeof(bayes_fmt));
      //world->bayes->count=0;
      //  world->bayes->hypercount=0;
      //  world->bayes->mu = FALSE;
      //  world->bayes->map = (longpair *) mycalloc(world->numpop2, sizeof(longpair));//ML
      //  setup_bayes_map(world->bayes->map, world, world->numpop2);//ML
      //}
      if(world->options->bayes_infer)
	{
	  npall =world->numpop2+world->bayes->mu+world->species_model_size*2 + world->grownum;
	  if (options->slice_sticksizes == NULL)
	    options->slice_sticksizes = (MYREAL *) mycalloc(npall, sizeof (MYREAL));
	  else
	    options->slice_sticksizes = (MYREAL *) myrealloc (options->slice_sticksizes,npall * sizeof (MYREAL));
	  if (world->options->slice_sticksizes == NULL)
	    world->options->slice_sticksizes = (MYREAL *) mycalloc (npall, sizeof (MYREAL));
	  else
	    world->options->slice_sticksizes = (MYREAL *) myrealloc (world->options->slice_sticksizes,npall * sizeof (MYREAL)); 
	  for(i=0;i<npall;i++)
	    {
	      options->slice_sticksizes[i] = (world->bayes->maxparam[i] - world->bayes->minparam[i])/10.;
	    }
	  memcpy(world->options->slice_sticksizes,options->slice_sticksizes,sizeof(MYREAL) * (size_t) (npall));
	}
      world->mig0list = (MYREAL *) mycalloc (world->numpop, sizeof (MYREAL));
      //for precalc_world
      world->migproblist = (MYREAL **) mycalloc ((2 * world->numpop), sizeof (MYREAL *));
      //for precalc_world
      for (i = 0; i < 2 * world->numpop; ++i)
	{
	  world->migproblist[i] =
	    (MYREAL *) mycalloc (world->numpop, sizeof (MYREAL));
	}
      world->speciesproblist = world->migproblist + world->numpop;
      world->design0list = (long *) mycalloc (world->numpop, sizeof (long));
      //for precalc_world
      rep = world->repstop =
        world->options->replicate ? (world->options->replicatenum >
                                     0 ? world->options->
                                     replicatenum : world->options->lchains) : 1;
      // ML:                 create time archive
      // Bayesian inference: Standard method, don't create time archive of trees
      //                     integrated likelihood, use time archive to calculate profiles
      if(options->bayes_infer)
	{
#ifdef INTEGRATEDLIKE
	  if(options->integrated_like)
	    {
	      create_timearchive (&(world->atl), world->loci, SAMPLETREE_GUESS,
				  world->numpop, rep);
	      if (world->numpop > 1)
		init_plotplane (world);
	      create_cov (world);
	    }
#endif
	}
      else
	{
	  error("Only Bayesian inference is allowed");
	}
      addition = world->species_model_size * 2 + world->bayes->mu + world->grownum;//#speciation*mu-rate+#growthmax
#ifdef LONGSUM
      
      addition = 1 + world->numpop * 3;
#endif /*LONGSUM*/
      world->alloclike = world->options->lsteps;
      world->likelihood = (MYREAL *) mycalloc (world->alloclike , sizeof (MYREAL)); //debug
      world->lineages = (long *) mycalloc (world->numpop, sizeof (long));
      world->param0 = (MYREAL *) mycalloc ((world->numpop2 + addition), sizeof (MYREAL));
      world->param00 = (MYREAL *) mycalloc ((world->numpop2 +addition), sizeof (MYREAL));
      // skyline parameters
      world->timeelements = options->timeelements;
      world->times = (MYREAL *) mycalloc ((world->timeelements +  world->timeelements * (world->numpop2)), sizeof (MYREAL));
      world->timek = world->times + world->timeelements;
      world->times[0] = 0.0;
      for(i=1;i<world->timeelements-1;i++)
	{
	  //TODO time range needs to be fixed
	  world->times[i] += world->times[i-1] + (double) options->eventbinsize;
	}
      world->times[i]=MYREAL_MAX;
      for(i=0;i<world->timeelements*world->numpop2;i++)
	{
	  world->timek[i]=1.0;
	}
      //-----------------------
      world->heat = 1.;
      world->heatid=0;
      world->proposal = (proposal_fmt *) calloc(1, sizeof(proposal_fmt));
      world->has_proposal_details=FALSE;
      world->has_proposal_first = TRUE;
      world->varheat = CHAINVARIANCEDELTA;
      world->essminimum = ESSMINIMUM;
      world->logprior = (double) -HUGE;
      world->mstart = (int *) mycalloc (world->numpop, sizeof (int));
      world->mend = (int *) mycalloc (world->numpop, sizeof (int));
      for (i = 0; i < world->numpop; i++)
	{
	  world->mstart[i] = (short) mstart (i, world->numpop);
	  world->mend[i] = (short) mend (i, world->numpop);
	}
      switch (options->datatype)
	{
        case 's':
        case 'n':
        case 'h':
        case 'u':
        case 'a':
        case 'f':
	case '@':
	  world->fstparam =
            (MYREAL **) mycalloc (1, sizeof (MYREAL *) * (size_t) (world->loci + sumloc));
	  for (locus = 0; locus < world->loci + sumloc; locus++)
	    world->fstparam[locus] =
	      (MYREAL *) mycalloc ((world->numpop2 * 2), sizeof (MYREAL));
	  break;
        case 'm':
	  world->options->steps =
	    (MYREAL ***) mycalloc (world->loci, sizeof (MYREAL **));
	  for(locus=0; locus < world->loci; locus++)
	    {
	      world->options->steps[locus]=(MYREAL **) mycalloc(world->options->micro_stepnum,sizeof(MYREAL *));
	      for (step = 0; step < world->options->micro_stepnum; step++)
		{
		  world->options->steps[locus][step] = (MYREAL *) mycalloc (world->options->micro_stepnum, sizeof (MYREAL) );
		}
	    }
	  /* do not break here */
	case 'b':
	  world->fstparam =
	    (MYREAL **) mycalloc ((world->loci + sumloc), sizeof (MYREAL *));
	  for (locus = 0; locus < world->loci + sumloc; locus++)
	    world->fstparam[locus] =
	      (MYREAL *) mycalloc ((world->numpop2 * 2), sizeof (MYREAL));
	  break;
	}
      /*
       * chainlikelihood - for over last chains or replicate combining
       * estimator
       */
      world->chainlikes = (MYREAL **) mycalloc ((world->loci + sumloc), sizeof (MYREAL *));
      world->chainlikes[0] =
        (MYREAL *) mycalloc (((world->loci + sumloc) *
			      (rep + world->options->pluschain)), sizeof (MYREAL));
      for (i = 1; i < world->loci + sumloc; i++)
        world->chainlikes[i] =
	  world->chainlikes[0] + (rep + world->options->pluschain) * i;
      //defined earlier:world->repstop = rep;
      world->rep = 0;
      world->lsteps = 0;

      maxreplicate = (options->replicate
		      && options->replicatenum > 0) ? options->replicatenum : 1;
      world->maxreplicate = maxreplicate;
      /* timer -- how long will the program run? */
      time (&world->starttime);
      /* replicates * loci * steps where
         steps = schains * (ssteps * sincr) + lchains * (lsteps*lincr) + (schains + lchains) * burnin*/

      world->treestotal = world->loci * (options->lincrement *
					 (options->lsteps + options->burn_in));
      if (world->options->replicatenum > 0)
        {
	  world->treestotal *= maxreplicate;
	}
#ifdef UEP
      //setup UEP
      setup_uep (world, options);
#endif
      
#ifdef LONGSUM
      // allows for three rates through time at GIVEN times [currently fixed at 1/3 of
      // the longest tree found [start with dbl_max if no start is given
      world->flucrates = (MYREAL *) mycalloc (world->numpop * 6, sizeof (MYREAL));
      world->lflucrates = (MYREAL *) mycalloc (world->numpop * 6, sizeof (MYREAL));
      setup_fluctuate(world,options);
#endif /*LONGSUM*/
      // only for cold chain
      if(world->cold)
	{
	  //migration histogram
	  setup_mighist (world, options);
	  //expected parameters (skyline plots) based on timelist events histogram
	  //needs setup_mighist()
	  setup_expected_events(world,options);
	  /*convergence material*/
	  numchains = options->lchains + options->schains;
	  world->burnin_stops_alloc = numchains*maxreplicate*world->loci+1;
	  world->burnin_stops = (burnin_record_fmt *) mycalloc(world->burnin_stops_alloc, sizeof(burnin_record_fmt));
	  //
	  if(options->bayes_infer)
	    convergence_len = world->numpop2 + 1;
	  else
	    convergence_len = world->numpop2 + world->numpop * 3;
	  convergence_len2 = maxreplicate;
	  convergence_len *= convergence_len2;
	  convergence_len2 *= convergence_len2;
	  world->convergence = (convergence_fmt *) mycalloc(1, sizeof(convergence_fmt));
	  world->convergence->gelmanmeanmaxR = (MYREAL *) mycalloc(convergence_len2, sizeof(MYREAL));
	  world->convergence->chain_s = (MYREAL *) mycalloc(convergence_len, sizeof(MYREAL));
	  world->convergence->chain_means = (MYREAL *) mycalloc(convergence_len, sizeof(MYREAL));
	  world->convergence->chain_counts = (long *) mycalloc(convergence_len, sizeof(long));
	}
      /*Bayes factor*/
      world->bfscale = (MYREAL *) mycalloc(((4 * options->heated_chains * world->loci + 5 * world->loci)), sizeof(MYREAL));
      world->hmscale = world->bfscale + world->loci;
      world->amscale = world->hmscale + world->loci;
      world->bf = world->amscale + world->loci;
      world->hm = world->bf + options->heated_chains * world->loci; 
      world->am = world->hm + world->loci;
      world->steppingstones = world->am + world->loci;
      world->steppingstone_scalars = world->steppingstones + options->heated_chains * world->loci;
      world->steppingstone_counters = world->steppingstone_scalars + options->heated_chains * world->loci; 
      for(locus=0;locus < world->loci; locus++)
	{
	  world->hmscale[locus] = 0.0;//HUGE
	  world->amscale[locus] = (double) -HUGE;
	  world->bfscale[locus] = (double) -HUGE;
#ifdef MPI
	  if(myID==MASTER)
	    {
	      world->hmscale[locus] = 0.0;//-HUGE
	      world->amscale[locus] = (double) HUGE;
	      world->bfscale[locus] = (double) HUGE;	    
	    }
#endif
	}
      /*Allocation of memory for effective sample size and autocorrelation 
        for all parameters (incl rate and speciation PLUS the tree*/
      size = 1 + world->numpop2+options->bayesmurates + world->species_model_size * 2 + world->grownum;
      world->autocorrelation =  (MYREAL *) mycalloc((2 * size),sizeof(MYREAL));
      world->effective_sample =  world->autocorrelation + size; 
      if(world->cold)
	{
	  world->auto_archive =  (MYREAL *) mycalloc((2 * size),sizeof(MYREAL));
	  world->ess_archive =  world->auto_archive + size; 
	  world->accept_archive = (long*) mycalloc((2 * size),sizeof(long));
	  world->trials_archive = world->accept_archive + size;
	  world->archive_n = 1;
	}      
      /* first attempt for garbage collection and recycle*/
#ifdef DISPENSER
      start_node_collection(world);
#endif
      /* tree and treetimes are not yet allocated */
#ifdef BEAGLE
      //init_mutationmodel(world, data, options);
      world->beagle = (beagle_fmt *) mycalloc(1,sizeof(beagle_fmt));
#endif
      world->data->haplotyping = options->haplotyping;      
      world->data->haplotyping_report = options->haplotyping_report;
      world->haplotypes    = (individualDB_fmt **) mycalloc(data->loci,sizeof(individualDB_fmt *));
      world->numhaplotypes    = (long *) mycalloc(data->loci,sizeof(long));
      // warning buffer init
      world->warningsize=0;
      world->warningallocsize=0;
      world->warning=NULL;
      //world->steppingstone = 0.0;
      //world->steppingstone_scalar = (double) -HUGE;
      //world->steppingstone_counter = 0;
}

#ifdef LONGSUM
/// setup the rates used for variable population sizes through time
void setup_fluctuate(world_fmt *world, option_fmt * options)
{
    long i,ii;
    long numpop = world->numpop;
    long numpop3 = 3 * numpop;
    long numpop6 = 6 * numpop;
    if (options->fluctuate)
    {
        world->options->fluctuate = options->fluctuate;
        for (i = 0, ii=0; i < numpop6; i+=2, ii++)
        {
            //even elements are rates, (odd are times)
            world->flucrates[ii] = options->flucrates[i];
            world->lflucrates[ii] = LOG (options->flucrates[i]);
            world->flucrates[ii+numpop3] = options->flucrates[i+1];
        }
    }
    else
    {
        for (i = 0; i < numpop3; i++)
        {
            world->flucrates[i] = 1.; //all have the same rate
            world->lflucrates[i] = 0.;
            world->flucrates[i+numpop3] = MYREAL_MAX/(numpop3-i); //all times are large
        }
    }
}
#endif /*LONGSUM*/


/// set start, end, and incremetn values for chains
void
set_bounds (long *increment, long *steps, long *chains,
            const option_fmt * options, const char type)
{
    switch (type)
    {
        case 's':
            *increment = options->sincrement;
            *steps = options->ssteps;
            *chains = options->schains;
            break;
        case 'l':
            *increment = options->lincrement;
            *steps = options->lsteps;
            *chains = options->lchains;
            break;
        default:
            error ("Wrong chain type\n");
            //break;
    }
}


/// print results
void
print_list (world_fmt ** universe, option_fmt * options, data_fmt * data)
{
#ifdef MPI
    long maxreplicate;
    maxreplicate = (options->replicate ?
                    ((options->replicatenum > 0) ?
                     options->replicatenum  : options->lchains ) : 1);
    if(!options->readsum)
      {
#ifdef DEBUG_MPI
	printf("%i> before unpacking results buffer in master",myID);
#endif
	mpi_results_master (MIGMPI_RESULT, EARTH, maxreplicate,
			    unpack_result_buffer);
      }
#endif
    
    print_results (universe, options, data);
#ifdef LONGSUM
    
    print_fluctuate(universe,options,data);
#endif /*LONGSUM*/
    
    if (!options->readsum)
    {
        if (options->printfst)
            print_fst (EARTH, options, data, EARTH->fstparam);
    }
    //if (EARTH->options->plot && EARTH->numpop>1)
    //{
    //    plot_surface (EARTH, options, data, EARTH->plane,
    //                  EARTH->options->plotintervals + 2);
    //}
    //#ifdef PRETTY
    //MLstuff	    pdf_print_results (universe, options, data);
	    //	    pdf_print_mcmctable(&page_height, page_width, EARTH, data, options);
	    //#endif
}

/// calculate location of tickmakrs
void set_ticks (MYREAL **ticks, MYREAL *plotrange, short type)
{
    long i;
    MYREAL diff;
    (*ticks)[0] = plotrange[0];
    (*ticks)[PLANEBIGTICKS - 1] = plotrange[1];
    if (type == PLOTSCALELOG)
    {
        diff =
        (log10 (plotrange[1]) - log10 (plotrange[0])) / (PLANEBIGTICKS - 1.);
        for (i = 1; i < PLANEBIGTICKS - 1; i++)
        {
            (*ticks)[i] = pow (10., (log10 ((*ticks)[0]) + i * diff));
        }
    }
    else
    {
        diff = (plotrange[1] - plotrange[0]) / (PLANEBIGTICKS - 1.);
        for (i = 1; i < PLANEBIGTICKS - 1; i++)
        {
            (*ticks)[i] = (*ticks)[0] + i * diff;
        }
    }
}

///// plot the likelihood surface
//void
//plot_surface (world_fmt * world, option_fmt * options,
//              data_fmt * data, char ****plane, long x)
//{
//    long locus;
//    long loci = world->loci;
//    FILE *outfile = world->outfile;
//    MYREAL *ticks;
//    MYREAL prangex[2];
//    MYREAL prangey[2];
//    ticks = (MYREAL *) mycalloc (1, sizeof (MYREAL) * PLANEBIGTICKS);
//    prangex[0] = world->options->plotrange[0];
//    prangex[1] = world->options->plotrange[1];
//    prangey[0] = world->options->plotrange[2];
//    prangey[1] = world->options->plotrange[3];
//    if (world->options->progress)
//        FPRINTF (stdout, "           Plotting the likelihood surfaces\n");
//    if (world->options->writelog)
//        FPRINTF (world->options->logfile,
//                 "           Plotting the likelihood surfaces\n");
//    PAGEFEED;
//    FPRINTF (outfile,
//             "Ln-Likelihood surfaces for each of the %3li populations\n",
//             world->numpop);
//    FPRINTF (outfile,
//             "-------------------------------------------------------\n\n");
//    FPRINTF (outfile, "Legend:\n");
//    FPRINTF (outfile, "   X = Maximum likelihood\n");
//    FPRINTF (outfile, "   * = in approximative 50%% confidence limit\n");
//    FPRINTF (outfile, "   + = in approximative 95%% confidence limit\n");
//    FPRINTF (outfile, "   - = in approximative 99%% confidence limit\n");
//    set_ticks (&ticks, prangex, world->options->plotscale);
//    FPRINTF (outfile, "X-tickmarks are (1) %f, (2) %f, (3) %f\n",
//             ticks[0], ticks[1], ticks[2]);
//    FPRINTF (outfile, "                (4) %f, (5) %f, (6) %f\n",
//             ticks[3], ticks[4], ticks[5]);
//    set_ticks (&ticks, prangey, world->options->plotscale);
//    FPRINTF (outfile, "Y-tickmarks are (1) %f, (2) %f, (3) %f\n",
//             ticks[0], ticks[1], ticks[2]);
//    FPRINTF (outfile, "                (4) %f, (5) %f, (6) %f\n",
//             ticks[3], ticks[4], ticks[5]);
//    /* change: show only over all */
//    if (!options->readsum)
//    {
//        if (world->loci == 1)
//        {
//            /*
//             * for (locus = 0; locus < loci; locus++) { if
//                 * (world->data->skiploci[locus]) { continue; }
//                 * FPRINTF(outfile, "\n\nLocus %li\n", locus + 1);
//                 */
//            
//            plot_surface_header2 (outfile, 0,
//                                  world->options->plotvar ==
//                                  PLOT4NM ? "xNm" : "M");
//            plot_surface2 (outfile, PLANESIZE + 2, 0, world->numpop, plane,
//                           world->plotmax, data->popnames,
//                           world->options->plotvar == PLOT4NM ? "xNm" : "M",
//                           world->options->migration_model);
//             }
//        }
//    if ((loci - world->skipped > 1) || (options->readsum))
//    {
//        FPRINTF (outfile, "\n%s\n\n",
//                 (options->readsum) ? ((loci - world->skipped >
//                                        1) ? "Over all loci" : "Over locus 1") :
//                 "Over all loci");
//        locus = (options->readsum) ? ((loci - world->skipped >
//                                       1) ? loci : 0) : loci;
//        
//        plot_surface_header2 (outfile, locus,
//                              world->options->plotvar == PLOT4NM ? "xNm" : "M");
//        plot_surface2 (outfile, PLANESIZE + 2, locus, world->numpop, plane,
//                       world->plotmax, data->popnames,
//                       world->options->plotvar == PLOT4NM ? "xNm" : "M",
//                       world->options->migration_model);
//    }
//    }
//
///// print plot header
//void
//plot_surface_header2 (FILE * outfile, long locus, char plotvar[])
//{
//    FPRINTF (outfile,
//             "    x-axis= %3.3s [xNm = effective population size * migration rate = Theta * M\n                  M   = migration rate / mutation rate = m/mu],\nx=1, 2, or 4 for mtDNA, haploid, or diploid data\n",
//             plotvar);
//    FPRINTF (outfile, "    y-axis = Theta,\n    units = see above\n");
//}
//
///// print plot surface
//void
//plot_surface2 (FILE * outfile, long x, long locus, long numpop,
//               char ****plane, plotmax_fmt ** plotmax,
//               char **popnames, char plotvar[], long migrmodel)
//{
//    long pop, i;
//    if (migrmodel == ISLAND)
//    {
//        FPRINTF (outfile, "N-ISLAND MODEL: For each Population\n");
//        FPRINTF (outfile,
//                 "  **Average** immigration: %s=%f, Theta=%f, log likelihood=%f\n",
//                 plotvar, plotmax[locus][0].x1, plotmax[locus][0].y1,
//                 plotmax[locus][0].l1);
//        FPRINTF (outfile, "[Remember: the maximum values are from a grid]\n");
//        FPRINTF (outfile, "\n            Mean(Immigration)  \n\n");
//        for (i = x; i >= 0; i--)
//            FPRINTF (outfile, "%46.46s\n", plane[locus][0][i]);
//        FPRINTF (outfile, "%46.46s\n", plane[locus][0][x]);
//        PAGEFEED;
//    }
//    else
//    {
//        for (pop = 0; pop < numpop; pop++)
//        {
//            if (popnames == NULL)
//                FPRINTF (outfile, "Population %li\n", pop + 1);
//            else
//                FPRINTF (outfile, "Population %li: %s\n", pop + 1, popnames[pop]);
//            FPRINTF (outfile,
//                     "  **Average** immigration: %s=%f, Theta=%f, log likelihood=%f\n",
//                     plotvar, plotmax[locus][pop].x1, plotmax[locus][pop].y1,
//                     plotmax[locus][pop].l1);
//            FPRINTF (outfile,
//                     "  **Average** emigration: %s=%f, Theta=%f, log likelihood=%f\n",
//                     plotvar, plotmax[locus][pop].x2, plotmax[locus][pop].y2,
//                     plotmax[locus][pop].l2);
//            FPRINTF (outfile,
//                     "[Remember: the maximum values are from a grid]\n");
//            FPRINTF (outfile,
//                     "\n           Mean(Immigration)                       Mean(Emigration)\n\n");
//            for (i = x; i >= 0; i--)
//                FPRINTF (outfile, "%s\n", plane[locus][pop][i]);
//            FPRINTF (outfile, "%s\n", plane[locus][pop][x]);
//            PAGEFEED;
//        }
//    }
//}


/// free all the world structures from the array universe
#ifdef MPI
void
free_universe (world_fmt ** worlds, long numworlds, option_fmt *options)
{
    long i;
    if (myID != MASTER)
      {
	for(i = numworlds-1; i>=0; i--)
	  free_world(worlds[i],options);
	myfree(worlds);
      }
    else
      free_world(worlds[0],options);
}
#else
void
free_universe (world_fmt ** worlds, long numworlds, option_fmt *options)
{
    long i;
    for(i = numworlds-1; i>=0; i--)
        free_world(worlds[i],options);
    myfree(worlds);
}
#endif
/// free the contents of world structure
void
free_world(world_fmt *world, option_fmt *options)
{
    long i,j,z, locus;
#ifdef MPI
    long maxnum;
#endif
    long sumloc = world->loci > 1 ? 1 : 0;
    // timevector is already freed in mcmc1.c

#ifdef MPI
    myfree(world->who);
    myfree(world->mpistack);
  if(myID==MASTER && world->cold)
    {
      maxnum = (options->replicate == TRUE) ? 
	((options->replicatenum > 0) ? options->replicatenum  : options->lchains ) : 1;
      maxnum = world->loci * maxnum;
      if(maxnum < numcpu)
      	maxnum = numcpu;
      for(i=0;i<maxnum;i++)
	myfree(world->mpistack_request[i].tempstr);
      myfree(world->mpistack_request);
    }
#ifdef SLOWNET
    myfree(world->profilewho);
#endif
#endif

//    if (world->options->plot)
//    {
//        for (i = 0; i < world->loci + sumloc; i++)
//        {
//            for (j = 0; j < world->numpop; j++)
//                myfree(world->plane[i][j]);
//            myfree(world->plane[i]);
//        }
//        myfree(world->plane);
//	myfree(world->options->plotxvalues);
//	myfree(world->options->plotyvalues);
//
//    }

    myfree(world->burnin_stops);
    
    destroy_mutationmodel(world);

    if(world->atl != NULL)
      {
	for (j = 0; j < world->repstop; j++)
	  {
	    for (i = 0; i < world->loci + 1; i++)
	      {
		for (z = 0; z < world->atl[j][i].allocT; z++)
		  {
		    myfree(world->atl[j][i].tl[z].data);
#ifdef ALTIVEC
		    myfree(world->atl[j][i].tl[z].vdata);
#endif
		  }
		myfree(world->atl[j][i].tl);
		myfree(world->atl[j][i].parameters);
	      }
	  }
	myfree(world->atl[0]);
	myfree(world->atl);
      }
    if (strchr (SEQUENCETYPES, world->options->datatype))
    {
        if (world->tbl != NULL)
        {
            for (i = 0; i < world->options->rcategs; i++)
            {
                for (j = 0; j < world->options->categs; j++)
                    myfree(world->tbl[i][j]);
                myfree(world->tbl[i]);
            }
	    myfree(world->tbl);
        }
	myfree(world->data->seq[0]->alias);
	myfree(world->data->seq[0]->ally);
	myfree(world->data->seq[0]->aliasweight);
	myfree(world->data->seq[0]->location);
	myfree(world->data->seq[0]->category);
	myfree(world->data->seq[0]->weight);
	myfree(world->data->seq[0]);
    }
    if (!options->readsum)
      {
	myfree(world->data->maxalleles); // frees also data->sites
      }
    myfree(world->data->seq);
    myfree(world->data->locusweight);
    if(world->cold)
      {
	if (options->recorddivtime)
	  myfree(world->divtime);
	// free skyline histograms
	destroy_expected_events(world);
	// free migrate histogram
	destroy_mighist (world);
	// free convergence section
	myfree(world->convergence->chain_s);
	myfree(world->convergence->chain_means);
	myfree(world->convergence->chain_counts);
	myfree(world->convergence->gelmanmeanmaxR);
	myfree(world->convergence);
      }
   
    empty_world_unassigned(world);
    if(world->has_unassigned)
      {
	myfree(world->unassigned[0]->key);
	myfree(world->unassigned[0]);
	myfree(world->unassigned);
      }
    
    // free bayes part if necessary
    if(world->options->bayes_infer)
        bayes_free(world);

    // free the marginal likelihood
    myfree(world->bfscale);
    // free the autocorrelation material
    myfree(world->autocorrelation);
    if(world->cold)
      myfree(world->auto_archive);
    // free the data section in world
    myfree(world->data->skiploci);
    myfree(world->data->geo);
    myfree(world->data->lgeo);
    myfree(world->data);

    // free the option section in world
    myfree(world->options->mutationrate_year);
    if(world->options->mu_rates!=NULL)
      {
	myfree(world->options->mu_rates); //lmu_rates are freed here, too
	//  printf("free murate in wopt");
      }
    if(world->options->inheritance_scalars!=NULL)
      {
	myfree(world->options->inheritance_scalars); 
      }
    myfree(world->options->micro_threshold);

    if(world->cold)
      {
	if(world->options->has_bayesmdimfile)
	  {
	    myfree(mdimfilecount);
	  }
      }
    myfree(world->options->custm); // custm2 is freed here, too
    myfree(world->options->thetag);
    myfree(world->options->mg);
    myfree(world->options->zeroparam);
    myfree(world->options->constparam);
    myfree(world->options->symparam);
    myfree(world->options->sym2param);
    myfree(world->options->mmparam);
    myfree(world->options->rrate);
    myfree(world->options->rate);
    myfree(world->options->probcat);
    myfree(world->options->meanmu);
    if(world->options->steps != NULL)
      {
	for (locus = 0; locus < world->loci; locus++)
	  {
	    for (i = 0; i < world->options->micro_stepnum; i++)
	      {
		myfree(world->options->steps[locus][i]);
	      }
	    myfree(world->options->steps[locus]);
	  }
	myfree(world->options->steps);
      }
    //
    // free the treefile repository if necessary
    if(world->options->treeinmemory)
      {
	for(i=0;i<world->loci;i++)
	  myfree(world->treespace[i]);
	myfree(world->treespace);
	myfree(world->treespacenum);
	myfree(world->treespacealloc);
	myfree(world->besttreelike);
      }
    myfree(world->options);
    if(world->fstparam != NULL)
      {
	for (i = 0; i < world->loci + sumloc; i++)
	  {
	    myfree(world->fstparam[i]);
	  }
	myfree(world->fstparam);
      }
#ifdef LONGSUM
    myfree(world->flucrates);
    myfree(world->lflucrates);
#endif /*LONGSUM*/
#ifdef UEP
    destroy_uep(world, world->data->uepsites);
#endif
    myfree(world->mig0list);
    // migration and speciation
    for (i = 0; i < 2*world->numpop; ++i)
      myfree(world->migproblist[i]);
    myfree(world->migproblist);
    myfree(world->design0list);

    myfree(world->contribution);
    //    printf("%i> temp=%f: world->contribution freed\n",myID, world->heat);
    myfree(world->buffer);
    // already down but it leaks some (valgrind) free_tree(world->root, world);
    //myfree(world->nodep);
    myfree(world->likelihood);
    myfree(world->lineages);
    myfree(world->param0);
    myfree(world->param00);
    myfree(world->chainlikes[0]);
    myfree(world->chainlikes);
    myfree(world->mstart);
    myfree(world->mend);
    // remove garbage/recycler  collector
    #ifndef MPI
// this fails with MPI, since this cleaning up step after everything
// is finished fails. with threading and standard runs it seems to work
#ifdef DISPENSER
    warning ("check me because with MPI this will not be done");
    stop_node_collection(world);
#endif
    #endif
    myfree(world->sublocistarts);
    myfree(world->numsubloci);
    myfree(world->maxnumpattern);
    // free world
    //    fprintf(stdout,"Freed: world with temperature %f\n",1./world->heat);
    myfree(world);
}

/// print progress report of a locus with replicate information
//OBSOLETE
void
print_menu_locus (FILE *file, world_fmt * world, long locus)
{
    // printing machinery
    char *buffer;
    char *bufptr;
    long bufsize=0;
    char nowstr[STRSIZE];
    buffer = (char *) mycalloc(STRSIZE,sizeof(char));
    bufptr = buffer;
#ifdef MPI
    bufsize = sprintf(buffer, "[%3i] ",myID);
#endif
    get_time (nowstr, "%H:%M:%S");
    if (world->options->replicate)
      sprintf(buffer + bufsize, "%s   Locus %li: Replicate %li\n", nowstr, locus + 1,
                               world->replicate + 1);
    else
      sprintf(bufptr + bufsize, "%s   Locus %li:\n", nowstr, locus + 1);
    FPRINTF(file,"%s",buffer); //with MPI it send this to the master
    myfree(buffer);
}


/// print progress information after a chain finished
void
print_menu_chain (char type, long chain, long steps,
                  world_fmt * world, option_fmt * options, long rep)
{
    // printing machinery
    char *buffer;
    char *tempbuffer;
    long allocbufsize=LINESIZE;
    long bufsize=0;
    //
    long i, j, k;
    char strllike[STRSIZE];
    char nowstr[STRSIZE];
    char spacer[] = "           ";
    MYREAL value;
    MYREAL mini, maxi;
    MYREAL summ = 0;
    boolean writelog = world->options->writelog;
    boolean progress = world->options->progress;
    
    if(options->bayes_infer)
        return;
    
    if(writelog || progress)
    {
        buffer = (char *) mycalloc(allocbufsize,sizeof(char));
        tempbuffer = (char *) mycalloc(allocbufsize,sizeof(char));
        print_llike (world->param_like, strllike);
        
        get_time (nowstr, "%H:%M:%S");
        if (chain == FIRSTCHAIN)
        {
	  print_to_buffer(&buffer, &allocbufsize, tempbuffer, &bufsize,"%s   [NODE:%i] Locus %li: Start conditions:\n", nowstr, myID, world->locus + 1);
	  print_param (&buffer, &bufsize, &allocbufsize, world->options->usem, world, world->numpop, spacer);
          print_to_buffer(&buffer, &allocbufsize, tempbuffer, &bufsize, "%sStart-tree-log(L)=%f\n", spacer,world->likelihood[0]);
	  prognose_time (nowstr, world,(steps == options->sincrement) ? options->sincrement : options->lincrement, 
			 steps, spacer, FALSE);
        }
        else
        {
	  if (world->repkind == SINGLECHAIN)
	    print_to_buffer(&buffer, &allocbufsize, tempbuffer, &bufsize,"%s   [NODE:%i] Locus %li: %*s chain %3li:\n", nowstr, myID,
			    world->locus + 1, type == 's' ? 5 : 4,
			    type == 's' ? "Short" : "Long", chain + 1);
            else
	      print_to_buffer(&buffer, &allocbufsize, tempbuffer, &bufsize, "%s   Multiple chain estimates:\n", nowstr);
            
            if(options!=NULL)
	      prognose_time (nowstr, world,(steps == options->sincrement) ? options->sincrement : options->lincrement, 
			     steps, spacer, FALSE);
            if(!world->options->bayes_infer)
            {
	      print_to_buffer(&buffer, &allocbufsize, tempbuffer, &bufsize,"%sln L(Param)=ln (P(G|Param)/P(G|Param0)):%15.15s\n", spacer, strllike);
                
	      if (world->repkind == SINGLECHAIN)
		print_to_buffer(&buffer, &allocbufsize, tempbuffer, &bufsize,"%sStart-tree-ln(L)=ln P(D|G):                %10.5f\n", spacer, world->likelihood[0]);
            }
	    
            print_param (&buffer, &bufsize, &allocbufsize, world->options->usem, world, world->numpop, spacer);
            
            if (world->options->verbose && world->repkind == SINGLECHAIN)
            {
                mini = MYREAL_MAX;
                maxi = -MYREAL_MAX;
                for (i = 0; i < steps; i++)
                {
                    value = world->likelihood[i];
                    if (mini > value)
                        mini = value;
                    if (maxi < value)
                        maxi = value;
                }
                print_to_buffer(&buffer, &allocbufsize, tempbuffer, &bufsize, "           Sampled tree-log(L)={%f .. %f},\n",
				mini, maxi);
                print_to_buffer(&buffer, &allocbufsize, tempbuffer, &bufsize, "           Best of all visited trees =%f\n",
                                   world->maxdatallike);
                /* report migration events */
                print_to_buffer(&buffer, &allocbufsize, tempbuffer, &bufsize,
                                   "           Average migration events per genealogy (x-axis=to, y=from)\n           ");
                for (i = 0; i < world->numpop; i++)
                {
                    for (j = 0; j < world->numpop; j++)
                    {
                        if (i != j)
                        {
                            summ = 0;
                            for (k = 0; k < world->atl[rep][world->locus].T;k++)
                                summ +=
                                    world->atl[rep][world->locus].tl[k].
                                    mindex[mm2m (j, i, world->numpop)] *
                                    world->atl[rep][world->locus].tl[k].copies;
                            print_to_buffer(&buffer, &allocbufsize, tempbuffer, &bufsize, "%6.2f ",
                                               ((MYREAL) summ) / ((MYREAL) steps));
                        }
                        else
                        {
			  print_to_buffer(&buffer, &allocbufsize, tempbuffer, &bufsize, "------- ");
                        }
                    }
                    print_to_buffer(&buffer, &allocbufsize, tempbuffer, &bufsize, "\n           ");
                }
                print_to_buffer(&buffer, &allocbufsize, tempbuffer, &bufsize, "\n");
            }
            
        }
        
        if(progress)
        {
	  LARGEFPRINTF(stdout,allocbufsize, "%s",buffer);
        }
        if(writelog)
        {
	  LARGEFPRINTF(world->options->logfile,allocbufsize, "%s",buffer);
        }
	myfree(tempbuffer);
        myfree(buffer);
    }
}

/// prognose how long the program will run
void
prognose_time (char *nowstr, world_fmt * world,
              long options_increment, long steps, char *spacer, boolean tobuffer)
{
  (void) options_increment;
  (void) steps;
#ifndef NOTIME_FUNC
  time_t starttime;
  MYREAL elapsedtime;
  MYREAL mytime;
    time_t nowbin;
    time_t proposetime;
    struct tm *nowstruct;
    //    long increment = options_increment;
    // printing machinery
    char *buffer;
    long bufsize=0;
    boolean writelog = world->options->writelog;
    boolean progress = world->options->progress;
    //
    if(writelog || progress)
    {
        buffer = (char *) mycalloc(STRSIZE,sizeof(char));
        if (world->treesdone == 0)
        {
            strcpy (nowstr, "indeterminate");
        }
        else
        {
	  mytime = (MYREAL) (time(&nowbin));
	    starttime = world->starttime;
	    elapsedtime = mytime - (MYREAL) starttime;
#ifdef MPI
            proposetime = (time_t) (elapsedtime / ((MYREAL) world->treesdone) * world->treestotal) / (numcpu - 1);
#else
	    proposetime = (time_t) (elapsedtime / ((MYREAL) world->treesdone) * world->treestotal);
#endif            
	    nowbin = ((time_t) starttime) + proposetime;
	    nowstruct = localtime (&nowbin);
	    strftime (nowstr, STRSIZE, "%H:%M %B %d %Y", nowstruct);
        }
        
        bufsize += sprintf(buffer+bufsize, "%sPrognosed end of sampling is    %s\n", spacer, nowstr);
        if(!tobuffer)
        {
            if (progress)
            {
	      LARGEFPRINTF(stdout,bufsize, "%s",buffer);
            }
            if(writelog)
            {
	      LARGEFPRINTF(world->options->logfile,bufsize, "%s",buffer);
            }
        }
        myfree(buffer);
    }
#endif /*time construct*/
}

boolean is_same(long i, long j)
{
  return (boolean) (i==j);
}

/// copy the times from the working genealogy to the minimal statistics database
void
copy_time (world_fmt * world, timelist_fmt * ltl, long from,
           long to, long np, long rep)
{
  long j;
    long T;
    MYREAL t;
    tarchive_fmt *tl;
    //    mighist_fmt *aa = NULL;
    mighistloci_fmt *bb = NULL;

    if (from == -1)
    {
        (*ltl).copies = 1;
        T = (*ltl).T - 1;
        increase_timearchive (world, world->locus, 1, np, world->rep);
        tl = world->atl[rep][world->locus].tl;
	//  @TEST DEBUG      tl[0].copies = (*ltl).copies;
        tl->copies = (*ltl).copies;
        archive_timelist (tl, (*ltl).tl, T, np, world);
        if (world->in_last_chain && world->options->mighist)
        {
	  //if(world->options->mighist_counter++ == world->options->mighist_increment)
	  //  {
	      world->options->mighist_counter = 0;
	      //    increase_mighist(&world->mighistloci[world->locus]);
	      //world->mighistloci[world->locus].mighist[world->mighistloci[world->locus].mighistnum].copies = 1;
	      //world->mighistloci[world->locus].mighistnum++;
	      //}
	}
        return;
    }
    if (from == to)
      {
        (*ltl).copies += 1;
        tl = world->atl[rep][world->locus].tl;
        tl[from].copies = (*ltl).copies;
        if (world->in_last_chain && world->options->mighist)
	  {
	    world->options->mighist_counter = 0;
	    bb = &(world->mighistloci[world->locus]);
	    if((*ltl).tl[0].eventnode->type != 't' && (*ltl).tl[0].visited)
	      {
		calculate_event_values(bb->migeventbins, bb->migeventbinnum, bb->eventbinsize, 
				       (*ltl).tl[0].age, 
				       (*ltl).tl[0].age, 
				       (*ltl).tl[0].eventnode->type,
				       (*ltl).tl[0].from, (*ltl).tl[0].to,
				       (*ltl).tl[0].lineages, world->numpop, FALSE, world);
		// calculate skyline values
		if(world->options->skyline)
		  {
		    calculate_expected_values(bb->eventbins, bb->eventbinnum, bb->eventbinsize, 
					      (*ltl).tl[0].age, 
					      (*ltl).tl[0].age, 
					      (*ltl).tl[0].eventnode->type,
					      (*ltl).tl[0].from, (*ltl).tl[0].to,
					      (*ltl).tl[0].lineages, world->numpop, world);
		  }
	      }
	    T = (*ltl).T - 1;
	    for (j = 1; j < T; j++)
	      {
		if((*ltl).tl[j].eventnode->type != 't' && (*ltl).tl[j].visited)
		  {
		    continue;
		  }
		t = (*ltl).tl[j].age - (*ltl).tl[j - 1].age;
		// calculate event values
		calculate_event_values(bb->migeventbins, bb->migeventbinnum, bb->eventbinsize, 
				       t, 
				       (*ltl).tl[j].age, 
				       (*ltl).tl[j].eventnode->type,
				       (*ltl).tl[j].from, (*ltl).tl[j].to,
				       (*ltl).tl[j].lineages, world->numpop, is_same(j,T-1), world);
		
		// calculate skyline values
		if(world->options->skyline)
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
    else
      {
        T = (*ltl).T - 1;
        increase_timearchive (world, world->locus, to + 1, np, world->rep);
        tl = world->atl[rep][world->locus].tl;
        tl[to].copies = (*ltl).copies = 1;
        archive_timelist (&(tl[to]), (*ltl).tl, T, np, world);
      }
}

/// precalculates values that are needed later during the MCMC chains and maximization process
/// added speciation processes in 2013 add 2x to the migproblist, first part mig second part spec
void
precalc_world (world_fmt * world)
{
    long numpop = world->numpop;
    long pop;
    long msta;
    long msto;
    long i, j;
    boolean usem = world->options->usem;//@@
    MYREAL px;
    MYREAL tmp = 0.;
    memset (world->mig0list, 0, sizeof (MYREAL) * (size_t) world->numpop);
    if (numpop > 1)
    {
        for (i = 0; i < 2*numpop; ++i)
            memset (world->migproblist[i], 0,
                    sizeof (MYREAL) * (size_t) (world->numpop)); //was -1 but I allocate numpop elements, and use -1
    }
    memset (world->design0list, 0, sizeof (long) * (size_t) world->numpop);
    for (pop = 0; pop < numpop; pop++)
    {
        msta = world->mstart[pop];
        msto = world->mend[pop];
        for (i = msta, j = 0; i < msto; i++, j++)
        {
	  //1229 deletion of 'd' section
	      px = usem ? world->param0[i] : world->param0[i]/world->param0[pop];//@@
	      tmp                        = world->data->geo[i] * px;
	      world->mig0list[pop]      += tmp;
	      world->migproblist[pop][j] = tmp;
	      world->design0list[pop]   += world->options->custm2[i] == '0' ? 1 : 0;
	      world->speciesproblist[pop][j] = 0.0;
	}
        for (j = 1; j < numpop - 1; ++j)
	  {
            world->migproblist[pop][j] += world->migproblist[pop][j - 1];
            world->speciesproblist[pop][j] += world->speciesproblist[pop][j - 1];
	  }
    }
    for (pop = 0; pop < numpop; pop++)
    {
        for (j = 0; j < numpop - 1; ++j)
        {
            if (world->migproblist[pop][numpop - 2] > 0.0)
                world->migproblist[pop][j] /= world->migproblist[pop][numpop - 2];
            else
                world->migproblist[pop][j] = 0.;
            if (world->speciesproblist[pop][numpop - 2] > 0.0)
                world->speciesproblist[pop][j] /= world->speciesproblist[pop][numpop - 2];
            else
                world->speciesproblist[pop][j] = 0.;
        }
    }
}

///
/// re-precalculates values that are needed later during the MCMC chains and maximization process
/// instead of recalculating everything only the values involved with parameter that are recalculated
void
reprecalc_world (world_fmt * world, long that)
{
    long numpop = world->numpop;
    long msta;
    long msto;
    long i, j;
    long thatpop = that;
    boolean usem = world->options->usem;//@@
    MYREAL px;
    MYREAL tmp = 0.;
    //    memset (world->mig0list, 0, sizeof (MYREAL) * world->numpop);
    if(that>=numpop && that < world->numpop2)
      {
	thatpop = (that - numpop) / (numpop-1);
	world->mig0list[thatpop] = 0;
	if (numpop > 1)
	  {
	    memset (world->migproblist[thatpop], 0,
                    sizeof (MYREAL) * (size_t) (world->numpop - 1));
	    memset (world->speciesproblist[thatpop], 0,
                    sizeof (MYREAL) * (size_t) (world->numpop - 1));
	  }
        //    memset (world->design0list, 0, sizeof (long) * world->numpop);
	world->design0list[thatpop] = 0;
	msta = world->mstart[thatpop];
	msto = world->mend[thatpop];
	for (i = msta, j = 0; i < msto; i++, j++)
	  {
	    //1229 deleted 'd' section
	      px = usem ? world->param0[i] : world->param0[i]/world->param0[thatpop];//@@
	      tmp                        = world->data->geo[i] * px;
	      world->mig0list[thatpop]      += tmp;
	      world->migproblist[thatpop][j] = tmp;
	      world->design0list[thatpop]   += world->options->custm2[i] == '0' ? 1 : 0;
	      world->speciesproblist[thatpop][j] = 0.0;
	  }
	for (j = 1; j < numpop - 1; ++j)
	  {
	    world->migproblist[thatpop][j] += world->migproblist[thatpop][j - 1];
            world->speciesproblist[thatpop][j] += world->speciesproblist[thatpop][j - 1];
	  }
	for (j = 0; j < numpop - 1; ++j)
	  {
	    if (world->migproblist[thatpop][numpop - 2] > 0.0)
	      world->migproblist[thatpop][j] /= world->migproblist[thatpop][numpop - 2];
	    else
	      world->migproblist[thatpop][j] = 0.;
            if (world->speciesproblist[thatpop][numpop - 2] > 0.0)
                world->speciesproblist[thatpop][j] /= world->speciesproblist[thatpop][numpop - 2];
            else
                world->speciesproblist[thatpop][j] = 0.;
	  }
      }
}

/// create the database that contains the minimal statistic for the maximization process
void
create_timearchive (timearchive_fmt *** atl, long loci,
                    long samples, long numpop, long replicates)
{
    long h, i, j;
    long numpop2 = numpop * numpop;
    // experiment with unrolling in probG4 padding so that data vector is multiplier of 3
    long size = numpop2 + numpop + numpop;
    long datasize = size + (size % 3);
#ifdef ALTIVEC
    size -= (size % 4) - 4;
    size /= 4;
#endif /*ALTIVEC*/
    
    (*atl) =
        (timearchive_fmt **) mycalloc (replicates + 1, sizeof (timearchive_fmt *));
    (*atl)[0] =
      (timearchive_fmt *) mycalloc (((replicates + 1) * (2 + loci)),
                                      sizeof (timearchive_fmt));
    for (h = 1; h < replicates+1; h++)
      (*atl)[h] = (*atl)[0] + h * (2 + loci);
    for (h = 0; h < replicates+1; h++)
    {
        for (i = 0; i < loci + 2; i++)
        {
#ifdef ALTIVEC
            //on ALTIVEC processors we want ot shift some data into streamlined
            // vectors so that they can be processed in calc_locus_like at
            // max speed
            //alllog(copies) padded so that they are allows multyplies of 4
            // because the calc_locus_function will stuff 4 vectors at the same time into the
            // altivec unit to fill the pipeline.
            //(*atl)[h][i].lcopiesvec = (FloatVec *) mycalloc(samples + samples % 4, sizeof(FloatVec));
            
            // contains all the compressed treesummaries
            //size = (2 * numpop + numpop2) * (samples + (samples % 4));
            //(*atl)[h][i].data = (FloatVec *) mycalloc(size, sizeof(FloatVec)); //all tl[].data
#endif /*ALTIVEC*/
            (*atl)[h][i].parameters =
	      (MYREAL *) mycalloc ((4 * (numpop2 + 1)), sizeof (MYREAL));
            //    (MYREAL *) mycalloc (3 * numpop2 + samples + 1, sizeof (MYREAL));
            (*atl)[h][i].param = (*atl)[h][i].parameters;
            (*atl)[h][i].param0 = (*atl)[h][i].parameters + numpop2 + 1;
            (*atl)[h][i].lparam0 = (*atl)[h][i].param0 + numpop2 + 1;
            (*atl)[h][i].likelihood = (*atl)[h][i].lparam0 + numpop2 + 1;
            (*atl)[h][i].tl =
                (tarchive_fmt *) mycalloc (samples, sizeof (tarchive_fmt));
            (*atl)[h][i].T = (*atl)[h][i].allocT = samples;
            for (j = 0; j < samples; j++)
            {
	      // experiment with unrolling in probG4 padding so that data vector is multiplier of 3
                (*atl)[h][i].tl[j].data =
                (MYREAL *) mycalloc (datasize, sizeof (MYREAL));
#ifdef ALTIVEC
                
                (*atl)[h][i].tl[j].vdata = (FloatVec *) mycalloc (size, sizeof (FloatVec));
#endif
                
                (*atl)[h][i].tl[j].point = (*atl)[h][i].tl[j].data;
                (*atl)[h][i].tl[j].wait = &(*atl)[h][i].tl[j].data[numpop2];
                (*atl)[h][i].tl[j].kt = (*atl)[h][i].tl[j].wait;
                (*atl)[h][i].tl[j].km = &(*atl)[h][i].tl[j].wait[numpop];
                (*atl)[h][i].tl[j].p = (*atl)[h][i].tl[j].point;
                (*atl)[h][i].tl[j].mindex = (*atl)[h][i].tl[j].point; // the first numpop element are not used
                                                                      // this is really ugly and should be changed
#ifdef LONGSUM
                
                (*atl)[h][i].tl[j].longsum = (longsum_fmt *) mycalloc(1,sizeof(longsum_fmt));
#endif /*LONGSUM*/
                
            }
        }
    }
}

/// increase the size of the data base of minimal statistics
void
increase_timearchive (world_fmt * world, long locus, long sample,
                      long numpop, long rep)
{
    long i = locus, j, oldT = 0, size;
    long numpop2 = world->numpop2;
    long datasize = numpop2 + world->numpop + world->numpop;
#ifdef ALTIVEC
    
    long vsize = datasize;
    vsize -= (vsize % 4) - 4;
    vsize /= 4;
#endif /*ALTIVEC*/
    // needs to be a multiple of 4 because of vdot_product;
    //datasize += datasize % 4;    
    if (sample >= world->atl[rep][i].allocT)
    {
        oldT = world->atl[rep][i].allocT;
        world->atl[rep][i].allocT = MAX (sample + 1,(long) (1.25 * oldT));
        size = world->atl[rep][i].allocT;
        world->atl[rep][i].tl = (tarchive_fmt *) myrealloc (world->atl[rep][i].tl, size * sizeof (tarchive_fmt));
        for (j = oldT; j < size; j++)
        {
            world->atl[rep][i].tl[j].data = (MYREAL *) mycalloc (datasize, sizeof (MYREAL));
#ifdef ALTIVEC
            world->atl[rep][i].tl[j].vdata = (FloatVec *) mycalloc (vsize, sizeof (FloatVec));
#endif
            
            world->atl[rep][i].tl[j].point = world->atl[rep][i].tl[j].data;
            world->atl[rep][i].tl[j].wait = &world->atl[rep][i].tl[j].data[numpop2];
            world->atl[rep][i].tl[j].kt = world->atl[rep][i].tl[j].wait;
            world->atl[rep][i].tl[j].km = &world->atl[rep][i].tl[j].wait[numpop];
            world->atl[rep][i].tl[j].p = world->atl[rep][i].tl[j].point;
            world->atl[rep][i].tl[j].mindex = world->atl[rep][i].tl[j].point; //ugly:the first numpop elem are not used
#ifdef LONGSUM            
            world->atl[rep][i].tl[j].longsum = (longsum_fmt *) mycalloc(1,sizeof(longsum_fmt));
#endif /*LONGSUM*/            
        }
        world->atl[rep][i].T = sample;
    }
    else
    {
        world->atl[rep][i].T = sample;
    }
}

///// initialize the structures for the plot planes that are then filled in create_plots()
//void
//init_plotplane (world_fmt * world)
//{
//    long pop, locus, i;
//    long plotsize = world->options->plotintervals;
//    world->plane =
//        (char ****) mycalloc (1, sizeof (char ***) * (world->loci + 1));
//    world->plotmax =
//        (plotmax_fmt **) mycalloc (1, sizeof (plotmax_fmt *) * (world->loci + 1));
//    world->plotmax[0] =
//        (plotmax_fmt *) mycalloc (1,
//                                  sizeof (plotmax_fmt) * (world->loci +
//                                                          1) * world->numpop);
//    for (pop = 1; pop < world->loci + 1; pop++)
//    {
//        world->plotmax[pop] = world->plotmax[0] + pop * (world->numpop);
//    }
//    for (locus = 0; locus < world->loci + 1; locus++)
//    {
//        world->plane[locus] =
//        (char ***) mycalloc (1, sizeof (char **) * world->numpop);
//        for (pop = 0; pop < world->numpop; pop++)
//        {
//            world->plane[locus][pop] =
//            (char **) mycalloc (1, sizeof (char *) * (plotsize + 3));
//            for (i = 0; i < plotsize + 3; i++)
//            {
//                world->plane[locus][pop][i] =
//                (char *) mycalloc (1, sizeof (char) * plotsize + plotsize + 20);
//            }
//        }
//    }
//}

void
create_cov (world_fmt * world)
{
    long locus, i;
    world->cov =
      (MYREAL ***) mycalloc ((world->loci+1), sizeof (MYREAL **));
    for (locus = 0; locus < world->loci + 1; locus++)
    {
        world->cov[locus] =
        (MYREAL **) mycalloc ((world->numpop2 + 1), sizeof (MYREAL *));
        world->cov[locus][0] =
	  (MYREAL *) mycalloc (((world->numpop2 + 1)*(world->numpop2 + 1)),
                                 sizeof (MYREAL));
        for (i = 1; i < world->numpop2 + 1; i++)
        {
            world->cov[locus][i] =
            world->cov[locus][0] + i * (world->numpop2 + 1);
        }
    }
}

void
print_menu_equilib (world_fmt * world)
{
    char nowstr[STRSIZE];
    
    // printing machinery
    char *buffer;
    long bufsize=0;
    boolean writelog = world->options->writelog;
    boolean progress = world->options->progress;
    
    //
    if(writelog || progress)
    {
      if(world->cold)
	{
	  buffer = (char *) mycalloc(STRSIZE,sizeof(char));
	  get_time (nowstr, "%H:%M:%S");
#ifdef MPI
	  bufsize = sprintf(buffer,
		  "[%3i] %8.8s   Burn-in of %li steps (Locus: %li/%li, Replicate: %li/%li) \n",
		  myID, nowstr, 
		  world->options->burn_in * world->increment, 
		  1+world->locus, world->loci, 
		  1+world->rep, world->maxreplicate);
#else
	  bufsize = sprintf(buffer,
		  "%8.8s   Burn-in of %li steps (Locus: %li/%li, Replicate: %li/%li) \n", 
		  nowstr, 
		  world->options->burn_in * world->increment, 
		  1+world->locus, world->loci, 
		  1+world->rep, world->maxreplicate);
#endif
	  if (progress)
	    {
	      FPRINTF(stdout,"%s",buffer);
	    }
	  if(writelog)
	    {
	      FPRINTF(world->options->logfile,"%s",buffer);
	    }
	  myfree(buffer);
	}
    }
}


void
print_finish (world_fmt * world, long filepos)
{
    char nowstr[LINESIZE];
    if (myID==MASTER && world->options->progress)
    {
        get_time (nowstr, "%H:%M:%S");
        FPRINTF (stdout, "%s   Program finished\n", nowstr);
        if (world->options->writelog)
            FPRINTF (world->options->logfile, "%s   Program finished\n", nowstr);
    }
    get_time (nowstr, "%c");
    if (nowstr[0] != '\0')
    {
        if (filepos > 0)
            fseek (world->outfile, filepos, SEEK_SET);
        FPRINTF (world->outfile, "         finished at %s\n", nowstr);
    }
}

///
/// Archives the current tree as a minimal statistic and also archives the events
/// in a list for histograms through time (skyline plot)
void
archive_timelist (tarchive_fmt * atl, vtlist * tl, long T, long np,
                  world_fmt * world)
{
    long         j;
    long         i;
    //    long         tt;
    MYREAL       t;
    MYREAL       line;
  //  mighist_fmt  *aa    = NULL;
    mighistloci_fmt *bb = NULL;

    for (i = 0; i < np; i++)
    {
        line = (MYREAL) tl[0].lineages[i];
        atl->km[i] = line * tl[0].age;
        atl->kt[i] = line * (line - 1) * tl[0].age;
	if(tl[0].eventnode->type != 't')
	  {
	    atl->p[i] = (MYREAL) (whichp (tl[0].from, tl[0].to, i));
	  }
	else
	  {
	    atl->p[i] = 0.;
	  }
    }
    memset (atl->mindex + np, 0, sizeof (MYREAL) * (size_t) (np * (np - 1)));
    if (tl[0].from != tl[0].to)
        atl->mindex[mm2m (tl[0].from, tl[0].to, np)] += 1;
    
    if (world->in_last_chain && world->options->mighist)
    {
      //aa = &(world->mighistloci[world->locus].mighist[0]);
      bb = &(world->mighistloci[world->locus]);
      //      aa->migeventsize=0;  
      // calculate event values
      calculate_event_values(bb->migeventbins, bb->migeventbinnum, bb->eventbinsize, 
			     tl[0].age, 
			     tl[0].age, 
			     tl[0].eventnode->type,
			     tl[0].from, tl[0].to, tl[0].lineages, world->numpop, FALSE, world);
      // calculate skyline values
      if(world->options->skyline)
	{
	  calculate_expected_values(bb->eventbins, bb->eventbinnum, bb->eventbinsize, tl[0].age, 
				    tl[0].age,
				    tl[0].eventnode->type, 
				    tl[0].from, tl[0].to, tl[0].lineages, 
				    world->numpop, world);
	}
    }
    for (j = 1; j < T; j++)
      {
	t = tl[j].age - tl[j - 1].age;
	for (i = 0; i < np; i++)
	  {
	    line = (MYREAL) tl[j].lineages[i];
	    atl->km[i] += line * t;
	    atl->kt[i] += line * (line - 1.) * t;
	    if(tl[j].eventnode->type != 't')
	      {
		atl->p[i] += (MYREAL) (whichp (tl[j].from, tl[j].to, i));
	      }
	    else
	      {
		atl->p[i] = 0.;
	      }
	  }
	if (tl[j].from != tl[j].to)
	  atl->mindex[mm2m (tl[j].from, tl[j].to, np)] += 1;
	if (world->in_last_chain && world->options->mighist)
	  {
	    // calculate event values
	    calculate_event_values(bb->migeventbins, bb->migeventbinnum, bb->eventbinsize, 
				   tl[j].age - tl[j-1].age, 
				   tl[j].age, 
				   tl[j].eventnode->type,
				   tl[j].from, tl[j].to, tl[j].lineages, world->numpop, is_same(j, T-1),world);
	    
	    // calculate skyline values
	    if(world->options->skyline)
	      {
		calculate_expected_values(bb->eventbins, bb->eventbinnum, bb->eventbinsize, tl[j].age - tl[j-1].age, 
					  tl[j].age, 
					  tl[j].eventnode->type,
					  tl[j].from, tl[j].to, tl[j].lineages, world->numpop, world);
	      }
	    
	  }
      }
    //
    if(!world->in_last_chain)
      {
        for (j = world->numpop; j < world->numpop2; j++)
	  if ((atl->mindex[j] < 1.0) && (world->param0[j] > 0.))
            {
	      atl->mindex[j] = world->options->minmigsumstat;
            }
      }
#ifdef ALTIVEC
    load_float_floatvec(atl->vdata,atl->data,world->numpop2+2*world->numpop);
#endif
}

/* this function is obsolete but contains some parts that should be rescued
///
/// old archive_timelist [keep because of longsum stuff and also in case
// the memory reduction does not work
void
archive_timelist_OLD (tarchive_fmt * atl, vtlist * tl, long T, long np,
                  world_fmt * world)
{
    long j, i;
    MYREAL t;
    MYREAL line;
    MYREAL sumlines=0;
    long msize = 0;

#ifdef LONGSUM    
    longsum_fmt *longsum;
#endif 

    
    mighist_fmt *aa = NULL;
    mighistloci_fmt *bb = NULL;

    for (i = 0; i < np; i++)
    {
        line = (MYREAL) tl[0].lineages[i];
        atl->km[i] = line * tl[0].age;
        atl->kt[i] = line * (line - 1) * tl[0].age;
        atl->p[i] = (MYREAL) whichp (tl[0].from, tl[0].to, i);
    }
    memset (atl->mindex + np, 0, sizeof (MYREAL) * np * (np - 1));
    if (tl[0].from != tl[0].to)
        atl->mindex[mm2m (tl[0].from, tl[0].to, np)] += 1;
#ifdef LONGSUM
    
    atl->longsumlen = T;
    atl->longsum = (longsum_fmt *) myrealloc(atl->longsum,sizeof(longsum_fmt)*T);
    atl->longsum[0].lineages = (long *) mycalloc(T,sizeof(long)*np);
    atl->longsum[0].lineages2 = (long *) mycalloc(T,sizeof(long)*np);
    longsum = atl->longsum;
    memcpy(longsum[0].lineages,tl[0].lineages,sizeof(long)*np);
    memcpy(longsum[0].lineages2,tl[0].lineages,sizeof(long)*np);
    for(i=0;i<np;i++)
        longsum[0].lineages2[i] *= longsum[0].lineages2[i] - 1;
    longsum[0].fromto = mm2m(tl[0].from, tl[0].to,np);
    longsum[0].to = tl[0].to;
    longsum[0].eventtime = tl[0].age;
    longsum[0].interval = tl[0].age;
    warning("needs updating to deal with 'd' world.c:2654");
    longsum[0].eventtype = tl[0].from != tl[0].to ? 'm' : 'c';
#endif 
    
    if (world->in_last_chain && world->options->mighist)
    {
      //        aa = &(world->mighistloci[world->locus].mighist[world->mighistloci[world->locus].mighistnum]);
      aa = &(world->mighistloci[world->locus].mighist[world->mighistloci[world->locus].mighistnum]);
      bb = &(world->mighistloci[world->locus]);
      //	msize=aa->migeventsize=0;  // allocation is not necessary anymore but check first
      //  if(aa->allocsize <= msize + 1)
      //  {
      //      aa->allocsize += DEFAULTALLOCSIZE;
      //      aa->migevents = myrealloc(aa->migevents, sizeof(migevent_fmt) * aa->allocsize);
      //  }
        //    printf("(%li) migeventsize=%li, allocsize=%li\n",world->mighistloci[world->locus].mighistnum,
        //    aa->migeventsize,aa->allocsize);
	if(world->options->mighist_all || (tl[0].from != tl[0].to))
	  {
	    sumlines =0.;
	    for(i=0; i< np; i++)
	      sumlines += tl[0].lineages[i];
	    aa->migevents[msize].age = tl[0].age;
	    aa->migevents[msize].from = (MYREAL) tl[0].from;
	    aa->migevents[msize].to = (MYREAL) tl[0].to;
	    aa->migevents[msize].sumlines = sumlines;

	    // calculate event values
	    calculate_event_values(bb->migeventbins, bb->migeventbinnum, bb->eventbinsize, 
				      tl[0].age, 
				      tl[0].age, tl[0].from, tl[0].to, tl[0].lineages, world->numpop);
     
	    // calculate skyline values
	    if(world->options->skyline)
	      {
		calculate_expected_values(bb->eventbins, bb->eventbinnum, bb->eventbinsize, tl[0].age, 
					  tl[0].age, tl[0].from, tl[0].to, tl[0].lineages, world->numpop, world);
	      }
	    aa->migeventsize++;
	    msize = aa->migeventsize;
	  }
        //        }
    }
    for (j = 1; j < T; j++)
    {
        t = tl[j].age - tl[j - 1].age;
        for (i = 0; i < np; i++)
        {
            line = (MYREAL) tl[j].lineages[i];
            atl->km[i] += line * t;
            atl->kt[i] += line * (line - 1.) * t;
            atl->p[i] += whichp (tl[j].from, tl[j].to, i);
        }
        if (tl[j].from != tl[j].to)
            atl->mindex[mm2m (tl[j].from, tl[j].to, np)] += 1;
#ifdef LONGSUM
        
        atl->longsum[j].lineages = (long *) mycalloc(1,sizeof(long)*np);
        atl->longsum[j].lineages2 = (long *) mycalloc(1,sizeof(long)*np);
        memcpy(longsum[j].lineages,tl[j].lineages,sizeof(long)*np);
        memcpy(longsum[j].lineages2,tl[j].lineages,sizeof(long)*np);
        for(i=0;i<np;i++)
            longsum[j].lineages2[i] *= longsum[j].lineages2[i] - 1;
        
        longsum[j].fromto = mm2m(tl[j].from,tl[j].to,np);
        longsum[j].to = tl[j].to;
        longsum[j].eventtime = tl[j].age;
        longsum[j].interval = t;
	warning("needs updating to deal with 'd' world.c:2654");
        longsum[j].eventtype = tl[j].from != tl[j].to ? 'm' : 'c';
#endif 
        
        if (world->in_last_chain && world->options->mighist)
        {
            //            if(tl[j].from != tl[j].to)
            //            {
            if(aa->allocsize <= aa->migeventsize + 1 )
            {
                aa->allocsize += DEFAULTALLOCSIZE;
                aa->migevents = myrealloc(aa->migevents, sizeof(migevent_fmt) * aa->allocsize);
            }
            // printf("(%li) %li> migeventsize=%li, allocsize=%li\n",world->mighistloci[world->locus].mighistnum, j,aa->migeventsize,aa->allocsize);
            sumlines =0.;
            for(i=0; i< np; i++)
                sumlines += tl[j].lineages[i];
            aa->migevents[msize].age = tl[j].age;
            aa->migevents[msize].from = (int) tl[j].from;
            aa->migevents[msize].to = (int) tl[j].to;
            aa->migevents[msize].sumlines = sumlines;
	    
	    // calculate event values
	    calculate_event_values(bb->migeventbins, bb->migeventbinnum, bb->eventbinsize, 
				      tl[j].age - tl[j-1].age, 
				      tl[j].age, tl[j].from, tl[j].to, tl[j].lineages, world->numpop);
     
	    // calculate skyline values
	    if(world->options->skyline)
	      {
		calculate_expected_values(bb->eventbins, bb->eventbinnum, bb->eventbinsize, tl[j].age - tl[j-1].age, 
					  tl[j].age, tl[j].from, tl[j].to, tl[j].lineages, world->numpop, world);
	      }
	    
            aa->migeventsize++;
            msize = aa->migeventsize;
            //            }
        }
    }
    if(!world->in_last_chain)
    {
        for (j = world->numpop; j < world->numpop2; j++)
            if ((atl->mindex[j] < 1.0) && (world->param0[j] > 0))
            {
                //DEBUG		printf("changing mindex[%li]=%f to %f\n",j,atl->mindex[j],world->options->minmigsumstat);
                atl->mindex[j] = world->options->minmigsumstat;
            }
    }
#ifdef ALTIVEC
        load_float_floatvec(atl->vdata,atl->data,world->numpop2+2*world->numpop);
#endif
    }

*/

void
copy_timelist (tarchive_fmt * from, tarchive_fmt * to, long np)
{
  memcpy (to->data, from->data, sizeof (MYREAL) * (size_t) (3 * np + np * (np - 1)));
    //memcpy(to->p, from->p, sizeof(long) * np);
    //memcpy(to->mindex, from->mindex, sizeof(MYREAL) * np * (np-1));
    //memcpy(to->km, from->km, sizeof(MYREAL) * np);
    //memcpy(to->kt, from->kt, sizeof(MYREAL) * np);
}

long
whichp (long from, long to, long pop)
{
    if (from == to)
    {
        if (from == pop)
            return 1;
    }
    return 0;
}
#ifdef LONGSUM
void print_fluctuate_header(world_fmt *world)
{
    FPRINTF(world->outfile,"\n\n\n");
    FPRINTF(world->outfile,"============================================================\n");
    FPRINTF (world->outfile, "Population size Rate changes in the past [today's rate is 1] \n");
    FPRINTF(world->outfile,"============================================================\n\n");
    FPRINTF(world->outfile,"Population      Loc.  Rate_2     Time_2       Rate_3     Time_3\n");
    FPRINTF(world->outfile,"--------------- ----  --------------------    -------------------- \n");
}

void print_fluctuate_results(world_fmt * world, long locus, long rep,
                             long pop)
{
    FPRINTF(world->outfile," %10.7f %10.7f   %10.7f %10.7f\n",
            world->atl[rep][locus].param[world->numpop2 + pop * 3+1],
            world->flucrates[world->numpop * 3 + pop * 3+1],
            world->atl[rep][locus].param[world->numpop2 + pop * 3+2],
            world->flucrates[world->numpop * 3 + pop * 3+2]);
}

void print_fluctuate(world_fmt **universe, option_fmt *options, data_fmt *data)
{
    long skipped = 0;
    long pop;
    long locus;
    world_fmt *world=EARTH;
    long maxrep = world->options->replicate ?
        (world->options->replicatenum > 0 ?
         world->options->replicatenum + 1 : world->options->lchains + 1) : 1;
    long rep;
    print_fluctuate_header(world);
    for(pop=0; pop < world->numpop; pop++)
    {
        print_popstring(pop, world, options, data);
        for (locus = 0; locus < world->loci; locus++)
        {
            if (world->data->skiploci[locus])
            {
                skipped++;
                continue;
            }
            for(rep=0; rep < maxrep; rep++)
            {
                print_replicate(world, maxrep, rep,locus);
                print_fluctuate_results(world, locus, rep, pop);
            }
        }
        if (world->loci - skipped > 1)
        {
            FPRINTF (world->outfile, "                All ");
            print_fluctuate_results(world, world->loci, 0, pop);
        }
        
    }
    FPRINTF(world->outfile,"\n");
}
#endif /*LONGSUM*/

void
print_results (world_fmt ** universe, option_fmt * options, data_fmt * data)
{
    long pop;
    FILE *outfile;
    world_fmt *world=EARTH;
    worldoption_fmt *wopt = world->options;
    char sch[10], lch[10], cva[50];
    long rep = world->loci > 1 ? 0 : (wopt->replicate ? world->repstop : 0);
    outfile = world->outfile;
    if (options->schains == 1)
        strcpy (sch, "chain");
    else
        strcpy (sch, "chains");
    if (options->lchains == 1)
        strcpy (lch, "chain");
    else
        strcpy (lch, "chains");
    FPRINTF (outfile, "\n\n");
    PAGEFEED;
    print_result_header ("MCMC estimates ", world);
    for (pop = 0; pop < world->numpop; pop++)
    {
        print_result_population (pop, world, options, data);
    }
    FPRINTF (outfile,
             "\nComments:\n The x is 1, 2, or 4 for mtDNA, haploid, or diploid data, respectively\nThere were %li short %s (%li used trees out of ",
             options->schains, sch, options->ssteps);
    FPRINTF (outfile, "sampled %li)\n", options->sincrement * options->ssteps);
    FPRINTF (outfile,
             "  and %li long %s (%li used trees out of sampled %li)\n",
             options->lchains, lch, options->lsteps,
             options->lincrement * options->lsteps);
    if (wopt->heating)
    {
        if(options->adaptiveheat)
        {
            FPRINTF (outfile,
                     "  Adaptive heating with %li chains was active\ncheck Log file (if present) for average temperatures\n",options->heated_chains);
            
            //Average last chain temp: 1.0",
            //                   options->heated_chains);
            //            for(pop=1;pop<options->heated_chains;pop++)
            //               FPRINTF(outfile, ", %f", universe[pop]->averageheat);
            //				FPRINTF (outfile,"\n\n");
        }
        else
            FPRINTF (outfile, "  Static heating with %li chains was active\n",options->heated_chains);
    }
    if (options->gamma)
    {
        if (world->atl[rep][world->loci].param[world->numpop2] < 10e-9)
            strcpy (cva, "0");
        else
            sprintf (cva, "%f",
                     sqrt (1. /
                           world->atl[rep][world->loci].param[world->numpop2]));
        FPRINTF (outfile,
                 "With shape parameter Alpha=%g ([1/CV(mu)]^2; CV(mu)=%s)\n",
                 world->atl[rep][world->loci].param[world->numpop2],
                 cva);
    }
    if (world->options->replicate)
    {
        if (world->repkind == MULTIPLECHAIN)
            FPRINTF (outfile, "  COMBINATION OF ALL LONG CHAINS\n");
        else
            FPRINTF (outfile, "  COMBINATION OF %li MULTIPLE RUNS)\n",
                     world->options->replicatenum);
    }
    if (world->atl[rep][world->loci].normd > LOCI_NORM)
        FPRINTF (outfile,
                 "  [Last maximization needed %li cycles of maximal %i,\n  Norm(first derivatives)=%f (Normal stopping criteria is < %f)]\n\n\n",
                 world->atl[rep][world->loci].trials,
                 NTRIALS, world->atl[rep][world->loci].normd, LOCI_NORM);
    FPRINTF (outfile,"\n\n");
}

void
print_fst (world_fmt * world, option_fmt * options, data_fmt * data,
           MYREAL **fstparam)
{
  (void) fstparam;
    long pop;
    long loci = world->loci;
    FILE *outfile = world->outfile;
    if (loci < 40)
    {
        PAGEFEED;
    }
    if (options->fsttype == 'T')
        print_result_header
            ("FST estimates (Thetas are variable)\n   [Only used as start values for MCMC run!!!!!]",
             world);
    else
        print_result_header
            ("FST estimates (Migration rates are variable)\n  [Only used as start values for MCMC run!!!]",
             world);
    for (pop = 0; pop < world->numpop; pop++)
    {
        print_result_fst (pop, world, data);
    }
    FPRINTF (world->outfile, "\nComments:\n");
    FPRINTF (world->outfile, "(-) can not be estimated\n");
    FPRINTF (world->outfile, "(0/0 or x/0) Divisions by zero\n");
    FPRINTF (world->outfile,
             "Maynard Smith, J. 1970. Population size, polymorphism, and the rate of\n");
    FPRINTF (world->outfile,
             "    non-Darwinian evolution. American Naturalist 104:231-237\n");
    FPRINTF (world->outfile,
             "Nei, M., and M. Feldman 1972. Identity of genes by descent within and\n");
    FPRINTF (world->outfile,
             "    between populations under mutation and migration pressures.\n");
    FPRINTF (world->outfile, "    Theoretical Population Biology 3: 460-465\n");
}

void
prepare_print_nu (MYREAL nu, char *str)
{
    if (nu <= -999)
        sprintf (str, "-");
    else
        sprintf (str, "% 12.6f", nu);
}

void
prepare_print_nm (MYREAL nm, MYREAL nmu, char *strllike)
{
    if ((fabs (nmu) > 10e-20) && (fabs (nm) > 10e-20))
        sprintf (strllike, "% 10.5f", nm * nmu);
    else
    {
        if ((fabs (nmu) < 10e-20) && (fabs (nm) < 10e-20))
        {
            sprintf (strllike, "0/0");
        }
        else
        {
            if ((fabs (nmu) < 10e-20) && (fabs (nm) > 10e-20))
                sprintf (strllike, "% 10.5f/0", nm);
            else
            {
                if ((fabs (nmu) > 10e-20) && (fabs (nm) < 10e-20))
                    sprintf (strllike, "0");
            }
        }
    }
}

void
print_menu_coalnodes (FILE * file, world_fmt * world, long G, long rep)
{
    long g, pop, minp = world->sumtips, maxp = 0;
    long maxp1;
    char ss[10];
    long **contribution;
    char *buffer;
    //char *bufptr;
    long allocbufsize=LINESIZE;
    long bufsize=0;
    long nodenum = world->sumtips;
    tarchive_fmt *tl = world->atl[rep][world->locus].tl;
    char * tmp = (char *) mycalloc(allocbufsize,sizeof(char));
    if (world->options->verbose)
    {
        buffer = (char *) mycalloc(allocbufsize,sizeof(char));
        contribution = (long **) mycalloc (1, sizeof (long *) * (size_t) nodenum);
        contribution[0] = (long *) mycalloc (1, sizeof (long) * (size_t) (nodenum * world->numpop));
        for (pop = 1; pop < nodenum; pop++)
        {
            contribution[pop] = contribution[0] + pop * world->numpop;
        }
        for (g = 0; g < G; g++)
        {
            for (pop = 0; pop < world->numpop; pop++)
            {
                contribution[pop][(long) tl[g].p[pop]] += tl[g].copies;
            }
        }
        for (g = 0; g < nodenum; g++)
        {
            for (pop = 0; pop < world->numpop; pop++)
            {
                if (maxp < g && contribution[pop][g] > 0)
                    maxp = g;
                if (minp > g && contribution[pop][g] > 0)
                    minp = g;
            }
        }
        sprintf (tmp, "           Coalescent nodes: ");
	add_to_buffer(tmp, &bufsize, &buffer, &allocbufsize);
	maxp1 = maxp + 1;
        for (g = minp; g < maxp1; g++)
        {
            sprintf (tmp, "%2li ", g);
	    add_to_buffer(tmp, &bufsize, &buffer, &allocbufsize);
		
        }
        sprintf (tmp, "\n");
	add_to_buffer(tmp, &bufsize, &buffer, &allocbufsize);
	for (pop = 0; pop < world->numpop; pop++)
        {
            sprintf (tmp, "             population %3li: ", pop);
	    add_to_buffer(tmp, &bufsize, &buffer, &allocbufsize);
            for (g = minp; g < maxp1; g++)
            {
                if (contribution[pop][g] == 0)
                {
                    strcpy (ss, "-");
                }
                else
                {
                    if (contribution[pop][g] >= 100)
                        strcpy (ss, "*");
                    else
                        sprintf (ss, "%-6li", contribution[pop][g]);
                }
		add_to_buffer(ss, &bufsize, &buffer, &allocbufsize);
            }
	    sprintf (tmp,"\n");
	    add_to_buffer(ss, &bufsize, &buffer, &allocbufsize);
        }
        LARGEFPRINTF(file, bufsize, "%s",buffer);
        myfree(contribution[0]);
        myfree(contribution);
        myfree(buffer);
	myfree(tmp);
    }
}


///
/// print the covariance matrix, inverse of the approximate second derivative matrix 
void
print_cov (world_fmt * world, long numpop, long loci, MYREAL ***cov)
{
  (void) numpop;
  (void) loci;
    FILE *outfile = world->outfile;
    long locus, skipped = 0;
    MYREAL *corr;
    corr = (MYREAL *) mycalloc ((world->numpop2 + 1), sizeof (MYREAL));
    FPRINTF (outfile, "\n\n");
    FPRINTF (outfile,
             "--------------------------------------------------------------------\n");
    FPRINTF (outfile, "MCMC estimation:\n");
    FPRINTF (outfile,
             "Covariance matrix(*)                         Correlation matrix\n");
    FPRINTF (outfile,
             "-------------------------------------------  -----------------------\n");
    for (locus = 0; locus < world->loci; locus++)
    {
        if (world->data->skiploci[locus])
        {
            skipped++;
            continue;
        }
        FPRINTF (outfile, "Locus %li:\n", locus + 1);
	print_cov_table (outfile, locus, world, cov, corr, 0);
	FPRINTF (outfile, "\n");
    }
    FPRINTF (outfile, "\n\n");
    if (world->loci - skipped > 1)
    {
        if (world->options->gamma)
        {
            FPRINTF (outfile, "Over all loci\n");
            FPRINTF (outfile,
                     "------------------------------------------------------  ");
            FPRINTF (outfile, "----------------------------\n");
            print_cov_table (outfile, locus, world, cov, corr, 1);
        }
        else
        {
            FPRINTF (outfile, "Over all loci\n");
            FPRINTF (outfile,
                     "-------------------------------------------  -----------------------\n");
            print_cov_table (outfile, locus, world, cov, corr, 0);
        }
    }
    myfree(corr);
}

/// print the table content of the covariance matrix
void
print_cov_table (FILE * outfile, long locus, world_fmt * world,
                 MYREAL ***cov, MYREAL *corr, long addvar)
{
  //#ifdef MYREAL == float
  //const MYREAL eps = FLT_EPSILON;
  //#else
  const MYREAL eps = DBL_EPSILON;
  //#endif
    long i, j;
    MYREAL denom, temp1, temp2;
    for (i = 0; i < world->numpop2 + addvar; i++)
    {
        for (j = 0; j < world->numpop2 + addvar; j++)
        {
            temp1 = fabs (cov[locus][i][i]);
            if (temp1 < eps)
                temp1 = eps;
            temp2 = fabs (cov[locus][j][j]);
            if (temp2 < eps)
                temp2 = eps;
            denom = 0.5 * (LOG (temp1) + LOG (temp2));
            if ((EXP (denom)) == 0.0)
                corr[j] = 9.99;
            else
                corr[j] = cov[locus][i][j] / EXP (denom);
        }
        for (j = 0; j < world->numpop2 + addvar; j++)
        {
            FPRINTF (outfile, "% 10.4f ", cov[locus][i][j]);
        }
        for (j = 0; j < world->numpop2 + addvar; j++)
        {
            if (corr[j] < 9.)
                FPRINTF (outfile, "% 4.2f ", corr[j]);
            else
                FPRINTF (outfile, "  --  ");
        }
        FPRINTF (outfile, "\n");
    }
}
///
/// print correlation matrix only using the covariance matrix
/*
void
print_cov2 (world_fmt * world, long numpop, long loci, MYREAL ***cov)
{
    FILE *outfile = world->outfile;
    long locus, skipped = 0;
    MYREAL *corr;
    long i;
    long frompop, topop;
    char head[20];
    corr = (MYREAL *) mycalloc (1, sizeof (MYREAL) * (world->numpop2 + 1));
    FPRINTF (outfile, "\n\n");
    FPRINTF (outfile,
             "--------------------------------------------------------------------\n");
    FPRINTF (outfile,
             "Correlation matrix of parameters\n");
    FPRINTF (outfile,
             "--------------------------------------------------------------------\n");
    for(i=0;i<world->numpop2;i++)
      {
	if (i<world->numpop)
	  {
	    sprintf(head,"Q%li",i+1);
	  }
	else
	  {
	    m2mm(i, numpop, &frompop, &topop);
	    sprintf(head,"M%li%li",frompop+1,topop+1);
	  }
	FPRINTF(outfile,"%5.5s ",head);
      }
    FPRINTF(outfile,"\n");
    for (locus = 0; locus < world->loci; locus++)
    {
        if (world->data->skiploci[locus])
        {
            skipped++;
            continue;
        }
        FPRINTF (outfile, "Locus %li:\n", locus + 1);
	print_cov_table2 (outfile, locus, world, cov, corr, world->bayes->mu);
	FPRINTF (outfile, "\n");
    }
    FPRINTF (outfile, "\n\n");
    if (world->loci - skipped > 1)
    {
      FPRINTF (outfile, "Over all loci\n");
      FPRINTF (outfile,
	       "--------------------------------------------------------------------\n");
      print_cov_table2 (outfile, locus, world, cov, corr, world->bayes->mu);
      FPRINTF (outfile,
	       "--------------------------------------------------------------------\n");
    }
    myfree(corr);
}

/// print the table content of the covariance matrix
void
print_cov_table2 (FILE * outfile, long locus, world_fmt * world,
                 MYREAL ***cov, MYREAL *corr, long addvar)
{
  const MYREAL eps = DBL_EPSILON;
  long i, j;
  MYREAL denom, temp1, temp2;
  for (i = 0; i < world->numpop2 + addvar; i++)
    {
      for (j = 0; j < world->numpop2 + addvar; j++)
        {
	  temp1 = fabs (cov[locus][i][i]);
	  if (temp1 < eps)
	    temp1 = eps;
	  temp2 = fabs (cov[locus][j][j]);
	  if (temp2 < eps)
	    temp2 = eps;
	  denom = 0.5 * (LOG (temp1) + LOG (temp2));
	  if ((EXP (denom)) == 0.0 || cov[locus][i][i] == 0.0 || cov[locus][j][j]==0.0)
	    corr[j] = 9.99;
	  else
	    corr[j] = cov[locus][i][j] / EXP (denom);
        }
        for (j = 0; j < world->numpop2 + addvar; j++)
        {
	  if (corr[j] < 9.)
	    FPRINTF (outfile, "% 4.2f ", corr[j]);
	  else
	    FPRINTF (outfile, "  --  ");
        }
        FPRINTF (outfile, "\n");
    }
}
*/

/// free the conditional likelihood vector in every node recursively
void
free_seqx (node * p, world_fmt * world)
{
    long j;
    if (p->type != 't')
    {
        if (p->next->back != NULL)
        {
            free_seqx (crawlback (p->next), world);
        }
        if (p->next->next->back != NULL)
        {
            free_seqx (crawlback (p->next->next), world);
        }
    }
    if (p->type == 'r')
        return;

    mutationmodel_fmt *s;
    long sublocus;
    long locus = world->locus;
    const long sublocistart = world->sublocistarts[locus];
    const long sublociend   = world->sublocistarts[locus+1];

    for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
      {
        s = &world->mutationmodels[sublocus];
	const long xs = sublocus - sublocistart;
	for (j = 0; j < s->numpatterns; j++)
	  {
	    myfree(p->x[xs].s[j]);
	  }
      }
}

/// free every cond likelihood for alleles in every node
void
free_x (node * p, world_fmt *world)
{
    if (p->type != 't')
    {
        if (p->next->back != NULL)
        {
	  free_x (crawlback (p->next), world);
        }
        if (p->next->next->back != NULL)
        {
	  free_x (crawlback (p->next->next), world);
        }
    }
    if (p->type == 'r')
        return;


 
   long sublocus;
    long locus = world->locus;
    const long sublocistart = world->sublocistarts[locus];
    const long sublociend   = world->sublocistarts[locus+1];

    for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
      {
	const long xs = sublocus - sublocistart;
	myfree(p->x[xs].a);
      }
}


/// print the coefficient of variation when using gamma distributed mutation rates among loci
void
print_CV (world_fmt * world)
{
    long i;
    long elem;
    char temp[100];
    FPRINTF (world->outfile, "\n\nVariance and coefficient of variance\n");
    if (world->loci - world->skipped > 1)
    {
        elem = world->numpop2 + (world->options->gamma ? 1 : 0);
        FPRINTF (world->outfile, "PARAM@ ");
        for (i = 0; i < elem; i++)
        {
            FPRINTF (world->outfile, "%20.20g ",
                     world->atl[0][world->loci].param[i]);
        }
        FPRINTF (world->outfile, "\nVAR@   ");
        for (i = 0; i < elem; i++)
        {
            FPRINTF (world->outfile, "%20.20g ", world->cov[world->loci][i][i]);
        }
        FPRINTF (world->outfile, "\nCV@    ");
        for (i = 0; i < elem; i++)
        {
            strcpy (temp, "-");
            if (world->cov[world->loci][i][i] >= 0)
            {
                sprintf (temp, "%20.20g",
                         (sqrt (world->cov[world->loci][i][i])) /
                         world->atl[0][world->loci].param[i]);
            }
            FPRINTF (world->outfile, "%-s ", temp);
        }
        FPRINTF (world->outfile, "\n");
    }
    else
    {    /* one locus */
        elem = world->numpop2;
        FPRINTF (world->outfile, "PARAM@ ");
        for (i = 0; i < elem; i++)
        {
            FPRINTF (world->outfile, "%20.20g ", world->atl[0][1].param[i]);
        }
        FPRINTF (world->outfile, "\nVAR@   ");
        for (i = 0; i < elem; i++)
        {
            FPRINTF (world->outfile, "%20.20g ", world->cov[0][i][i]);
        }
        FPRINTF (world->outfile, "\nCV@    ");
        for (i = 0; i < elem; i++)
        {
            strcpy (temp, "-");
            if (world->cov[0][i][i] >= 0)
            {
                sprintf (temp, "%20.20g",
                         (sqrt (world->cov[0][i][i])) /
                         world->atl[0][1].param[i]);
            }
            FPRINTF (world->outfile, "%-s ", temp);
        }
        FPRINTF (world->outfile, "\n");
    }
}

///
/// print progress report during run
void
print_progress(worldoption_fmt * options, world_fmt * world, long rep,
               long visited, long accepted)
{
    if (options->progress && !options->bayes_infer)
    {
        print_menu_coalnodes(stdout, world, world->G + 1, rep);
        print_menu_accratio(accepted, visited, world);
        if (options->writelog)
        {
            print_menu_coalnodes(options->logfile, world, world->G + 1, rep);
        }
        if (!world->start && world->chains > 1 && world->options->gelman)
            print_gelmanr(world->convergence->gelmanmeanRall, world->convergence->gelmanmaxRall);
    }
}



/// print the MLE parameter table
void
print_param (char ** file, long *bufsize, long *allocbufsize, boolean usem, world_fmt *world, long nn, char spacer[])
{
    long i, j;   //, fpos;
    long tt = nn;
    long counter;
    long spacersize = (long) strlen(spacer)+1;
    long newbufsize =  *bufsize + spacersize + 25 + nn * ((spacersize + 4 + 13) + nn * 33); 
    MYREAL *param = world->param0;
    if(newbufsize > *allocbufsize)
      {
	*file = (char *) myrealloc(*file, sizeof(char)* (size_t) newbufsize);
	*allocbufsize = newbufsize;
      }
    *bufsize += sprintf (*file + *bufsize, "%sPop. Theta    %s\n", spacer, usem ? "M" : "Theta*M");
    for (i = 0; i < nn; i++)
    {
        //FPRINTF (file, " ");
        counter = 0;
        *bufsize += sprintf (*file + *bufsize, "%s%3li % 7.5f", spacer, i + 1, param[i]);
        for (j = 0; j < nn; j++)
        {
            if (i != j)
            {
                if (usem)
                    *bufsize += sprintf (*file + *bufsize, "% 7.1f", param[tt++]);
                else
                    *bufsize += sprintf (*file + *bufsize, "% 7.5f", param[i] * param[tt++]);
            }
            else
                *bufsize += sprintf (*file + *bufsize, " ------");
            if (counter++ > 10)
            {
                counter = 0;
                *bufsize += sprintf (*file + *bufsize, "\n%s         ", spacer);
            }
        }
        *bufsize += sprintf (*file + *bufsize, " \n"); //, spacer);
    }
    *bufsize += sprintf (*file + *bufsize, " \n");
#ifdef LONGSUM
    // this needs to be fixed when I revisit LONGSUM -- not on top ten list
    bufsize += sprintf (file + bufsize,"\nRates (at specific times):\n");
    for(i=0;i<nn; i++ )
    {
        j = i*3; //because we use only 3 different rates per pop
        bufsize += sprintf (file + bufsize,"Pop %li: %g (%g) %g (%g) %g (%g)\n", i+1,
                            world->flucrates[j],world->flucrates[j+nn*3],
                            world->flucrates[j+1], world->flucrates[j+1+nn*3],
                            world->flucrates[j+2], world->flucrates[j+2+nn*3]);
    }
#endif
}


/*
 * print header
 * ===========================================================================
 * === Titletext
 * ===========================================================================
 * === Population     Loci  Ln(L)   Theta    xNm [xNe mu] xx     xx     xx
 * xx     xx     xx xx.... -------------- ---- -------- --------
 * ----------------------------------------
 */
/// print the MLE table header
void
print_result_header (char *titletext, world_fmt * world)
{
    long p1, zz;
    char dline[] =
        "==============================================================================";
    char sline[] =
        "-------------- ---- -------- -------- ----------------------------------------";
    FPRINTF (world->outfile, "%s\n", dline);
    FPRINTF (world->outfile, "%s\n", titletext);
    FPRINTF (world->outfile, "%s\n", dline);
    if (!world->options->usem)
    {
        FPRINTF (world->outfile,
                 "Population [x] Loc.  Ln(L)   Theta    xNm [+=receiving population]\n");
        FPRINTF (world->outfile, "                             [xN mu] ");
    }
    else
    {
        FPRINTF (world->outfile,
                 "Population [x] Loc.  Ln(L)   Theta    M [m/mu] [+=receiving population]  \n");
        FPRINTF (world->outfile, "                             [xN mu]   ");
    }
    zz = 38;
    for (p1 = 0; p1 < world->numpop; p1++)
    {
        zz = zz + 7;
        if (zz > 79)
        {
            zz = 38;
            FPRINTF (world->outfile,
                     "\n                                      ");
        }
        FPRINTF (world->outfile, "%2li,+    ", p1 + 1);
    }
    FPRINTF (world->outfile, "\n%s\n", sline);
}


/// print the string for a population
void print_popstring(long pop, world_fmt *world, option_fmt *options, data_fmt *data)
{
    char popstring[LINESIZE];
    if (options->readsum)
    {
        sprintf (popstring, "%2li: ", pop + 1);
    }
    else
    {
        sprintf (popstring, "%2li: %s", pop + 1, data->popnames[options->newpops[pop]-1]);
    }
    FPRINTF (world->outfile, "%-14.14s ", popstring);
}


/// print the replicate number
void print_replicate(world_fmt *world, long maxrep, long rep, long locus)
{
    char repstring[LINESIZE];
    sprintf (repstring, "%2li", rep + 1);
    FPRINTF (world->outfile, "%s%2li%2s ", locus == 0
             && rep == 0 ? "" : "               ", locus + 1,
             maxrep > 1 ? (rep ==
                           maxrep - 1 ? " A" : repstring) : "  ");
}

/// print the MLE table content for each population
void
print_result_population (long pop, world_fmt * world,
                         option_fmt * options, data_fmt * data)
{
    long skipped = 0, locus;
    long maxrep = world->options->replicate ?
        (world->options->replicatenum > 0 ?
         world->options->replicatenum + 1 : world->options->lchains + 1) : 1;
    long rep;
    print_popstring(pop, world, options, data);
    for (locus = 0; locus < world->loci; locus++)
    {
        if (world->data->skiploci[locus])
        {
            skipped++;
            continue;
        }
        for (rep = 0; rep < maxrep; rep++)
        {
            print_replicate(world, maxrep, rep, locus);
            FPRINTF (world->outfile, "% 8.3f ",
                     world->atl[rep][locus].param_like);
            print_result_param (world->outfile, world->atl[rep][locus].param,
                                world->numpop, pop, world->options->usem);
        }
    }
    if (world->loci - skipped > 1)
    {
        FPRINTF (world->outfile, "                All ");
        //locus is exactly world->loci
        // re is always one because we have only replication of single locus chains
        FPRINTF (world->outfile, "% 8.3f ", world->atl[0][locus].param_like);
        print_result_param (world->outfile, world->atl[0][locus].param,
                            world->numpop, pop, world->options->usem);
    }
    /* FPRINTF(world->outfile,"%s\n",sline);     */
}


/// print the parameter for the MLE parameter printout
void
print_result_param (FILE * file, MYREAL *param, long numpop, long pop,
                    boolean usem)
{
    long i;
    long linelen = 0;
    long msta = mstart (pop, numpop);
    long msto = mend (pop, numpop);
    MYREAL tmp = 0;
    if (param[pop] <= SICK_VALUE)
        FPRINTF (file, "     -    ");
    else
    {
        if (param[pop] < 0.0001 && param[pop] > 0)
            FPRINTF (file, "%3.2e ", param[pop]);
        else
            FPRINTF (file, "%8.5f ", param[pop]);
    }
    for (i = msta; i < msto; i++)
    {
        if (pop == i - msta)
        {
            FPRINTF (file, "------- ");
            linelen++;
        }
        if (linelen > 4)
        {
            FPRINTF (file, "\n                                      ");
            linelen = 0;
        }
        linelen++;
        if ((param[i] <= SICK_VALUE) || (param[pop] <= SICK_VALUE))
            FPRINTF (file, "    -    ");
        else
        {
            if (usem)
            {
                tmp = param[i];
                if (tmp < 0.00001)
                    FPRINTF (file, " 0.0000 ");
                else
                    FPRINTF (file, "%7.4f ", tmp);
            }
            else
            {
                tmp = param[pop] * param[i];
                if (tmp < 0.00001)
                    FPRINTF (file, " 0.0000 ");
                else
                    FPRINTF (file, "%7.4f ", tmp);
            }
        }
    }
    if (pop == numpop - 1)
        FPRINTF (file, "-------");
    FPRINTF (file, "\n");
}


/// print the results for the FST calculations (used to start program)
void
print_result_fst (long pop, world_fmt * world, data_fmt * data)
{
    
    char popstring[LINESIZE];
    long skipped = 0, locus;
    sprintf (popstring, "%2li: %s", pop + 1, data->popnames[pop]);
    FPRINTF (world->outfile, "%14.14s ", popstring);
    for (locus = 0; locus < world->loci; locus++)
    {
        if (world->data->skiploci[locus])
        {
            skipped++;
            continue;
        }
        FPRINTF (world->outfile, "%s%4li ",
                 locus == 1 ? "" : "               ", locus + 1);
        FPRINTF (world->outfile, "   -    ");
        print_result_param (world->outfile, world->fstparam[locus],
                            world->numpop, pop, world->options->usem);
    }
    if (world->loci - skipped > 1)
    {
        FPRINTF (world->outfile, "                All ");
        FPRINTF (world->outfile, "   -    ");
        print_result_param (world->outfile, world->fstparam[locus],
                            world->numpop, pop, world->options->usem);
    }
    /* FPRINTF(world->outfile,"%s\n",sline);     */
}


/// klone cleanup -- tree related
void klone_tree_setup(world_fmt *kopie, option_fmt *options)
{
  create_treetimelist (kopie, &kopie->treetimes);
  fix_times (kopie, options);
  set_tree_dirty (kopie->root);
  first_smooth(kopie,kopie->locus);
  kopie->likelihood[0] = treelikelihood (kopie);
}

///
/// Clones world structure, used when heating is used
void
klone (world_fmt * original, world_fmt * kopie,
       option_fmt * options, data_fmt * data, 
       MYREAL temperature)
{
  long i; 
#ifdef UEP
  long j;
#endif
  long locus=0;
  const long np = original->numpop;
  const long npp = original->numpop2 + original->bayes->mu + original->species_model_size * 2 + original->grownum;
  kopie->repkind = SINGLECHAIN;
  kopie->loci = original->loci;
  kopie->skipped = original->skipped;
  kopie->locus = original->locus;
  kopie->replicate = original->replicate;
  kopie->numpop = original->numpop;
  kopie->numpop2 = original->numpop2;
  kopie->sumtips = original->sumtips;
  memcpy (kopie->mstart, original->mstart, sizeof (int) * (size_t) np);
  memcpy (kopie->mend, original->mend, sizeof (int) * (size_t) np);
  kopie->atl = NULL;

  // copy options  
  fill_worldoptions (kopie->options, options, np);
  
  kopie->has_unassigned = options->has_unassigned;
  kopie->has_unassignedfreq = options->has_unassignedfreq;
  kopie->has_unassignedpoplist = options->has_unassignedpoplist;
  if (options->has_unassignedpoplist)
    {
      kopie->unassignedpoplistnum = options->unassignedpoplistnum;
      kopie->unassignedpoplist = (long *) mycalloc(options->unassignedpoplistnum,sizeof(long));
      memcpy(kopie->unassignedpoplist,options->unassignedpoplist,sizeof(long)*kopie->unassignedpoplistnum);
    }
  if (kopie->has_unassigned)
    {
      if(kopie->unassigned == NULL)
	//empty_world_unassigned(kopie);
	//else
	fill_world_unassigned(kopie);
    }
  switch (options->datatype)
    {
        case 'm':
            for (locus = 0; locus < original->loci; locus++)
            {
                for (i = 0; i < options->micro_stepnum; i++)
                {
		  memcpy(kopie->options->steps[locus][i],
			 original->options->steps[locus][i],
			 sizeof(MYREAL)* (size_t) options->micro_stepnum);
                }
            }
            break;
    case 'a':
      kopie->data->freq = original->data->freq;
      kopie->data->freqlast = original->data->freqlast;
    }

  // copy data parts
  fill_worlddata (kopie->data, original, data, np, options->readsum);
  /* mighistloci not copied */
  kopie->options->datatype = original->options->datatype;
  // the tree was copied in earlier versions [<2.3] using:   copy_tree (original, kopie);
  // now this changed to a complete rebuild of the tree from scratch to guarantee that the
  // the starting points are random (if random was chosen)
  memcpy (kopie->param0, original->param0,
	  sizeof (MYREAL) * (size_t) npp);
  memcpy (kopie->param00, original->param00,
	  sizeof (MYREAL) * (size_t) npp);
  /* fstparam not copied */
  //the tree is generated from scratch therefore no: memcpy (kopie->lineages, original->lineages,
  //                                                         sizeof (long) * original->numpop);
  //old  create_treetimelist (kopie, &kopie->treetimes);
  //old fix_times (kopie, options);
  memcpy (kopie->mig0list, original->mig0list,
	  sizeof (MYREAL) * (size_t) np);
  for (i = 0; i < 2 * np; ++i)
    {
      memcpy (kopie->migproblist[i], original->migproblist[i],
	      sizeof (MYREAL) * (size_t) original->numpop);
    }
  kopie->speciesproblist = kopie->migproblist + np;
  memcpy (kopie->design0list, original->design0list,
	  sizeof (long) * (size_t) np);
#ifdef UEP
  
  if (original->options->uep)
    {
      for (j = 0; j < original->data->uepsites; ++j)
        {
	  memcpy (kopie->ueplike[j], original->ueplike[j],
		  sizeof (MYREAL) * original->numpop);
	  memcpy (kopie->ueptime[j].populations,
		  original->ueptime[j].populations,
		  sizeof (long) * original->ueptime[j].size);
	  memcpy (kopie->ueptime[j].ueptime,
		  original->ueptime[j].ueptime,
		  sizeof (MYREAL) * original->ueptime[j].size);
	  kopie->ueptime[j].size = original->ueptime[j].size;
        }
    }
  if (original->options->uep)
    for (j = 0; j < original->loci; j++)
      memcpy (kopie->uepanc[j],
	      original->uepanc[j],
	      sizeof (long) * original->data->uepsites *
	      np * 2);
#endif
  /* here would come some things that are filled later */
  kopie->heat = 1. / temperature;
  kopie->varheat = pow(CHAINVARIANCEDELTA,kopie->heat);
  kopie->essminimum = pow(ESSMINIMUM,kopie->heat);
  kopie->starttime = original->starttime;
  kopie->treesdone = original->treesdone;
  kopie->treestotal = original->treestotal;
  kopie->migration_counts = original->migration_counts;
  kopie->chains = original->chains;
  kopie->start = original->start;

  //old  set_tree_dirty (kopie->root);
  //old first_smooth(kopie,kopie->locus);
  //old kopie->likelihood[0] = treelikelihood (kopie);

  if (original->has_estimateseqerror)
    {
      kopie->has_estimateseqerror = original->has_estimateseqerror;
      fill_world_seqerror(kopie,options);
    }
  kopie->has_speciation = original->has_speciation;
  kopie->has_migration = original->has_migration;
  //init_speciesvector(kopie);
  //fill_speciesvector(kopie,options);
  reset_growth(kopie);
}


/// clones only parts of the world structure during the heating runs
void
klone_part (world_fmt * original, world_fmt * kopie,
            option_fmt * options, data_fmt * data, MYREAL temperature)
{
  (void) options;
  (void) data;
  (void) temperature;
  long npp = original->numpop2 + original->bayes->mu + original->species_model_size * 2 + original->grownum;
    long i;
#ifdef UEP
    
    long j;
#endif
    
    kopie->repkind = SINGLECHAIN;
    kopie->loci = original->loci;
    kopie->skipped = original->skipped;
    kopie->locus = original->locus;
    kopie->numpop = original->numpop;
    kopie->numpop2 = original->numpop2;
    kopie->sumtips = original->sumtips;
    memcpy (kopie->param0, original->param0,
            sizeof (MYREAL) * (size_t) npp);
    memcpy (kopie->param00, original->param00,
            sizeof (MYREAL) * (size_t) npp);
    memcpy (kopie->mig0list, original->mig0list,
            sizeof (MYREAL) * (size_t) original->numpop);
    for (i = 0; i < 2 * original->numpop; ++i)
    {
        memcpy (kopie->migproblist[i], original->migproblist[i],
                sizeof (MYREAL) * (size_t) (original->numpop - 1));
    }
    memcpy (kopie->design0list, original->design0list,
            sizeof (long) * (size_t) original->numpop);
#ifdef UEP
    
    if (original->options->uep)
    {
        for (j = 0; j < original->data->uepsites; ++j)
        {
            memcpy (kopie->ueplike[j], original->ueplike[j],
                    sizeof (MYREAL) * original->numpop);
            memcpy (kopie->ueptime[j].populations,
                    original->ueptime[j].populations,
                    sizeof (long) * original->ueptime[j].size);
            memcpy (kopie->ueptime[j].ueptime,
                    original->ueptime[j].ueptime,
                    sizeof (MYREAL) * original->ueptime[j].size);
            kopie->ueptime[j].size = original->ueptime[j].size;
        }
    }
    if (original->options->uep)
        for (j = 0; j < original->loci; j++)
            memcpy (kopie->uepanc[j],
                    original->uepanc[j],
                    sizeof (long) * original->data->uepsites *
                    original->numpop * 2);
#endif
}

/// finish up the cloning of world structure
void
clone_polish (world_fmt * original, world_fmt * kopie)
{
    long i;
    kopie->increment = original->increment;
    kopie->chains = original->chains;
    kopie->lsteps = original->lsteps;
    kopie->options->lincr = original->options->lincr;
    kopie->options->lsteps = original->options->lsteps;
    kopie->numlike = original->numlike;
    memcpy (kopie->mig0list, original->mig0list, sizeof (MYREAL) * (size_t) original->numpop);
    for (i = 0; i < 2 * original->numpop; ++i)
    {
        memcpy (kopie->migproblist[i], original->migproblist[i],
                sizeof (MYREAL) * (size_t) (original->numpop - 1));
    }
    memcpy (kopie->design0list, original->design0list,
            sizeof (long) * (size_t) original->numpop);
    polish_world (kopie);
    if (kopie->likelihood[0]<SMALL_VALUE && kopie->likelihood[0]> -SMALL_VALUE)
        kopie->likelihood[0] = original->likelihood[0];
}


/// calculate the probability to swap two heated chains and swap them
/// temperatures are only used for the tree and not for the parameters
/// therefore no swapping on parameters.
long
chance_swap_tree (world_fmt * tthis, world_fmt * that)
{
  MYREAL templike;
  MYREAL a, b, rr, quot;
  MYREAL ha = tthis->heat;
  MYREAL hb = that->heat;
#ifdef UEP
  MYREAL **tempuep;
  ueptime_fmt *tempueptime;
  long *tempanc;
  MYREAL treelen;
#endif
  timelist_fmt *templist;

  a = tthis->likelihood[tthis->numlike - 1];
  b = that->likelihood[that->numlike - 1];

  if(fabs(a-b) < DBL_EPSILON)
    return 0;

  if (tthis == that)
    error ("mcmcmc swap error: original and target chain are the same");

  rr = LOG (RANDUM ());
  // ------------ begin check section -----------------------------------------
  //to check whether thermodynamic integration with swapping makes a difference
  //quot = rr - 1.;
  //this prohibits the swapping.
  // ------------ end check section -------------------------------------------
  quot =  ((a * hb) + (b * ha)) -  ((a * ha) +  (b * hb));
  if (rr < quot)
    {
      swap_tree (tthis, that);
#ifdef DISPENSER
      swap_node_collection(tthis,that);
#endif
      if (tthis->has_unassigned)
	swap_unassigned_nodecollection(tthis,that);
      //swap likelihood value that goes with tree
      templike                              = tthis->likelihood[tthis->numlike - 1];
      tthis->likelihood[tthis->numlike - 1] = that->likelihood[that->numlike - 1];
      that->likelihood[that->numlike - 1]   = templike;
      // swap timelist that goes with tree
      templist                              = tthis->treetimes;
      tthis->treetimes                      = that->treetimes;
      that->treetimes                        = templist;
#ifdef UEP        
      if (tthis->options->uep)
        {
	  tempuep = tthis->ueplike;
	  tthis->ueplike = that->ueplike;
	  that->ueplike = tempuep;
	  tempueptime = tthis->ueptime;
	  tthis->ueptime = that->ueptime;
	  that->ueptime = tempueptime;
	  templike = tthis->ueplikelihood;
	  tthis->ueplikelihood = that->ueplikelihood;
	  tthis->ueplikelihood = templike;
	  tempanc = tthis->oldrootuep;
	  tthis->oldrootuep = that->oldrootuep;
	  tthis->oldrootuep = tempanc;
	  treelen = tthis->treelen;
	  tthis->treelen = that->treelen;
	  that->treelen = treelen;
        }
#endif
      that->treeswapcount++;
      return 1;
    }
  //  if(that->heat==1.0)
  //  printf("    : %5.3f %5.3f || %f %f %f %f %f %f\n",tthis->heat, that->heat, rr, quot, a, b, ha,hb);
  return 0;
}

/// prepare heated chains for next round of MCMC updates
void
advance_clone_like (world_fmt * world, long accepted, long *j)
{
    if (accepted > 0)
    {
        world->numlike = 1 + *j + (accepted > 0 ? 1 : 0);
        *j += 1;
        if(world->alloclike <= world->numlike)
        {
            world->alloclike = world->numlike + 100;
            world->likelihood =
                (MYREAL *) myrealloc (world->likelihood,
                                      sizeof (MYREAL) * (size_t) world->alloclike);
        }
        world->likelihood[*j] = world->likelihood[*j - 1];
    }
}

void
polish_world (world_fmt * world)
{
    MYREAL tmp = world->likelihood[0];
    world->numlike = 1;
    world->likelihood[0] = tmp;
        world->treelen = 0.0;
        calc_treelength (world->root->next->back, &world->treelen);
#ifdef UEP    
    if (world->options->uep)
    {
        world->treelen = 0.0;
        calc_treelength (world->root->next->back, &world->treelen);
        update_uep (world->root->next->back, world);
        check_uep_root (world->root->next->back, world);
    }
#endif
}

/// print Gelman's convergence criterium
void
print_gelmanr (MYREAL average, MYREAL biggest)
{
    if (!MYISNAN ((float)average) && average < 100000)
        FPRINTF (stdout, "           Average Gelman's R = %f\n", average);
    else
        FPRINTF (stdout, "           Average Gelman's R = -failed-\n");
    if (!MYISNAN ((float)biggest) && biggest < 10000000)
        FPRINTF (stdout,
                 "           Largest Gelman's R = %f [Value < 1.2 show convergence]\n",
                 biggest);
    else
        FPRINTF (stdout,
                 "           Largest Gelman's R = -failed- [Value < 1.2 show convergence]\n");
}

/// calculates the start and end of a replication step, sets repstop = repstart+1
void set_replicates (world_fmt * world, long repkind, long rep, long *repstart,
		     long *repstop)
{
    if (world->options->replicate)
    {
        *repstart = (repkind == SINGLECHAIN) ? rep : 0;
        *repstop = (repkind == SINGLECHAIN) ? rep + 1 : world->repstop;
    }
    else
    {
        *repstart = 0;
        *repstop = 1;
    }
}


boolean updating(world_fmt *world)
{
  long i;
  boolean success = FALSE;
  long count = 0;
  double *choices = world->options->choices;
  long choice=0;
  const boolean has_mu = world->bayes->mu;
  const long np = world->numpop2 + has_mu + world->species_model_size * 2 + world->grownum;
  double r = RANDUM();
  while(r>choices[choice])
    {
      choice++;
    }
  //update ratio for trees versus parameters
  //  printf("%li-->{%f %f %f}\n", choice, choices[0],choices[1],choices[2]);
  switch(choice)
    {
    case 0:/*TREEUPDATE*/
      count = tree_update (world, world->G, NOASSIGNING);
      world->accept += count;
      world->bayes->accept[np] += count;	//permanent recorder for the tree acceptance, 
      world->bayes->trials[np] += 1;
      //printf("@ updatin: [%i,%li] tree accept=%li (%li)\n",world->cold,np,world->bayes->accept[np], count);
      world->param_like = world->bayes->oldval;
      break; 
    case 1:/*PARAMETERUPDATE*/
      world->bayesaccept = bayes_update (world);
      world->param_like = world->bayes->oldval;
      break;
    case 2: /*HAPLOTYPEUPDATE*/
      success = (boolean) swap_haplotypes(world);
      break;
    case 3: /*SKYLINETIMEUPDATE*/
      // time update
      world->param_like = bayes_update_timeparam(world, &success);
      break;
    case 4: /* ASSIGNMENTUPDATE*/
      // assignment update
      count = tree_update (world, world->G, world->has_unassigned);
      world->assigncount += count;
      world->accept += count;
      world->bayes->accept[np] += count;	//permanent recorder for the tree acceptance, 
      world->bayes->trials[np] += 1;

      // reporting associated population
      // input order is that equal to nodep order?
      if (world->options->verbose && world->cold)
	{
	  for (i=0;i<world->sumtips;i++)
	    {
	      fprintf(stdout,"%ld",world->nodep[i]->actualpop);
	    }
	  fprintf(stdout,"\n");
	}
      break;
    case 5: /*SEQUENCEERRORUPDATE*/
      //seqerror-rate update
      //assignment update
      change_freq(world);
      break;
    case 6: /*MITTAGLEFFLERUPDATE*/
      change_mittag_leffler(world);
      break;
    default:
      printf("%li-->{%f %f %f %f %f}\n", choice, choices[0],choices[1],choices[2],choices[3],choices[4]);
      error("failure in choosing among updates -- ML");
    }
  return success;
}

long get_numparam(world_fmt* world)
{
  //returns the number of total parameters
  return world->numpop2 + world->bayes->mu + 2 * world->species_model_size + world->grownum;
}


/*boolean autostopping(world_fmt *world, 
		     boolean autostop, 
		     boolean essstop,
		     long *z,
		     long step, 
		     long delta,
		     MYREAL *var,
		     MYREAL *oldvar,
		     MYREAL *oldvar2,
		     MYREAL *autocorrelation,
		     MYREAL *effective_sample)
{
  boolean done=FALSE;
  !!!!!const long nn = world->numpop2 + world->bayes->mu * world->loci + 1;
  // stopping automatically for Bayes and ML but parts copied to reduce calls on if statements
  if(autostop)
    {
      if((step % delta) == 0)
	{
	  *oldvar2= *var;
	  single_chain_var (world, step, var, autocorrelation, effective_sample);
	  if(essstop)
	    {
	      done = max_ess(effective_sample,nn,world->essminimum);
	    }
	  else
	    {
	      done = (fabs(*var/ *oldvar - 1.0) < world->varheat);
	    }
	  if(done)
	    {
	      if(world->cold)
		{
#ifdef MPI
		  char p[LINESIZE];
		  char p1[LINESIZE];
		  long bufsize;
		  if(myID!=MASTER)
		    {
		      bufsize = sprintf(p,"%li %f %f %li\n",world->locus, *var, *oldvar,step);
		      sprintf(p1,"B%li",bufsize);
		      MYMPISEND (p1, SMALLBUFSIZE, MPI_CHAR, (MYINT) MASTER, (MYINT) myID+BURNTAG, comm_world);
		      MYMPISEND (p, bufsize, MPI_CHAR, (MYINT) MASTER, (MYINT) myID+BURNTAG, comm_world);
		    }
#else
                  if(world->burnin_stops_alloc <= *z)
		  {
		  world->burnin_stops_alloc += *z;
		  world->burnin_stops = (burnin_record_fmt *) myrealloc(world->burnin_stops, world->burnin_stops_alloc * sizeof(burnin_record_fmt));  
		  }

		  world->burnin_stops[*z].locus        =  world->locus; 
		  world->burnin_stops[*z].replicate =  world->rep;
		  world->burnin_stops[*z].stopstep = step; 
		  world->burnin_stops[*z].variance = *var; 
		  world->burnin_stops[*z].oldvariance = *oldvar;
		  world->burnin_stops[*z].worker = myID;
		  *z += 1;
		  world->burnin_counter = *z;
#endif
		}
	    }
	}
    }
  return done;
}
*/
