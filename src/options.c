/* \File options.c */
/*------------------------------------------------------
 Bayesian estimation 
 of migration rate, divergence,  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm   
 ------------------------------------------------------- 
 O P T I O N S   R O U T I N E S 
 
 creates options structures,
 reads options from parmfile if present
 
 prints options,
 and finally helps to destroy itself.
 
 Peter Beerli 1996-2018
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2018 Peter Beerli, Tallahassee FL
 
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
 
 
$Id: options.c 2158 2013-04-29 01:56:20Z beerli $
-------------------------------------------------------*/
#pragma clang diagnostic ignored "-Wformat-nonliteral"

//#include <stdio.h>
#include <time.h>
#include "migration.h"
#include "sighandler.h"
#include <stdarg.h>

#include "tools.h"
#include "migrate_mpi.h"
#include "bayes.h"
#include "priors.h"
#include "mutationmodel.h"
#include "menu.h"
#include "random.h"
#include "options.h"

#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif
#include "options.h"

extern char * generator;

/* parmfile parameter specifications and keywords */
//#define LONGLINESIZE 10000
#define NUMBOOL 36
#define BOOLTOKENS {\
"menu", /*0*/					\
"recover",\
"print-data",\
"mixfile",\
  "moving-steps",\
"freqs-from-data",\
"useroldtree", \
  "autocorrelation",\
 "simulation",\
"bayes-all-posteriors", \
"weights",\
  "read-summary" /*11*/,				\
"write-summary",\
"NODATA",\
 "include-unknown",\
 "print-fst",\
  "distfile",\
"geo",\
"gelman-convergence",\
 "randomtree",\
      "fast-likelihood",\
  "pdf-terse" /*21*/,	\
 "use-M",\
 "use-Nm",\
"bayes-update",\
 "bayes-allfile",\
 "bayes-file",\
 "divergence-times",\
"tipdate-file",\
"heated-swap",\
 "haplotyping",\
  "auto-tune" /*31*/,				\
  "assign",\
  "bayes-hyperpriors",\
  "inheritance-scalars", "mittag-leffler-alpha"}
#define NUMNUMBER 69
#define NUMBERTOKENS {"ttratio","rate",\
 "split","splitstd","long-chains",\
 "long-steps", "long-inc", "theta", \
 "nmlength","random-seed","migration","mutation",\
 "datatype", "categories","rates","prob-rates", \
 "micro-max", "micro-threshold", "delimiter","burn-in",\
 "infile", "outfile", "smoothing", "title", \
 "long-chain-epsilon","print-tree","progress","l-ratio",\
 "fst-type","profile","custom-migration","population-growth","unused",\
 "long-sample", "replicate","datamodel","logfile", "seqerror-rate","uep", \
 "uep-rates","uep-bases", "mu-rates","heating","fluctuate", "resistance",\
 "bayes-updatefreq", "bayesfile","bayes-prior", "usertree", "bayes-posteriorbins",\
 "mig-histogram", "bayes-posteriormaxtype", "pdf-outfile",\
 "bayes-allfileinterval", "bayes-priors","skyline","rates-gamma", "bayes-proposals", \
      "generation-per-year", "mutationrate-per-year","inheritance-scalars", "micro-submodel","random-subset", "population-relabel", "sequence-submodel", "updatefreq", "analyze-loci","divergence-distrib","mittag-leffler-alpha"};

// myID is a definition for the executing node (master or worker)
extern int myID;

static fpos_t thePos;   // used to track pposition in the parmfile

/* prototypes ------------------------------------------- */
//void set_usem_related (option_fmt * options);
//void create_options (option_fmt ** options);
//void init_options (option_fmt * options);
//void set_param (world_fmt * world, data_fmt * data, option_fmt * options,long locus);
//void set_profile_options (option_fmt * options);
//void print_menu_options (world_fmt * world, option_fmt * options,
//                         data_fmt * data);
//void print_options (FILE * file, world_fmt * world, option_fmt * options,
//                    data_fmt * data);
//void decide_plot (worldoption_fmt * options, long chain, long chains,
//                  char type);
//void destroy_options (option_fmt * options);

//void read_options_master (option_fmt * options);
//void read_options_worker (char **buffer, option_fmt * options);
///* private functions */
boolean booleancheck (option_fmt * options, char *var, char *value);
//long boolcheck (char ch);
boolean numbercheck (option_fmt * options, char *var, char *value);
//void reset_oneline (option_fmt * options, long position);
void reset_oneline (option_fmt * options, fpos_t * position);
//void read_theta (option_fmt * options, char *parmvar, char * varvalue, char **buffer);
void read_startparameter (long parampos, option_fmt * options, char *parmvar, char *varvalue, char **buffer);
//void add_thetaguess(option_fmt *options, char *tmp);
//void read_mig (option_fmt * options, char *parmvar, char *varvalue, char **buffer);
//void add_mguess(option_fmt *options, char *tmp);
#ifdef MPI
//void read_theta_worker (char **buffer, option_fmt * options, char *parmvar, char *varvalue);
//void read_mig_worker (char **buffer, option_fmt * options, char *parmvar, char *varvalue);
char skip_sspace (char **buffer);
void read_custom_migration_worker (char **buffer, option_fmt * options,
                               char *value, long customnumpop);
#endif
char skip_space (option_fmt * options);
//void read_custom_migration (FILE * file, option_fmt * options, char *value,
//                            long customnumpop);
//void synchronize_param (world_fmt * world, option_fmt * options);
//void resynchronize_param (world_fmt * world);
void specify_migration_type (option_fmt * options);
//void print_options_nomcmc (FILE * file, option_fmt * options,
//                               world_fmt * world);
//void destroy_options (option_fmt * options);
//void set_partmean_mig (long **mmparam, MYREAL *param, char *custm2, long migm,
//                       long numpop2);
//long scan_connect (char *custm2, long start, long stop, int check);
//void set_plot (option_fmt * options);
//void set_plot_values (MYREAL **values, MYREAL plotrange[], long intervals, int type);
//void set_grid_param (world_fmt * world, long gridpoints);
void print_arbitrary_migration_table (FILE * file, world_fmt * world, option_fmt *options,
                                      data_fmt * data);
void print_distance_table (FILE * file, world_fmt * world,
                           option_fmt * options, data_fmt * data);

void fillup_custm (long len, world_fmt * world, option_fmt * options);
//

//long save_options_buffer (char **buffer, long *allocbufsize, option_fmt * options, data_fmt *data);
//void print_parm_delimiter(long *bufsize, char **buffer, long *allocbufsize);
//void print_parm_br(long *bufsize, char **buffer, long *allocbufsize);
//
void prior_consistency(option_fmt *options, int paramtype, char priortype);

void init_filename (char **filename, char initstring[]);
void show_pretty_datamodel(int datamodel, char *buffer);
char * show_prioralpha(char *tmp, prior_fmt *prior);
char * show_priorupdatefreq(char *tmp, prior_fmt *prior);
void add_startparam(option_fmt *options, short key, float value);
void set_startparam_randomstart(world_fmt *world, option_fmt *options);
void set_theta_nrandomstart(world_fmt *world, option_fmt *options);
void set_mystartparams(long i, long numx,  long guess, float *ppp, world_fmt * world, option_fmt * options, prior_fmt * priors);
void set_param_fromstartparam(world_fmt *world, option_fmt *options);
MYREAL set_paramvalue_data_fst(char datatype,MYREAL val1, MYREAL val2);
void set_theta_fststart(world_fmt *world, option_fmt *options, long locus);
void set_mig_urandomstart(world_fmt *world, option_fmt *options);
void set_mig_nrandomstart(world_fmt *world, option_fmt *options);
void set_mig_ownstart(world_fmt *world, option_fmt *options, data_fmt *data);
void set_mig_fststart(world_fmt *world, option_fmt *options, long locus);
void set_partmean_mig (long **mmparam, MYREAL *param, char *custm2, long migm, long numpop2);
void free_options_filenames(option_fmt * options);
void destroy_startparameters(option_fmt* options);
void print_parm_delimiter(long *bufsize, char **buffer, long *allocbufsize);
void print_parm_smalldelimiter(long *bufsize, char **buffer, long *allocbufsize);
void print_parm_br(long *bufsize, char **buffer, long *allocbufsize);
void print_parm_comment(long *bufsize, char **buffer, long *allocbufsize, char message[]);

/// prints parmfile mutable comment line
void print_parm_mutable_comment(long *bufsize, char **buffer, long *allocbufsize, char string[], ...);
void print_parm(long *bufsize, char **buffer, long *allocbufsize, char string[]);
void print_parm_title(long *bufsize, char **buffer, long *allocbufsize, char message[]);
void print_parm_ttratio(long *bufsize, char **buffer,  long *allocbufsize, option_fmt *options);
void print_parm_freqfrom(long *bufsize, char **buffer,  long *allocbufsize, option_fmt * options);
void print_parm_categs(long *bufsize, char **buffer,  long *allocbufsize, option_fmt * options);
void print_parm_weights(long *bufsize, char **buffer, long *allocbufsize, option_fmt * options);
void print_parm_rates(long *bufsize, char **buffer, long *allocbufsize, option_fmt * options);
void print_parm_datatype(long *bufsize, char **buffer, long * allocbufsize, option_fmt *options);
void print_parm_tipdate(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options, data_fmt *data);
void print_parm_inheritance(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options, data_fmt *data);
void print_parm_newpops(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options, data_fmt *data);
void print_parm_growpops(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options, data_fmt *data);
void print_parm_randomsubset(long * bufsize, char **buffer, long *allocbufsize, option_fmt *options);
void print_parm_usertree(long * bufsize, char **buffer, long *allocbufsize, option_fmt *options);
void print_parm_theta(long *bufsize, char ** buffer, long *allocbufsize, option_fmt * options);
void print_parm_m(long *bufsize, char **buffer, long *allocbufsize, option_fmt * options);
void print_parm_split(long *bufsize, char **buffer, long *allocbufsize, option_fmt * options);
void print_parm_splitstd(long *bufsize, char **buffer, long *allocbufsize, option_fmt * options);
void print_parm_heating(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options);
void print_parm_haplotyping(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options);
boolean  set_filename(char *value, char comparison[], char ** filename);
void  set_filename_only(boolean check, char *value, char ** filename);
void   set_parm_prior_values(prior_fmt * prior, char * mytext);
char * show_proposaltype(boolean priorset);
char * show_parmpriortype(int priorset);
void print_parm_proposal(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options);
char * getpriortype(int kind);
void print_parm_prior(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options);
void print_parm_hyperprior(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options);
long save_inheritance_scalars_buffer (char **buffer, long *allocbufsize, option_fmt * options);
long save_newpops_buffer (char **buffer, long *allocbufsize, option_fmt * options);
void set_growth(char **value, char **tmp, option_fmt *options);
long save_growpops_buffer (char **buffer, long *allocbufsize, option_fmt * options);
void set_bayes_options(char *value, option_fmt *options);
long boolcheck (char ch);
void read_random_theta(option_fmt *options, char ** buffer);
void read_random_mig(option_fmt *options, char ** buffer);
void read_theta (option_fmt * options, char *parmvar, char *varvalue, char **buffer);
void read_mig (option_fmt * options, char *parmvar, char * varvalue, char **buffer);
void add_custom_mig(option_fmt *options, long *z, char ch);
int sequence_submodeltype(char *mutmodel);
void print_parm_mutable(long *bufsize, char **buffer, long *allocbufsize, char string[], ...);
void print_parm_smoothing(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options);
//##
//
///*======================================================*/
/// initialize filename
void init_filename (char **filename, char initstring[])
{
  int len=LINESIZE;
  //  len = strlen(initstring);
    *filename = (char *) mycalloc(len+2, sizeof(char));
    sprintf(*filename, "%-*s", len, initstring);
    remove_trailing_blanks(filename);
}


void
create_options (option_fmt ** options)
{
    (*options) = (option_fmt *) mycalloc (1, sizeof (option_fmt));
}

/// initialize options
void init_options (option_fmt * options)
{
    long i;
    unsigned long timeseed;
    /* General options --------------------------------------- */
    //
    options->smoothing[0]=KERNEL;
    options->smoothing[1]=SAVITZKY;
    options->smoothing_window[0] = DEFAULTKERNELWINDOW;
    options->smoothing_window[1] = DEFAULTKERNELWINDOW;
    // name length
    options->nmlength = DEFAULT_NMLENGTH;
    options->popnmlength = DEFAULT_POPNMLENGTH;
    options->allelenmlength = DEFAULT_ALLELENMLENGTH;
    //
    // custom migration matrix setup
    options->custm = (char *) mycalloc (1, sizeof (char) * 1000); /// \todo needs
    //options->custm2 is allocated later
    //options->symn = 0;
    //options->sym2n = 0;
    //options->zeron = 0;
    //options->constn = 0;
    //options->mmn = 0;
    //options->tmn = 0;
    
    /* input/output options ---------------------------------- */
    options->menu = TRUE;
    options->progress = TRUE;
    options->verbose = FALSE;
    options->writelog = FALSE;
    options->geo = FALSE;
    options->div = FALSE;
    options->recorddivtime = FALSE;
#ifdef UEP
    options->uep = FALSE;
    options->ueprate = 1.0;
    options->uepmu = 0.9999;
    options->uepnu = 0.0001;
    options->uepfreq0=0.5;
    options->uepfreq1=0.5;
    //  options->uep_last = FALSE;
#endif
    options->printdata = FALSE;
    options->printcov = FALSE;
    options->usertree = FALSE;
    options->usertreewithmig = FALSE;
    options->randomtree = TRUE;  //default changed May 29 2006 [convergence statistic]
    options->fastlike = FALSE;
    options->treeprint = myNONE;
    options->treeinc = 1;
    options->printfst = FALSE;
    options->fsttype = THETAVARIABLE;
    options->usem = TRUE;
    options->migvar = PLOT4NM;

    options->plot = FALSE;
    //options->plotmethod = PLOTALL; /* outfile and mathematica file */
    //options->plotvar = PLOT4NM;
    //options->plotscale = PLOTSCALELOG;
    //options->plotrange[0] = PLANESTART; /* start x axis */
    //options->plotrange[1] = PLANEEND; /*end x axis */
    //options->plotrange[2] = PLANESTART; /*start y axis */
    //options->plotrange[3] = PLANEEND; /* end y axis */
    //options->plotintervals = PLANEINTERVALS;
    options->simulation = FALSE;
    options->movingsteps = FALSE;
    options->acceptfreq = 0.0;
    options->mighist = FALSE;
    options->mighist_all = FALSE;
    options->mighist_counter = 0;
    options->mighist_increment = 1;
    options->skyline = FALSE;
    options->skyline_param = FALSE;
    options->timeelements = 2;
    options->eventbinsize = (float) 0.001;
    options->mixplot = FALSE;
    init_filename( &options->parmfilename, PARMFILE);
    init_filename( &options->infilename, INFILE);
    init_filename( &options->outfilename, OUTFILE);
#ifdef PRETTY
    init_filename( &options->pdfoutfilename, PDFOUTFILE);
#endif
#ifdef THERMOCHECK
    init_filename( &options->mixfilename, MIXFILE);
#endif
    init_filename( &options->logfilename, LOGFILE);
    //init_filename( &options->mathfilename, MATHFILE);
    init_filename( &options->sumfilename, SUMFILE);
    init_filename( &options->treefilename, TREEFILE);
    init_filename( &options->utreefilename, UTREEFILE);
    init_filename( &options->catfilename, CATFILE);
    init_filename( &options->weightfilename, WEIGHTFILE);

    init_filename( &options->mighistfilename, MIGHISTFILE);
    init_filename( &options->skylinefilename, SKYLINEFILE);

    init_filename( &options->distfilename, DISTFILE);
    init_filename( &options->geofilename, GEOFILE);
    init_filename( &options->divfilename, DIVFILE);
    init_filename( &options->bootfilename, BOOTFILE);
    init_filename( &options->divtimefilename, DIVTIMEFILE);
    init_filename( &options->seedfilename, SEEDFILE);    
#ifdef UEP
    init_filename( &options->uepfilename, UEPFILE);
#endif
    init_filename( &options->bayesfilename, BAYESFILE);
#ifdef ZNZ
    init_filename( &options->bayesmdimfilename, BAYESMDIMFILE);
    options->use_compressed=1;
#else
    init_filename( &options->bayesmdimfilename, BAYESMDIMFILE2);
    options->use_compressed=0;
#endif
    init_filename( &options->datefilename, TIPDATEFILE);
    // allocation of prior parameter settings
    options->bayes_priors_num = 0; // priors will be set for all parameters
    options->bayes_priors = NULL;// once they are read from parmfile or then assigned from default in check_bayes_priors()
    options->automatic_bins = TRUE;

    options->hyperprior=FALSE;

    for(i=0;i<PRIOR_SIZE;i++)
      options->slice_sampling[i] = FALSE;// option for proposal setting
#ifdef PRETTY
    options->bayespretty = PRETTY_P100;
#endif
    strcpy (options->title, " ");
    options->lratio = (lratio_fmt *) mycalloc (1, sizeof (lratio_fmt));
    options->lratio->alloccounter = 1;
    options->lratio->data =
      (lr_data_fmt *) mycalloc (options->lratio->alloccounter, sizeof (lr_data_fmt));
    options->lratio->data[0].value1 =
        (char *) mycalloc (1, sizeof (char) * LONGLINESIZE);
    options->lratio->data[0].value2 =
        (char *) mycalloc (1, sizeof (char) * LONGLINESIZE);
    options->lratio->data[0].connect =
        (char *) mycalloc (1, sizeof (char) * LONGLINESIZE);
    options->tersepdf = FALSE;
    options->profile = ALL;
    //options->profilemethod = 'p';
    options->df = 1;
    //options->qdprofile = FALSE;
    //options->printprofsummary = TRUE;
    //options->printprofile = TRUE;
    //if (options->usem)
    //    options->profileparamtype = options->plotvar = options->migvar = PLOTM;
    //else
    //    options->profileparamtype = options->plotvar = options->migvar = PLOT4NM;
    /* data options ------------------------------------------ */
    options->datatype = 's';
    options->include_unknown = FALSE;
    options->migration_model = MATRIX;
    //options->thetag = (MYREAL *) mycalloc (1, sizeof (MYREAL) * NUMPOP);
    //options->mg = (MYREAL *) mycalloc (1, sizeof (MYREAL) * NUMPOP);
    allocate_startparam(options,NUMPOP);
    options->gamma = FALSE;
    options->alphavalue = START_ALPHA;
    options->murates = FALSE;
    options->bayesmurates = FALSE;
#ifdef LONGSUM
    options->fluctuate=FALSE;
    options->flucrates = (MYREAL *) mycalloc(NUMPOP * 3, sizeof(MYREAL));
#endif
    options->mu_rates = NULL;
    options->inheritance_scalars = NULL;
    options->newpops = NULL;
    options->growpops = NULL;
    options->prioralone=FALSE;
    /* EP data */
    options->dlm = '\0';
    /* microsat data */
    options->micro_threshold = MICRO_THRESHOLD;
    options->micro_stepnum = MAX_MICROSTEPNUM;
    options->msat_option = SINGLESTEP;
    options->msat_tuning[0] = 0.0;
    options->msat_tuning[1] = 0.5;
    /*sequence data */
    options->sequence_model = F84;
    options->sequence_model_parameters = (MYREAL *) mycalloc(16, sizeof(MYREAL));
    options->sequence_model_numparam = 1;
    options->sequence_model_parameters[0] = 2.0;
    options->seqerror[0]=options->seqerror[1]=options->seqerror[2]=options->seqerror[3]=0.0;
    options->has_estimateseqerror = FALSE;
    options->seqerrorcombined = FALSE;
    options->interleaved = FALSE;
    options->ttratio = (MYREAL *) mycalloc (1, sizeof (MYREAL) * 3);
    options->ttratio[0] = 2.0;
    options->freqsfrom = TRUE;
    options->categs = ONECATEG;
    options->rate = (MYREAL *) mycalloc (1, sizeof (MYREAL));
    options->rate[0] = 1.0;
    options->rcategs = 1;
    options->rrate = (MYREAL *) mycalloc (1, sizeof (MYREAL));
    options->probcat = (MYREAL *) mycalloc (1, sizeof (MYREAL));
    options->autocorr = FALSE;
    options->rrate[0] = 1.0;
    options->probcat[0] = 1.0;

    options->probsum = 0.0;

    options->lambda = 1.0;
    options->weights = FALSE;
    /* random number options --------------------------------- */
    options->autoseed = AUTO;
    options->autoseed = AUTO;
    timeseed = (unsigned long) time (NULL) / 4;
    options->inseed = (unsigned long) timeseed + 1;
    /* mcmc options ------------------------------------------ */
    //options->thetaguess = OWN;
    //options->migrguess = OWN;
    for (i = 0; i < NUMPOP; i++)
    {
    //    options->thetag[i] = 1.0;
    //    options->mg[i] = 1.0;
        options->custm[i] = '*';
    }
    options->custm[i] = '\0';
    options->custm2 = (char * ) strdup(options->custm);
    //options->numthetag = options->nummg = 0;
    options->bayes_infer = TRUE; //DEFAULT TO BAYESIAN INFERERENCE July 2012
#ifdef INTEGRATEDLIKE
    options->integrated_like = TRUE;
#else
    options->integrated_like = FALSE;
#endif
    options->updateratio = HALF;
    options->tree_updatefreq = 0.2;
    options->parameter_updatefreq =0.2;
    options->haplotype_updatefreq = 0.2;
    options->timeparam_updatefreq = 0.2;
    options->unassigned_updatefreq = 0.2;
    options->seqerror_updatefreq = 0.1;
#ifdef NEWVERSION
    options->mlalpha_updatefreq = 0.1;
#endif
    options->bayesmdiminterval = 1;
    options->schains = 0;
    options->sincrement = 100;
    options->ssteps = 500;
    options->lchains = 1;
    options->lincrement = INCREMENT;
    options->lsteps = SAMPLESTAKEN;
    options->burn_in = BURNINPERIOD;
    options->burnin_autostop = ' ';
    //
    // heating options
    options->heating = 1;  /* default was changed to heating January 23*/
    options->adaptiveheat=NOTADAPTIVE;
    options->heating_interval = 1; //set to 1 on January 23
    options->heatedswap_off=FALSE;
    options->heat[0] = COLD;
    options->heat[1] = WARM;
    options->heat[2] = HOT;
    options->heat[3] = VERYHOT;
    options->heated_chains = 4; // default is 4 chains
    //
    // obscure options to accept more different trees
    options->lcepsilon = LONGCHAINEPSILON;
    options->gelman = FALSE; 
    options->gelmanpairs = FALSE;
    options->pluschain = PLUSCHAIN;
    //
    // replication options
    options->replicate = FALSE;
    options->replicatenum = 0;
    options->gridpoints = 0;
    /* genealogy summary options-[this relates to ML and options have beeen
       excised -- these settings remain to make sure that downstream nothing
       breaks */
    options->readsum = FALSE;
    options->checkpointing = FALSE;
    options->writesum = FALSE;
    /*threading over loci */
    options->cpu = 1;
    //
    // fatal attraction to zero resistance
    options->minmigsumstat = MINMIGSUMSTAT;
    //
    options->segregs=NULL;
    options->wattersons=NULL;

    options->onlyvariable=FALSE; //analyzes only variable sequence loci:default: analyze all
    options->has_variableandone = FALSE; // analyze only one invariant locus (not more) and reweight
    options->firstinvariant = -1; // locus ID of invariant locus analyzed

    options->mutationrate_year = (MYREAL*) mycalloc(1, sizeof(MYREAL));
    options->mutationrate_year_numalloc = 1;
    options->generation_year = 1.0;
    options->randomsubset = 0;
    options->randomsubsetseed = (unsigned long) -1;
    options->inheritance_scalars = (MYREAL*) mycalloc(1, sizeof(MYREAL));
    options->inheritance_scalars[0] = 1.0;
    options->inheritance_scalars_numalloc  = 1;
    options->newpops = (long*) mycalloc(1, sizeof(long));
    options->newpops[0] = 1;
    options->newpops_numalloc  = 1;
    options->growpops=(long*) mycalloc(1, sizeof(long));
    options->growpops[0] = 0;
    options->growpops_numalloc  = 1;

    options->slice_sticksizes = (MYREAL*) mycalloc(1, sizeof(MYREAL));
    options->slice_sticksizes[0] = 1.0;
    options->heatedswap_off=FALSE;
    //haplotyping
    options->haplotyping=FALSE;
    options->haplotyping_report=FALSE;
    options->has_autotune = TRUE;
    options->autotune=AUTOTUNEDEFAULT;
    options->totalsites = 0;
    // speciation
    options->species_model_dist = NORMAL_DIST; //NORMALSHORTCUT_DIST;
    options->has_speciation=FALSE;
    options->bayes_posterior_bins[0]=BAYESNUMBIN;
    options->bayes_posterior_bins[1]=BAYESNUMBIN;
    options->bayes_posterior_bins[2]=BAYESNUMBIN;
    options->bayes_posterior_bins[3]=BAYESNUMBIN;
    options->bayes_posterior_bins[4]=BAYESNUMBIN;
    options->bayes_posterior_bins[5]=BAYESNUMBIN;
    options->mlalpha = 1.0;
    options->mlinheritance = 4.0;
}

void
read_options_master (option_fmt * options)
{
  const unsigned long linecpmax = LINESIZE -1; 
#ifdef DEBUG
  printf("In read option master\n");
#endif
    long counter = 0;
    char varvalue[LONGLINESIZE];
    char parmvar[LONGLINESIZE];
    char input[LONGLINESIZE];
    char *p, *tmp;
  
    options->parmfile = fopen(options->parmfilename, "r");
    if (options->parmfile)
    {
        counter = 0;
	fgetpos(options->parmfile, &thePos);
	FGETS (input, LINESIZE, options->parmfile);// read first set of ###
	FGETS (input, LINESIZE, options->parmfile);// read "# Parmfile for Migrate"
	if(strncmp(input,"# Parmfile for Migrate", 20L))
	  {
	    usererror("This file is not a parameter file (parmfile) for migrate\nSyntax for calling from the commandline is\nmigrate-n parmfile\n and not (!!!!)\nmigrate-n datafile\n");
	    //exit(-1);
	  }
	fprintf(stdout,"Reading parmfile \"%s\"....\n",options->parmfilename);
	//fix this after all the addition of the things: 
	options->bayes_priors_num=0;

        while (FGETS (input, LINESIZE, options->parmfile) != EOF)
        {
            counter++;
#ifdef DEBUG
	    //fprintf(stdout,"%s: [%s]\n", options->parmfilename, input);
#endif
            if ((input[0] == '#') || isspace ((int) input[0])
                    || input[0] == ';')
            {
	      //position = ftell (options->parmfile);
	      fgetpos(options->parmfile, &thePos);
	      continue;
            }
            else
            {
                if (!(isalnum ((int) input[0]) || strchr ("{}", input[0])))
                {
                    usererror ("The parmfile contains an error on line %li\n",
                               counter);
                }
            }
            if ((p = strchr (input, '#')) != NULL)
                *p = '\n';
            if (!strncmp (input, "end", 3))
                break;
            tmp = strtok (input, "=");
            if (tmp != NULL)
                strncpy (parmvar, tmp, linecpmax);
            else
            {
                if(input[0]!='\0')
                    fprintf (stderr,
                             "WARNING: error in parmfile on line %li with %s\n",
                             counter, input);
                continue;
            }

#ifdef DEBUG
	    //fprintf(stdout,"%s [%s]\n", parmvar,input);
#endif
            tmp = strtok (NULL, "\n");
            if (tmp != NULL)
	      {
		if(strlen(tmp)<LONGLINESIZE)
		  strcpy (varvalue, tmp);
		else
		  strncpy (varvalue, tmp,LONGLINESIZE-1);
	      }
#ifdef DEBUG
	    fprintf(stdout,"parmvar=%s [varval=%s]\n", parmvar,varvalue);
#endif

            if (!booleancheck (options, parmvar, varvalue))
            {
                if (!numbercheck (options, parmvar, varvalue))
                {
                    warning ("Inappropiate entry in parmfile: %s [%s %s] ignored\n",
                             input, parmvar, varvalue);
                }
            }
            //position = ftell (options->parmfile);
	    fgetpos(options->parmfile, &thePos);
        }
    }
}

#ifdef MPI
void
read_options_worker (char **buffer, option_fmt * options)
{
  const unsigned long linecpmax = LINESIZE - 1;

    long counter = 0;
    char *position;
    char varvalue[LONGLINESIZE];
    char parmvar[LONGLINESIZE];
    char input[LONGLINESIZE];
    char *p, *tmp;
    prior_fmt *ptr;
    options->buffer = buffer;
    options->parmfile = NULL;
    //options->bayespriorthetanum=0;
    //options->bayespriormnum=0;
    //options->bayespriorsnum=0;
    if (strlen (*buffer) > 0)
    {
        counter = 0;
        position = *buffer;
        while (sgets (input, LINESIZE, buffer) != NULL)
        {
            counter++;
#ifdef DEBUG
	    // fprintf(stdout,"%i: [%s]\n", myID, input);
#endif

            if ((input[0] == '#') || isspace ((int) input[0])
                    || input[0] == ';')
            {
                position = *buffer;
                continue;
            }
            else
            {
                if (!(isalnum ((int) input[0]) || strchr ("{}", input[0])))
                {
                    usererror ("The parmfile contains an error on line %li\n",
                               counter);
                }
            }
            if ((p = strchr (input, '#')) != NULL)
                *p = '\n';
            if (!strncmp (input, "end", 3))
                break;
            tmp = strtok (input, "=");
            if (tmp != NULL)
                strncpy (parmvar, tmp, linecpmax);
            else
            {
                if(input[0]!='\0')
                    fprintf (stderr,
                             "WARNING: error in parmfile on line %li with %s\n",
                             counter, input);
                continue;
            }
            tmp = strtok (NULL, "\n");
            if (tmp != NULL)
                strncpy (varvalue, tmp, LINESIZE - 1);
            if (!strncmp (parmvar, "theta", 5))
            {
	      //*buffer = position;
	      //read_theta (options, parmvar, varvalue, buffer);
	      read_startparameter (0, options, parmvar, varvalue, buffer);
                //position = *buffer;
                continue;
            }
            if (!strncmp (parmvar, "migration", 5))
            {
	      //*buffer = position;
	      //read_mig(options, parmvar, varvalue, buffer);
	      read_startparameter (1, options, parmvar, varvalue, buffer);
                //position = *buffer;
                continue;
	    }
            if (!strcmp (parmvar, "rate"))//there is also 'rates'
            {
	      read_startparameter (2, options, parmvar, varvalue, buffer);
	      continue;
	    } 
            if (!strcmp (parmvar, "split"))
            {
	      read_startparameter (3, options, parmvar, varvalue, buffer);
	      continue;
	    }
            if (!strncmp (parmvar, "splitstd", 6))
            {
	      read_startparameter (4, options, parmvar, varvalue, buffer);
	      continue;
	    }
           if (!booleancheck (options, parmvar, varvalue))
            {
                if (!numbercheck (options, parmvar, varvalue))
                {
                    warning ("@Inappropiate entry in parmfile: %s [%s %s] ignored\n",
                             input,parmvar, varvalue);
                }
            }
            position = *buffer;
        }
    }
    if(options->bayes_infer)
      {
	options->schains=0;
      }
#ifdef DEBUG_MPI
    printf("%i> DEBUG_MPI custm=%s custm2=%s\n",myID,options->custm, options->custm2);
    printf("%i> DEBUG_MPI generation/year=%f mutationrate/year=%f\n",myID, options->generation_year, options->mutationrate_year[0]);
#endif

}
#endif

void
print_menu_options (world_fmt * world, option_fmt * options, data_fmt * data)
{ 
  //if (options->numpop < data->numpop)
  //    usererror ("Inconsistency between your Menu/Parmfile and your datafile\n 
  //                 Check the number of populations!\n");
    if (options->progress)
    {
        print_options (stdout, world, options, data);
        if (options->writelog)
            print_options (options->logfile, world, options, data);
    }
}

///
/// prints the datamodel for the parmfile into a buffe
void show_pretty_datamodel(int datamodel, char *buffer)
{
  switch(datamodel)
    {
    case JC69: sprintf(buffer,"JC69"); break;//0
    case K2P: sprintf(buffer,"K2P"); break;//   1
    case F81: sprintf(buffer,"F81"); break;//  2
    case F84: sprintf(buffer,"F84"); break;//  3
    case HKY: sprintf(buffer,"HKY"); break;//  4
    case TN: sprintf(buffer,"TN"); break;//   5
    case GTR: sprintf(buffer,"GTR"); break;//  6
    case SSM: sprintf(buffer,"SSM"); break;//  7
    case MSM: sprintf(buffer,"MSM"); break;//  8
    case BM: sprintf(buffer,"BM"); break;//   9
    case IAM: sprintf(buffer,"IAM"); break;//  10
    default:
       sprintf(buffer,"OTHER"); break;//  -1
    }
}

///
/// prints the value of minimum value of the prior distribution
char * show_priormin(char *tmp, prior_fmt *prior)
{
  if (prior==NULL)
    error("showpriormin got a null pointer");
  else
    {
      sprintf(tmp, "%f ", prior->min);
    }
  return tmp;
}

///
/// prints the value of mean value of the prior distribution
/// or for the gamma proposal the alpha value, and the ,multipler proposal the multiplicator
char * show_priormean(char *tmp, prior_fmt *prior)
{
  sprintf(tmp, "%f ", prior->mean);
  return tmp;
}

/// prints the value of maximum value of the prior distribution
char * show_priormax(char *tmp, prior_fmt *prior)
{
  sprintf(tmp, "%f ", prior->max);
  return tmp;
}

/// prints the alpha value of the gamma prior distribution
char * show_prioralpha(char *tmp, prior_fmt *prior)
{
  sprintf(tmp, "%f ", prior->alpha);
  return tmp;
}

///
/// prints the value of delta value of the prior distribution
/// for some distribution this has no meaning
char * show_priordelta(char *tmp, prior_fmt *prior)
{
  switch(prior->kind)
    {
    case EXPPRIOR:
      sprintf(tmp, "-"); 
      break;
    case GAMMAPRIOR:
      sprintf(tmp, "(a=%1.3f)",prior->alpha);
      break;
    default:
      sprintf(tmp, "%f", prior->delta);
      break;
    }
  return tmp;
}

///
/// prints the bins used for the prior/posterior distribution
char * show_priorbins(char *tmp, prior_fmt *prior)
{
  sprintf(tmp, "%li ", prior->bins);
  return tmp;
}

///
/// prints the updatefreq used for the proposal distribution
char * show_priorupdatefreq(char *tmp, prior_fmt *prior)
{
  sprintf(tmp, "%.4f ", prior->updatefreq);
  return tmp;
}


void fill_printvar_startparam(option_fmt *options, char *(paramtgen[]), char *(parammgen[]))
{
  switch (startguessvalue(options,THETAPRIOR))
    {
    case OWN:
      strcpy (*paramtgen, "from guessed values");
      break;
    case PRIOR:
      strcpy (*paramtgen, "Using a percent value of the prior");
      break;
    case RANDOMPRIOR:
      strcpy (*paramtgen, "RANDOM start value from the prior");
      break;
    default:
      strcpy (*paramtgen, "ERROR");
      break;
    }
  switch (startguessvalue(options,MIGPRIOR))
    {
    case OWN:
      strcpy (*parammgen, "from guessed values");
      break;
    case PRIOR:
      strcpy (*parammgen, "Using a percent value of the prior");
      break;
    case RANDOMPRIOR:
      strcpy (*parammgen, "RANDOM start value from the prior");
      break;
    default:
      strcpy (*parammgen, "ERROR");
      break;
    }
}

void
print_options (FILE * file, world_fmt * world, option_fmt * options,
               data_fmt * data)
{
  //  const char text[8][20]={"All", "Multiplier", "Exponential", "Exp window", "Gamma", "Uniform", "Truncated Normal", "-"};
  boolean is_same = FALSE;
  long i;
  long z = 0;
  int count;
  char s[LINESIZE];
  char mytext[LINESIZE];
  char mytext1[LINESIZE];
  char mytext2[LINESIZE];
  char mytext3[LINESIZE];
  char mytext4[LINESIZE];
  char mytext5[LINESIZE];
  char mytext6[LINESIZE];
  char seedgen[LINESIZE],
    spacer[LINESIZE];
  char *paramtgen, *parammgen;
  long from;
  long to;
  const long numpop = world->numpop;
  const long npp = world->numpop2 + world->bayes->mu + 2 * world->species_model_size;
  char * priorkind = (char *) mycalloc(LINESIZE, sizeof(char));
  paramtgen = (char *) mycalloc(2*LINESIZE,sizeof(char));
  parammgen = paramtgen + LINESIZE;
  if (options->datatype != 'g')
    {
        switch ((short) options->autoseed)
        {
        case AUTO:
            strcpy (seedgen, "with internal timer");
            strcpy (spacer, "  ");
            break;
        case NOAUTOSELF:
            strcpy (seedgen, "from parmfile");
            strcpy (spacer, "      ");
            break;
        case NOAUTO:
            strcpy (seedgen, "from seedfile");
            strcpy (spacer, "      ");
            break;
        default:
            strcpy (seedgen, "ERROR");
            strcpy (spacer, " ");
            break;
        }
	fill_printvar_startparam(options,&paramtgen, &parammgen);
    }
    fprintf (file, "Options in use:\n");
    fprintf (file, "---------------\n\n");
    fprintf (file, "Analysis strategy is BAYESIAN INFERENCE\n");

    switch(options->tri_mlalpha)
      {
      case FIXED:
	sprintf(mytext6,"Mittag-Leffler with alpha=%.2f",options->mlalpha);
	break;
      case NO:
	sprintf(mytext6,"Exponential Distribution");
	break;
      case ESTIMATE:
	sprintf(mytext6,"Mittag-Leffler with alpha=Estimate");
	break;
      }
    fprintf (file, "    - Population size estimation: Theta [%s]\n",mytext6);

    if(world->has_growth)
      fprintf (file, "    - Population growth estimation: Growth [Exponential]\n");
    
    if (world->numpop>1)
      {
	if(world->has_migration)
	  fprintf (file, "    - Geneflow estimation: Migration [%s]\n",mytext6);
	if(world->has_speciation)
	  {
	    switch(world->species_model_dist)
	      {
	      case WEIBULL_DIST:
		strcpy (mytext5, "Weibull Distribution");
		fprintf (file, "    - Divergence estimation: Divergence time [Weibull Distribution with mean and]\n");
		fprintf (file, "                                             [and spread parameter sigma        ]\n");
		break;
	      case NORMAL_DIST:
		strcpy (mytext5, "Normal Distribution");
		fprintf (file, "    - Divergence estimation: Divergence time [Normal Distribution with mean and]\n");
		fprintf (file, "                                             [and standard deviation sigma     ]\n");
		break;
	      case NORMALSHORTCUT_DIST:
		strcpy (mytext5, "Normal Distribution (Shortcut)");
		fprintf (file, "    - Divergence estimation: Divergence time [Normal Distribution [Shortcut] with mean and]\n");
		fprintf (file, "                                             [and standard deviation sigma     ]\n");
		break;
	      case 9:
		error("what is case 9");
		//break;
	      case EXP_DIST:
		strcpy (mytext5, "Exponential Distribution");
		fprintf (file, "    - Divergence estimation: Divergence time [Exponential Distribution with mean]\n");
		break;
	      }
	  }
      }
    fprintf (file, "\nProposal distribution:\n");
    fprintf (file, "Parameter group          Proposal type\n");
    fprintf (file, "-----------------------  -------------------\n");
    fprintf (file, "Population size (Theta) %20s\n", 
	     is_proposaltype(options->slice_sampling[THETAPRIOR]));
    if(numpop > 1)
      {
	fprintf (file, "Migration rate  %7.7s %20s\n",
		 options->usem ? "(M)" : "(xNm)",
		 is_proposaltype(options->slice_sampling[MIGPRIOR]));
	fprintf (file, "Divergence Time (D) %20s\n", 
		 is_proposaltype(options->slice_sampling[SPECIESTIMEPRIOR]));
    fprintf (file, "Divergence time spread (STD) %14s\n", 
	     is_proposaltype(options->slice_sampling[SPECIESSTDPRIOR]));
      }
    if(world->bayes->mu)
      {
	fprintf (file, "Mutation rate modifier  %20s\n",
		 is_proposaltype(options->slice_sampling[RATEPRIOR]));
      }
    if(world->has_growth)
      {
	fprintf (file, "Population growth  %20s\n",
		 is_proposaltype(options->slice_sampling[GROWTHPRIOR]));
      }
    
    fprintf (file, "Genealogy               %20s\n", 
	     "Metropolis-Hastings");    
    fprintf(file,"\n\n");
    if(world->options->has_autotune)
      fprintf (file, "Prior distribution (Proposal-delta will be tuned to acceptance frequency %f):\n",world->options->autotune);
    else
      {
	fprintf (file, "Prior distribution:\n");
      }
    fprintf (file, "Parameter group            Prior type   Minimum    Mean(*)    Maximum    Delta      Bins   Updatefreq\n");
    fprintf (file, "-------------------------  ------------ ---------- ---------- ---------- ---------- ------ -------\n");
    prior_fmt *p = options->bayes_priors;
    long pnum = options->bayes_priors_num;
    i = 0;
    z = 0;
    long d;
    long pa=0;
    long numparam = 0;
    for(i=0;i<pnum;i++)
      {
	if(shortcut(i,world,&pa))
	  {
	    continue;
	  }
	numparam++;
      }
    for(i=0;i<pnum;i++)
      {
	prior_fmt *ptr = &p[i];
	switch(ptr->type)
	  {
	  case THETAPRIOR:
	    sprintf(s, "Population size (Theta_%li)",i+1);
	    break;
	  case MIGPRIOR:
	    pa=0;
	    if(shortcut(i,world,&pa))
	      {
		continue;
	      }
	    m2mm(pa,numpop,&from,&to);
	    sprintf (s, "Migration %li to %li %s", from+1, to+1, options->usem ? "  (M)   " : " (xNm)  ");
	    break;
	  case RATEPRIOR:
	    sprintf (s, "Mutation rate modifier");
	    break;
	  case SPECIESTIMEPRIOR:
	  case SPECIESSTDPRIOR:
	    if (world->species_model==NULL)
	      {
		continue;
	      }
	    from = world->species_model[z/2].from;
	    to = world->species_model[z/2].to;
	    sprintf (s, "Ancestor %li to %li %s",
		     from+1, to+1, (z % 2) == 0 ? "(D_time)" : "(S_time)");
	    z++;
	    break;
	  case GROWTHPRIOR:
	    d=i-npp;
	    while (world->options->growpops[d]==0 && d < world->options->growpops_numalloc)
	      d++;		       
	    sprintf(s, "Population growth (Growth_%li)",d+1);
	    break;

	  }
	is_priorkind(ptr, priorkind);
	fprintf (file, "%s %15s %10.10s %10.10s %10.10s %10.10s %7.7s %5.5f\n",
		 s,
		 priorkind,
		 show_priormin(mytext1, ptr),
		 show_priormean(mytext2, ptr),
		 show_priormax(mytext3, ptr),
		 show_priordelta(mytext4, ptr),
		 show_priorbins(mytext5,ptr),
		 (world->options->choices[1]-world->options->choices[0])/numparam
		 //show_priorupdatefreq(mytext6,ptr)
		 );
      }
    char smoothing1[STRSIZE] = {"--"};
    char smoothing2[STRSIZE] = {"--"};
    switch(options->smoothing[0])
      {
      case BINNING: sprintf(smoothing1,"BINNING"); break;
      case KERNEL: sprintf(smoothing1,"KERNEL SMOOTHING (window=%li)",options->smoothing_window[0]*2+1); break;
      }
    switch(options->smoothing[1])
      {
      case NONE: sprintf(smoothing2,"NO SMOOTHING"); break;
      case AVERAGE: sprintf(smoothing2,"MOVING AVERAGE SMOOTHING (window=%li)",options->smoothing_window[1]*2+1); break;
      case SAVITZKY: sprintf(smoothing2,"SAVITZKY-GOLAY SMOOTHING (window=%li)",options->smoothing_window[1]*2+1); break;
      }
    /* TODO currently only bayes_smooth() is used 

    fprintf (file, "\n\nPosterior distribution:\nParameter values were collected using MCMC, these values\n");
    fprintf (file, "were then used to generate the posterior histograms using %s\n",smoothing1);
    fprintf (file, "and subsequent %s for combination over loci\n", smoothing2);
    //fprintf (file, "Divergence random variable assumes %s\n", mytext7);
    fprintf(file,"\n\n\n");	 
    */
    
    switch (options->datatype)
    {
    case 'a':
      fprintf (file, "Datatype: Allelic data\n");
      fprintf (file, "Missing data is %s\n",options->include_unknown ? "included" : "not included");
      break;
    case 'b':
      fprintf (file, "Datatype: Microsatellite data [Brownian motion]\n");
      fprintf (file, "Missing data is %s\n",options->include_unknown ? "included" : "not included");
      break;
    case 'm':
      if(options->msat_option==SINGLESTEP)
        fprintf (file, "Datatype: Microsatellite data [Singlestep model]\n");
      else
        fprintf (file, "Datatype: Microsatellite data [Multistep model (Tune=%f, P_increase=%f)]\n",options->msat_tuning[0], options->msat_tuning[1]);
      fprintf (file, "Missing data is %s\n",options->include_unknown ? "included" : "not included");
      break;
    case 's':
      fprintf (file, "Datatype: DNA sequence data\n");
      break;
    case 'n':
      fprintf (file, "Datatype: Single nucleotide polymorphism data\n");
      break;
    case 'h':
      fprintf (file, "Datatype: Hapmap Single nucleotide polymorphism data\n");
      break;
      //case 'u':
      //fprintf (file,
      //         "Datatype: Single nucleotide polymorphism data (PANEL)\n");
      //break;
      //case 'f':
      //fprintf (file, "Datatype: Ancestral state method\n");
      //break;
      //case 'g':
      //fprintf (file, "Datatype: Genealogy summary of an older run\n");
      //break;
    }

    count=0;
    is_same = TRUE;
    fprintf (file, "\nInheritance multipliers in use for Thetas (specified # %li)\n",options->inheritance_scalars_numalloc);
    double first = options->inheritance_scalars[0];
    for(i=1;i<options->inheritance_scalars_numalloc;i++)
      {
	//fprintf (file, "@ [%f] %f\n",first,options->inheritance_scalars[i]);
	if (fabs(options->inheritance_scalars[i] -  first) < (double) FLT_EPSILON)
	  is_same = TRUE;
	else
	  {
	    is_same = FALSE;
	    break;
	  }
      }
    if (is_same)
      fprintf (file, "All inheritance multipliers are the same [%f]\n",first);
    else
      {
	fprintf (file, "The inheritance multipliers are:\n");
	for(i=0;i<data->loci;i++)
	  {
	    if(i<options->inheritance_scalars_numalloc)
	      fprintf (file, "%2.2f ", options->inheritance_scalars[i]);
	    else
	      fprintf (file, "%2.2f ", options->inheritance_scalars[options->inheritance_scalars_numalloc-1]);
	    if(count++ > 7)
	      {
		fprintf(file,"\n");
		count=0;
	      }
	    //else 
	    //  {
	    //		count++;
	    //  }
	  }
      }
    //fprintf (file, "\n[Each Theta uses the (true) inheritance scalar of the first locus as a reference]\n\n");

    if(options->randomsubset > 0)
      {
	if(options->randomsubsetseed>0)
	  fprintf (file, "\nData set was subsampled: used a random sample of size: %li\nwith private random number seed %li\n\n", options->randomsubset, options->randomsubsetseed);
	else
	  fprintf (file, "\nData set was subsampled: used a random sample of size: %li\n(no specific random number stream for subset was specified)\n", options->randomsubset);
      }

    if (options->datatype != 'g')
    {
      fprintf (file, "\n%-80s\n", generator);
#ifndef QUASIRANDOM      
      fprintf (file, "Random number seed (%s)%s%20li\n", seedgen, " ",
	       options->saveseed);
#endif
      fprintf (file, "\nStart parameters:\n");
      fprintf(file,  "   First genealogy was started using a %s\n",
	      options->usertree ? "user tree" : (options->randomtree ? "random tree" : 
						   (options->dist ? "tree from userdistances" : 
						    "UPGMA-tree")));
      fprintf(file,  "   Start parameter values were generated\n");
      show_allstartparam(file, options);
      print_arbitrary_migration_table (file, world, options, data);
      print_distance_table (file, world, options, data);
      fprintf(file,"\n");
      // mutation related material
      if (options->gamma)
	{
	  fprintf (file, "Mutation rate among loci is Gamma-distributed\n");
	  fprintf (file, "Initial scale parameter alpha = %f\n",
		   options->alphavalue);
	  if (options->custm[world->numpop2] == 'c')
            fprintf (file, "and is constant [will not be estimated]\n");
	}
      else /*Gamma*/
	{
	  if (options->murates && world->loci > 1)
	    {
	      fprintf (file, "Mutations rate among loci is varying with\n");
	      fprintf (file, "   Rates per locus: ");
	      for (i = 0; i < world->loci - 1; i++)
		{
		  fprintf (file, "%.3f, ", options->mu_rates[i]);
		  if ((i + 1) % 6 == 0)
		    fprintf (file, "\n                    ");
		}
	      fprintf (file, "%.3f\n", options->mu_rates[i]);
	      if(options->murates_fromdata)
		fprintf (file, "[Estimated from the data using the Watterson estimator (ignoring migration)]");
	      else
		fprintf (file, "[User defined]");
	    }
	  else
	    fprintf (file, "Mutation rate is constant %s\n", world->loci > 1 ?
		     "for all loci" : "");
	}
#ifdef UEP
      if (options->uep)
	{
	  fprintf (file, "0/1 polymorphism analysis, with 0/1 data in file:%s\n",
		   options->uepfilename);
	  fprintf (file, "     with forward mutation rate %f*mu\n",
		   options->uepmu);
	  fprintf (file, "     with back    mutation rate %f*mu\n",
		   options->uepnu);
	  fprintf (file, "     with base frequencies \"0\"=%f and \"1\"=%f\n",
		   options->uepfreq0,options->uepfreq1);
	}
#endif /*UEP*/
      if (options->datatype != 'g')
      {
        fprintf (file, "\nMarkov chain settings:\n");
	if(!options->bayes_infer)
	  {
	    fprintf (file, "   Short chains (short-chains):         %20li\n",
		     options->schains);
	    fprintf (file, "      Trees sampled (short-inc*samples):%20li\n",
		     options->sincrement * options->ssteps);
	    fprintf (file, "      Trees recorded (short-sample):    %20li\n",
		     options->ssteps);
	  }
        fprintf (file, "   Long chains (long-chains):           %20li\n",
                 options->lchains);
	if(options->bayes_infer)
	  {
	    if(options->replicate)
	      {
		fprintf (file, "      Steps sampled (inc*samples*rep):  %20li\n",
			 options->lincrement * options->lsteps * options->replicatenum);
		fprintf (file, "      Steps recorded (sample*rep):      %20li\n",
			 options->lsteps*options->replicatenum);
	      }
	    else
	      {
		fprintf (file, "      Steps sampled (long-inc*samples): %20li\n",
			 options->lincrement * options->lsteps);
		fprintf (file, "      Steps recorded (long-sample):     %20li\n",
			 options->lsteps);
	      }
	  }
	else
	  {
	    fprintf (file, "      Trees sampled (long-inc*samples): %20li\n",
		     options->lincrement * options->lsteps);
	    fprintf (file, "      Trees recorded (long-sample):     %20li\n",
		     options->lsteps);
	  }
        if (options->replicate)
        {
            if (options->replicatenum == 0)
                fprintf (file, "   Averaging over long chains\n");
            else
	      {
		if(!options->bayes_infer)
		  fprintf (file, "   Averaging over replicates:           %20li\n",
			   options->replicatenum);
		else
		  fprintf (file, "   Combining over replicates:           %20li\n",
			   options->replicatenum);
	      }
        }
        if (options->heating > 0)
        {
            fprintf (file,
                     "   %s heating scheme\n      %li chains with %s temperatures\n      ",
                     options->adaptiveheat!=NOTADAPTIVE ? (options->adaptiveheat==STANDARD ? "Adaptive_Standard" : "Bounded_adaptive") : "Static", options->heated_chains, options->adaptiveheat ? "start" : "" );
            for (i = 0; i < options->heated_chains - 1; i++)
                fprintf (file, "%5.2f,", options->heat[i]);
	    if(!options->heatedswap_off)
	      fprintf (file, "%5.2f\n      Swapping interval is %li\n",
		       options->heat[i], options->heating_interval);
	    else
	      fprintf (file, "%5.2f\n      No swapping\n",
		       options->heat[i]);
        }
        if (options->movingsteps)
        {
            fprintf (file, "   Forcing at least this percentage of new genealogies:%6.2f\n",
                     (MYREAL) options->acceptfreq);
        }
        if (options->burn_in > 0)
        {
	  switch (options->burnin_autostop)
	    {
	    case 'a':
	      sprintf(mytext,"(var)%10li", options->burn_in);
	      break;
	    case 't':
	      sprintf(mytext,"(acc)%10li", options->burn_in);
	      break;
	    case 'e':
	      sprintf(mytext,"(ESS)%10li", options->burn_in);
	      break;
	    case ' ':
	    default:
	      sprintf(mytext,"%10li", options->burn_in * options->lincrement);
	    }
	  if (options->bayes_infer)
	    fprintf (file, "   Burn-in per replicate (samples*inc): %20.20s\n", mytext);
	  else
	    fprintf (file, "   Burn-in per chain: %35.35s\n", mytext);
        }
        if (options->lcepsilon < LONGCHAINEPSILON)
        {
            fprintf (file, "   Parameter-likelihood epsilon:        %20.5f\n",
                     options->lcepsilon);
        }
    }
    fprintf (file, "\nPrint options:\n");
    if (options->datatype != 'g')
    {
      fprintf (file, "   Data file: %46.46s\n", options->infilename);
      fprintf (file, "   Parameter file: %41.41s\n", options->parmfilename);
      sprintf(mytext,"%s", options->haplotyping ? (options->haplotyping_report ? "YES: report of haplotype probabilities" : "YES: NO report of haplotype probabilities") : "NO");
      fprintf (file, "   Haplotyping is turned on: %31.31s\n", mytext);
#ifdef PRETTY
        fprintf (file, "   Output file (ASCII text): %31.31s\n", options->outfilename);
        fprintf (file, "   Output file (PDF):        %31.31s\n", options->pdfoutfilename);
#else
        fprintf (file, "   Output file: %44.44s\n", options->outfilename);
#endif
        if(options->bayes_infer)
        {
	  if (options->has_bayesfile)
            fprintf (file,   "   Posterior distribution: %33.33s\n", options->bayesfilename);
	  if(options->has_bayesmdimfile)
	    fprintf (file, "   All values of Post.Dist:%33.33s\n", options->bayesmdimfilename);
        }
        fprintf (file, "   Print data: %45.45s\n",
                 options->printdata ? "Yes" : "No");
        switch (options->treeprint)
        {
        case myNONE:
            fprintf (file, "   Print genealogies: %38.38s\n", "No");
            break;
        case ALL:
            fprintf (file, "   Print genealogies: %38.38s\n", "Yes, all");
            break;
        case LASTCHAIN:
            fprintf (file, "   Print genealogies: %32.32s%li\n",
                     "Yes, only those in last chain, every ", options->treeinc);
            break;
        case BEST:
            fprintf (file, "   Print genealogies: %38.38s\n",
                     "Yes, only the best");
            break;
        }
    }

    if (options->mighist)
    {
      if(options->mighist_all)
	sprintf(mytext,"Yes: All events");
      else
	sprintf(mytext,"Yes: Migration events");
      fprintf (file,
	       "   Frequency histogram of events  %26.26s\n", mytext);
      fprintf (file,
	       "   Time of events are saved in file %24.24s\n",
	       options->mighistfilename);
    }
    if (options->mighist && options->skyline)
    {
      sprintf(mytext,"%3s","Yes");
      fprintf (file,
	       "   Histogram of the parameter values through time %10s\n", mytext);
      fprintf (file,
	       "   Parameters through time are saved in file %15.15s\n",
	       options->skylinefilename);
    }
    if(!options->bayes_infer)
      {
	switch (options->profile)
	  {
	  case myNONE:
	    strcpy (mytext, "No");
	    break;
	  case ALL:
	    strcpy (mytext, "Yes, tables and summary");
	    break;
	  case TABLES:
	    strcpy (mytext, "Yes, tables");
	    break;
	  case SUMMARY:
	    strcpy (mytext, "Yes, summary");
	    break;
	  }
	fprintf (file, "   Profile likelihood: %-36.36s\n", mytext);
	if (options->profile != myNONE)
	  {
	    switch (options->profilemethod)
	      {
	      case 'p':
		fprintf (file, "             Percentile method\n");
		break;
	      case 'q':
		fprintf (file, "             Quick method\n");
		break;
	      case 'f':
		fprintf (file, "             Fast method\n");
		break;
	      case 'd':
		fprintf (file, "             Discrete method\n");
		break;
	      case 's':
		fprintf (file, "             Spline method\n");
		break;
	      default:
		fprintf (file, "             UNKOWN method????\n");
		break;
	      }
	    fprintf (file, "             with df=%li and for Theta and %s\n\n\n\n",
		     options->df, options->profileparamtype ? "M=m/mu" : "4Nm");
	  }
      }
    }
    myfree(paramtgen);
    myfree(priorkind);
}
///
/// new start parameter procedures
/// the options->startparam structure is filled and maintained using
/// allocate_startparam()
/// realloc_startparam()
/// fill_startparam(options,key, index) key = /*prior type*/
void allocate_startparam(option_fmt *options, long numpop)
{
  startparam_fmt *startparam = &options->startparam;
  startparam->allocated=TRUE;
  startparam->numpop = numpop;
  startparam->theta = (float *) mycalloc(numpop, sizeof(float));
  startparam->mig = (float *) mycalloc((numpop*(numpop-1)), sizeof(float)); 
  startparam->rate = (float *) mycalloc(numpop, sizeof(float)); 
  startparam->split = (float *) mycalloc(numpop, sizeof(float));
  startparam->splitstd = (float *) mycalloc(numpop, sizeof(float)); 
}
void realloc_startparam(option_fmt *options, long numpop)
{
  startparam_fmt *startparam = &options->startparam;
  if(startparam->allocated)
    {
      startparam->numpop = numpop;
      startparam->theta = (float *) myrealloc(startparam->theta,startparam->numpop * sizeof(float));
      startparam->mig = (float *) myrealloc(startparam->mig,(size_t) (numpop*(numpop-1)) * sizeof(float)); 
      startparam->rate = (float *) myrealloc(startparam->rate, (size_t) numpop * sizeof(float)); 
      startparam->split = (float *) myrealloc(startparam->split, (size_t) numpop * sizeof(float));
      startparam->splitstd = (float *) myrealloc(startparam->splitstd, (size_t) numpop * sizeof(float)); 
    }
}

void fill_startparam(option_fmt *options, short key, long index, float value)
{
  startparam_fmt *startparam = &options->startparam;
  switch(key)
    {
    case THETAPRIOR:
      if (index >= startparam->numtheta)
	startparam->numtheta = index+1;
      startparam->theta[index] = value;
      break;
    case MIGPRIOR:
      if (index >= startparam->nummig)
	startparam->nummig = index+1;
      startparam->mig[index] = value;
      break;
    case RATEPRIOR:
      if (index >= startparam->numrate)
	startparam->numrate = index+1;
      startparam->rate[index] = value;
      break;
    case SPECIESTIMEPRIOR:
      if (index >= startparam->numsplit)
	startparam->numsplit = index+1;
      startparam->split[index] = value;
      break;
    case SPECIESSTDPRIOR:
      if (index >= startparam->numsplitstd)
	startparam->numsplitstd = index+1;
      startparam->splitstd[index] = value;
      break;
    }
}

void add_startparam(option_fmt *options, short key, float value)
{
  startparam_fmt *startparam = &options->startparam;
  switch(key)
    {
    case THETAPRIOR:
      fill_startparam(options,THETAPRIOR,startparam->numtheta,value);
      break;
    case MIGPRIOR:
      fill_startparam(options,MIGPRIOR,startparam->nummig,value);
      break;
    case RATEPRIOR:
      fill_startparam(options,RATEPRIOR,startparam->numrate,value);
      break;
    case SPECIESTIMEPRIOR:
      fill_startparam(options,SPLITPRIOR,startparam->numsplit,value);
      break;
    case SPECIESSTDPRIOR:
      fill_startparam(options,SPLITSTDPRIOR,startparam->numsplitstd,value);
      break;
    }
}


long startguessvalue(option_fmt *options, short key)
{

  switch(key)
    {
    case THETAPRIOR:
      return options->startguess[0][0];
    case MIGPRIOR:
      return options->startguess[1][0];
    case RATEPRIOR:
      return options->startguess[2][0];
    case SPLITPRIOR:
      return options->startguess[3][0];
    case SPLITSTDPRIOR:
      return options->startguess[4][0];
    default:
      usererror("startguessvalue() show not reach this point, no method defined!");
    }
  //return -1;
}

/// for theta
long show_thetaownparam(FILE *file, option_fmt *options, char **temp)
{
  (void) file;
  long pop;
  long tempsize = 0;
  startparam_fmt *startparam = &options->startparam;
  for(pop=0;pop<startparam->numtheta;pop++)
    {
      tempsize += sprintf(*temp+tempsize, "%f ",(double) startparam->theta[pop]);
    }
  return tempsize;
}

/// for mig
long show_migownparam(FILE *file, option_fmt *options, char **temp)
{
  (void) file;
  long pop1, pop2;
  long tempsize = 0;
  long z=0;
  startparam_fmt *startparam = &options->startparam;
  for(pop1=0;pop1<startparam->numpop;pop1++)
    {
      for(pop2=0;pop2<startparam->numpop;pop2++)
	{
	  if (pop1==pop2)
	    tempsize += sprintf(*temp+tempsize, " -- ");
	  else
	    {if (z<startparam->nummig)
		tempsize += sprintf(*temp+tempsize, "%f ",(double) startparam->mig[z]);
	      else
		tempsize += sprintf(*temp+tempsize, "  ? ");
	      z++;
	    }
	}
    }
  return tempsize;
}

/// for rate
long show_rateownparam(FILE *file, option_fmt *options, char **temp)
{
  (void) file;
  long pop;
  long tempsize = 0;
  startparam_fmt *startparam = &options->startparam;
  for(pop=0;pop<startparam->numrate;pop++)
    {
      tempsize += sprintf(*temp+tempsize, "%ff ",(double) startparam->rate[pop]);
    }
  return tempsize;
}

/// for split
long show_splitownparam(FILE *file, option_fmt *options, char **temp)
{
  (void) file;
  long pop;
  long tempsize = 0;
  startparam_fmt *startparam = &options->startparam;
  for(pop=0;pop<startparam->numsplit;pop++)
    {
      tempsize += sprintf(*temp+tempsize, "%ff ",(double) startparam->split[pop]);
    }
  return tempsize;
}
/// for split std
long show_splitstdownparam(FILE *file, option_fmt *options, char **temp)
{
  (void) file;
  long pop;
  long tempsize = 0;
  startparam_fmt *startparam = &options->startparam;
  for(pop=0;pop<startparam->numsplitstd;pop++)
    {
      tempsize += sprintf(*temp+tempsize, "%ff ",(double) startparam->splitstd[pop]);
    }
  return tempsize;
}

void show_startparamtype(twin_fmt *guess, long index, char **temp)
{
  switch(guess[index][0])
    {
    case OWN:
      sprintf(*temp,"User specified values");
      break;
    case PRIOR:
      sprintf(*temp,"Values at %li%% Prior CDF",guess[index][1]);
      break;
    case RANDOMPRIOR:
      sprintf(*temp,"Random values from  Prior PDF");
      break;
    }
}

 void show_allstartparam(FILE *file, option_fmt *options)
 {
   show_startparam(file,options,THETAPRIOR,VERBOSE);
   show_startparam(file,options, MIGPRIOR,VERBOSE);
   show_startparam(file,options,RATEPRIOR,VERBOSE);
   show_startparam(file,options,SPLITPRIOR,VERBOSE);
   show_startparam(file,options,SPLITSTDPRIOR,VERBOSE);
 }

 void show_startparam(FILE *file, option_fmt *options, short key, boolean verbose)
 {
   //startparam_fmt *startparam = &options->startparam;
  twin_fmt * guess = options->startguess;
  //long tempsize=0;
  char *temp = (char * ) mycalloc(LINESIZE,sizeof(char));
  
  switch(key)
    {
    case THETAPRIOR:
      if (verbose)
	show_thetaownparam(file,options, &temp);
      else
	show_startparamtype(guess,0,&temp);
      break;
    case MIGPRIOR:
      if (verbose)
	show_migownparam(file,options, &temp);
      else
	show_startparamtype(guess,1,&temp);
      break;
    case RATEPRIOR:
      if (verbose)
	show_rateownparam(file,options, &temp);
      else
	show_startparamtype(guess,2,&temp);
      break;
    case SPLITPRIOR:
      if (verbose)
	show_splitownparam(file,options, &temp);
      else
	show_startparamtype(guess,3,&temp);
      break;
    case SPLITSTDPRIOR:
      if (verbose)
	show_splitstdownparam(file,options, &temp);
      else
	show_startparamtype(guess,4,&temp);
      break;
    }   
  myfree(temp);
 }


/// \brief sets the start parameters for prior %,  random start parameter setting
/// 
/// Sets the startparameter to a random number derived from a prior distribution
void set_startparam_randomstart(world_fmt *world, option_fmt *options)
{
    long i;
    long pos=0;
    char temp[SUPERLINESIZE];
    if(world->options->progress)
      pos = sprintf(temp,"Start parameter values:");
    for (i = 0; i < world->numpop; i++)
    {
        world->param0[i] = 
options->thetag[0] + RANDUM() * (options->thetag[1]-options->thetag[0]);
	if(world->options->progress)
	  pos += sprintf(temp+pos," %f",world->param0[i]);
    }
    if(world->options->progress)
      FPRINTF(stdout,"[%3i] %s\n",myID, temp);
    if(world->options->writelog)
      FPRINTF(world->options->logfile,"[%3i] %s\n",myID, temp);
}

/// \brief sets the theta parameters for random start parameter setting
/// 
/// Sets the startparameter to a random number derived from a normal distribution
/// with men options->thetag[0] and standard deviation options-.thetag[1]
void set_theta_nrandomstart(world_fmt *world, option_fmt *options)
{
    long i;
    long pos=0;
    char temp[SUPERLINESIZE];
    if(world->options->progress)
      pos = sprintf(temp,"Random start Theta values:");
    for (i = 0; i < world->numpop; i++)
    {
        world->param0[i] = rannor (options->thetag[0], options->thetag[1]);
        while (world->param0[i] < 0)
	  {
            world->param0[i] =
                rannor (options->thetag[0], options->thetag[1]);
	  }
	if(world->options->progress)
	  pos += sprintf(temp+pos," %f",world->param0[i]);
    }
    if(world->options->progress)
      FPRINTF(stdout,"%i> %s\n",myID, temp);
    if(world->options->writelog)
      FPRINTF(world->options->logfile,"%i> %s\n",myID, temp);
}


//set parameters to start values
void set_mystartparams(long i, long numx,  long guess, float *ppp, world_fmt * world, option_fmt * options, prior_fmt * priors)
{
  long ii;
  switch(options->startguess[guess][0])
    {
    case PRIOR:
      world->param0[i] = (double) priors[i].cdf((float) (options->startguess[guess][1]/100.),priors[i].v);
      break;
    case RANDOMPRIOR:
      world->param0[i] = (double) priors[i].random(priors[i].v);
      break;
    case OWN:
      if (i < numx - 1)
	ii = i;
      else
	ii = numx - 1;
      world->param0[i] = (double) ppp[ii];
      if(i < world->numpop)
	{
	  if (world->param0[i] < SMALLEST_THETA)
	    world->param0[i] = SMALLEST_THETA;
	}
      break;
    default:
      world->param0[i] = (double) priors[i].cdf(0.5,priors[i].v);
    }
}

/// \brief sets the start parameters for OWN start parameter setting
/// 
/// Sets the startparameter to a user-defined value
void set_param_fromstartparam(world_fmt *world, option_fmt *options)
{
    long i;
    prior_fmt *priors = options->bayes_priors;
    float *ppp;
    short guess=0;
    long numx=0;
    long b1 = world->numpop;
    long b2 = world->numpop2;
    long b3 = world->bayes->mu ? world->numpop2 : -1;
    long xxx = b3 > 0 ? 1 : 0;
    long b4 = world->numpop2 + xxx;
    long b5 = world->numpop2 + xxx + 2 * world->species_model_size;
    for (i = 0; i < b4; i++)
      {
	if (i < b1)
	  {
	    guess = 0;
	    numx = options->startparam.numtheta;
	    ppp = options->startparam.theta;
	  }
	else if (i < b2)
	  {
	    guess = 1;
	    numx = options->startparam.nummig;
	    ppp = options->startparam.mig;
	  }
	else if (i == b3)
	  {
	    guess = 2;
	    numx = options->startparam.numrate;
	    ppp = options->startparam.rate;
	  }
	else
	  {
	    continue;
	  } 
	set_mystartparams(i, numx, guess, ppp, world, options, priors);
      }
    for (i=b4; i<b5; i+=2)
      {
	numx = options->startparam.numsplit;
	ppp = options->startparam.split;
	set_mystartparams(i, numx, 3, ppp, world, options, priors);
	numx = options->startparam.numsplitstd;
	ppp = options->startparam.splitstd;
	set_mystartparams(i+1, numx, 4, ppp, world, options, priors);
      }
    //if(options->automatic_bins)
    //	{
    //	  priors[i]->bins = MIN(sqrt(options->lsteps), priors[i]->bins);
    //	}
#ifdef DEBUG
    printf("%i> b1=%li\n    b2=numpop=%li\n    b3=numpop2=%li\n    b4=mu=%li\n     b5=%li",myID,b1,b2,b3,b4,b5); 
    for (i=0;i<b5;i++)
      {
	printf("%i> startparam %li: %f\n",myID, i, world->param0[i]);  
      }
#endif
}

/// \brief sets the theta parameters to a start parameter using an FST value
/// 
/// Sets the startparameter to a value from the FST calculation
MYREAL set_paramvalue_data_fst(char datatype,MYREAL val1, MYREAL val2)
{
  if (strchr (SEQUENCETYPES, datatype))
    {
      return val1;
    }
  else
    {
      return val2;
    }
}

void set_theta_fststart(world_fmt *world, option_fmt *options, long locus)
{
    long i;
    MYREAL theta=0.0, mtheta=0.0;
    for (i = 0; i < world->numpop; i++)
    {
      if (world->options->bayes_infer)
	{
	  theta = options->bayes_priors[i].mean;
	  mtheta = options->bayes_priors[i].max;
	}
      else
	{
	  theta = set_paramvalue_data_fst(options->datatype,0.1,10.0);
	  mtheta = 10.0 * theta;
	}
      if (world->fstparam[locus][i] < SMALLEST_THETA)
	{
	  world->param0[i] = theta;
	}
      else
        {
	  if (world->fstparam[locus][i] > mtheta)
	    world->param0[i] = mtheta;
	  else
	    world->param0[i] = world->fstparam[locus][i];
        }
    }
}


/// \brief sets the migration parameters for random start parameter setting
/// 
/// Sets the startparameter to a random number derived from a uniform distribution
/// with minimum options->mg[0] and maximum in options->mg[1]
void set_mig_urandomstart(world_fmt *world, option_fmt *options)
{
    long i;
    long pos=0;
    char temp[SUPERLINESIZE];
    if(world->options->progress)
      pos = sprintf(temp,"Random start M values:");

    for (i = world->numpop; i < world->numpop2; i++)
    {
        world->param0[i] =  options->mg[0] + RANDUM() * (options->mg[1]-options->mg[0]);
	if(world->options->progress)
	  pos += sprintf(temp+pos," %f",world->param0[i]);
    }
    if(world->options->progress)
      FPRINTF(stdout,"%i> %s\n",myID, temp);
    if(world->options->writelog)
      FPRINTF(world->options->logfile,"%i> %s\n",myID, temp);
}

/// \brief sets the migration parameters for random start parameter setting
/// 
/// Sets the startparameter to a random number derived from a normal distribution
/// with men options->mg[0] and standard deviation options->mg[1]
void set_mig_nrandomstart(world_fmt *world, option_fmt *options)
{
    long i;
    long pos=0;
    char temp[SUPERLINESIZE];
    if(world->options->progress)
      pos = sprintf(temp,"Random start M values:");

    for (i = world->numpop; i < world->numpop2; i++)
    {
        world->param0[i] = rannor (options->mg[0], options->mg[1]);
        while (world->param0[i] <= 0.0)
	  {
            world->param0[i] = rannor (options->mg[0], options->mg[1]);
	  }
	if(world->options->progress)
	  pos += sprintf(temp+pos," %f",world->param0[i]);
    }
    if(world->options->progress)
      FPRINTF(stdout,"%i> %s\n",myID, temp);
    if(world->options->writelog)
      FPRINTF(world->options->logfile,"%i> %s\n",myID, temp);
}

/// \brief sets the migration parameters for OWN start parameter setting
/// 
/// Sets the startparameter to a user-defined value
void set_mig_ownstart(world_fmt *world, option_fmt *options, data_fmt *data)
{
    long i;
    long j;
    long ii;
    long iitest;
    long numpop = world->numpop;
    MYREAL tmp;
    for (i = 0; i < numpop; i++)
    {
        for (j = 0; j < numpop; j++)
        {
            if (i == j)
                continue;
            if ((iitest = mm2m (j, i, numpop) - numpop) < options->nummg)
                ii = iitest;
            else
	      {
                ii = options->nummg - 1;
	      }
            //@@if (options->usem)
                tmp = options->mg[ii];
		//@@else
                //@@tmp = options->mg[ii] / world->param0[i];
	    iitest += numpop;
            if (options->geo)
            {
                world->param0[iitest] = data->geo[iitest] * tmp;
            }
            else
            {
                world->param0[iitest] = tmp;
            }
        }
    }
}

/// \brief sets the migration parameters to a start parameter using an FST value
/// 
/// Sets the startparameter to a value from the FST calculation
void set_mig_fststart(world_fmt *world, option_fmt *options, long locus)
{
    long i;
    MYREAL mig=0.0, mmig=0.0;
    for (i = world->numpop; i < world->numpop2; i++)
    {
      if (world->options->bayes_infer)
	{
	  mig = options->bayes_priors[i].mean;
	  mmig = options->bayes_priors[i].max;
	}
      else
	{
	  mig = set_paramvalue_data_fst(options->datatype,100.0,1.0);
	  mmig = 10.0 * mig;
	}
      if (world->fstparam[locus][i] <= 0.0)
	{
	  world->param0[i] = mig;
	}
      else
        {
	  if (world->fstparam[locus][i] > mmig)
	    {
	      world->param0[i] =
                1.0 / world->param0[(i - world->numpop) /
				    (world->numpop)];
	      if (world->param0[i] > mmig)
		world->param0[i] = mig;
            }
            else
	      world->param0[i] = world->fstparam[locus][i];
        }
    }
}


/// \brief set the starting parameters in main structure world from options
///
/// Set the starting parameters in main structure world from options
// - Random starting parameters:'\n'
///   will use an average and a standard deviation to generate a starting parameter
/// - Own starting parameters
/// - Starting parameters derived from FST
/// - check that values from FST are not ridicoulously off and that Bayes start is in min/max bounds
/// Once all parameters are set they are synchronized using the custom-migration matrix
/// that will set some values to zero or to the same value. The synchronization will
/// override user-error, but does provide little control about the outcome
/// error control is only done later when user can check the values in the logfile
void
set_param (world_fmt * world, data_fmt * data, option_fmt * options,
           long locus)
{
  (void) data;
  (void) locus;
  set_param_fromstartparam(world, options);
  // set mutation rate to middle
  if(world->bayes->mu)
      {
	if (!(options->murates && options->murates_fromdata && options->gamma))
	  world->options->mu_rates[world->locus] = 1.0;
      }    
  if(options->bayes_infer)
    {
      bayes_check_and_fix_param(world);
    }
  synchronize_param (world, options);
  if (options->gamma)
    world->options->alphavalue = options->alphavalue;
}


/// \brief synchronize parameters
void
synchronize_param (world_fmt * world, option_fmt * options)
{
    char    type;
    boolean found       = FALSE;
    boolean allsymmig   = FALSE;
    boolean allsamesize = FALSE;
    boolean partmig     = FALSE;
    long    i, ii, jj, len;
    long    z=0, zz=0, zzz=0, ns = 0, ss = 0, ss2 = 0, xs = 0, migm = 0;
    MYREAL  summ;
    MYREAL  nsum = 0;
    len = (long) strlen (world->options->custm2);
    if(world->options->bayes_infer)
      {
	world->bayes->custm2 = world->options->custm2;
      }
    if (len < world->numpop2)
      {
        fillup_custm (len, world, options);
      }
    migm = scan_connect (world->options->custm2, 0L, world->numpop, 'm');
    world->options->tmn = migm; // check for averaging of thetas
    migm = scan_connect (world->options->custm2, world->numpop, world->numpop2, 'm');
    world->options->mmn = migm; // check for averaging of migration rates
    for (i = 0; i < world->numpop2; i++)
    {
        if (!(allsymmig && allsamesize))
        {
            type = world->options->custm2[i];
            switch (type)
            {
	    case 'T':
	    case 'D':
	    case 'x':
            case '*':
                xs++;
                break;
            case 's':  // M is symmetric
                if (i >= world->numpop)
                {
                    z = i;
                    m2mm (z, world->numpop, &ii, &jj);
                    zz = mm2m (jj, ii, world->numpop);
                    world->options->symparam = (twin_fmt *) myrealloc
		      (world->options->symparam, sizeof (twin_fmt) * (size_t) (ss + 2));
                    world->options->symparam[ss][0] = zz;
                    world->options->symparam[ss++][1] = z;
                    world->options->symn = ss;
                    summ = (world->param0[z] + world->param0[zz]) / 2.;
                    world->param0[zz] = world->param0[z] = summ;
                }
                break;
            case 'S': 
                if (i >= world->numpop)
                {
                    z = i;
                    m2mm (z, world->numpop, &ii, &jj);
                    zz = mm2m (jj, ii, world->numpop);
                    zzz = 0;
                    found = FALSE;
                    while (zzz < ss2)
                    {
                        if (world->options->sym2param[zzz][1] == zz)
                            found = TRUE;
                        zzz++;
                    }
                    if (found)
                        break;
                    world->options->sym2param = (quad_fmt *)
                                                 myrealloc (world->options->sym2param,
							    sizeof (quad_fmt) * (size_t) (ss2 + 2));
                    world->options->sym2param[ss2][0] = zz;
                    world->options->sym2param[ss2][1] = z;
                    world->options->sym2param[ss2][2] = ii;
                    world->options->sym2param[ss2++][3] = jj;
                    world->options->sym2n = ss2;
                    summ = (world->param0[z] * world->param0[jj] +
                            world->param0[zz] * world->param0[ii]) / 2.;
                    world->param0[z] = summ / world->param0[jj];
                    world->param0[zz] = summ / world->param0[ii];
                }
                break;
            case 'C':
            case 'c':
                world->options->constparam = (long *) myrealloc
		  (world->options->constparam, sizeof (long) * (size_t) (ns + 2));
                world->options->constparam[ns++] = i;
                world->options->constn = ns;
                break;
            case '0':
	    case 'd':
	    case 't':
                z = i;
                world->param0[z] = 0;
                world->options->zeroparam = (long *) myrealloc
		  (world->options->zeroparam, sizeof (long) * (size_t) (ns + 2));
                world->options->zeroparam[ns++] = i;
                world->options->zeron = ns;
                break;
            case 'm': /*average*/
                summ = 0;
                if (i < world->numpop) /*theta */
                {
                    if (!allsamesize)
                    {
                        nsum = 0;
                        allsamesize = TRUE;
                        for (z = 0; z < world->numpop; z++)
                        {
                            if (world->options->custm2[z] == 'm')
                            {
                                nsum++;
                                summ += world->param0[z];
                            }
                        }
                        summ /= nsum;
                        for (z = 0; z < world->numpop; z++)
                            if (world->options->custm2[z] == 'm')
                                world->param0[z] = summ;
                    }
                }
                else  /* migration */
                {
                    summ = 0;
                    if (!partmig)
                    {
                        nsum = 0;
                        partmig = TRUE;
                        for (z = world->numpop; z < world->numpop2; z++)
                        {
                            if (world->options->custm2[z] == 'm')
                            {
                                nsum++;
                                summ += world->param0[z];
                            }
                        }
                        summ /= nsum;
                        for (z = world->numpop; z < world->numpop2; z++)
                            if (world->options->custm2[z] == 'm')
                                world->param0[z] = summ;
                    }
                }
                break;
            case 'M':
                summ = 0;
                if (i < world->numpop) /*theta */
                {
                    if (!allsamesize)
                    {
                        nsum = 0;
                        allsamesize = TRUE;
                        for (z = 0; z < world->numpop; z++)
                        {
                            if (world->options->custm2[z] == 'm')
                            {
                                nsum++;
                                summ += world->param0[z];
                            }
                        }
                        summ /= nsum;
                        for (z = 0; z < world->numpop; z++)
                            if (world->options->custm2[z] == 'm')
                                world->param0[z] = summ;
                    }
                }
                else  /* migration */
                {
                    summ = 0;
                    if (!partmig)
                    {
                        nsum = 0;
                        partmig = TRUE;
                        for (z = world->numpop; z < world->numpop2; z++)
                        {
                            if (world->options->custm2[z] == 'M')
                            {
                                nsum++;
				m2mm (z, world->numpop, &ii, &jj);
                                summ += world->param0[z] * world->param0[jj];
                            }
                        }
                        summ /= nsum;
                        for (z = world->numpop; z < world->numpop2; z++)
			  {
                            if (world->options->custm2[z] == 'M')
			      {
				m2mm (z, world->numpop, &ii, &jj);
                                world->param0[z] = summ/world->param0[jj];
			      }
			  }
                    }
                }
                break;
            default:
	      break;
	      // do nothing , this assumes that the user has specified multiple groups using
	      // the a,b,c,.... syntax
	      //                warning ("The migration connection matrix is misspecified\nOnly these are allowed\n* s S m M 0 (=zero)\nSupplied value was: %c\ncustom-migration=%s\n", type, world->options->custm);
	      //error("[Program exits]");
            }
        }
    }
    //--------gamma stuff
    if (world->options->gamma)
    {
        if (world->options->custm2[world->numpop2] == 'c')
        {
            world->options->constparam = (long *) myrealloc
	      (world->options->constparam, sizeof (long) * (size_t)(ns + 2));
            world->options->constparam[ns++] = i;
            world->options->constn = ns;
        }
        else
        {
            world->options->custm2[world->numpop2] = '*';
            world->options->custm2[world->numpop2 + 1] = '\0';
        }
    }
}


void
resynchronize_param (world_fmt * world)
{
    char type;
    long i, j = 0, z = 0, zz = 0;
    //long len;
    long ii, jj;
    long ns = 0, ss = 0, ss2 = 0, xs = 0, migm = 0;
    boolean allsymmig = FALSE;
    boolean allsamesize = FALSE;
    boolean partmig = FALSE;
    MYREAL summ;
    MYREAL nsum = 0;
 	//xcode      len = (long) strlen (world->options->custm2);
    world->options->symn = 0;
    world->options->sym2n = 0;
    world->options->zeron = 0;
    world->options->constn = 0;
    world->options->mmn = 0;
    migm = scan_connect (world->options->custm2, 0L, world->numpop, 'm');
    world->options->tmn = migm; // are there mean theta values?
    migm = scan_connect (world->options->custm2, world->numpop, world->numpop2, 'm');
    world->options->mmn = migm;
    for (i = 0; i < world->numpop2; i++)
    {
        if (!(allsymmig && allsamesize))
        {
            type = world->options->custm2[i];
            switch (type)
            {
	    case 'x':
            case '*':
                xs++;
                break;
            case 's':  // M is symmetric
                if (i >= world->numpop)
                {
                    z = i;
                    m2mm (z, world->numpop, &ii, &jj);
                    zz = mm2m (jj, ii, world->numpop);
                    world->options->symparam = (twin_fmt *) myrealloc
		      (world->options->symparam, sizeof (twin_fmt) * (size_t) (ss + 2));
                    world->options->symparam[ss][0] = zz;
                    world->options->symparam[ss++][1] = z;
                    world->options->symn = ss;
                    summ = (world->param0[z] + world->param0[zz]) / 2.;
                    world->param0[zz] = world->param0[z] = summ;
                }
                break;
            case 'S':  // 4Nm is symmetric, not completely
                // implemented yet, -> derivatives.c
                if (i >= world->numpop)
                {
                    z = i;
                    m2mm (z, world->numpop, &ii, &jj);
                    zz = mm2m (jj, ii, world->numpop);
                    world->options->sym2param = (quad_fmt *)
                                                 myrealloc (world->options->sym2param,
							    sizeof (quad_fmt) * (size_t) (ss2 + 2));
                    world->options->sym2param[ss2][0] = zz;
                    world->options->sym2param[ss2][1] = z;
                    world->options->sym2param[ss2][2] = i;
                    world->options->sym2param[ss2++][3] = j;
                    world->options->sym2n = ss2;
                    summ = (world->param0[z] * world->param0[i] +
                            world->param0[zz] * world->param0[j]) / 2.;
                    world->param0[z] = summ / world->param0[i];
                    world->param0[zz] = summ / world->param0[j];
                }
                break;
            case 'C':
            case 'c':
                world->options->constparam = (long *) myrealloc
		  (world->options->constparam, sizeof (long) * (size_t) (ns + 2));
                world->options->constparam[ns++] = i;
                world->options->constn = ns;
                break;
            case '0':
                z = i;
                world->param0[z] = 0;
                world->options->zeroparam = (long *) myrealloc
		  (world->options->zeroparam, sizeof (long) * (size_t) (ns + 2));
                world->options->zeroparam[ns++] = i;
                world->options->zeron = ns;
                break;
            case 'm':
                summ = 0;
                if (i < world->numpop) /*theta */
                {
                    if (!allsamesize)
                    {
                        allsamesize = TRUE;
                        nsum = 0;
                        for (z = 0; z < world->numpop; z++)
                        {
                            if (world->options->custm2[z] == 'm')
                            {
                                nsum++;
                                summ += world->param0[z];
                            }
                        }
                        summ /= nsum;
                        for (z = 0; z < world->numpop; z++)
                            if (world->options->custm2[z] == 'm')
                                world->param0[z] = summ;
                    }
                }
                else  /* migration */
                {
                    summ = 0;
                    if (!partmig)
                    {
                        nsum = 0;
                        partmig = TRUE;
                        for (z = world->numpop; z < world->numpop2; z++)
                        {
                            if (world->options->custm2[z] == 'm')
                            {
                                nsum++;
                                summ += world->param0[z];
                            }
                        }
                        summ /= nsum;
                        for (z = world->numpop; z < world->numpop2; z++)
                            if (world->options->custm2[z] == 'm')
                                world->param0[z] = summ;
                    }
                }
                break;
            default:
                error ("no defaults allowed in resynchronize_param()\n");
            }
        }
    }
    //--------gamma stuff
    if (world->options->gamma)
    {
        if (world->options->custm2[world->numpop2] == 'c')
        {
            world->options->constparam = (long *) myrealloc
	      (world->options->constparam, sizeof (long) * (size_t) (ns + 2));
            world->options->constparam[ns++] = i;
            world->options->constn = ns;
        }
        else
        {
            world->options->custm2[world->numpop2] = '*';
            world->options->custm2[world->numpop2 + 1] = '\0';
        }
    }
}


long
scan_connect (char *custm2, long start, long stop, int check)
{
    long i, summ = 0;
    for (i = start; i < stop; i++)
    {
        
        if (check == custm2[i])
            summ++;
    }
    return summ;
}

void set_partmean_mig (long **mmparam, MYREAL *param, char *custm2, long migm,
                  long numpop2)
{
    long i, z = 0;
    MYREAL summ = 0;
    long start = (long) (sqrt ((MYREAL) numpop2));
    (*mmparam) = (long *) myrealloc ((*mmparam), sizeof (long) * (size_t) (migm + 2));

    for (i = start; i < numpop2; i++)
    {
        if (custm2[i] == 'm')
        {
            summ += param[i];
            (*mmparam)[z++] = i;
        }
    }
    summ /= migm;
    for (i = start; i < numpop2; i++)
    {
        if (custm2[i] == 'm')
            param[i] = summ;
    }
}

void free_options_filenames(option_fmt * options)
{
    myfree( options->parmfilename);
    myfree( options->infilename);
    myfree( options->outfilename);
#ifdef PRETTY
    myfree( options->pdfoutfilename);
#endif
    myfree( options->logfilename);
    myfree( options->sumfilename);
    myfree( options->treefilename);
    myfree( options->utreefilename);
    myfree( options->catfilename);
    myfree( options->weightfilename);
    myfree( options->mighistfilename);
    myfree( options->skylinefilename);
    myfree( options->distfilename);
    myfree( options->divtimefilename);
    myfree( options->geofilename);
    myfree( options->divfilename);
    myfree( options->bootfilename);
    myfree( options->seedfilename);
    
#ifdef UEP
    myfree( options->uepfilename);
#endif
    myfree( options->bayesfilename);
    myfree( options->bayesmdimfilename);
    myfree( options->datefilename);
}

void destroy_startparameters(option_fmt* options)
{
  myfree(options->startparam.theta);
  myfree(options->startparam.mig);
  myfree(options->startparam.rate);
  myfree(options->startparam.split);
  myfree(options->startparam.splitstd);
}

void
destroy_options (option_fmt * options)
{
    myfree(options->thetag);
    myfree(options->mg);
    destroy_startparameters(options);
    myfree(options->mu_rates);
    myfree(options->inheritance_scalars);
    myfree(options->newpops);
    myfree(options->growpops);
    myfree(options->ttratio);
    myfree(options->rate);
    myfree(options->probcat);
    myfree(options->rrate);

    while (--options->lratio->alloccounter >= 0)
    {
        myfree(options->lratio->data[options->lratio->alloccounter].value1);
        myfree(options->lratio->data[options->lratio->alloccounter].value2);
        myfree(options->lratio->data[options->lratio->alloccounter].connect);
    }
    myfree(options->lratio->data);
    myfree(options->lratio);
    myfree(options->custm);
    myfree(options->custm2);
    myfree(options->bayes_priors);
    myfree(options->mutationrate_year);
    free_options_filenames(options);
    //myfree(options);
}

//void
//decide_plot (worldoption_fmt * options, long chain, long chains, char type)
//{
//    if (options->plot && (chain >= chains - 1) && (type == 'l'))
//        options->plotnow = TRUE;
//    else
//        options->plotnow = FALSE;
//}

//void
//set_plot (option_fmt * options)
//{
//    long intervals = options->plotintervals;
//    MYREAL prangex[2];
//    MYREAL prangey[2];
//    if (!options->plot)
//        return;
//    prangex[0] = options->plotrange[0];
//    prangex[1] = options->plotrange[1];
//    prangey[0] = options->plotrange[2];
//    prangey[1] = options->plotrange[3];
//
//
//    options->plotxvalues = (MYREAL *) mycalloc (1, sizeof (MYREAL) * intervals);
//    options->plotyvalues = (MYREAL *) mycalloc (1, sizeof (MYREAL) * intervals);
//    set_plot_values (&options->plotxvalues, prangex, options->plotintervals,
//                     options->plotscale);
//    set_plot_values (&options->plotyvalues, prangey, options->plotintervals,
//                     options->plotscale);
//}
//
//void
//set_plot_values (MYREAL **values, MYREAL plotrange[], long intervals,
//                 int type)
//{
//    long i;
//    MYREAL diff = 0;
//    MYREAL logstart = 0;
//    (*values)[0] = plotrange[0];
//    (*values)[intervals - 1] = plotrange[1];
//    if (type == PLOTSCALELOG)
//    {
//        logstart = log10 ((*values)[0]);
//        diff =
//            (log10 ((*values)[intervals - 1]) - logstart) / (MYREAL) (intervals -
//                    1);
//        for (i = 1; i < intervals - 1; i++)
//        {
//            (*values)[i] = pow (10., (logstart + i * diff));
//        }
//    }
//    else
//    {
//        diff =
//            ((*values)[intervals - 1] - (*values)[0]) / (MYREAL) (intervals - 1);
//        for (i = 1; i < intervals - 1; i++)
//        {
//            (*values)[i] = (*values)[i - 1] + diff;
//        }
//    }
//}

///
/// save all the options into a file (default) called parmfile;
/// uses save_options_buffer()
long save_parmfile (option_fmt * options, world_fmt * world, data_fmt *data)
{
  (void) world;
    FILE *fp; 
    long bufsize;
    long allocbufsize = LONGLINESIZE;
    char *sbuffer;
    char *buffer;
    sbuffer = (char *) mycalloc (allocbufsize, sizeof (char));
    fp = options->parmfile;
    if (fp)
    {
        fclose (fp);
        openfile (&fp, options->parmfilename, "w", NULL);
    }
    else
    {
        openfile (&fp, options->parmfilename, "w", NULL);
    }
    options->parmfile = fp;
    bufsize = save_options_buffer (&sbuffer, &allocbufsize, options, data);
    buffer = sbuffer;
    fprintf (fp, "%s", buffer);
    fflush (fp);
    printf ("\n\n+++ Options were written to file %s in current directory +++\n\n", options->parmfilename);
    myfree(sbuffer);
    return bufsize;
}



/// prints delimiters between sections in the parmfile
void print_parm_delimiter(long *bufsize, char **buffer, long *allocbufsize)
{
    add_to_buffer("################################################################################\n",bufsize,buffer, allocbufsize);
}

/// prints delimiters between sections in the parmfile
void print_parm_smalldelimiter(long *bufsize, char **buffer, long *allocbufsize)
{
	add_to_buffer("#-------------------------------------------------------------------------------\n", bufsize,buffer, allocbufsize);
}

/// prints delimiters between sections in the parmfile
void print_parm_br(long *bufsize, char **buffer, long *allocbufsize)
{
	add_to_buffer("#\n", bufsize,buffer, allocbufsize);
}
/// prints parmfile single comment line
void print_parm_comment(long *bufsize, char **buffer, long *allocbufsize, char message[])
{
    char fp[LINESIZE];
	sprintf(fp,"# %s\n",message);
	add_to_buffer(fp, bufsize,buffer, allocbufsize);
}

/// prints parmfile mutable comment line
void print_parm_mutable_comment(long *bufsize, char **buffer, long *allocbufsize, char string[], ...)
{
    char message[LINESIZE];
	char fp[LINESIZE];
	va_list args;
    va_start (args, string);
    vsprintf (message, string, args);
    va_end (args);	
	sprintf(fp,"# %s\n",message);
	add_to_buffer(fp, bufsize,buffer, allocbufsize);
}
/// prints parmfile mutable option line
void print_parm_mutable(long *bufsize, char **buffer, long *allocbufsize, char string[], ...)
{
    char message[LINESIZE];
	char fp[LINESIZE];
	va_list args;
    va_start (args, string);
    vsprintf (message, string, args);
    va_end (args);	
	sprintf(fp,"%s\n",message);
	add_to_buffer(fp, bufsize,buffer, allocbufsize);
}

/// prints parmfile fixed option line
void print_parm(long *bufsize, char **buffer, long *allocbufsize, char string[])
{
	char fp[LINESIZE];
	sprintf(fp,"%s\n",string);
	add_to_buffer(fp, bufsize,buffer, allocbufsize);
}
/// prints a title of section in the parmfile
void print_parm_title(long *bufsize, char **buffer, long *allocbufsize, char message[])
{
	print_parm_delimiter(bufsize, buffer, allocbufsize);
	print_parm_comment(bufsize, buffer, allocbufsize, message);
	print_parm_delimiter(bufsize, buffer, allocbufsize);
}



/// print the ttratio into the parmfile buffer
void print_parm_ttratio(long *bufsize, char **buffer,  long *allocbufsize, option_fmt *options)
{
  if(options->datamodel != MSM && options->datamodel != SSM && options->datamodel != BM && options->datamodel != IAM)
    {
      if (options->datamodel == TN)
	print_parm_mutable(bufsize, buffer,  allocbufsize, "ttratio=%f %f %f\n", options->ttratio[0],options->ttratio[1],options->ttratio[1]);
      else
	print_parm_mutable(bufsize, buffer,  allocbufsize, "ttratio=%f\n", options->ttratio[0]); 
    }
}

/// print base frequency into parmfile buffer
void print_parm_freqfrom(long *bufsize, char **buffer,  long *allocbufsize, option_fmt * options)
{
  char totalsites[LINESIZE];
  if (options->totalsites>0)
    {
      sprintf(totalsites,"[%li]",options->totalsites);
    }
  else
    {
      sprintf(totalsites,"");
    }
  if (options->freqsfrom)
    print_parm_mutable(bufsize, buffer, allocbufsize, "freqs-from-data=YES %s\n",totalsites);
  else
    print_parm_mutable(bufsize,buffer, allocbufsize, "freqs-from-data=NO:%f,%f, %f, %f %s\n", options->freqa,
		       options->freqc, options->freqg, options->freqt,totalsites);
}

/// print whether one uses categories and in what file they are into parmfile buffer
void print_parm_categs(long *bufsize, char **buffer,  long *allocbufsize, option_fmt * options)
{
    if (options->categs>1)
        print_parm_mutable(bufsize, buffer,  allocbufsize, "categories=%li:%s",options->categs, options->catfilename);
    else
        print_parm(bufsize, buffer, allocbufsize,  "categories=1 #no categories file specified");
}

/// print whether one uses weights and in what file they are into parmfile buffer
void print_parm_weights(long *bufsize, char **buffer, long *allocbufsize, option_fmt * options)
{
    if (options->weights)
      print_parm_mutable(bufsize, buffer, allocbufsize, "weights=YES:%s", options->weightfilename);
    else
        print_parm(bufsize, buffer,  allocbufsize, "weights=NO");
}

/// print whether one uses rate categories, autocorrelation and in what file they are into parmfile buffer
void print_parm_rates(long *bufsize, char **buffer, long *allocbufsize, option_fmt * options)
{
    long i;
    char fp[LINESIZE];
    long lbufsize = 0L;
    char *lbuffer;
    long alloclbufsize;
    // rates
    lbuffer = (char *) mycalloc(LINESIZE,sizeof(char));
    alloclbufsize = LINESIZE;
    for (i = 0; i < options->rcategs; i++)
    {
        sprintf (fp, "%f ", options->rrate[i]);
        add_to_buffer(fp,&lbufsize, &lbuffer, &alloclbufsize);
    }    
    print_parm_mutable(bufsize, buffer, allocbufsize, "rates=%li: %s",options->rcategs, lbuffer);
    lbufsize = 0L;
    //reset the lbuffer
    lbuffer[0] = '\0';
    // probablities
    for (i = 0; i < options->rcategs; i++)
    {
        sprintf (fp, "%f ", options->probcat[i]);
        add_to_buffer(fp,&lbufsize,&lbuffer, &alloclbufsize);
    }
    print_parm_mutable(bufsize, buffer, allocbufsize, "prob-rates=%li: %s",options->rcategs, lbuffer);
    myfree(lbuffer);
    // autocorrelation
    if (!options->autocorr)
      print_parm(bufsize, buffer, allocbufsize, "autocorrelation=NO");
    else
      print_parm_mutable(bufsize, buffer, allocbufsize, "autocorrelation=YES:%f", 1. / options->lambda);
}


/// print data-option into parmfile buffer
void print_parm_datatype(long *bufsize, char **buffer, long * allocbufsize, option_fmt *options)
{
  char datamodeltext[LINESIZE];
    switch (options->datatype)
    {
        case 'a':
	  print_parm(bufsize, buffer, allocbufsize, "datatype=AllelicData");
            print_parm_mutable(bufsize, buffer, allocbufsize,"include-unknown=%s", options->include_unknown ? "YES" : "NO"); 
	  print_parm(bufsize, buffer, allocbufsize, "datamodel=IAM");
            break;
        case 'b':
            print_parm(bufsize, buffer, allocbufsize, "datatype=BrownianMicrosatelliteData");
            print_parm_mutable(bufsize, buffer, allocbufsize, "include-unknown=%s", options->include_unknown ? "YES" : "NO"); 
	  print_parm(bufsize, buffer, allocbufsize, "datamodel=BM (BROWNIAN MOTION)");
            break;
        case 'm':
            print_parm(bufsize, buffer, allocbufsize, "datatype=MicrosatelliteData");
	    if(options->msat_option == SINGLESTEP)
	      {
		print_parm(bufsize, buffer, allocbufsize, "datamodel=SSM (singlestep)");
		print_parm_mutable(bufsize, buffer, allocbufsize, "micro-submodel=%li", options->msat_option);
	      }
	    else
	      {
		print_parm(bufsize, buffer, allocbufsize, "datamodel=MSM (multistep)");
		print_parm_mutable(bufsize, buffer, allocbufsize, "micro-submodel=%li:{%f, %f}", options->msat_option, options->msat_tuning[0], options->msat_tuning[1]);
	      }
	    print_parm_mutable(bufsize, buffer, allocbufsize, "micro-threshold=%li", options->micro_threshold);
            print_parm_mutable(bufsize, buffer, allocbufsize, "include-unknown=%s", options->include_unknown ? "YES" : "NO"); 
            break;
        case 's':
            print_parm(bufsize, buffer, allocbufsize, "datatype=SequenceData");
	    show_pretty_datamodel(options->datamodel, &datamodeltext[0]);
	    //printf("debug: %s\n",datamodeltext);
	    print_parm_mutable(bufsize, buffer, allocbufsize, "datamodel=%s",datamodeltext);
            print_parm_ttratio(bufsize, buffer,  allocbufsize, options);
            print_parm_freqfrom(bufsize, buffer, allocbufsize, options);
	    if(options->has_estimateseqerror)
	      print_parm_mutable(bufsize, buffer, allocbufsize, "seqerror-rate=Estimate:%li",options->seqerrorcombined ? 1 : 4);
	    else
	      print_parm_mutable(bufsize, buffer, allocbufsize, "seqerror-rate={%f,%f,%f,%f}", 
				 options->seqerror[0],options->seqerror[1],options->seqerror[2],options->seqerror[3]);
            print_parm_categs(bufsize, buffer, allocbufsize, options);
            print_parm_rates(bufsize, buffer, allocbufsize, options);
            print_parm_weights(bufsize, buffer, allocbufsize, options);
            print_parm_mutable(bufsize, buffer, allocbufsize, "recover=%s", options->checkpointing ? "YES" : "NO"); 
            print_parm_mutable(bufsize, buffer, allocbufsize, "fast-likelihood=%s", options->fastlike ? "YES" : "NO"); 
            break;
        case 'h':
	  print_parm(bufsize, buffer, allocbufsize, "datatype=HapmapSNPfrequencydata\n");
	  show_pretty_datamodel(options->datamodel, &datamodeltext[0]);
	  //printf("debug: %s\n",datamodeltext);
	  print_parm_mutable(bufsize, buffer, allocbufsize, "datamodel=%s",datamodeltext);
	  
	  print_parm_ttratio(bufsize, buffer, allocbufsize, options);
	  print_parm_freqfrom(bufsize, buffer, allocbufsize, options);
	  if(options->has_estimateseqerror)
	    print_parm_mutable(bufsize, buffer, allocbufsize, "seqerror-rate=Estimate:%li",options->seqerrorcombined ? 1 : 4);
	  else
	    print_parm_mutable(bufsize, buffer, allocbufsize, "seqerror-rate={%f,%f,%f,%f}",
			       options->seqerror[0],options->seqerror[1],options->seqerror[2],options->seqerror[3]);
	  print_parm_categs(bufsize, buffer, allocbufsize, options);
	  print_parm_rates(bufsize, buffer, allocbufsize, options);
	  print_parm_weights(bufsize, buffer, allocbufsize, options);
	  print_parm_mutable(bufsize, buffer, allocbufsize, "recover=%s", options->checkpointing ? "YES" : "NO"); 
	  print_parm_mutable(bufsize, buffer, allocbufsize, "fast-likelihood=%s", options->fastlike ? "YES" : "NO"); 
	  break;
        case 'n':
            print_parm(bufsize, buffer, allocbufsize,"datatype=NucleotidePolymorphismData");
	    show_pretty_datamodel(options->datamodel, &datamodeltext[0]);
	    //printf("debug: %s\n",datamodeltext);
	    print_parm_mutable(bufsize, buffer, allocbufsize, "datamodel=%s",datamodeltext);
	    print_parm_ttratio(bufsize, buffer, allocbufsize, options);
            print_parm_freqfrom(bufsize, buffer, allocbufsize, options);
            if(options->has_estimateseqerror)
              print_parm_mutable(bufsize, buffer, allocbufsize, "seqerror-rate=Estimate:%li",options->seqerrorcombined ? 1 : 4);
            else
              print_parm_mutable(bufsize, buffer, allocbufsize, "seqerror-rate={%f,%f,%f,%f}",
                                 options->seqerror[0],options->seqerror[1],options->seqerror[2],options->seqerror[3]);
            print_parm_categs(bufsize, buffer, allocbufsize, options);
            print_parm_rates(bufsize, buffer, allocbufsize, options);
            print_parm_weights(bufsize, buffer, allocbufsize, options);
            print_parm_mutable(bufsize, buffer, allocbufsize, "recover=%s", options->checkpointing ? "YES" : "NO"); 
            print_parm_mutable(bufsize, buffer, allocbufsize, "fast-likelihood=%s", options->fastlike ? "YES" : "NO"); 
            break;
	    //case 'u':
            //print_parm(bufsize, buffer, allocbufsize,"datatype=UnlinkedSNPData");
            //print_parm_ttratio(bufsize, buffer, allocbufsize, options);
            //print_parm_freqfrom(bufsize, buffer, allocbufsize, options);
            //if(options->has_estimateseqerror)
            //  print_parm_mutable(bufsize, buffer, allocbufsize, "seqerror-rate=Estimate:%li",options->seqerrorcombined ? 1 : 4);
            //else
            //  print_parm_mutable(bufsize, buffer, allocbufsize, "seqerror-rate={%f,%f,%f,%f}",
            //                     options->seqerror[0],options->seqerror[1],options->seqerror[2],options->seqerror[3]);
            //print_parm_categs(bufsize, buffer, allocbufsize, options);
            //print_parm_rates(bufsize, buffer, allocbufsize, options);
            //print_parm_weights(bufsize, buffer, allocbufsize, options);
            //print_parm_mutable(bufsize, buffer, allocbufsize, "recover=%s", options->checkpointing ? "YES" : "NO"); 
            //print_parm_mutable(bufsize, buffer, allocbufsize, "fast-likelihood=%s", options->fastlike ? "YES" : "NO"); 
            //break;
	    //case 'f':
            //print_parm(bufsize, buffer, allocbufsize,"datatype=F-Ancestral state method");
            //print_parm_ttratio(bufsize, buffer, allocbufsize, options);
            //print_parm_freqfrom(bufsize, buffer, allocbufsize, options);
            //if(options->has_estimateseqerror)
            //  print_parm_mutable(bufsize, buffer, allocbufsize, "seqerror-rate=Estimate:%li",options->seqerrorcombined ? 1 : 4);
            //else
            //  print_parm_mutable(bufsize, buffer, allocbufsize, "seqerror-rate={%f,%f,%f,%f}",
            //                     options->seqerror[0],options->seqerror[1],options->seqerror[2],options->seqerror[3]);
            //print_parm_categs(bufsize, buffer, allocbufsize, options);
            //print_parm_rates(bufsize, buffer, allocbufsize, options);
            //print_parm_weights(bufsize, buffer, allocbufsize, options);
            //print_parm_mutable(bufsize, buffer, allocbufsize, "recover=%s", options->checkpointing ? "YES" : "NO"); 
            //print_parm_mutable(bufsize, buffer, allocbufsize, "fast-likelihood=%s", options->fastlike ? "YES" : "NO"); 
            //break;
	    //case 'g':
            //print_parm(bufsize, buffer, allocbufsize,"datatype=GenealogySummaryOlderRun");
            //break;
        default:
            error ("the parmfile-writer contains an error");
    }
}
 
void print_parm_tipdate(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options, data_fmt *data)
{
  (void) data;
  long pos;
  long locus;
  char *input;
  if(options->has_datefile)
  {
    print_parm_mutable(bufsize, buffer, allocbufsize, "tipdate-file=YES:%s", options->datefilename); 
    print_parm_mutable(bufsize, buffer, allocbufsize, "generation-per-year=%f", options->generation_year); 
    input = (char *) mycalloc(options->mutationrate_year_numalloc * 50, sizeof(char));
    pos = sprintf(input, "{%20.20f", options->mutationrate_year[0]);
    for(locus=1; locus < options->mutationrate_year_numalloc; locus++)
      {
	pos += sprintf(input + pos,", %20.20f", options->mutationrate_year[locus]);
      }
    	//xcode   pos += 
    sprintf(input + pos,"}");
    print_parm_mutable(bufsize, buffer, allocbufsize, "mutationrate-per-year=%s", input);
    myfree(input);
  }
}

///
/// print the parmfile entry for the inheritance scalar
void print_parm_inheritance(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options, data_fmt *data)
{
  (void) data;
  long pos;
  long locus;
  char *input;
  input = (char *) mycalloc(options->inheritance_scalars_numalloc * 50, sizeof(char));
  print_parm_br(bufsize, buffer, allocbufsize);
  print_parm_comment(bufsize, buffer, allocbufsize, " inheritance-scalars=<NO|YES:{values for each locus}>");
  print_parm_comment(bufsize, buffer, allocbufsize, "       these values are multiplied with Theta, for example having");
  print_parm_comment(bufsize, buffer, allocbufsize, "       two autosomal and a locus on X- and one on Y-chromosome we would give ");
  print_parm_comment(bufsize, buffer, allocbufsize, "       inheritance-scalars=YES:{1 1 1.33 4.0}");
  print_parm_comment(bufsize, buffer, allocbufsize, "       # this will scale Theta = 4*N*mu");
  print_parm_comment(bufsize, buffer, allocbufsize, "       # We can scale differently:");
  print_parm_comment(bufsize, buffer, allocbufsize, "       inheritance-scalars=YES:{0.25 0.25 0.33 1.00}");
  print_parm_comment(bufsize, buffer, allocbufsize, "       # this will scale Theta = N*mu all is in scale for Y or mtDNA");
  print_parm_comment(bufsize, buffer, allocbufsize, "       the first locus is the reference and the combined estimate has its inheritance-scalar");
  print_parm_comment(bufsize, buffer, allocbufsize, "       if loci are {mtDNA,nucDNA} then {1.0 0.25} --> Theta=N*mu for all");
  if (options->has_inheritance)
    {
      pos = sprintf(input, "inheritance-scalars=YES:{%4.2f", options->inheritance_scalars[0]);
      for(locus=1; locus < options->inheritance_scalars_numalloc; locus++)
	{
	  pos += sprintf(input + pos,", %4.2f", options->inheritance_scalars[locus]);
	}
      sprintf(input + pos,"}");
    }
  else
    pos = sprintf(input, "inheritance-scalars=NO");
  print_parm_mutable(bufsize, buffer, allocbufsize, "%s", input);
  print_parm_br(bufsize, buffer, allocbufsize);
  myfree(input);
}

///
/// print the parmfile entry for the population relabeling
void print_parm_newpops(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options, data_fmt *data)
{
  (void) data;
  long pos;
  long i;
  char *input;
  print_parm_br(bufsize, buffer, allocbufsize);
  print_parm_comment(bufsize, buffer, allocbufsize, "      population-relabel={assignment for each location in the infile}");
  print_parm_comment(bufsize, buffer, allocbufsize, "            example is population-relabel={1 2 2}");
  
  input = (char *) mycalloc(options->newpops_numalloc * 50, sizeof(char));
  pos = sprintf(input, "population-relabel={%li", options->newpops[0]);
  for(i=1; i < options->newpops_numalloc; i++)
    {
      pos += sprintf(input + pos,", %li", options->newpops[i]);
    }
  sprintf(input + pos,"}");
  print_parm_mutable(bufsize, buffer, allocbufsize, "%s", input);
  print_parm_br(bufsize, buffer, allocbufsize);
  myfree(input);
}

///
/// print the parmfile entry for the population growth parameter setting and labeling
/// 0=no growth, other labels are marking growth either individual populations or combinations
/// for example:
/// population-growth={0} # all populations are constant
/// population-growth={1} # all populations are exponentially growing wit hthe same growth rate
/// population-growth={1 2 3} # all 3 populations are growing with individual rates, if there are more than
///                           # than 3 populations, pop 4 etc will grow in lockstep with 3
/// population-growth={1 0 1} # population 1 and 3 grow in lockstep, population 2 is constant in size
void print_parm_growpops(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options, data_fmt *data)
{
  (void) data;
  long pos;
  long i;
  char *input;
  print_parm_br(bufsize, buffer, allocbufsize);
  print_parm_comment(bufsize, buffer, allocbufsize, "      print the parmfile entry for the population growth parameter setting and labeling");
  print_parm_comment(bufsize, buffer, allocbufsize, "      0=no growth, other labels are marking growth either individual populations or combinations");
  print_parm_comment(bufsize, buffer, allocbufsize, "      for example:");
  print_parm_comment(bufsize, buffer, allocbufsize, "      population-growth={0} # all populations are constant");
  print_parm_comment(bufsize, buffer, allocbufsize, "      population-growth={1} # all populations are exponentially growing with the same growth rate");
  print_parm_comment(bufsize, buffer, allocbufsize, "      population-growth={1 2 3} # all 3 populations are growing with individual rates, if there are more than");
  print_parm_comment(bufsize, buffer, allocbufsize, "                                # than 3 populations, pop 4 etc will grow in lockstep with 3");
  print_parm_comment(bufsize, buffer, allocbufsize, "      population-growth={1 0 1} # population 1 and 3 grow in lockstep, population 2 is constant in size.");
  input = (char *) mycalloc(options->growpops_numalloc * 50, sizeof(char));
  pos = sprintf(input, "population-growth={%li", options->growpops[0]);
  for(i=1; i < options->growpops_numalloc; i++)
    {
      pos += sprintf(input + pos,", %li", options->growpops[i]);
    }
  sprintf(input + pos,"}");
  print_parm_mutable(bufsize, buffer, allocbufsize, "%s", input);
  print_parm_br(bufsize, buffer, allocbufsize);
  myfree(input);
}

void print_parm_randomsubset(long * bufsize, char **buffer, long *allocbufsize, option_fmt *options)
{
    print_parm_comment(bufsize, buffer, allocbufsize, " random-subset=number:random_numberseed");
    print_parm_comment(bufsize, buffer, allocbufsize, "       allows to subset the dataset randomly, if number > sample in population");
    print_parm_comment(bufsize, buffer, allocbufsize, "       all samples are taken, if number is smaller then the pop sample is shuffled and");
    print_parm_comment(bufsize, buffer, allocbufsize, "       and the first number samples are taken");
    if(options->randomsubset>0)
    {
	if(options->randomsubsetseed>0)
	  print_parm_mutable(bufsize, buffer, allocbufsize, "random-subset=%li:%li", options->randomsubset,options->randomsubsetseed);
        else
	  print_parm_mutable(bufsize, buffer, allocbufsize, "random-subset=%li", options->randomsubset);
	return;
    }
}
   
void print_parm_usertree(long * bufsize, char **buffer, long *allocbufsize, option_fmt *options)
{
  (void) options;
    print_parm_br(bufsize, buffer, allocbufsize);
    print_parm_comment(bufsize, buffer, allocbufsize, "         usertree=<RANDOM>");
    print_parm_comment(bufsize, buffer, allocbufsize, "               Default is RANDOM,");
    print_parm_comment(bufsize, buffer, allocbufsize, "               currently no other start trees are allowed");
    //if (options->randomtree)
    //{
    print_parm(bufsize, buffer, allocbufsize, "usertree=RANDOMTREE");
    print_parm_br(bufsize, buffer, allocbufsize);
}    

/// print the theta starting parameters
void print_parm_theta(long *bufsize, char ** buffer, long *allocbufsize, option_fmt * options)
{
    long i;
    char fp[LINESIZE];
    char *lbuffer;
    long lbufsize = 0L;
    long alloclbufsize = LINESIZE;
    lbuffer = (char *) mycalloc(alloclbufsize,sizeof(char));
    
    switch (options->numthetag)
      {
      case 0:
        if (strchr ("snupf", options->datatype))
	  {
	    sprintf (fp, "theta=PRIOR:50");
            add_to_buffer(fp, &lbufsize, &lbuffer, &alloclbufsize);
	  }
        else
	  { 
	    sprintf (fp, "theta=PRIOR:50");
            add_to_buffer(fp, &lbufsize, &lbuffer,  &alloclbufsize);
	  }
	break;
      case 1:
            sprintf (fp, "theta=Own:%f",options->thetag[0]);
            add_to_buffer(fp, &lbufsize, &lbuffer, &alloclbufsize);
	    break;
      default:
	sprintf (fp, "theta=Own:{");
	add_to_buffer(fp,&lbufsize, &lbuffer, &alloclbufsize);
	for (i = 0; i < options->numthetag - 1; i++)
	  {
	    sprintf (fp, "%f ", options->thetag[i]);
	    add_to_buffer(fp,&lbufsize,&lbuffer, &alloclbufsize);
	  }
	sprintf (fp, "%f}", options->thetag[i]);
	add_to_buffer(fp,&lbufsize, &lbuffer, &alloclbufsize);
	break;
      }
#ifdef DEBUG
    printf("%i>print_parm_theta>>> [%s]\n",myID,lbuffer);
#endif
    print_parm_mutable(bufsize, buffer, allocbufsize, "%s", lbuffer);
}

/// print M starting parameters to the parmfile buffer
void print_parm_m(long *bufsize, char **buffer, long *allocbufsize, option_fmt * options)
{
    char fp[LINESIZE];
    long i, j, z, num;
    switch (options->nummg)
    {
        case 0:
            sprintf (fp, "migration=PRIOR:10\n");
            add_to_buffer(fp,bufsize,buffer, allocbufsize);
            break;
        case 1:
            sprintf (fp, "migration=Own:%f\n", options->mg[0]);
            add_to_buffer(fp,bufsize,buffer, allocbufsize);
            break;
        default:
            sprintf (fp, "migration=Own:{ ");
            add_to_buffer(fp,bufsize,buffer, allocbufsize);
            z = 0;
            num = (long) (1. + sqrt (4. * (MYREAL) options->nummg + 1.) / 2.);
            for (i = 0; i < num; i++)
            {
                for (j = 0; j < num; j++)
                {
                    if (i == j)
                    {
                        sprintf (fp, "- ");
                        add_to_buffer(fp,bufsize,buffer, allocbufsize);
                    }
                    else
                    {
                        sprintf (fp, "%f ", options->mg[z++]);
                        add_to_buffer(fp,bufsize,buffer, allocbufsize);
                    }
                }
            }
	    sprintf (fp, "}\n");
            add_to_buffer(fp,bufsize,buffer, allocbufsize);
    }
}

/// print split starting parameters to the parmfile buffer
void print_parm_split(long *bufsize, char **buffer, long *allocbufsize, option_fmt * options)
{
    char fp[LINESIZE];
    long i, j, z, num;
    switch (options->numsplitg)
    {
        case 0:
            sprintf (fp, "split=PRIOR:10\n");
            add_to_buffer(fp,bufsize,buffer, allocbufsize);
            break;
        case 1:
            sprintf (fp, "split=Own:%f\n", options->splitg[0]);
            add_to_buffer(fp,bufsize,buffer, allocbufsize);
            break;
        default:
            sprintf (fp, "split=Own:{ ");
            add_to_buffer(fp,bufsize,buffer, allocbufsize);
            z = 0;
            num = (long) (1. + sqrt (4. * (MYREAL) options->numsplitg + 1.) / 2.);
            for (i = 0; i < num; i++)
            {
                for (j = 0; j < num; j++)
		  {
		    if ((uppercase(options->custm[i*num + j])!='D') && (uppercase(options->custm[i*num + j])!='T'))
                    {
                        sprintf (fp, "- ");
                        add_to_buffer(fp,bufsize,buffer, allocbufsize);
                    }
                    else
                    {
                        sprintf (fp, "%f ", options->splitg[z++]);
                        add_to_buffer(fp,bufsize,buffer, allocbufsize);
                    }
                }
            }
	    sprintf (fp, "}\n");
            add_to_buffer(fp,bufsize,buffer, allocbufsize);
    }
}
/// print split starting parameters to the parmfile buffer
void print_parm_splitstd(long *bufsize, char **buffer, long *allocbufsize, option_fmt * options)
{
    char fp[LINESIZE];
    long i, j, z, num;
    switch (options->numsplitg)
    {
        case 0:
            sprintf (fp, "splitstd=PRIOR:10\n");
            add_to_buffer(fp,bufsize,buffer, allocbufsize);
            break;
        case 1:
            sprintf (fp, "splitstd=Own:%f\n", options->splitg[0]);
            add_to_buffer(fp,bufsize,buffer, allocbufsize);
            break;
        default:
            sprintf (fp, "splitstd=Own:{ ");
            add_to_buffer(fp,bufsize,buffer, allocbufsize);
            z = 0;
            num = (long) (1. + sqrt (4. * (MYREAL) options->numsplitg + 1.) / 2.);
            for (i = 0; i < num; i++)
            {
                for (j = 0; j < num; j++)
                {
		  if ((uppercase(options->custm[i*num + j])!='D') && (uppercase(options->custm[i*num + j])!='T'))
                    {
                        sprintf (fp, "- ");
                        add_to_buffer(fp,bufsize,buffer, allocbufsize);
                    }
                    else
                    {
                        sprintf (fp, "%f ", options->splitg[z++]);
                        add_to_buffer(fp,bufsize,buffer, allocbufsize);
                    }
                }
            }
	    sprintf (fp, "}\n");
            add_to_buffer(fp,bufsize,buffer, allocbufsize);
    }
}

void print_parm_heating(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options)
{
    long i;
    char fp[LINESIZE];
    sprintf (fp, "heating=%s", options->heating ? (options->adaptiveheat!=NOTADAPTIVE  ? (options->adaptiveheat==STANDARD ? "ADAPTIVE_standard" : "Bounded_adaptive") : "YES") : "NO\n");
    add_to_buffer(fp,bufsize,buffer, allocbufsize);

    if (options->heating)
    {
        sprintf (fp, ":%li:{%f,", options->heating_interval, options->heat[0]);
        add_to_buffer(fp,bufsize,buffer, allocbufsize);
        
        for (i = 1; i < options->heated_chains - 1; i++)
        {
            sprintf (fp, "%f,", options->heat[i]);
            add_to_buffer(fp,bufsize,buffer, allocbufsize);
        }
        sprintf (fp, "%f}\nheated-swap=%s\n", options->heat[i], (options->heatedswap_off ? "NO" : "YES") );
        add_to_buffer(fp,bufsize,buffer, allocbufsize);
    }
}

void print_parm_haplotyping(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options)
{
    print_parm_br(bufsize, buffer, allocbufsize);
    print_parm_comment(bufsize, buffer, allocbufsize,"   Haplotyping: ");	
    print_parm_comment(bufsize, buffer, allocbufsize,"   Syntax: haplotyping=<YES:<report|no-report> | NO > ");	
    print_parm_comment(bufsize, buffer, allocbufsize,"   If you have diplotype data input two lines per individual marking");
    print_parm_comment(bufsize, buffer, allocbufsize,"   individual names name:other where name must be on both lines, it will");
    print_parm_comment(bufsize, buffer, allocbufsize,"   will be used to build the indiviudal database for the haplotyping.");
    print_parm_comment(bufsize, buffer, allocbufsize,"   If you are not interested in the haplotypes, refrain from reprting the");
    print_parm_comment(bufsize, buffer, allocbufsize,"   the haplotypes because this will be very time consuming and slows the program");
    print_parm_mutable(bufsize, buffer, allocbufsize,"haplotyping=%s", options->haplotyping ? (options->haplotyping_report ? "YES:reporting\n" : "YES:not-reporting") : "NO");
    print_parm_br(bufsize, buffer, allocbufsize);
}

///
/// returns a TRUE when the option is set and sets a filename, returns FALSE when the
/// option is not used and then of course does not set the filename
boolean  set_filename(char *value, char comparison[], char ** filename)
{
    unsigned long len=1;
    char *newvalue;
    char *temp;
    temp = (char *) mycalloc(LINESIZE,sizeof(char));
    upper(value, &temp);
    len = strlen(comparison);
    if (strncmp (temp, comparison, len)==0)
    {
        newvalue = strchr(value,':');
        if(newvalue!=NULL)
	  get_filename (filename, newvalue+1);
        myfree(temp);
        return TRUE;
    }
    myfree(temp);
    return FALSE;
}

/// assumes that an options will be used and takes the filename after the ":"
void  set_filename_only(boolean check, char *value, char ** filename)
{
    char *newvalue;
    if(check)
    {
        newvalue = strchr(value,':');
        if(newvalue != NULL)
	  get_filename(filename, newvalue + 1);
    }
}
    

void   set_parm_prior_values(prior_fmt * prior, char * mytext)
{
  int priortype;
  if (prior !=NULL)
    priortype = prior->kind;
  else
    {
      error("prior is not allocated");
    }
  char tmp1[LINESIZE];
  char tmp2[LINESIZE];
  char tmp3[LINESIZE];
  char tmp4[LINESIZE];
  mytext[0]='\0';
  show_priormin(tmp1, prior);
  strcat(mytext,tmp1);
  switch(priortype)
    {
      //    case SLICE:
      //show_priormax(tmp3, prior, priortype);
      //strcat(mytext,tmp3);
      //break;
      //    case MULTPRIOR:  
      //	  show_priormax(tmp3, prior, priortype);
      //	  show_priordelta(tmp4, prior,priortype);
      //	  strcat(mytext,tmp3);
      //	  strcat(mytext,tmp4);
      //	  break; 
	case EXPPRIOR: 
	  show_priormean(tmp2,prior);
	  show_priormax(tmp3, prior);
	  strcat(mytext,tmp2);
	  strcat(mytext,tmp3);
	  break;
	case WEXPPRIOR:
	  show_priormean(tmp2, prior);
	  show_priormax(tmp3, prior);
	  show_priordelta(tmp4, prior);
	  strcat(mytext,tmp2);
	  strcat(mytext,tmp3);
	  strcat(mytext,tmp4);
	  break;
        case GAMMAPRIOR:
	  show_priormean(tmp2, prior);
	  show_priormax(tmp3, prior);
	  show_prioralpha(tmp4, prior);
	  strcat(mytext,tmp2);	  
	  strcat(mytext,tmp3);
	  strcat(mytext,tmp4);
	  break;
        case NORMALPRIOR:
	  show_priormean(tmp2, prior);
	  show_priormax(tmp3, prior);
	  show_prioralpha(tmp4, prior);
	  strcat(mytext,tmp2);	  
	  strcat(mytext,tmp3);
	  strcat(mytext,tmp4);
	  break;
	case UNIFORMPRIOR:
	default:	  
	  show_priormax(tmp3, prior);
	  show_priordelta(tmp4, prior);
	  strcat(mytext,tmp3);
	  strcat(mytext,tmp4);
	  break;
    }
}


/// \brief returns proposal type string
/// returns a string that shows what proposal type is set
char * show_proposaltype(boolean priorset)
{
  if(priorset)
    {
      return  "SLICE Sampler";
    }
  else
    {
      return "METROPOLIS-HASTINGS Sampler";
    }
}

/// \brief returns priortype sting
/// returns a string that shows what prior distribution is set
char * show_parmpriortype(int priorset)
{
  switch(priorset)
    {
      //    case MULTPRIOR: return  "MULTPRIOR"; 
    case EXPPRIOR: return "EXPPRIOR"; 
    case WEXPPRIOR: return "WEXPPRIOR";
    case GAMMAPRIOR: return "GAMMAPRIOR";
    case UNIFORMPRIOR: return "UNIFORMPRIOR";
    case NORMALPRIOR: return "NORMALPRIOR";
    default: return "UNIFORMPRIOR";
    }
}

///
/// print prior values to buffer for parmfile
void print_parm_proposal(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options)
{
  print_parm_mutable(bufsize, buffer, allocbufsize, "bayes-proposals= THETA %s",
		     show_proposaltype(options->slice_sampling[THETAPRIOR]));

  print_parm_mutable(bufsize, buffer, allocbufsize, "bayes-proposals= MIG %s",
		     show_proposaltype(options->slice_sampling[MIGPRIOR]));
  print_parm_mutable(bufsize, buffer, allocbufsize, "bayes-proposals= DIVERGENCE %s",
		     show_proposaltype(options->slice_sampling[SPECIESTIMEPRIOR]));

  print_parm_mutable(bufsize, buffer, allocbufsize, "bayes-proposals= DIVERGENCESTD %s",
		     show_proposaltype(options->slice_sampling[SPECIESSTDPRIOR]));

  if(options->bayesmurates)
    {
      print_parm_mutable(bufsize, buffer, allocbufsize, "bayes-proposals= RATE %s",
			 show_proposaltype(options->slice_sampling[RATEPRIOR]));
    }
  if(options->growpops_numalloc>0)
    {
      print_parm_mutable(bufsize, buffer, allocbufsize, "bayes-proposals= GROWTH %s",
			 show_proposaltype(options->slice_sampling[GROWTHPRIOR]));
    }
}

char * getpriortype(int kind)
{
  switch(kind)
    {
    case EXPPRIOR:
      return "EXPPRIOR";
    case WEXPPRIOR:
      return "WEXPPRIOR";
    case GAMMAPRIOR:
      return "GAMMAPRIOR";
    case UNIFORMPRIOR:
      return "UNIFORMPRIOR";
    case NORMALPRIOR:  
      return "NORMALPRIOR";
    }
  return "";
}

///
/// print prior values to buffer for parmfile
void print_parm_prior(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options)
{
  char mytext[LINESIZE]="";
  prior_fmt * p = options->bayes_priors;
  long counter = 0;
  long pnum = options->bayes_priors_num;
#ifdef DEBUG
  printf("%i> print_parm_prior(): pnum=%li ",myID,pnum);
#endif
  if (pnum==0)
    {      
      pnum = set_all_default_priors(options);
      p = options->bayes_priors;
    }
#ifdef DEBUG
  printf(" pnum2=%li\n",pnum);
#endif
  for (counter=0; counter < pnum; counter++)
    {
      prior_fmt *ptr = &p[counter];
      set_parm_prior_values(ptr, mytext);
      if (ptr->from == -1 && ptr->to == -1)
	print_parm_mutable(bufsize, buffer, allocbufsize, "bayes-priors= %s * * %s: %s",
			   ptr->ptypename, getpriortype(ptr->kind), mytext);
      else
	print_parm_mutable(bufsize, buffer, allocbufsize, "bayes-priors= %s %li %li %s: %s",
			   ptr->ptypename, ptr->from+1, ptr->to+1, getpriortype(ptr->kind), mytext);
    }
#ifdef DEBUG
  printf("%i>print_parm_prior: %s\n",myID, *buffer);
#endif
}
///
/// print hyperprior values to buffer for parmfile
void print_parm_hyperprior(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options)
{
  char tmp[LINESIZE];
  if (options->hyperprior)
    sprintf(tmp,"YES:%li:%f:%f",options->hyperinterval,options->hyperfactormean,options->hyperfactoralpha);
  else
    sprintf(tmp,"NO");
  print_parm_comment(bufsize, buffer, allocbufsize, " Hyper-prior for all parameters");
  print_parm_comment(bufsize, buffer, allocbufsize, " The parameter of the prior is drawn from a Gamma distribution with mean and alpha");
  print_parm_comment(bufsize, buffer, allocbufsize, " for example:");
  print_parm_comment(bufsize, buffer, allocbufsize, "   bayes-hyperprior=YES:10000:1.0:5.0");
  print_parm_comment(bufsize, buffer, allocbufsize, " uses a hyper prior with the mean of the specified prior");
  print_parm_comment(bufsize, buffer, allocbufsize, " and and alpha so that this specifies ~Normal");
  print_parm_br(bufsize, buffer, allocbufsize);

  print_parm_mutable(bufsize, buffer, allocbufsize, "bayes-hyperpriors=%s",tmp);
}

//print parmfile smoothing
void print_parm_smoothing(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options)
{
  char tmp[LINESIZE];
  //size_t tmpsize=0;
  if (options->smoothing[0]==BINNING)
    //tmpsize =
    sprintf(tmp,"smoothing={%s,","BIN");
  else
    //tmpsize =
    sprintf(tmp,"smoothing={%s:%li,","KERNEL",options->smoothing_window[0]);
  if (options->smoothing[1]==NONE)
    print_parm_mutable(bufsize, buffer, allocbufsize, "%s%s}",tmp,"NONE");
  else
    {
      if (options->smoothing[1]==AVERAGE)
	print_parm_mutable(bufsize, buffer, allocbufsize, "%s%s:%li}",tmp, "AVERAGE",options->smoothing_window[1]);
      else
	print_parm_mutable(bufsize, buffer, allocbufsize, "%s%s:%li}",tmp, "SAVITZKY",options->smoothing_window[1]);
    }
}

///copy mu_rates into world->options
void set_meanmu(worldoption_fmt * wopt, option_fmt * options, long loci)
{
  long i;
  long n = 0;
  MYREAL last;

  if(options->mutationrate_year_numalloc == 0)
    {
      last = 1.0;
    }
  else
    {
      last = 0.0;
    }

  wopt->meanmu = (MYREAL *) mycalloc(loci, sizeof(MYREAL));
  for(i=0; i < loci; i++)
    {
      wopt->meanmu[i] = 1.0;
    }
  for(i=0;i< options->mutationrate_year_numalloc;i++)
    {
      n++;
      if (options->mutationrate_year[i]>0.0)
	wopt->meanmu[i] = options->mutationrate_year[i];
      last += (wopt->meanmu[i] - last ) / n;
    }

  for(i = options->mutationrate_year_numalloc; i < loci;i++)
      {
	wopt->meanmu[i] = last;
      }
}

///save option->mu_rates into a buffer
long save_mu_rates_buffer (char **buffer, long *allocbufsize, option_fmt * options)
{
    long i;
    long bufsize = 0;
    print_parm_mutable(&bufsize, buffer, allocbufsize, "%li ",options->muloci);
    for(i=0; i < options->muloci; i++)
      {
	print_parm_mutable(&bufsize, buffer, allocbufsize, "%f ",options->mu_rates[i]);
      }
    //    fprintf(stdout,"muloci=%li, muratebuffer:>%s<\n",options->muloci, *buffer, allocbufsize);
    return bufsize;
}

///save option->inheritance_scalars into a buffer
long save_inheritance_scalars_buffer (char **buffer, long *allocbufsize, option_fmt * options)
{
    long i;
    long bufsize = 0;
    print_parm_mutable(&bufsize, buffer, allocbufsize, "%li ",options->inheritance_scalars_numalloc);
    for(i=0; i < options->inheritance_scalars_numalloc; i++)
      {
	print_parm_mutable(&bufsize, buffer, allocbufsize, "%f ",options->inheritance_scalars[i]);
      }
    return bufsize;
}
///save option->newpops into a buffer
long save_newpops_buffer (char **buffer, long *allocbufsize, option_fmt * options)
{
    long i;
    long bufsize = 0;
    print_parm_mutable(&bufsize, buffer, allocbufsize, "%li ",options->newpops_numalloc);
    for(i=0; i < options->newpops_numalloc; i++)
      {
	print_parm_mutable(&bufsize, buffer, allocbufsize, "%f ",options->newpops[i]);
      }
    return bufsize;
}

///save option->growpops into a buffer
long save_growpops_buffer (char **buffer, long *allocbufsize, option_fmt * options)
{
    long i;
    long bufsize = 0;
    print_parm_mutable(&bufsize, buffer, allocbufsize, "%li ",options->growpops_numalloc);
    for(i=0; i < options->growpops_numalloc; i++)
      {
	print_parm_mutable(&bufsize, buffer, allocbufsize, "%f ",options->growpops[i]);
      }
    return bufsize;
}

///save options into a buffer
long save_options_buffer (char **buffer, long *allocbufsize, option_fmt * options, data_fmt *data)
{
  long i;
  char mytext[LINESIZE]="";
  char nowstr[LINESIZE] = "----";
  char fp[LINESIZE];
  long bufsize = 0;
  get_time (nowstr, "%c");
  if(options->bayes_infer)
    options->schains=0;
  // header for parmfile
  print_parm_delimiter(&bufsize, buffer, allocbufsize);	
  print_parm_mutable_comment(&bufsize, buffer, allocbufsize,  "Parmfile for Migrate %s-%s", MIGRATEVERSION,MIGRATESUBVERSION);
  print_parm_comment(&bufsize, buffer, allocbufsize,  "DO NOT REMOVE THE FIRST TWO LINES!");
  print_parm_br(&bufsize, buffer, allocbufsize);
  print_parm_comment(&bufsize, buffer, allocbufsize, "generated automatically on");
  print_parm_mutable_comment(&bufsize, buffer, allocbufsize, "%s", nowstr);
  print_parm_br(&bufsize, buffer, allocbufsize);
  print_parm_comment(&bufsize, buffer, allocbufsize, "please report problems to Peter Beerli");
  print_parm_comment(&bufsize, buffer, allocbufsize, " email: beerli@fsu.edu");
  print_parm_comment(&bufsize, buffer, allocbufsize," http://popgen.sc.fsu.edu/migrate.html");
  print_parm_delimiter(&bufsize, buffer, allocbufsize);
  // general options
  print_parm_br(&bufsize, buffer, allocbufsize);
  print_parm_title(&bufsize, buffer, allocbufsize, "General options");
  print_parm_br(&bufsize, buffer, allocbufsize);
  print_parm_comment(&bufsize, buffer, allocbufsize,"Interactive or batch job usage");
  print_parm_comment(&bufsize, buffer, allocbufsize,"  Syntax: menu= < YES | NO > ");
  print_parm_comment(&bufsize, buffer, allocbufsize,"For batch runs it needs to be set to NO");
  print_parm_mutable(&bufsize, buffer, allocbufsize,"menu=%s", options->menu ? "YES" : "NO ");
  print_parm_br(&bufsize, buffer, allocbufsize);
  print_parm_comment(&bufsize, buffer, allocbufsize,"Specification of length of names of individuals (default=10)");	
  print_parm_comment(&bufsize, buffer, allocbufsize,"   Syntax: nmlength=<INTEGER between 0 .. 30>");	
  print_parm_mutable(&bufsize, buffer, allocbufsize,"nmlength=%li", options->nmlength);
  print_parm_br(&bufsize, buffer, allocbufsize);
  print_parm_br(&bufsize, buffer, allocbufsize);
  print_parm_comment(&bufsize, buffer, allocbufsize,"If you long runs fails because of timelimits");
  print_parm_comment(&bufsize, buffer, allocbufsize,"AND you have specified bayes-allfile=YES....");
  print_parm_comment(&bufsize, buffer, allocbufsize,"then you can recover and continue");
  print_parm_comment(&bufsize, buffer, allocbufsize,"by setting recover to yes, default is NO");
  print_parm_comment(&bufsize, buffer, allocbufsize,"This option fails if your bayesallfile is complete!");
  print_parm_comment(&bufsize, buffer, allocbufsize,"At the moment, I have no experience about failure rate etc. [March 30 2015]");
  print_parm_comment(&bufsize, buffer, allocbufsize, "         recover=<YES | NO> ");
  print_parm_br(&bufsize, buffer, allocbufsize);
  print_parm_title(&bufsize, buffer, allocbufsize, "Data options");
  print_parm_br(&bufsize, buffer, allocbufsize);
  print_parm_comment(&bufsize, buffer, allocbufsize,"Several different main datatypes are possible:");
  //NOT YETprint_parm_comment(&bufsize, buffer, allocbufsize,"* : this datatype is defined in the infile, this can allow a mixture of datatypes");
  //print_parm_comment(&bufsize, buffer, allocbufsize," ");
  print_parm_comment(&bufsize, buffer, allocbufsize,"INFINITE ALLELE: usable for electrophoretic markers,");
  print_parm_comment(&bufsize, buffer, allocbufsize,"                 other markers with unknown mutation model");
  print_parm_comment(&bufsize, buffer, allocbufsize,"STEPWISE MUTATION: usable for microsatellite data or");
  print_parm_comment(&bufsize, buffer, allocbufsize,"                 other markers with stepwise change");
  print_parm_comment(&bufsize, buffer, allocbufsize,"                 from one allele to another");
  print_parm_comment(&bufsize, buffer, allocbufsize,"                 [singlestep versus multistep model, see micro-submodel option]");
  print_parm_comment(&bufsize, buffer, allocbufsize,"FINITE SITES MUTATION: standard DNA/RNA sequence mutation");
  print_parm_comment(&bufsize, buffer, allocbufsize,"                 model, usable for DNA or RNA contiguous");
  print_parm_comment(&bufsize, buffer, allocbufsize,"                 sequences or variable sites only (SNP)");
  //print_parm_comment(&bufsize, buffer, allocbufsize,"GENEALOGY SUMMARY: reanalyzing an old migrate run");
  //print_parm_br(&bufsize, buffer, allocbufsize);
  print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
  print_parm_comment(&bufsize, buffer, allocbufsize,"INFINITE ALLELE");
  print_parm_comment(&bufsize, buffer, allocbufsize," Syntax: datatype=ALLELICDATA ");
  print_parm_comment(&bufsize, buffer, allocbufsize,"         include-unknown=<YES | NO> with YES unknown alleles");
  print_parm_comment(&bufsize, buffer, allocbufsize,"               are included into analysis, NO is the default");
  print_parm_br(&bufsize, buffer, allocbufsize);
  print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
  print_parm_br(&bufsize, buffer, allocbufsize);
  print_parm_comment(&bufsize, buffer, allocbufsize,"STEPWISE MUTATION");
  print_parm_comment(&bufsize, buffer, allocbufsize," Syntax: datatype=<MICROSATELLITEDATA | BROWNIANDATA");
  print_parm_comment(&bufsize, buffer, allocbufsize,"               MICRO specifies the standard stepwise mutation");
  print_parm_comment(&bufsize, buffer, allocbufsize,"               model, the BROWNIAN is an approximation to this");
  print_parm_comment(&bufsize, buffer, allocbufsize,"         FOR 99% of all cases BROWNIAN WILL BE FASTEST AND 'BEST'");
	print_parm_comment(&bufsize, buffer, allocbufsize,"         for MICRO several suboptions are available, again use BROWNIAN! ");
	print_parm_comment(&bufsize, buffer, allocbufsize,"         micro-submodel=<1|2:{tune,pinc}>");
	print_parm_comment(&bufsize, buffer, allocbufsize,"                1 means singlestep mutation model (this is the default and the standard");
	print_parm_comment(&bufsize, buffer, allocbufsize,"                2 is the Multistep model (see Watkins 2007 TPB, section 4.2) it needs");
	print_parm_comment(&bufsize, buffer, allocbufsize,"                  two parameters: tune specifies how close the model is to a singlestep model");
	print_parm_comment(&bufsize, buffer, allocbufsize,"                  so tune=0 --> singlestep, tune=1 --> infinite allele model;");
	print_parm_comment(&bufsize, buffer, allocbufsize,"                  the second parameter defines the probability that the repeat number");
	print_parm_comment(&bufsize, buffer, allocbufsize,"                  is increasing, this value cannot be larger than 0.666, I suggest 0.5.");
	print_parm_comment(&bufsize, buffer, allocbufsize,"                  Example: micro-submodel=2:{0.5,0.5}");
	print_parm_comment(&bufsize, buffer, allocbufsize,"         micro-threshold=<INTEGER> Default is 10 [MICRO only, NEEDS TO BE EVEN!],");
	print_parm_comment(&bufsize, buffer, allocbufsize,"               smaller values speed up analysis, but might also");
	print_parm_comment(&bufsize, buffer, allocbufsize,"               crash, large values slow down analysis considerably.");
	print_parm_comment(&bufsize, buffer, allocbufsize,"               Change this value only when you suspect that your");
	print_parm_comment(&bufsize, buffer, allocbufsize,"               data has huge gaps in repeat length.");
	print_parm_comment(&bufsize, buffer, allocbufsize,"         include-unknown=<YES | NO> with YES unknown alleles");
	print_parm_comment(&bufsize, buffer, allocbufsize,"               are included into analysis, NO is the default");
	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_comment(&bufsize, buffer, allocbufsize,"FINITE SITES MUTATION");
	print_parm_comment(&bufsize, buffer, allocbufsize," Syntax: datatype=<SEQUENCEDATA | NUCLEOTIDE");
	print_parm_comment(&bufsize, buffer, allocbufsize,"        SEQENCEDATA: typical linked stretches of DNA, for example mtDNA");
    print_parm_comment(&bufsize, buffer, allocbufsize,"        NUCLEOTIDE: linked DNA stretches, all invariable sites removed");
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize,"         freqs-from-data=<YES | NO: freq(A), freq(C), freq(G), freq(T)>");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               calculate the prior base frequencies from the data,");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               or specify the frequencies");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         ttratio=<RATIO1 RATIO2 ....> Default is 2.0,");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               ratio between transitions and transversions.");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         seqerror-rate=<{VALUE,VALUE,VALUE,VALUE}|Estimate:1|4> Default is 0.0, typical values for ABI 3700 ");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               sequencers after base calling are around 0.001 (1/650)");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         categories=<VALUE:CATFILE> The categories are integers or letters");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               specified in file called CATFILE, this assumes that all");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               sites belong to known categories, this can be used to");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               weight third positions etc.");
	print_parm_comment(&bufsize, buffer, allocbufsize, "         rates=<VALUE1 VALUE2 ...> the rates are specified arbitrarily or");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               then are from a Gamma distribution with alpha=x, currently");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                the alpha value gets lost and is not recorded in the parmfile");
	print_parm_comment(&bufsize, buffer, allocbufsize, "         prob-rates=<RATE2 RATE1 ... > These rates can be arbitrary or ");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               generated with gamma-deviated rates and then are derived");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               using Laguerre's quadrature, this should get better");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               results than equal probability methods.");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         autocorrelation=<NO | YES:VALUE> Default is NO");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               autocorrelation makes only sense with rates,");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               VALUE should be >1.0");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         weights=<NO | YES:WEIGHTFILE> The weights are specified");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               in file called WEIGHTFILE, this assumes that all sites");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               belong to known weights, this can be used to weight");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               portions of the sequence etc.");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         fast-likelihood=<YES | NO> Default is NO");
    //	print_parm_comment(&bufsize, buffer, allocbufsize, "               have many hundred individuals and get strange errors");
    //	print_parm_comment(&bufsize, buffer, allocbufsize, "               during a run, NO is scaling the conditional likelihood");
    //print_parm_comment(&bufsize, buffer, allocbufsize, "               so that very small values are >0.00000");
    print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
    print_parm_comment(&bufsize, buffer, allocbufsize, "         population-relabel={assignment for each location in the infile}");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               example is population-relabel={1 2 2}");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         random-subset=number<:seed>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               allows to subset the dataset randomly, if number > sample in population");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               all samples are taken, if number is smaller then the pop sample is shuffled and");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               and the first number samples are taken.");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               the random number seed guarantees that the");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               same subset is chosen in different runs");    
    print_parm_datatype(&bufsize, buffer, allocbufsize, options);
    print_parm_tipdate(&bufsize, buffer, allocbufsize, options, data);
    print_parm_inheritance(&bufsize, buffer, allocbufsize, options, data);
    print_parm_haplotyping(&bufsize, buffer, allocbufsize, options);
    print_parm_newpops(&bufsize, buffer, allocbufsize, options, data);
    print_parm_randomsubset(&bufsize, buffer, allocbufsize, options);
    //print_parm_usertree(&bufsize, buffer, allocbufsize, options);

#ifdef UEP
    // unique event polymorphisms
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_title(&bufsize, buffer, allocbufsize,  "Unique event polymorphism options");
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax uep=<NO | YES:UEPFILE >");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               Default is NO, with YES the user needs to ");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               give a file for each individual (same order as in datafile");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               the value: indiviudal<10 characters> uep-state<0|1>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         uep-rates= mu[0->1] : nu[0->1] mutation rate between the two UEP states");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         uep-bases= FREQ1 : FREQ2   prior frequency of the two alleles");
    print_parm_br(&bufsize, buffer, allocbufsize);

    if (options->uep)
      {
	print_parm_mutable(&bufsize, buffer, allocbufsize, "uep=YES:%s", options->uepfilename);
	print_parm_mutable(&bufsize, buffer, allocbufsize, "uep-rates=%f:%f", options->uepmu, options->uepnu);
        print_parm_mutable(&bufsize, buffer, allocbufsize, "uep-bases=%f:%f", options->uepfreq0, options->uepfreq1);
      }
    else
      {
	print_parm(&bufsize, buffer, allocbufsize, "uep=NO");
      }
#endif
    //input and output options
	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_title(&bufsize, buffer, allocbufsize, "Input options");
	print_parm_br(&bufsize, buffer, allocbufsize);
    
    print_parm_comment(&bufsize, buffer, allocbufsize, "input file location");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax infile=FILEPATH");
    print_parm_mutable(&bufsize, buffer, allocbufsize, "infile=%s", options->infilename);
    print_parm_br(&bufsize, buffer, allocbufsize);
    if(options->prioralone)
      {
	print_parm_mutable(&bufsize, buffer, allocbufsize, "NODATA=Yes");
	print_parm_br(&bufsize, buffer, allocbufsize);
      }
    print_parm_comment(&bufsize, buffer, allocbufsize, "Random number seed specification");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax random-seed=<AUTO | OWN:< seedfile | value >");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     AUTO           uses computer system clock to generate seed");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     OWN:seedfile   uses file seedfile with random number seed");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     OWN:value      uses number value for seed");
    switch (options->autoseed)
      {
      case NOAUTO:
        print_parm_mutable(&bufsize, buffer, allocbufsize, "random-seed=OWN:%s",options->seedfilename);
        break;
      case AUTO:
        print_parm_mutable(&bufsize, buffer, allocbufsize, "random-seed=AUTO #OWN:%li", options->inseed);
        break;
      case NOAUTOSELF:
        print_parm_mutable(&bufsize, buffer, allocbufsize, "random-seed=OWN:%li", options->inseed);
        break;
      default:
        error ("error in writing parmfile start seed method unknown");
      }
    print_parm_br(&bufsize, buffer, allocbufsize);
	
    print_parm_comment(&bufsize, buffer, allocbufsize, "Specify the title of the run, will be overridden by title in datafile");
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax: title=title text [up to 80 characters]");
    print_parm_mutable(&bufsize, buffer, allocbufsize, "title=%s", options->title);
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_title(&bufsize, buffer, allocbufsize, "Output options");
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "Progress report to the window where the program was started");
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax: progress=<NO | YES>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         NO       nothing is printed to the console");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         YES      some messages about progress are reported [default]");
    //    print_parm_comment(&bufsize, buffer, allocbufsize, "         VERBOSE  more messages are reported to console");    
    print_parm_mutable(&bufsize, buffer, allocbufsize, "progress=%s",
             options->progress ? (options->verbose ? "YES" : "YES") : "NO ");
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
    print_parm_br(&bufsize, buffer, allocbufsize);
    
    print_parm_comment(&bufsize, buffer, allocbufsize, "Recording messages to screen and into logfile");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax logfile=<NO | YES:logfilename>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "      NONE     no recording of progress");
    print_parm_comment(&bufsize, buffer, allocbufsize, "      logfilename  path to logfile");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  [this is usually a bad choice! Because in batch jobs]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "   the system logs the screen, on macs you can copy the screen,");
    print_parm_comment(&bufsize, buffer, allocbufsize, "   only on windows you may need this; this slows down the run a little");
    if(options->writelog)
      print_parm_mutable(&bufsize, buffer, allocbufsize, "logfile=YES:%s",  options->logfilename);
    else
      print_parm(&bufsize, buffer, allocbufsize, "logfile=NO");
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
    print_parm_br(&bufsize, buffer, allocbufsize);
    
    print_parm_comment(&bufsize, buffer, allocbufsize, "Print the data as read into the program");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax print-data=<NO | YES>");
    print_parm_mutable(&bufsize, buffer, allocbufsize, "print-data=%s", options->printdata ? "YES" : "NO");
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
    print_parm_br(&bufsize, buffer, allocbufsize);
    
    print_parm_comment(&bufsize, buffer, allocbufsize, "Print output to file [default is outfile]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax outfile=outfilename");
    print_parm_mutable(&bufsize, buffer, allocbufsize, "outfile=%s", options->outfilename);
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
    print_parm_br(&bufsize, buffer, allocbufsize);

#ifdef PRETTY
    print_parm_comment(&bufsize, buffer, allocbufsize, "Print output to a PDF file [default is outfile.pdf]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax pdf-outfile=outfilename.pdf");
    print_parm_mutable(&bufsize, buffer, allocbufsize, "pdf-outfile=%s", options->pdfoutfilename);
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
    print_parm_br(&bufsize, buffer, allocbufsize);

    print_parm_comment(&bufsize, buffer, allocbufsize, "Reduce size of PDF -- print/plot only summaries [default is NO]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax pdf-terse=<NO|YES>");
    print_parm_mutable(&bufsize, buffer, allocbufsize, "pdf-terse=%s", options->tersepdf ? "YES" : "NO");
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
    print_parm_br(&bufsize, buffer, allocbufsize);

#endif
#ifdef NEWVERSION
    print_parm_comment(&bufsize, buffer, allocbufsize, "Use an alternative to exponential distribution [mittag-leffler]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax mittag-leffler-alpha=<NO|YES|YES:ESTIMATE|YES:number>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  where numbers can have the range of 0.01 to 1.0, (NO=1.0=default=exp distrib)");
    switch(options->tri_mlalpha)
      {
      case FIXED:	
	print_parm_mutable(&bufsize, buffer, allocbufsize, "mittag-leffler-alpha=YES:%.2f", options->mlalpha);
	break;
      case NO:	
	print_parm_mutable(&bufsize, buffer, allocbufsize, "mittag-leffler-alpha=NO");
	break;
      case ESTIMATE:	
	print_parm_mutable(&bufsize, buffer, allocbufsize, "mittag-leffler-alpha=YES:ESTIMATE");
	break;
      }
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
    print_parm_br(&bufsize, buffer, allocbufsize);
#endif
    //
print_parm_comment(&bufsize, buffer, allocbufsize, "Report M (=migration rate/mutation rate) instead of 4Nm or 2 Nm or Nm");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax use-M=<NO | YES> Default is YES, the name 4Nm is ambiguous");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     for non-diploid data");
    print_parm_mutable(&bufsize, buffer, allocbufsize, "use-M=%s", options->usem ? "YES" : "NO");
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
    print_parm_br(&bufsize, buffer, allocbufsize);
        
    print_parm_comment(&bufsize, buffer, allocbufsize, "Print tree into treefile");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax print-tree=< NONE | <ALL | BEST | LASTCHAIN:Increment>:treefile >");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        NONE no tree printed [Default, and only choice using parallel");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        ALL  print all visited genealogies [careful this will be huge]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        BEST print only the best tree visited");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        LASTCHAIN print all trees in last chain");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        with increment INCREMENT");
    switch (options->treeprint)
      {
      case myNONE:
        print_parm(&bufsize, buffer, allocbufsize, "print-tree=NONE");
        break;
      case ALL:
        print_parm_mutable(&bufsize, buffer, allocbufsize, "print-tree=ALL:%s",options->treefilename);
        break;
      case BEST:
        print_parm_mutable(&bufsize, buffer, allocbufsize, "print-tree=BEST:%s",options->treefilename);
        break;
      case LASTCHAIN:
        print_parm_mutable(&bufsize, buffer, allocbufsize, "print-tree=LASTCHAIN:%li:%s", options->treeinc, options->treefilename);
        break;
      default:
        print_parm(&bufsize, buffer, allocbufsize, "print-tree=NONE");
        break;
      }
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
    print_parm_br(&bufsize, buffer, allocbufsize);
    
    
    print_parm_comment(&bufsize, buffer, allocbufsize, "Print a histogram of the time of migration events for each M(i->j)");
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax  mig-histogram=<NO | <ALL | MIGRATIONEVENTSONLY>:binsize:mighistfile >");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        NO            do not record any events");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        ALL           record migration and coalescence event");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        MIGRATIONEVENTSONLY record only migration events");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        binsize has to be in mutation units, with an average Theta=0.01 try 0.001");
    print_parm_comment(&bufsize, buffer, allocbufsize, "Print a histogram of the parameters through time (skyline plot)");
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax  skyline=<NO | YES | PARAM:n >:binsize:skylinefile >");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        NO            do not calculate parameter estimates through time");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        YES           calculate parameters through time");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        binsize has to be in mutation units, with an average Theta=0.01 try 0.001");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        If the interval is too fine the output will be very noisy");
    
    if(options->mighist)
      {
	if(options->mighist_all)
	  print_parm_mutable(&bufsize, buffer, allocbufsize, "mig-histogram=ALL:%f:%s", (double) options->eventbinsize,options->mighistfilename);
	else
	  print_parm_mutable(&bufsize, buffer, allocbufsize, "mig-histogram=YES:%f:%s", (double) options->eventbinsize, options->mighistfilename);
	if(options->skyline)
	  {
	    if(options->skyline_param)
	      print_parm_mutable(&bufsize, buffer, allocbufsize, "skyline=PARAM:%li:%f:%s  #needs mig-histogram", options->timeelements, (double) options->eventbinsize, options->skylinefilename);
	    else
	      print_parm_mutable(&bufsize, buffer, allocbufsize, "skyline=YES:%f:%s  #needs mig-histogram", (double) options->eventbinsize, options->skylinefilename);
	  }
      }
    else
      {
        print_parm(&bufsize, buffer, allocbufsize, "mig-histogram=NO");
        print_parm(&bufsize, buffer, allocbufsize, "skyline=NO #needs mig-histogram=ALL:...");
      }
    print_parm_br(&bufsize, buffer, allocbufsize);


    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_title(&bufsize, buffer, allocbufsize, "Parameter start settings");
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax: theta=<RANDOM | OWN:<{value} | {value1, value2, ...., valuen} | PRIOR:percentage >");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     migration==<RANDOM | OWN:<{value} | {value1, value2, ...., valuen} | PRIOR:percentage >");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       OWN     starting values are supplied by user");          
    print_parm_comment(&bufsize, buffer, allocbufsize, "          {value}   if only one value is supplied then all population");     
    print_parm_comment(&bufsize, buffer, allocbufsize, "                    have the same starting value");   
    print_parm_comment(&bufsize, buffer, allocbufsize, "          {value1, value2, ..., valuen} each population has its");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                    own starting value, if the number of values is");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                    insuffient, then the last value is the template");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                    for the remaining populations");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       RANDOM  starting parameter is drawn randomely from prior distribution");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       PRIOR   starting parameter is the value at the CDF(percentage) of the prior distribution");
    print_parm_comment(&bufsize, buffer, allocbufsize, "          percentage is a value between 0 and 100");
    //theta
    switch(options->startguess[0][0])
      {
      case PRIOR:
	print_parm_mutable(&bufsize, buffer, allocbufsize, "theta=PRIOR:%li",options->startguess[0][1]);
	break;
      case RANDOMPRIOR:
	print_parm(&bufsize, buffer, allocbufsize,"theta=RANDOM");
	break;
      default:
	print_parm_theta(&bufsize, buffer, allocbufsize, options);
	break;
      }
    // migration
    switch(options->startguess[1][0])
      {
      case PRIOR:
	print_parm_mutable(&bufsize, buffer, allocbufsize, "migration=PRIOR:%li",options->startguess[1][1]);
	break;
      case RANDOMPRIOR:
	print_parm(&bufsize, buffer, allocbufsize,"migration=RANDOM");
	break;
      default:
	print_parm_m(&bufsize, buffer, allocbufsize, options);
	break;
      }
    // RATE
    switch(options->startguess[2][0])
      {
      case PRIOR:
	print_parm_mutable(&bufsize, buffer, allocbufsize, "rate=PRIOR:%li",options->startguess[2][1]);
	break;
      case RANDOMPRIOR:
	print_parm(&bufsize, buffer, allocbufsize,"rate=RANDOM");
	break;
      default:
	print_parm_mutable(&bufsize, buffer, allocbufsize, "rate=PRIOR:50");// insert arbitrary menu value here TODO
	break;
      }
    // split time
    switch(options->startguess[3][0])
      {
      case PRIOR:
	print_parm_mutable(&bufsize, buffer, allocbufsize, "split=PRIOR:%li",options->startguess[3][1]);
	break;
      case RANDOMPRIOR:
	print_parm(&bufsize, buffer, allocbufsize,"split=RANDOM");
	break;
      default:
	print_parm_split(&bufsize, buffer, allocbufsize, options);
	break;
      }
    // split std
    switch(options->startguess[4][0])
      {
      case PRIOR:
	print_parm_mutable(&bufsize, buffer, allocbufsize, "splitstd=PRIOR:%li",options->startguess[4][1]);
	break;
      case RANDOMPRIOR:
	print_parm(&bufsize, buffer, allocbufsize,"splitstd=RANDOM");
	break;
      default:
	print_parm_splitstd(&bufsize, buffer, allocbufsize, options);
	break;
      }
    print_parm_br(&bufsize, buffer, allocbufsize);

    print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "Mutation rate modifiers");
    print_parm_br(&bufsize, buffer, allocbufsize);
    //    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax: mutation=<NOGAMMA | CONSTANT | ESTIMATE | GAMMA:alpha | OWN:loci: rate1 rate2 ... rate_loci>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax: mutation=<CONSTANT | ESTIMATE | OWN:loci: rate1 rate2 ... rate_loci>");
    //   print_parm_comment(&bufsize, buffer, allocbufsize, "     NOGAMMA      all loci have same mutation rate");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     CONSTANT     all loci have same mutation rate");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     ESTIMATE     BAYESIAN estimate: mutation rate is drawn from prior");
    //    print_parm_comment(&bufsize, buffer, allocbufsize, "     GAMMA:alpha  ML estimate: mutation rate has Gamma distribution with alpha");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     OWN          mutation rate is different for every locus, but fixed");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        :loci: rate1, ...     number of loci, rate of locus 1, locus 2 etc.");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     DATA         mutation rate modifier is deducted from loci in the data");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                  using Watterson's Theta and then scaling all rates Theta_locus/mean(Theta_loci");

    if (options->gamma)
      {
        print_parm_mutable(&bufsize, buffer, allocbufsize, "mutation=%s:%f", "GAMMA", options->alphavalue);
      }
    else
      {
        if (options->murates)
	  {
	    if(options->murates_fromdata)
	      {
		sprintf (fp, "mutation=DATA\n");
		add_to_buffer(fp,&bufsize,buffer, allocbufsize);	      
	      }
	    else
	      {
		sprintf (fp, "mutation=OWN:%li: ", options->muloci);
		add_to_buffer(fp,&bufsize,buffer, allocbufsize);
		for (i = 0; i < options->muloci; i++)
		  {
		    sprintf (fp, "%f ", options->mu_rates[i]);
		    add_to_buffer(fp,&bufsize,buffer, allocbufsize);
		  }
		sprintf (fp, " \n");
		add_to_buffer(fp,&bufsize,buffer, allocbufsize);
	      }
	  }
        else
	  {
	    if(options->bayesmurates)
	      {
		print_parm(&bufsize, buffer, allocbufsize, "mutation=ESTIMATE");
	      }
	    else
	      {
		print_parm(&bufsize, buffer, allocbufsize, "mutation=CONSTANT");
	      }
	  }
      }
    print_parm_br(&bufsize, buffer, allocbufsize);
    /*
    print_parm_comment(&bufsize, buffer, allocbufsize, "Smoothing of the posterior distribution");
    print_parm_comment(&bufsize, buffer, allocbufsize, "Syntax: smoothing={<BIN|KERNEL<:el>>,<NONE|AVERAGE|SAVITZKY<:el>>}");
    print_parm_comment(&bufsize, buffer, allocbufsize, "Default is smoothing={KERNEL:20,SAVITZKY:20}");
    print_parm_comment(&bufsize, buffer, allocbufsize, "el is the windowsize=2*el+1");
    print_parm_smoothing(&bufsize, buffer, allocbufsize, options);
    */
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "Treatment of invariant sequence loci");
    print_parm_comment(&bufsize, buffer, allocbufsize, "Syntax: analyze-loci=<A | F | V>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        A = analyze all loci (Default!)");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        F = analyze all variable loci and ONE invariant and extrapolate");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        V = analyze only variable loci");
    print_parm_mutable(&bufsize, buffer, allocbufsize,  "analyze-loci=%c",( options->onlyvariable ? 'V' :(options->has_variableandone ? 'F' : 'A')));
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "Distribution of divergence paramemter");
    //print_parm_comment(&bufsize, buffer, allocbufsize, "Syntax: divergence-distrib=<E | W | N | S>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "Syntax: divergence-distrib=<E | N>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        E = Exponential = use exponential distribution (Default!)");
    //    print_parm_comment(&bufsize, buffer, allocbufsize, "        W = Weibull = use Weibull distribution");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        N = Normal = use Normal distribution");
    //print_parm_comment(&bufsize, buffer, allocbufsize, "        S = Normal Shortcut = use Normal distribution using a shortcut");
    // print_parm_mutable(&bufsize, buffer, allocbufsize,  "divergence-distrib=%c",( options->species_model_dist == 2 ? 'E' :(options->species_model_dist==0 ? 'W' : (options->species_model_dist == NORMAL_DIST ? 'N' : 'S'))));
    print_parm_mutable(&bufsize, buffer, allocbufsize,  "divergence-distrib=%c",( options->species_model_dist == 2 ? 'E' :(options->species_model_dist==0 ? 'N' : (options->species_model_dist == NORMAL_DIST ? 'N' : 'N'))));
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "Custom migration model");
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax: custom-migration={ab..bbab..ba ... a}");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       the {} is a square matrix with values for the population sizes");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       on the diagonal and migration rates or divergences off-diagonal");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       the values _a_ for the diagonal can be any of these:");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       c constant, define the values in the theta option");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       * [or x] free to vary, the default is * for every parameter");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       m mean of theta, this can be a subgroup");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         for example: theta 1-3 are averaged, thetas 4,5 are estimated");
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "       the values _b_ for the migration/divergence rates can be any of these:");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       c       constant, the value needs to be defined in the migration option");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       *       migration rate free to vary, this is the default");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       d       row population is offspring of column population");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       D       row population is offspring of column population with migration");
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "       m       mean of M_ij, this can be a subgroup of migration rates");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               for example the M_1-3i are averaged and M_4,5i are estimated");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       M       means of 4Nm (diploid), 2Nm (haploid), Nm (mtDNA, Y-chromosome)");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       s       symmetric migration rates M");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       S       symmetric migrants 4Nm");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       a,b,c,d,...  any parameter group label for migration and theta");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               this allows using averages of multiple groups [se example]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       ");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       an example for 5 populations could look like this:");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       custom-migration={*s00s"); 
    print_parm_comment(&bufsize, buffer, allocbufsize, "                         s*s00");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                         0s*s0");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                         00s*s");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                         s00s*}");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       this describes a circular stepping stone model with 5 symmetric rates");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        and independent sizes, a stepping stone with 3 parameters would");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       look like this custom-migration={ma00b bma00 0bma0 00bma a00bm}");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       showing two different migration rates: clockwise and counter-clockwise");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       population splitting example with 3 populations:");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       custom-migration={*D0 **d 00*}"); 
    print_parm_comment(&bufsize, buffer, allocbufsize, "                    *-----1     1 receives migrants from 2");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                    |MMMMM        and split from 2 ");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                 *--------2     2 receives migrants from 1");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                 |                and split from 3");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        ---------*--------3     3 is the persisting ancestor");
    print_parm_mutable(&bufsize, buffer, allocbufsize, "custom-migration={%s}", options->custm);
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_growpops(&bufsize, buffer, allocbufsize, options, data);
    print_parm_comment(&bufsize, buffer, allocbufsize, "Influence of geography on migration rate");
    print_parm_comment(&bufsize, buffer, allocbufsize, "a distance matrix between populations changes the migration rate matrix so that");
    print_parm_comment(&bufsize, buffer, allocbufsize, "(genetic?) migration rates =  inferred migration rate / distance ~ a dispersion coefficient");
    print_parm_comment(&bufsize, buffer, allocbufsize, "the geofile contains a number of populations, names for populations (10 characters), they");
    print_parm_comment(&bufsize, buffer, allocbufsize, "need to be in order of the dataset. And the distances between the populations, they do not");
    print_parm_comment(&bufsize, buffer, allocbufsize, "need to be symmetric; useful distances are relative to 1.0");
    print_parm_comment(&bufsize, buffer, allocbufsize, "for example using km with large distances may lead to very large");
    print_parm_comment(&bufsize, buffer, allocbufsize, "values since they are per distance unit.");
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax: geo:<NO | YES:filename>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "            NO       distances among populations are considered to be 1 [all equal]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "            YES      distances are read from a file");

    if (options->geo)
      {
        print_parm_mutable(&bufsize, buffer, allocbufsize, "geo=YES:%s", options->geofilename);
      }
    else
      {
        print_parm(&bufsize, buffer, allocbufsize, "geo=NO");
      }
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_br(&bufsize, buffer, allocbufsize);
    // SEARCH STRATEGIES
    print_parm_title(&bufsize, buffer, allocbufsize, "Search strategies");
    print_parm_br(&bufsize, buffer, allocbufsize);
    
    print_parm_comment(&bufsize, buffer, allocbufsize, "Bayesian MCMC Strategy method");
    print_parm_comment(&bufsize, buffer, allocbufsize, "      updatefreq=VALUE VALUE VALUE VALUE VALUE VALUE");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        VALUE is a ratio between 0 and 1");
    print_parm_comment(&bufsize, buffer, allocbufsize, "              ratio of how many times the genealogy is updated compared to the parameters");
    print_parm_comment(&bufsize, buffer, allocbufsize, "              If the value is 0.4 in a 2-population scenario and with 1000000 steps");
    print_parm_comment(&bufsize, buffer, allocbufsize, "              The tree will be evaluated 400000 times, Theta_1, Theta_2, M_21, and M_12");
    print_parm_comment(&bufsize, buffer, allocbufsize, "              will be each evaluated 125000 times. The second value is the ratio for parameter");
    print_parm_comment(&bufsize, buffer, allocbufsize, "              updates, and the third value is the frequency of hapltype updates.");
    print_parm_comment(&bufsize, buffer, allocbufsize, "              the fourth values is for assignment of individual updates");
    print_parm_comment(&bufsize, buffer, allocbufsize, "              The values do not need add up to 1.0 but will be recalculated to do so");
    print_parm_comment(&bufsize, buffer, allocbufsize, "              For example: 1.0 2.0 0.1 0.0 results in 0.32 treeupdates 0.65 parameter ");
    print_parm_comment(&bufsize, buffer, allocbufsize, "              updates and 0.03 haplotype updates");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       bayes-posteriorbins=VALUE VALUE");
    print_parm_comment(&bufsize, buffer, allocbufsize, "           VALUE      is the number of bins in the posterior distribution histogram for Theta or M");
#ifdef PRETTY
    print_parm_comment(&bufsize, buffer, allocbufsize, "       bayes-posteriormaxtype=< ALL | P99 | MAXP99 | TOTAL >");
    print_parm_comment(&bufsize, buffer, allocbufsize, "           ALL        plots the WHOLE prior-parameter range");
    print_parm_comment(&bufsize, buffer, allocbufsize, "           P99        plots from the minimum prior range value to");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                      the 99% percentile value of EACH parameter");
    print_parm_comment(&bufsize, buffer, allocbufsize, "           MAXP99     sets all axes from minimum to the maximal");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                      99% percentile value of ALL parameter");
    print_parm_comment(&bufsize, buffer, allocbufsize, "           TOTAL      plots from the minimum prior range value to");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                      the 100% percentile value of EACH parameter");
#endif
    print_parm_comment(&bufsize, buffer, allocbufsize, "       bayes-file=<YES:FILENAME|NO>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "           FILENAME is the name of the file that will contain");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                   the results for the posterior distribution");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       bayes-allfile=<<YES|TEMP>:INTERVAL:FILENAME|NO>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "           FILENAME is the name of the file that will contain");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                   all parameters of the posterior distribution [HUGE]");
    //OLD    print_parm_comment(&bufsize, buffer, allocbufsize, "       bayes-allfileinterval=INTERVAL");
    print_parm_comment(&bufsize, buffer, allocbufsize, "           INTERVAL is the interval at which all parameters are written to file\n");
    print_parm_comment(&bufsize, buffer, allocbufsize, "    PROPOSAL:");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       bayes-proposals= THETA < SLICE | METROPOLIS >");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       bayes-proposals= MIG < SLICE | METROPOLIS >");
    //    print_parm_comment(&bufsize, buffer, allocbufsize, "       bayes-proposal= RATE < SLICE | METROPOLIS >");
    print_parm_comment(&bufsize, buffer, allocbufsize, "              SLICE uses the slice sampler to propose new parameter values");
    print_parm_comment(&bufsize, buffer, allocbufsize, "              METROPOLIS uses the Metropolis-Hastings sampler");
    print_parm_comment(&bufsize, buffer, allocbufsize, "              (this is done for each parameter group: THETA or MIGration)");
    print_parm_comment(&bufsize, buffer, allocbufsize, "    PRIORS:");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       Priors can be set for each parameter");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       There are several ways to set: a. old format, b. new format for all, c. individual");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       bayes-priors= FORCE <PRIORdistribution priorvalues> # a. all in FORCE");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       bayes-priors= FORCE * * <PRIORdistribution priorvalues> # b. all in FORCE");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       bayes-priors= FORCE from to <PRIORdistribution priorvalues> # c. individual,");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                                   from and to are population numbers");    
    print_parm_comment(&bufsize, buffer, allocbufsize, "       FORCE is one of THETA, MIG, RATE, SPLIT, or SPLITSTD");
    print_parm_comment(&bufsize, buffer, allocbufsize, "              THETA is used for population size parameters");    
    print_parm_comment(&bufsize, buffer, allocbufsize, "              MIG is used for migration rate parameters");
    print_parm_comment(&bufsize, buffer, allocbufsize, "              RATE is used for evolutionary rate differences (use only with date samples)");    
    print_parm_comment(&bufsize, buffer, allocbufsize, "              SPLIT is used for mean of the normal distributed population divergence");    
    print_parm_comment(&bufsize, buffer, allocbufsize, "              SPLITSTD is used for the standard deviation of the population divergence");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       PRIORdistribution is one of UNIFORMPRIOR, EXPPRIOR, WEXPPRIOR, GAMMAPRIOR, NORMALPRIOR");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               unipriorvalues: min max delta");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               exppriorvalues: min mean max");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               wexppriorvalues: min mean max delta");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               gammapriorvalues: min mean max alpha");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               normalpriorvalues: min mean max std");
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "Search OPTIONS");    
    print_parm_comment(&bufsize, buffer, allocbufsize, "       long-inc=VALUE      VALUE is the number of updates that are not recorded");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       long-sample=VALUE   VALUE is the number of sampled updates");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       burn-in=VALUE       VALUE is the number of updates to discard at the beginning");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       auto-tune=<NO | YES:VALUE>  VALUE the the target acceptance ratio");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                            if value is missing, it is set to 0.44");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       assign=<YES:<FREQ|UNIFORM<:{pop1,pop2,...}>> | NO>    YES will assign individuals to populations");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                            with '?' as first character");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                            [FREQ uses past assignment during the MCMC");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                            UNFORM picks a random population");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                            optional list if populations to assign to");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                            can be supplied using :{pop1,pop2,...}");
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "updatefreq= tree parameter haplotype timeparam assignment seqerror, ml-alpha");
    print_parm_mutable(&bufsize, buffer, allocbufsize, "updatefreq=%f %f %f %f %f %f %f ",
		       options->tree_updatefreq,
		       options->parameter_updatefreq,
		       options->haplotype_updatefreq,
		       options->timeparam_updatefreq,
		       options->unassigned_updatefreq,
		       options->seqerror_updatefreq,
		       options->mlalpha_updatefreq);
    long count=0;
    long cc;
    //prior_fmt * p = options->bayes_priors;
    //while (p!= NULL)
    for(cc=0;cc<5;cc++)
      {
	count += sprintf(mytext+count,"%li ",options->bayes_posterior_bins[cc]);
	//p = p->next;
	//count += sprintf(mytext+count,"%li ",p[count]);
	//p = p->next;
      }
    //printf("%i> \n\n\n\n\n\n\n\n\n",myID);
    print_parm_mutable(&bufsize, buffer, allocbufsize, "bayes-posteriorbins=%s", mytext);
#ifdef DEBUG
    printf("@@@@@@@@@@@@@ bins in save_options_buffer: %s\n",mytext);
#endif
    print_parm_mutable(&bufsize, buffer, allocbufsize, "bayes-posteriormaxtype=%s", 
		       (options->bayespretty == PRETTY_P99 ? "P99" :
			(options->bayespretty == PRETTY_MAX ? "ALL" :
			 (options->bayespretty == PRETTY_P100 ? "TOTAL" : "TOTAL"))));
    
    if(options->has_bayesfile)
      print_parm_mutable(&bufsize, buffer, allocbufsize, "bayes-file=YES:%s", options->bayesfilename);
    else
      print_parm_mutable(&bufsize, buffer, allocbufsize, "bayes-file=NO");
    
    if(options->has_bayesmdimfile)
      {
	print_parm_mutable(&bufsize, buffer, allocbufsize, "bayes-allfile=%s:%li:%s",
			   (options->mdimdelete ? "TEMP" : "YES"),	    
			   options->bayesmdiminterval, options->bayesmdimfilename);
      }
    else
      print_parm_mutable(&bufsize, buffer, allocbufsize, "bayes-allfile=NO");
    if(options->allposteriors)
      print_parm_mutable(&bufsize, buffer, allocbufsize, "bayes-all-posteriors=YES /*this can results in very many (param*loci) plots */");
    else
      print_parm_mutable(&bufsize, buffer, allocbufsize, "bayes-all-posteriors=NO");
    print_parm_proposal(&bufsize, buffer, allocbufsize, options);
    //
    print_parm_prior(&bufsize, buffer, allocbufsize, options);  //<-- print priors into buffer
    print_parm_br(&bufsize, buffer, allocbufsize);
    //
    print_parm_hyperprior(&bufsize,buffer,allocbufsize,options);
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_mutable(&bufsize, buffer, allocbufsize, "long-chains=%li", options->lchains);
    print_parm_mutable(&bufsize, buffer, allocbufsize, "long-inc=%li", options->lincrement);
    print_parm_mutable(&bufsize, buffer, allocbufsize, "long-sample=%li", options->lsteps);
    print_parm_mutable(&bufsize, buffer, allocbufsize, "burn-in=%li %c", options->burn_in, options->burnin_autostop);
    if (options->has_autotune)
      print_parm_mutable(&bufsize, buffer, allocbufsize, "auto-tune=YES:%f", options->autotune);
    else
      print_parm_mutable(&bufsize, buffer, allocbufsize, "auto-tune=NO");
    if (options->has_unassigned)
      {
	char assignopt[LINESIZE];
	char *asptr = assignopt;
	long mysize=0;
	if (options->has_unassignedfreq)
	  mysize += sprintf(asptr,"YES:FREQ");
	else
	  mysize += sprintf(asptr,"YES:UNIFORM");

	if (options->has_unassignedpoplist)
	  {
	    size_t ai;
	    mysize += sprintf(asptr+mysize,":{");
	    for (ai=0; ai < options->unassignedpoplistnum;ai++)
	      mysize += sprintf(asptr+mysize,"%li,",options->unassignedpoplist[ai]);
	    mysize = sprintf(asptr+mysize-1,"}");
	  }
	print_parm_mutable(&bufsize, buffer, allocbufsize, "assign=%s\n",assignopt);
      }
    else
      print_parm_mutable(&bufsize, buffer, allocbufsize, "assign=NO");

    print_parm_br(&bufsize, buffer, allocbufsize);    
    print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "Schemes to improve MCMC searching and/or thermodynamic integration");
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "Heating schemes {MCMCMC = MC cubed}");
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax: heating=< NO | <YES | ADAPTIVE>:SKIP:TEMPERATURES");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       NO    No heating");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       YES   heating using TEMPERATURES");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       ADAPTIVE adaptive heating using start TEMPERATURES [fails sometimes!!]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       SKIP skip that many comparisons, this lengthens the run by SKIP");
    print_parm_comment(&bufsize, buffer, allocbufsize, "           TEMPERATURES    { 1.0, temp1, temp2, temp3 .. tempn}");
    print_parm_comment(&bufsize, buffer, allocbufsize, "    Example: heating=YES:1:{1.0, 1.2, 3.0,1000000.0}");
    print_parm_comment(&bufsize, buffer, allocbufsize, "Heating:  swapping chains");
    print_parm_comment(&bufsize, buffer, allocbufsize, "    Syntax: heated-swap=< YES | NO >");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        YES  swapping of chains enabled [DEFAULT]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        NO   swapping of chains disabled");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     Example: heated-swap=YES");
    print_parm_heating(&bufsize, buffer, allocbufsize, options);
    
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "Lengthening chain schemes");
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax: moving-steps=< NO | YES:VALUE>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "      VALUE   frequency is between 0..1");

    if (options->movingsteps)
      {
        print_parm_mutable(&bufsize, buffer, allocbufsize, "moving-steps=YES:%f", options->acceptfreq);
      }
    else
      {
        print_parm(&bufsize, buffer, allocbufsize, "moving-steps=NO");
      }
    print_parm_br(&bufsize, buffer, allocbufsize);
    //print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax: long-chain-epsilon=VALUE");
    //print_parm_comment(&bufsize, buffer, allocbufsize, "      VALUE    is between 0..INFINITY");
    //print_parm_comment(&bufsize, buffer, allocbufsize, "               the VALUE is the likelihood ratio between the old and thew chain");
    //print_parm_comment(&bufsize, buffer, allocbufsize, "               the VALUE depends on the number of parameters: with 1 values of 0.5 are great");
    //print_parm_comment(&bufsize, buffer, allocbufsize, "               but with many parameters values and bad data >20 is more reasonable");
    //if (options->lcepsilon < LONGCHAINEPSILON)
    //   print_parm_mutable(&bufsize, buffer, allocbufsize, "long-chain-epsilon=%f", options->lcepsilon);
    //else
    //    print_parm_mutable(&bufsize, buffer, allocbufsize, "long-chain-epsilon=INFINITY");

    
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Convergence statistic [Gelman and Rubin]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax: gelman-convergence=< YES:Pairs|Summary | NO >");    
    print_parm_comment(&bufsize, buffer, allocbufsize, "      NO      do not use Gelman's convergence criterium");
    print_parm_comment(&bufsize, buffer, allocbufsize, "      YES     use Gelman's convergence criteria between chain i, and i-1");
    print_parm_comment(&bufsize, buffer, allocbufsize, "              PAIRS reports all replicate pairs");
    print_parm_comment(&bufsize, buffer, allocbufsize, "              SUM   reports only mean and maxima");
    print_parm_mutable(&bufsize, buffer, allocbufsize, "gelman-convergence=%s", options->gelman ? (options->gelmanpairs ? "Yes:Pairs" : "Yes:Sum" ) : "No");
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "   REPLICATON:");
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax: replicate=< NO | YES:<VALUE> >");
    print_parm_comment(&bufsize, buffer, allocbufsize, "      NO     no replication of run");
    print_parm_comment(&bufsize, buffer, allocbufsize, "      YES    replicate run");
    print_parm_comment(&bufsize, buffer, allocbufsize, "          VALUE     number between 2 and many, complete replicates");
    if(!options->replicate)
      {
        print_parm(&bufsize, buffer, allocbufsize, "replicate=NO");
      }
    else
      {
        if(options->replicatenum!=0)
	  print_parm_mutable(&bufsize, buffer, allocbufsize, "replicate=YES:%li", options->replicatenum);
      }
    print_parm_br(&bufsize, buffer, allocbufsize);
    
    //    print_parm_mutable(&bufsize, buffer, allocbufsize, "resistance=%f", options->minmigsumstat);
    //    print_parm_br(&bufsize, buffer, allocbufsize);
    
    print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm(&bufsize, buffer, allocbufsize, "end");
    return bufsize;
}


/*private functions============================================= */
///
/// fill the prior information of bayes theta, migrations, and rates 
/// this fills the first theta prior, migration prior, rate prior.
void set_bayes_options(char *value, option_fmt *options)
{
  //prior_fmt *last;
  float mini;
  float maxi;
  float meani;
  float alpha;
  float std;
  float delta;
  int ptype;
  //long key = 0;
  long from = -1;
  long to = -1;
  char sfrom[LINESIZE], sto[LINESIZE];
  char paramtype[LINESIZE], priortype[LINESIZE];
  prior_fmt *prior=NULL;
  char ptypename[PRIOR_SIZE][20];
  strcpy(ptypename[THETAPRIOR],"THETA");
  strcpy(ptypename[MIGPRIOR],"MIG");
  strcpy(ptypename[SPECIESTIMEPRIOR],"SPLIT");
  strcpy(ptypename[SPECIESSTDPRIOR],"SPLITSTD");
  strcpy(ptypename[RATEPRIOR],"RATE");
  strcpy(ptypename[GROWTHPRIOR],"GROWTH");
  
  sscanf(value,"%s%s", paramtype, priortype);
  char *valueptr = strstr(value,priortype);
  if (strchr("-*1234567890", uppercase(priortype[0])))
    {
      sscanf(value,"%s%s%s", paramtype, sfrom, priortype);
      valueptr = strstr(value,priortype);
      if (strchr("-*1234567890", uppercase(priortype[0])))
	{
	  sscanf(value,"%s%s%s%s", paramtype, sfrom, sto, priortype);
	  valueptr = strstr(value,priortype);
	  if (strchr("TMSRG", uppercase(paramtype[0])))
	    {
	      if (sfrom[0]=='*' || sfrom[0]=='-')
		{
		  from = -1;
		}
	      else
		{
		  from = atol(sfrom) - 1;
		}
	      if (sto[0]=='*' || sto[0] == '-')
		{
		  to = -1;
		}
	      else
		{
		  to = atol(sto) - 1;
		}
	    }
	  else
	    {
	      warning("Problem with prior definition: %s\n",value);
	    }
	}
      else
	{
	  if (strchr("TR", uppercase(paramtype[0])))
	    {
	      if (sfrom[0]=='*')
		{
		  from = -1;
		}
	      else
		{
		  from = atol(sfrom);
		}
	      to = from;
	    }
	  else
	    {
	      warning("Problem with prior definition: %s\n",value);
	    }
	}
    }
  else
    {
      from = -1;
      to   = -1;
    }
  switch(uppercase(paramtype[0]))
    {
    case 'T'/* THETA*/: 
      ptype = THETAPRIOR;
      
      break;
    case 'M'/* MIG  */: 
      ptype = MIGPRIOR; 
      break;
    case 'S'/* SPLIT or SPECIES  */:
      if (strstr(paramtype,"STD"))
	{
	  ptype = SPECIESSTDPRIOR; 
	}
      else
	{
	  ptype = SPECIESTIMEPRIOR;
	} 
      break;
    case 'R'/* RATE */: 
      ptype = RATEPRIOR; 
      break;
    case 'G'/* THETA*/: 
      ptype = GROWTHPRIOR;      
      break;
    default:
      return;
    }
  if (options->bayes_priors == NULL)
    {
      options->bayes_priors = (prior_fmt *) mycalloc(1,sizeof(prior_fmt));
      options->bayes_priors_num = 1;
    }
  else
    {
      options->bayes_priors_num += 1;
      options->bayes_priors = (prior_fmt *) myrealloc(options->bayes_priors,
				      (unsigned long) options->bayes_priors_num*sizeof(prior_fmt));
    }
  prior = &(options->bayes_priors[ options->bayes_priors_num-1]);
  prior->from = from;
  prior->to = to;
  prior->type = ptype;
  strcpy(prior->ptypename, ptypename[ptype]);  
  prior->number = options->bayes_priors_num - 1;
  switch(uppercase(priortype[0]))
    {
    case 'Z'/*slice sampler with uniform prior*/:
      sscanf(valueptr,"%s%f%f", priortype, &mini, &maxi);
      prior->min = (MYREAL) mini;
      prior->max = (MYREAL) maxi;
      prior->kind = SLICE;
      options->slice_sampling[ptype] = TRUE;
      break;
    case 'M'/*multprior   */:
      sscanf(valueptr,"%s%f%f%f", priortype, &mini, &maxi, &delta); 
      prior->min = (MYREAL) mini;
      prior->max = (MYREAL) maxi;
      prior->delta = (MYREAL) delta;
      prior->kind = MULTPRIOR;
      break;
      case 'E'/*expprior    */:   
	sscanf(valueptr,"%s%f%f%f", priortype, &mini,&meani, &maxi);
	prior->min = (MYREAL) mini;
	prior->mean = (MYREAL) meani;
	prior->max = (MYREAL) maxi;
	prior->delta = (MYREAL) (prior->max + prior->min) / 10.;
	prior->kind = EXPPRIOR;
	break;
      case 'W'/*wexpprior   */:   
	sscanf(valueptr,"%s%f%f%f%f", priortype, &mini,&meani, &maxi, &delta);
	prior->min = (MYREAL) mini;
	prior->mean = (MYREAL) meani;
	prior->max = (MYREAL) maxi;
	prior->delta = (MYREAL) delta;
	prior->kind = WEXPPRIOR;
	break;
      case 'G'/*gammaprior  */:   
	sscanf(valueptr,"%s%f%f%f%f", priortype, &mini,&meani, &maxi, &alpha);
	prior->min = (MYREAL) mini;
	prior->mean = (MYREAL) meani;
	prior->max = (MYREAL) maxi;
	prior->alpha = (MYREAL) alpha;
	prior->kind = GAMMAPRIOR;
	prior->delta = (prior->max + prior->min) / 10.;
	break;
      case 'N'/*normalprior  */:
	//warning("not supported yet!");
	sscanf(valueptr,"%s%f%f%f%f", priortype, &mini,&meani, &maxi, &std);
	prior->min = (MYREAL) mini;
	prior->mean = (MYREAL) meani;
	prior->max = (MYREAL) maxi;
	prior->alpha = (MYREAL) std;
	prior->kind = NORMALPRIOR;
	prior->delta = (fabs(prior->max) + fabs(prior->min)) / 10.;
	break;
      case 'U'/*uniformprior*/:   
	sscanf(valueptr,"%s%f%f%f", priortype, &mini, &maxi, &delta);
	prior->min = (MYREAL) mini;
	prior->max = (MYREAL) maxi;
	prior->mean = (prior->max + prior->min) / 2.;
	prior->delta = (MYREAL) delta;
	prior->kind = UNIFORMPRIOR;
	break;
    }
}

///
/// checks all elements in the options file (parmfile) that involve boolean yes no type options
/// typical format is option=<NO | YES<: filename or parameters etc>
long boolcheck (char ch)
{
    char c = uppercase (ch);
    if ((c == 'F') || (c == 'N'))
        return 0;
    else if ((c == 'T') || (c == 'Y'))
        return 1;
    else
        return -1;
}    /* boolcheck */

boolean
booleancheck (option_fmt * options, char *var, char *value)
{
    long i, check;
    char *booltokens[NUMBOOL] = BOOLTOKENS;
    char *tmp;
    //long ltemp;
    char *extension;
    char val;
    if (value != NULL)
      val = value[0];
    else
      val = ' ';
    check = boolcheck (val);
    if (check == -1)
      {
        return FALSE;
      }
    i = 0;
    while (i < NUMBOOL && strcmp (var, booltokens[i]))
        i++;
#ifdef DEBUG
    printf("%i> booleancheck: %s\n",myID, var);
#endif
    boolean has_mlalpha = FALSE;
    switch ((short) i)
    {
    case 0:   /*menu = <yes | no> */
        options->menu = (boolean) (check);
        break;
    case 1:   /*checkpointing =<yes | no> */
        options->checkpointing = (boolean) (check);
        break;
    case 2:   /*print-data = <yes | no> */
        options->printdata = (boolean) (check);
        break;
    case 3:   /* mixplot=<yes:mixfilename | no> */
      options->mixplot = set_filename(value, "YES", &options->mixfilename);
        break;
    case 4:   /* moving-steps = <yes | no> */
        options->movingsteps = (boolean) (check);
        if (options->movingsteps)
        {
            strtok (value, ":");
            tmp = strtok (NULL, " ,\n");
            if (tmp != NULL)
                options->acceptfreq = atof ((char *) tmp);
            else
                options->acceptfreq = 0.1;
        }
        break;
    case 5:   /* freqs-from-data =  <yes | no> */
        options->freqsfrom = (boolean) (check);
        if (options->freqsfrom)
	  {
	    strtok (value, "[");
            tmp = strtok (NULL, "]");
            if (tmp != NULL)
                options->totalsites = atol ((char *) tmp);          
	  }
	else
	  {
            strtok (value, ":");
            tmp = strtok (NULL, " ,");
            if (tmp != NULL)
                options->freqa = atof ((char *) tmp);
            tmp = strtok (NULL, " ,");
            if (tmp != NULL)
                options->freqc = atof ((char *) tmp);
            tmp = strtok (NULL, " ,");
            if (tmp != NULL)
                options->freqg = atof ((char *) tmp);
            tmp = strtok (NULL, "[,\n");
            if (tmp != NULL)
	      {
                options->freqt = atof ((char *) tmp);
		// if we encounter a bracket here then we asssume the user
		// gave us a number of total sites examined S, this will then
		// be used to calculate weight for the invariant sites per locus (averaged),
		// freq * (S-N)
		tmp = strtok (NULL, "\n");
		if (tmp != NULL)
		  {
		    options->totalsites = atol ((char *) tmp);
		  }
	      }
        }
        break;
    case 6:   /* useroldtree =  < NO | YES OBSOLETE use usetree*/
        options->usertree = set_filename(value, "YES", &options->utreefilename); //USERTREE
        break;
    case 7:
        {    /* autocorrelation=<YES:value | NO> */
            options->autocorr = (boolean) (check);
            if (options->autocorr)
            {
                strtok (value, ":");
                tmp = strtok (NULL, ";\n");
                if (tmp != NULL)
                    options->lambda = 1.0 / atof ((char *) tmp);
            }
            break;
        }
    case 8:   /* simulation =  <yes | no> */
        options->simulation = (boolean) check;
        break;
    case 9:   /* changed to bayes-all-posterior=YES|NO*/
        options->allposteriors = (boolean) check;
        break;
    case 10:   /* weights =  <yes | no> */
        options->weights = (boolean) check;
        set_filename(value, "YES", &options->weightfilename);
        break;
    case 11:   /* read-summary  <yes | no> */
        options->readsum = (boolean) check;
        options->datatype = 'g';
        break;
	//case 12:   /* write-summary =  <yes | no> */
        //options->writesum = (boolean) check;
        //set_filename(value, "YES", &options->sumfilename);
        //break;
    case 13: /* NODATA = <yes | no >  */
        options->prioralone = (boolean) check;
        break;
    case 14:   /* include-unknown=<yes | no> */
        options->include_unknown = (boolean) check;
        break;
        // old case 14(heating) moved to numbercheck
    case 15:   /* print-fst =  <yes | no> */
        options->printfst = (boolean) check;
        break;
    case 16:   /* distfile =  <yes | no> */
        warning("OBSOLETE option distance=YES|NO ===> use usertree=DISTANCE:distfile or usertree=NO\n");
        options->dist = (boolean) check;
        break;
    case 17:   /* geofile =  <yes:geofile | no> */
        options->geo = (boolean) check;
        set_filename(value, "YES", &options->geofilename);
        break;
    case 18:   /* gelman-convergence =  <yes:pairs|sum | no> */
        options->gelman = (boolean) check;
	options->gelmanpairs = FALSE;
        if (options->gelman)
        {
            get_next_word(&value, ":\n", &tmp);
            if (value != NULL)
            {
                if (uppercase (value[0]) == 'P')
		  {
		    options->gelmanpairs = TRUE;
		  }
	    }
	}
        break;
    case 19:   /* randomtree start =  <yes | no> */
        warning("OBSOLETE option randomtree=YES|NO ===> use usertree=RANDOM or usertree=NO\n");
        options->randomtree = (boolean) check;
        break;
    case 20:   /* fast-likelihood calculator =  <yes | no> */
        options->fastlike = (boolean) check;
        break;
    case 21:   /*pdf-terse=yes/no */
        options->tersepdf = (boolean) check;
        break;
    case 22:   //use-M=true/false
        options->usem = (boolean) check;
        //set_usem_related (options);
        break;
    case 23:   //alternative to above is use-4Nm=false/true
        options->usem = !(boolean) check;
        //set_usem_related (options);
        break;
    case 24: //bayes-update
        options->bayes_infer = TRUE;
	if(options->bayes_infer)
	  options->lchains = 1;
        break;
    case 25: // bayes-allfile: write all bayes parameter estimates to a file
        options->has_bayesmdimfile = (boolean) check;
	if(options->has_bayesmdimfile)
	  {
	    get_next_word(&value,":",&tmp);// extract the word YES, value is now interval:filename
	    options->mdimdelete=FALSE;
	    if (!strcmp(tmp,"TEMP"))
		       {
			 options->mdimdelete=TRUE;
		       }			 
	    get_next_word(&value,":",&tmp);// extract the interval, value is now filename
	    if(value==NULL)
	      {
		strncpy (options->bayesmdimfilename, tmp, 255);
	      }
	    else
	      {
		options->bayesmdiminterval = atol(tmp);
		strncpy (options->bayesmdimfilename, value, 255);
	      }
	    unpad(options->bayesmdimfilename," ");
	    extension = strrchr(options->bayesmdimfilename,'.');
#ifdef ZNZ
	    if(extension!=NULL && !strncmp(extension,".gz",3))
	      {
		options->use_compressed = 1;
	      }
	    else
	      {
		options->use_compressed = 0;
	      }  
#else
	    options->use_compressed = 0;
#endif
	  }
        break;
    case 26: // bayes-file: write marginal parameter posterior distribution
        options->has_bayesfile = (boolean) check;
        set_filename(value, "YES", &options->bayesfilename);
        break;
    /* divfile datefile */
    case 27:   /* divfile =  <yes:divfile | no> */
        options->div = (boolean) check;
        set_filename(value, "YES", &options->divfilename);
        break;
    case 28: /* tipdate-file: <yes:datefile | no >  */
        options->has_datefile = (boolean) check;
        set_filename(value, "YES", &options->datefilename);
        break;
    case 29: /* heated-swap: <yes | no >  */
      options->heatedswap_off = !((boolean) check);
      break;
    case 30: /* haplotyping: <yes | no >  */
#ifdef NEWVERSION
      options->haplotyping = (boolean) check;
      if(options->haplotyping)
	  {
	    get_next_word(&value,":",&tmp);// extract the word YES, value is now "reporting" or "non-reporting"
	    if(value[0]=='r')
	      options->haplotyping_report = TRUE;
	    else
	      options->haplotyping_report = FALSE;
	  }
#endif
      break;
    case 31: // auto-tune=<YES:acceptance-ratio | NO >
      options->has_autotune = (boolean) check;
      if(options->has_autotune)
	{
	  get_next_word(&value,":",&tmp);// extract the word YES, value is now interval:filename
	  get_next_word(&value,":",&tmp);// extract the interval, value is now filename
	  if(value==NULL)
	    {
	      options->autotune = 0.44;// using a value suggested by Roberts and Rosenthal (2009)
	    }
	  else
	    {
	      options->autotune = atof(value);
	    }
	}
      break;
    case 32: // assign=<YES | NO > or assign=YES:UNIFORM:{1,2}
      options->has_unassigned = (boolean) check;
      options->has_unassignedpoplist = FALSE;
      if (options->has_unassigned)
	{
	  printf("%i> option reading: %s\n",myID,value);
	  get_next_word(&value,":",&tmp);// extract the word YES, value is now FREQ or UNIFORM
	  if (value != NULL && (value[0]=='f' || value[0]=='F'))
	    options->has_unassignedfreq = TRUE;
	  else
	    options->has_unassignedfreq = FALSE;	  
	  get_next_word(&value,":",&tmp);// extract FREQ or UNIFORM, value can be now {1,2,...}
	  if (value != NULL)
	    {
	      get_next_word(&value,"{}:,; ",&tmp);
	      long z=0;
	      options->has_unassignedpoplist = TRUE;
	      options->unassignedpoplistnum = 1;
	      options->unassignedpoplist =  (long *) mycalloc(options->unassignedpoplistnum, sizeof(long));
	      while(tmp != NULL)
		{
		  if(z >= options->unassignedpoplistnum)
		    {
		      options->unassignedpoplistnum = z+1;
		      options->unassignedpoplist = (long *) myrealloc(options->unassignedpoplist, options->unassignedpoplistnum * sizeof(long));
		    }
		  options->unassignedpoplist[z] = atol(tmp);
		  printf("%i> option reading: unassignedpoplist[%li]=%li %s\n",myID,z, options->unassignedpoplist[z],tmp);
		  z++;
		  get_next_word(&value,"{}:,; ",&tmp);
		}	      
	    }
	}
      break;
    case 33: // bayes-hyperpriors=<YES | NO >
      options->hyperprior = (boolean) check;
      if (options->hyperprior)
	{
	  get_next_word(&value,":",&tmp);// extract the word YES, value is now interval:filename
	  get_next_word(&value,":",&tmp);// extract the interval, value is now filename
	  if(tmp==NULL)
	    {
	      options->hyperinterval = 1000;
	      options->hyperfactormean = 1.0;
	      options->hyperfactoralpha = 0.0;
	    }
	  else
	    {
	      options->hyperinterval = atol(tmp);
	      if(value==NULL)
		{
		  options->hyperfactormean=1.0;
		  options->hyperfactoralpha=0.0;
		}
	      else
		{
		  get_next_word(&value,":",&tmp);		 
		  options->hyperfactormean = atof(tmp);
		  if(value==NULL)
		    options->hyperfactoralpha=0.0;
		  else
		    options->hyperfactoralpha=atof(value);
		}
	    }
	}
      break;
    case 34:/*inheritance-scalars defines the scalars so that all reference*/
      /* is made to 4 Ne mu when the scalar is used, otherwise the scalar is*/
      /* assume to be 1*/
      options->has_inheritance = (boolean) check;
      if(options->has_inheritance)
	{
	  get_next_word(&value,":",&tmp);// extract the word YES, value is now "{x,x,x}"
	  get_next_word(&value,"{}:,; ",&tmp);
	  long z=0;
	  while(tmp != NULL)
	    {
	      if(z >= options->inheritance_scalars_numalloc)
		{
		  options->inheritance_scalars_numalloc = z+1;
		  options->inheritance_scalars = (MYREAL*) myrealloc(options->inheritance_scalars, options->inheritance_scalars_numalloc * sizeof(MYREAL));
		}
	      options->inheritance_scalars[z] = atof(tmp);
	      //printf("%i> option reading: inheritance scalar[%li]=%g %s\n",myID,z, options->inheritance_scalars[z],tmp);
	      z++;
	      get_next_word(&value,"{}:,; ",&tmp);
	    }
	}
      else
	{
	  options->has_inheritance = FALSE;
	  options->inheritance_scalars[0] = 1.0;
	}
      break;
    case 35: //mittag-leffler-alpha
      has_mlalpha = (boolean) check;
      if (!has_mlalpha)
	{
	  options->tri_mlalpha = NO;
	  options->mlalpha = 1.0;
	}
      else
	{
	  if (value!=NULL)
	    {
#ifdef NEWVERSION
	      get_next_word(&value,":",&tmp);
	      if (value==NULL)
		{
		  options->tri_mlalpha = ESTIMATE;
		  options->mlalpha = 1.0;
		  options->mlinheritance = 2.0;
		}
	      else
		{
		  if (value[0] == 'E' || value[0] == 'e')
		    {
		      options->tri_mlalpha = ESTIMATE;
		      options->mlalpha = 1.0;;
		      options->mlinheritance = 2.0; 
		    }
		  else
		    {
		      options->tri_mlalpha = FIXED;
		      options->mlalpha = atof(value);
		      options->mlinheritance = 2.0; //atof(value);
		    }
		}
	    }
#else
          options->tri_mlalpha = NO;
	  options->mlalpha = 1.0;
#endif
	}
      //else
      //	{
      //  if(options->tri_mlalpha)
      //    options->mlalpha = 1.0;
      //}
      break;
    default:
        return FALSE;
    }
    return TRUE;
}    /* booleancheck */


boolean
numbercheck (option_fmt * options, char *var, char *value)
{
  //int retval;
  long    tbins = BAYESNUMBIN;
  long    mbins = BAYESNUMBIN;
  long    sbins = BAYESNUMBIN;
  long    rbins = BAYESNUMBIN;
  
    MYREAL musum = 0., lastrate = 1.;
    long i = 0, z=0, cc = 0;
    double scale;
    boolean decided=FALSE;
    char *tmp, *tmp2, *temp, *temp2,  *keeptmp;
    char *numbertokens[NUMNUMBER] = NUMBERTOKENS;
    tmp2 = (char *) mycalloc (LINESIZE, sizeof (char));
    keeptmp = tmp2;
    while (i < NUMNUMBER && strcmp (var, numbertokens[i]))
        i++;
#ifdef DEBUG
    printf("%i> numbercheck: %s\n",myID, var);
#endif
    switch ((short) i)
    {
    case 0:   /*ttratio = value */
        z = 0;
        temp = strtok (value, " :,;\n\0");
        while (temp != NULL)
        {
            options->ttratio[z] = atof (temp);
            //options->ttratio[z] = 0.0;
	    options->sequence_model_parameters[z] = options->ttratio[z];
	    z++;
            options->ttratio =
	      (MYREAL *) myrealloc (options->ttratio, sizeof (MYREAL) * (size_t) z);
            temp = strtok (NULL, " ,;\n\0");
        }
        break;
    case 1: /*rate*/
      if(myID==MASTER)
	read_startparameter (2, options, var, value, NULL);
      break;  

    case 2: /*split*/
      if(myID==MASTER)
	read_startparameter (3, options, var, value, NULL);
      break;  
    case 3: /*splitstd*/
      if(myID==MASTER)
       	read_startparameter (4, options, var, value, NULL);
      break;  
      //case 1:   /*short-chains = value */
      //  options->schains = atol (value);
      //  break;
      //case 2:   /*short-steps = value */
      //case 32:   /* short-sample = value */
      //  options->ssteps = atol (value);
      //  break;
      //case 3:   /*short-increment = value */
      //  options->sincrement = atol (value);
      //  break;
    case 4:   /*long-chains = value */
        options->lchains = atol (value);
        break;
    case 5:   /*long-steps = value */
    case 33:   /*long-sample = value */
        options->lsteps = atol (value);
        break;
    case 6:   /*long-increment = value */
        options->lincrement = atol (value);
        break;
    case 7:  /* theta=; worker nodes read this outside of numbercheck, using a buffer */
      if(myID==MASTER)
	read_startparameter (0, options, var, value, NULL);
      break;  
    case 8:   /*nmlength = value */
        options->nmlength = atoi (value); //strtol (value, (char **) NULL, 10); 
        break;
    case 9:   /* seed = <Auto | seedfile | Own:value> */
        switch (value[0])
        {
        case 'A':
        case 'a':
        case '0':
            options->autoseed = AUTO;
            options->inseed = (unsigned long) (time (0) / 4 + 1);
            break;
        case 'S':
        case 's':
        case '1':
            options->autoseed = NOAUTO;
            openfile(&options->seedfile, options->seedfilename, "r", NULL);
            if (options->seedfile == NULL)
            {
                usererror ("cannot find seedfile\n");
            }
            fscanf (options->seedfile, "%ld%*[^\n]", &options->inseed);
            fclose (options->seedfile);
            break;
        case 'O':
        case 'o':
        case '2':
            options->autoseed = NOAUTOSELF;
            strtok (value, ":");
            tmp = strtok (NULL, ";\n");
            if (tmp != NULL)
	      options->inseed = (unsigned long) atol ((char *) tmp);
            if (options->inseed > 0)
                break;
        default:
            options->autoseed = AUTO;
            options->inseed = (unsigned long) time (0) / 4 + 1;
            usererror ("Failure to read seed method, should be\n \
                       random-seed=auto or random-seed=seedfile or random-seed=own:value\nwhere value is a positive integer\nUsing AUTOMATIC seed=%li\n", options->inseed);
            //break;
        }
        break;
    case 10:  /*"migration=", worker read this outside of numbercheck() */
      if(myID==MASTER)
	read_startparameter (1, options, var, value, NULL);
      break;  
    case 11:   /*mutation= <auto=gamma | nogamma | constant | estimate | own | data> */
        switch (value[0])
        {
        case 'A':  /*automatic */
        case 'a':
        case 'E':  /*automatic */
        case 'e':
	    options->bayesmurates=TRUE;
            options->murates = FALSE;
            options->murates_fromdata = FALSE;
	    options->gamma = FALSE;
	    break;
        case 'G':
        case 'g':
            options->gamma = TRUE;
	    options->bayesmurates=TRUE;
            options->murates = FALSE;
            options->murates_fromdata = FALSE;
            (void) strtok (value, " :");
            temp = strtok (NULL, " ,;\n");
            if (temp != NULL)
                options->alphavalue = atof (temp);
            else
                options->alphavalue = START_ALPHA;
            break;
        case 'O':
        case 'o':
            options->murates = TRUE;
            options->gamma = FALSE;
            options->murates_fromdata = FALSE;
	    options->bayesmurates = FALSE;
            (void) strtok (value, " :");
            temp = strtok (NULL, " ,;\n");
            if (temp != NULL)
            {
                options->muloci = atol (temp);
                options->mu_rates =
                    (MYREAL *) mycalloc (options->muloci, sizeof (MYREAL));
                musum = 0.;
                for (i = 0; i < options->muloci; i++)
                {
                    temp = strtok (NULL, " ,;\n");
                    if (temp == NULL)
                    {
                        while (i < options->muloci)
                        {
                            options->mu_rates[i] = lastrate;
                            musum += options->mu_rates[i];
                            i++;
                        }
                    }
                    lastrate = options->mu_rates[i] = atof (temp);
                    musum += options->mu_rates[i];
                }
                // mean must be 1.
                musum /= options->muloci;
                for (i = 0; i < options->muloci; i++)
                {
                    options->mu_rates[i] /= musum;
                }
            }
            break;
	case 'D': /* mutaton is varying and calculated from the data*/ 
	case 'd':
            options->gamma = FALSE;
            options->murates = TRUE;
            options->murates_fromdata = TRUE;
	    options->bayesmurates = FALSE;
	    break;
        case 'N':  /*nogamma, none, all loci have same mu */
        case 'n':
	case 'C':
	case 'c':
        default:
            options->murates = FALSE;
            options->gamma = FALSE;
            options->murates_fromdata = FALSE;
	    options->bayesmurates = FALSE;
            break;
        }
        break;
    case 12:   /*datatype=<allele|microsatellite|brownian|sequence|f-ancestral states|genealogies> */
        switch (value[0])
        {
        case 'a':
        case 'A':
            options->datatype = 'a';
            break;
        case 'm':
        case 'M':
            options->datatype = 'm';
            break;
        case 'b':
        case 'B':
            options->datatype = 'b';
            break;
        case 's':
        case 'S':
            options->datatype = 's';
            break;
        case 'n':
        case 'N':
            options->datatype = 'n';
            break;
        case 'h':
        case 'H':
            options->datatype = 'h';
	    options->fastlike=FALSE;
            break;
        case 'u':
        case 'U':
            options->datatype = 'u';
            break;
        case 'f':
        case 'F':
            options->datatype = 'f';
            break;
        case 'g':
        case 'G':
            options->datatype = 'g';
            options->readsum = TRUE;
            break;
        default:
            options->datatype = 's';
            break;
        }
        break;
    case 13:   /* categories=<None | value> */
        if (uppercase (value[0] == 'N'))
        {
            options->categs = ONECATEG;
            break;
        }
        else
        {
            options->categs = strtol (value, (char **) NULL, 10);
            /* needs to read auxilliary file catfile */
            sprintf(tmp2,"%li",options->categs);
            set_filename(value, tmp2, &options->catfilename);
        }
        break;
    case 14:   /*create rates=value:list of rates */
        strncpy (tmp2, value, strcspn (value, ":"));
        if (strtol (tmp2, (char **) NULL, 10) /*;atoi (tmp) */  > 1)
        {   /* rate categories */
            options->rcategs = strtol (tmp2, (char **) NULL, 10);
            options->rrate =
                (MYREAL *) myrealloc (options->rrate,
				      sizeof (MYREAL) * (size_t) (options->rcategs + 1));
            (void) strtok (value, " :");
            temp = strtok (NULL, " ,;\n");
            z = 0;
            while (temp != NULL)
            {
                if (z > options->rcategs)
                    usererror ("check parmfile-option  rates, missing rate\n");
                options->rrate[z++] = atof (temp);
                temp = strtok (NULL, " ,;\n");
            }
        }
        break;
    case 15:   /* probabilities for each rate category */
        strncpy (tmp2, value, strcspn (value, ":"));
        if (strtol (tmp2, (char **) NULL, 10) > 1)
        {   /* probabilities for each rate category */
            options->rcategs = strtol (tmp2, (char **) NULL, 10);
            options->probcat =
                (MYREAL *) myrealloc (options->probcat,
				      sizeof (MYREAL) * (size_t) (options->rcategs + 1));
            (void) strtok (value, " :");
            temp = strtok (NULL, " ,;\n");
            z = 0;
            while (temp != NULL)
            {
                if (z > options->rcategs)
                    usererror
                    ("check parmfile prob-rates, missing rate probability\n");
                options->probcat[z++] = atof (temp);
                temp = strtok (NULL, " ,;\n");
            }
        }
        break;
    case 16:   /*micro-stepmax */
        // options->micro_stepnum = strtol (value, (char **) NULL, 10);

        break;
    case 17:   /*micro-threshold */
      options->micro_threshold = atol(value);
      if(options->micro_threshold % 2 != 0)
	options->micro_threshold += 1;
      break;
    case 18:   /*delimiter */
      options->dlm = value[0];
      break;
    case 19:   /*burn-in */
      options->burn_in = atol (value);
      if(strchr(value,'a')!=NULL)
	options->burnin_autostop = 'a';
      else
	{

	  if(strchr(value,'t')!=NULL)
	    {
	      options->burnin_autostop = 't';
	    }
	  else

	  if(strchr(value,'e')!=NULL)
	    {
	      options->burnin_autostop = 'e';
	    }
	  else
	    options->burnin_autostop = ' ';
	}
      break;
    case 20:   /*infilename */
      get_filename (&options->infilename, value);
      break;
    case 21:   /*outfilename */
      get_filename(&options->outfilename,value);
        break;
    case 22:   /*smoothing*/
      // smoothing= {<BIN | KERNEL<:el>> , < NONE | SAVITZKY<:el> | AVERAGE<:el> >}
      // el is the window size --> 2*el+1
      // example:
      // smoothing={BIN,AVERAGE:20}
      // smoothing={BIN,NONE}
      // smoothing={KERNEL:20,SAVITZKY:20}  //default
      // smoothing={KERNEL,SAVITZKY}  //default using el=20 --> window=41
      if (value[0]=='{')
	temp = strtok (value+1, ",};\n");
      else
	temp = strtok (value, ",};\n");
      if (temp != NULL)
	{
	  char choice = uppercase (temp[0]);
	  switch(choice)
	    {
	    case 'B': // using the oldstyle binning system
	      options->smoothing[0]=BINNING;
	      break;
	    case 'K': // using kernel-smoothing
	    default:
	      options->smoothing[0]=KERNEL;
	      temp2 = strchr(temp,':');
	      if (temp2!=NULL)
		options->smoothing_window[0] = atol(temp2+1);
	      break;
	    }
	  temp = strtok (NULL, ",:;}\n");
	  if (temp != NULL)
	    {
	      choice = uppercase (temp[0]);
	      switch(choice)
		{
		case 'A': // using bayes_smooth, whic his  moving average
		  options->smoothing[1]=AVERAGE;
		  temp2 = strtok (NULL, ",:;}\n");
		  if (temp2!=NULL)
		    options->smoothing_window[1] = atol(temp2);
		  break;
		case 'N':
		  options->smoothing[1]=NONE;
		  options->smoothing_window[0] = DEFAULTKERNELWINDOW;
		  break;
		case 'S': // using kernel-smoothing
		default:
		  options->smoothing[1]=SAVITZKY;
		  temp2 = strtok (NULL, " ,:;}\n");
		  if (temp2!=NULL)
		    options->smoothing_window[1] = atol(temp2);
		  break;
		}
	    }
	}
      else
	{
	  options->smoothing[0]=KERNEL;
	  options->smoothing[1]=SAVITZKY;
	  options->smoothing_window[0] = DEFAULTKERNELWINDOW;
	  options->smoothing_window[1] = DEFAULTKERNELWINDOW;
	}
      break;
    case 23:   /*title */
        strncpy (options->title, value, 80);
        break;
    case 24:   /*long-chain-epsilon */
        options->lcepsilon = atof (value);
        if (options->lcepsilon <= 0)
            options->lcepsilon = LONGCHAINEPSILON;
        break;
    case 25:   /* print-tree options */
        switch (uppercase (value[0]))
        {
        case 'N':
            options->treeprint = myNONE;
            break;
        case 'A':
            options->treeprint = ALL;
	    options->treeinc = 1;
	    set_filename(value, "ALL", &options->treefilename);
            break;
        case 'B':
            options->treeprint = BEST;
	    set_filename(value, "BEST", &options->treefilename);
            break;
        case 'L':
            options->treeprint = LASTCHAIN;
	    get_next_word(&value,":",&tmp);// the word LASTCHAIN
	    get_next_word(&value,":",&tmp);
	    options->treeinc = atol(tmp);
	    set_filename(value, "LASTCHAIN", &options->treefilename);
            break;
        default:
            options->treeprint = myNONE;
            break;
        }
        break;
    case 26:   /* progress: No, Yes, Verbose */
        switch (uppercase (value[0]))
        {
        case 'F':
        case 'N':
            options->progress = FALSE;
            options->verbose = FALSE;
            break;
        case 'T':
        case 'Y':
            options->progress = TRUE;
            options->verbose = FALSE;
            break;
        case 'V':
            options->progress = TRUE;
            options->verbose = TRUE;
            break;
        }
        break;
    case 27:   /* l-ratio: <NO | YES>:val1,val2,val3,val4,val5 */
        cc = options->lratio->counter;
        switch (uppercase (value[0]))
        {
        case 'Y':
            options->lratio->data[cc].type = MLE;
            break;
        case 'N':
        default:
	  myfree(keeptmp);
            return FALSE;
        }
        (void) strtok (value, ":");
        temp = strtok (NULL, "\n");
        if (temp != NULL)
        {
            temp2 = strchr (temp, ':');
            if (temp2 != NULL)
            {
                strcpy (options->lratio->data[cc].value2, temp2);
                *temp2 = '\0';
                strcpy (options->lratio->data[cc].value1, temp);
            }
            else
                strcpy (options->lratio->data[cc].value1, temp);
        }
        if (cc + 1 == options->lratio->alloccounter)
        {
            options->lratio->alloccounter += 2;
            options->lratio->data =
                (lr_data_fmt *) myrealloc (options->lratio->data,
                                         sizeof (lr_data_fmt) *
					   (size_t) (options->lratio->alloccounter+1));
            for (i = cc + 1; i < options->lratio->alloccounter; i++)
            {
                options->lratio->data[i].elem = 0;
                options->lratio->data[i].value1 =
                    (char *) mycalloc (1, sizeof (char) * LONGLINESIZE);
                options->lratio->data[i].elem = 0;
                options->lratio->data[i].value2 =
                    (char *) mycalloc (1, sizeof (char) * LONGLINESIZE);
                options->lratio->data[i].connect =
                    (char *) mycalloc (1, sizeof (char) * LONGLINESIZE);

            }
        }
      	options->lratio->counter += 1;
        break;
    case 28:   /* fst-type: <Theta | Migration> */
        switch (uppercase (value[0]))
        {
        case 'T':
            options->fsttype = 'T';
            break;
        case 'M':
        default:
            options->fsttype = 'M';
            break;
        }

        break;
    case 29:   /*profile=<NO| NONE | YES | ALL | TABLES | SUMMARY>><: <FAST |  */
      warning("profile option is defunct");
      break;
    case 30:   /* custom-migration:<{> migration matrix and theta on
                           diagonal:
                           0 means not estimated,
                           x means estimated, s means symmetrically
                           estimated, m means all are the same  <}> */
        if (myID == MASTER)
            read_custom_migration (options->parmfile, options, value,
                                   options->numpop,0);
#ifdef MPI

        else
            read_custom_migration_worker (options->buffer, options, value,
                                          options->numpop);
#endif

        break;
    case 31:   /* population-growth */
      set_growth(&value, &tmp, options);
        break;
        /*case 32 and case 33 are fallthroughs to 2 and 3 */
    case 34:   /*replicate */
        switch (uppercase (value[0]))
        {
        case 'T':
        case 'Y':
            options->replicate = TRUE;
            temp = strtok (value, ":;\n");
            if (temp != NULL)
            {
                temp = strtok (NULL, ":;\n");
                if (uppercase (temp[0]) == 'L')
                    options->replicatenum = 0;
                else
                    options->replicatenum = strtol (temp, (char **) NULL, 10);
            }
            else
                options->replicatenum = 0;
            break;
        default:
            options->replicate = FALSE;
            options->replicatenum = 0;
        }
        break;
    case 35:   /* datamodel: JC, K2P, F84, F81 (D), HKY, TN  */
      options->datamodel = get_mutationmodel(value[0]);
      options->sequence_model = options->datamodel;
        break;
    case 36:   /* logfile=<YES:logfile | NO> do we write a logfile or not */
        options->writelog = set_filename(value,"YES",&options->logfilename);
        break;
    case 37:   /* sequencing error  seqerror-rate=Estimate::1|4> , seqerror-rate=0.01, seqerror-rate={0.01,0.02,003,0.01}*/
      if (value[0]=='E' || value[0]=='e')
	{
	  options->seqerror[0] = 0.0001;
	  options->seqerror[1] = 0.0001;
	  options->seqerror[2] = 0.0001;
	  options->seqerror[3] = 0.0001;
	  options->has_estimateseqerror = TRUE;
	  (void) strtok (value, ",:;\n");//extract ESTIMATE
	  temp = strtok (NULL, ",:;\n");
	  if (temp!=NULL)
	    {
	      if (atol(temp)==1)
		options->seqerrorcombined=TRUE;
	      else
		options->seqerrorcombined=FALSE;
	    }
	  else
	    options->seqerrorcombined=FALSE;
	}
      else
	{
	  options->has_estimateseqerror = FALSE;
	  if(value[0]=='{')
	    {
	      temp = strtok (value, ",:;\n");
	      if (temp != NULL)
		{
		  options->seqerror[0] = atof (temp);
		  temp = strtok (NULL, ",:;\n");
		  if (temp != NULL)
		    {
		      options->seqerror[1] = atof (temp);
		      temp = strtok (NULL, ",:;\n");
		      if (temp != NULL)
			{
			  options->seqerror[2] = atof (temp);
			  temp = strtok (value, ",:;\n");
			  if (temp != NULL)
			    {
			      options->seqerror[3] = atof (temp);
			    }
			}
		    }
		}    
	    }
	  else
	    {
	      options->seqerror[0] = atof (value);
	      options->seqerror[1]=options->seqerror[2]=options->seqerror[3]=options->seqerror[0];
	    }
	  if ((options->seqerror[0] < 0.0) || (options->seqerror[1] < 0.0) || (options->seqerror[2] < 0.0) || (options->seqerror[3] < 0.0))
	    {
	      warning
		("Sequencing error was misspecified in parmfile, reset to 0.0\n");
	      options->seqerror[0] = 0.0;
	      options->seqerror[1]=options->seqerror[2]=options->seqerror[3]=0.0;
	    }
	}
        break;
#ifdef UEP

    case 38:   /* do we have a uep file or not, function returns TRUE if uep=YES:filename*/
        options->uep = set_filename(value, "YES", &options->uepfilename);
        break;
    case 39:   /* do we have uep-rates */
        temp = strtok (value, ":;\n");
        if (temp != NULL)
        {
            options->uepmu = atof (temp);
            temp = strtok (NULL, ":;\n");
            if (temp != NULL)
            {
                options->uepnu = atof (temp);
            }
            options->ueprate = options->uepmu;
        }
        break;
    case 40:   /* do we have uep-bases */
        temp = strtok (value, ":; \n");
        if (temp != NULL)
        {
            options->uepfreq1 = atof (temp);
            temp = strtok (NULL, " :;\n");
            if (temp != NULL)
            {
                options->uepfreq0 = atof (temp);
            }
            else
                options->uepfreq0=1. - options->uepfreq1;
        }
        break;
#endif

    case 41: // mu-rates??????
        break;
    case 42: /* heating=<no | <yes | adaptive | bounded>:numintervals:{temperatures}> */
        switch (uppercase (value[0]))
        {
        case 'A': //adaptive heating on
            options->heating = 1;
            options->adaptiveheat = STANDARD;
            break;
        case 'B': //adaptive heating on
            options->heating = 1;
            options->adaptiveheat = BOUNDED;
            break;
        case 'Y':
        case 'P':
            options->heating = 1;
            options->adaptiveheat = NOTADAPTIVE;
            break;
        case 'N':
            options->heating=0;
            options->adaptiveheat=NOTADAPTIVE;
	    options->heated_chains = 1;
            break;
        default:
            options->heating = 1;
            options->adaptiveheat = NOTADAPTIVE;
            break;
        }
        if (options->heating == 1)
        {
            strtok (value, ":\n");
            tmp = strtok (NULL, ": ");
            if (tmp != NULL)
            {
                options->heating_interval = atol (tmp);
                tmp = strtok (NULL, "{, ");
                if (tmp != NULL)
                {
                    z = 0;
                    while (1)
                    {
                        options->heat[z++] = atof (tmp);
                        tmp = strtok (NULL, ", :\n");
                        if (tmp == NULL || z >= 1000)
			  {
			    if(z>=1000)
			      warning("Ignored the all temperature above 1000 heated chains\n");
                            break;
			  }
                    }
                    options->heated_chains = z;
                }
            }
        }
        break;
#ifdef LONGSUM

    case 43: /* fluctuate=<no | <yes>:{rate_1today, time1today,rate_1middle,time1middle, rate_1past, time1past,rate2_today,time2_today,...}> */
        switch (uppercase (value[0]))
        {
        case 'Y':
        case 'P':
            options->fluctuate = TRUE;
            break;
        case 'N':
        default:
            options->fluctuate = FALSE;
            break;
        }
        if (options->fluctuate)
        {
            strtok (value, ":\n");
            tmp = strtok (NULL, ":{,");
            if (tmp != NULL)
            {
                z = 0;
                while (1)
                {
                    options->flucrates = myrealloc(options->flucrates, sizeof(MYREAL) * (z+1));
                    options->flucrates[z++] = atof (tmp);
                    tmp = strtok (NULL, ", {}:\n");
                    if (tmp == NULL)
                        break;
                }
            }
        }
        if(z%3!=0)
            error("Fluctuating rates and times need to be 3 rates and 3 times per population");
        break;
#endif /*LONGSUM*/
    case 44:   /*resistance = value  [to fatal attraction to zero*/
    	//xcode       z = 0;
        temp = strtok (value, " ,;\n\0");
        if (temp != NULL)
            options->minmigsumstat = atof (temp);
        else
            options->minmigsumstat = MINMIGSUMSTAT;
        break;
	//    case 45: /* bayes-updatefreq=[0..1] */
        //temp = strtok (value, " ,;\n\0");
        //if (temp != NULL)
        //    options->updateratio = atof (temp);
        //else
        //    options->updateratio= HALF;
        //break;
    case 65: /* bayes-updatefreq from old parmfile will go here too*/
    case 45: /* updatefreq= treeupdate parameteterupdate haplotypeupdate timeparamupdate unassignedupdate seqerrorupdate mlalphaupdate*/
        temp = strtok (value, ";\n\r\0");
	// temp is a string containing maximal 4 values
        if (temp != NULL)
	  {
	    int elements =sscanf(temp,"%lf%lf%lf%lf%lf%lf%lf", 
				 &options->tree_updatefreq,
				 &options->parameter_updatefreq,
				 &options->haplotype_updatefreq,
			         &options->timeparam_updatefreq,
				 &options->unassigned_updatefreq,
				 &options->seqerror_updatefreq,
				 &options->mlalpha_updatefreq
				 );

	  }
        break;
    case 49: /*bayes-posteriorbins*/    
      tbins = BAYESNUMBIN;
      mbins = BAYESNUMBIN;
      rbins = BAYESNUMBIN;
      sbins = BAYESNUMBIN;
      temp = strtok (value, " ,;\n\0");
      if (temp != NULL)
        {
	  tbins = atol (temp);
	  temp = strtok (NULL, " ,;\n\0");
	  if (temp != NULL)
	      {
		mbins = atol (temp);
		temp = strtok (NULL, " ,;\n\0");
		if (temp != NULL)
		  {
		    rbins = atol (temp);
		    temp = strtok (NULL, " ,;\n\0");
		    if (temp != NULL)
		      sbins = atol (temp);
		  }
	      }
	}
      if (tbins>0)
	options->bayes_posterior_bins[0]=tbins;
      if (mbins>0)
	options->bayes_posterior_bins[1]=mbins;
      if (rbins>0)
	options->bayes_posterior_bins[2]=rbins;
      if (sbins>0)
	{
	  options->bayes_posterior_bins[3]=sbins;
	  options->bayes_posterior_bins[4]=sbins;
	}
#ifdef DEBUG
      printf("%i> read from parmfile: bayesposteriobins %li %li %li %li %li\n\n\n\n\n\n",
	     myID,
	     options->bayes_posterior_bins[0],
	     options->bayes_posterior_bins[1],
	     options->bayes_posterior_bins[2],
	     options->bayes_posterior_bins[3],
	     options->bayes_posterior_bins[4]);
#endif 
      break;
    case 46:   /*bayesfile=FILENAME  set bayesfile and bayes analysis? */
      strncpy (options->bayesfilename, value, 255); //this is legacy code and got replaced
      // by bayes-file=<YES | NO>:FILENAME
      options->has_bayesfile = TRUE;
	    break;
		
    case 47:   /*set bayes-prior*/
      set_bayes_options(value, options);//same as number 54
      break;
    case 48:   /* usertree =  < NO |  UPGMA | AUTOMATIC | TREE:intreefilename | RANDOM | DISTANCE:distfilename */
        options->usertree = set_filename(value, "T", &options->utreefilename); //USERTREE
        options->randomtree = (value[0] =='R') ? TRUE : FALSE; //RANDOMTREE
        options->dist = set_filename(value, "D", &options->distfilename); //DISTANCE
        break;
	/* case 49 is further back*/
    case 50: /* mig-histogram=<<yes|all>:histogram-binsize:filename | no> */
      get_next_word(&value,":",&tmp);
      options->mighist_all = ((uppercase(tmp[0]) == 'A') ? TRUE : FALSE);
      options->mighist     = ((uppercase(tmp[0]) == 'Y') ? TRUE : FALSE);
      if(options->mighist_all || options->mighist)
	{
	  options->mighist = TRUE;
	  get_next_word(&value,":",&tmp);
	  if(value!=NULL)
	    {
	      options->eventbinsize = (float) atof(tmp);
	      strncpy (options->mighistfilename, value, 255);
	    }
	  else
	    {
	      strncpy (options->mighistfilename, tmp, 255);
	    }
	}
      break;
    case 51:
        switch(toupper(value[0]))
           {
	   case 'A':
	     options->bayespretty = PRETTY_MAX;
	     break;
	   case 'P':
	     options->bayespretty = PRETTY_P99;
	     break;
	   case 'M':
	     options->bayespretty = PRETTY_P99MAX;
	     break;
	   case 'T':
	   default:
	     options->bayespretty = PRETTY_P100;
	     break;
	   }
	break;
#ifdef PRETTY
    case 52:   /*PDF-outfile name */
        strcpy (options->pdfoutfilename, value);
        break;
#endif
    case 53:   /*bayes interval for writing all parameters to file */
        options->bayesmdiminterval = atol(value);
        break;
    case 54:   /*set bayes-priorS*/
      // new parmfile variable so that the old stype setting still work
      set_bayes_options(value, options);
      break;
    case 55: /*set skyline-histogram*/
      get_next_word(&value,":",&tmp); 
      options->skyline = (tmp[0] == 'Y' ? TRUE : FALSE);
      if (tmp[0]=='P')
	{
	  options->skyline = TRUE;
	  options->skyline_param = TRUE;
	  get_next_word(&value,":",&tmp);
	  options->timeelements = atol(tmp);
	}
      if(options->skyline)
	{
	  get_next_word(&value,":",&tmp);
	  options->eventbinsize = (float) atof(tmp);
	  strncpy (options->skylinefilename, value, 255);
	}
      break;
    case 56: /* rates-gamma */
      get_next_word(&value,":,; ",&tmp);
      options->seqrate_gamma_num = 0;
      while(tmp!=NULL)
	{
	  options->seqrate_gamma[options->seqrate_gamma_num++] = atof(tmp);
	}
      break;
    case 57: /* Bayes-proposals*/
      get_next_word(&value,":,; ",&tmp);
      if(tmp != NULL)
	{
	  switch(uppercase(tmp[0]))
	    {
	    case 'T':       get_next_word(&value,":,; ",&tmp);
	      if(tmp != NULL)
		{
		  if(uppercase(tmp[0])=='S')
		    options->slice_sampling[THETAPRIOR] = TRUE;
		  else
		    options->slice_sampling[THETAPRIOR] = FALSE;
		}
	      break;
	    case 'M':        get_next_word(&value,":,; ",&tmp);
	      if(tmp != NULL)
		{
		  if(uppercase(tmp[0])=='S')
		    options->slice_sampling[MIGPRIOR] = TRUE;
		  else
		    options->slice_sampling[MIGPRIOR] = FALSE;
		}
	      break;
	    case 'D':
	      if (strstr(tmp,"STD")==NULL)
		{
		  get_next_word(&value,":,; ",&tmp);
		  if(tmp != NULL)
		    {
		      if(uppercase(tmp[0])=='S')
			options->slice_sampling[SPECIESTIMEPRIOR] = TRUE;
		      else
			options->slice_sampling[SPECIESTIMEPRIOR] = FALSE;
		    }
		}
	      else
		{
		  get_next_word(&value,":,; ",&tmp);
		  if(tmp != NULL)
		    {
		      if(uppercase(tmp[0])=='S')
			options->slice_sampling[SPECIESSTDPRIOR] = TRUE;
		      else
			options->slice_sampling[SPECIESSTDPRIOR] = FALSE;
		    }
		}	      		  
	      break;
	    case 'R':       get_next_word(&value,":,; ",&tmp);
	      if(tmp != NULL)
		{
		  if(uppercase(tmp[0])=='S')
		    options->slice_sampling[RATEPRIOR] = TRUE;
		  else
		    options->slice_sampling[RATEPRIOR] = FALSE;
		}
	      break;
	    case 'G':       get_next_word(&value,":,; ",&tmp);
	      if(tmp != NULL)
		{
		  if(uppercase(tmp[0])=='S')
		    options->slice_sampling[GROWTHPRIOR] = TRUE;
		  else
		    options->slice_sampling[GROWTHPRIOR] = FALSE;
		}
	      break;
	    }
	}
      break;
    case 58: /*generation-per-year*/
      get_next_word(&value,":,; ",&tmp);
      options->generation_year = atof(tmp);
      break;
    case 59: /*mutationrate-per-year absolute mutation rate per year for each locus*/
      options->mutationrate_year[0] = 1.0;
      get_next_word(&value,"{}:,; ",&tmp);
      z=0;
      while(tmp != NULL)
	{
	  if(z >= options->mutationrate_year_numalloc)
	    {
	      options->mutationrate_year_numalloc = z+1;
	      options->mutationrate_year = (MYREAL*) myrealloc(options->mutationrate_year, options->mutationrate_year_numalloc * sizeof(MYREAL));
	    }
	  options->mutationrate_year[z] = atof(tmp);
	  printf("%i> mutationrate_year[%li]=%g %s\n",myID,z, options->mutationrate_year[z],tmp);
	  z++;
	  get_next_word(&value,"{}:,; ",&tmp);
	}
      break;
    case 60: /*botched inheritance-scalars without YES/NO*/
      get_next_word(&value,"{}:,; ",&tmp); /*this reads the {} structure*/
      z=0;
      decided = FALSE;
      options->has_inheritance = FALSE;
      while(tmp != NULL)
	{
	  if(z >= options->inheritance_scalars_numalloc)
	    {
	      options->inheritance_scalars_numalloc = z+1;
	      options->inheritance_scalars = (MYREAL*) myrealloc(options->inheritance_scalars, options->inheritance_scalars_numalloc * sizeof(MYREAL));
	    }
	  scale = atof(tmp);
	  if (!decided)
	    if (fabs(scale-1.0) > EPSILON)
	      {
		options->has_inheritance = TRUE;
		decided = TRUE;
	      }
	  options->inheritance_scalars[z] = scale;
	  //printf("%i> second version option reading: inheritance scalar[%li]=%g %s\n",myID,z, options->inheritance_scalars[z],tmp);
	  z++;
	  get_next_word(&value,"{}:,; ",&tmp);
	}      
      break;
    case 61: /* micro-submodel */
      get_next_word(&value,":",&tmp); 
      if(tmp==NULL)
	options->msat_option = (int) atol(value);
      else
	{
	  options->msat_option = (int) atol(tmp);
	  if(value==NULL)
	    break;
	    get_next_word(&value,", ",&tmp);
	  if(tmp[0]=='{')
	    options->msat_tuning[0] = atof(tmp+1);
	  else
	    options->msat_tuning[0] = atof(tmp);
	  get_next_word(&value," }",&tmp);
	  options->msat_tuning[1] = atof(tmp);
	}
      break;

    case 62: /*random-subset*/
      get_next_word(&value,":,; ",&tmp);
      options->randomsubset = atol(tmp);
      get_next_word(&value,":,; ",&tmp);
      if (tmp != NULL)
	options->randomsubsetseed = (unsigned long) atol(tmp);
      break;

    case 63: /* pooling of populations an array of numbers indicates the pooling*/
      set_localities(&value, &tmp, options);
      break;
    case 64: /* sequence-submodel */
      get_next_word(&value,":",&tmp); 
      if(tmp==NULL)
	options->sequence_option = (int) atol(value);
      else
	{
	  options->sequence_option = (int) atol(tmp);
	  if(value==NULL)
	    break;
	  get_next_word(&value,", ",&tmp);
	  if(tmp[0]=='{')
	    options->sequence_model_parameters[0] = atof(tmp+1);
	  else
	    options->sequence_model_parameters[0] = atof(tmp);
	  if(value==NULL)
	    break;
	  get_next_word(&value," ,}",&tmp);
	  options->sequence_model_parameters[1] = atof(tmp);
	  if(value==NULL)
	    break;
	  get_next_word(&value," }",&tmp);
	  options->sequence_model_parameters[2] = atof(tmp);
	}
      break;
    case 66:/* analyze-loci=<All | Variable | Fast > */
      if (value!=NULL)
	{
	  switch(value[0])
	    {
	    case 'V': // analyze only variable loci
	      options->onlyvariable = TRUE;	  
	      break;
	    case 'F': // analyze all variable and one invariant weighted to speed up run
	      options->onlyvariable = FALSE;
	      options->has_variableandone = TRUE;
	      break;
	    case 'A': //analyze all, default
	    default:
	      options->onlyvariable = FALSE;
	      break;
	    }
	}
      break;
    case 67:/* regress on April 2021 to two choices E | S 
	       S is labeled N for users, the real N will we resurrected some time later, 
	       divergence-distrib=<E | W | N | S > */
      if (value!=NULL)
	{
	  switch(value[0])
	    {
	      //case 'W': 
	      //options->species_model_dist = WEIBULL_DIST;	  
	      //break;
	    case 'N': 
	      options->species_model_dist = NORMAL_DIST;	  
	      break;
	    case 'E':	      
	      options->species_model_dist = EXP_DIST;	  
	      break;
	      //case 'N':
	    case 'W':
	    case 'S':
	      //needs fixing it seems to call the lambda' instead of lambda
	      //and does not deliver good results in simulation test!
	    default:
	      options->species_model_dist = NORMAL_DIST;	  
	      //options->species_model_dist = NORMALSHORTCUT_DIST;	  
	      break;
	    }
	}
      break;
    case 68: /*mittag-leffler-alpha=value*/
      //mittag-leffler
      if (value!=NULL)
	{
#ifdef NEWVERSION
	  get_next_word(&value,":",&tmp);
	  if (value==NULL)
	    {
	      //if (tmp[0] === 'E' || tmp == 'e')
	      options->tri_mlalpha = FIXED;
	      options->mlalpha = atof(tmp);
	      options->mlinheritance = 2.0;
	    }
	  else
	    {
	      options->tri_mlalpha = FIXED;
	      options->mlalpha = atof(tmp);
	      options->mlinheritance = atof(value);
	    }
#else
	  options->tri_mlalpha = NO;
	  options->mlalpha = 1.0;
#endif
	}
      else
	{
	  options->tri_mlalpha = NO;
	  options->mlalpha = 1.0;
	}
      break;
    default:
      myfree(keeptmp);
      return FALSE;      
    }
    myfree(keeptmp);
    
    return TRUE;
}    /* numbercheck */

/*void
reset_oneline (option_fmt * options, long position)
{
    fseek (options->parmfile, position, SEEK_SET);
    }*/

void
reset_oneline (option_fmt * options, fpos_t * position)
{
    fsetpos (options->parmfile, position);
}

void read_random_theta(option_fmt *options, char ** buffer)
{
   char *tempstr;
   char *tmp;
   char *otempstr;
   //char *otmp;
   //char *retval;

   tempstr = (char *) mycalloc(LINESIZE,sizeof(char));
   //tmp = (char *) mycalloc(LINESIZE,sizeof(char));
   otempstr = tempstr;
   //otmp = tmp;

  options->numthetag = 2;
  options->thetag =
    (MYREAL *) myrealloc (options->thetag, sizeof (MYREAL) * 3);
  if( myID == MASTER)
    {
      fgets (tempstr, LINESIZE, options->parmfile);
    }
  else
    {
      sgets (tempstr, LINESIZE, buffer);
    }
  get_next_word(&tempstr,",{} \n",&tmp);
  while(strchr("{ ", (tmp)[0]))
    {
      get_next_word(&tempstr,",{} \n",&tmp);
    }
  options->thetag[0] = atof(tmp);
  get_next_word(&tempstr,",{} \n",&tmp);
  options->thetag[1] = atof(tmp);
  //  myfree(otmp);
  myfree(otempstr);
}



void read_random_mig(option_fmt *options, char ** buffer)
{
   char *tempstr;
   char *tmp;
   char *otempstr;
   char *otmp;
   //char *retval;
   tempstr = (char *) mycalloc(LINESIZE,sizeof(char));
   tmp = (char *) mycalloc(LINESIZE,sizeof(char));
   otempstr = tempstr;
   otmp = tmp;
  options->nummg = 2;
  options->mg = (MYREAL *) myrealloc (options->mg, sizeof (MYREAL) * 3);
  if(myID == MASTER)
    {
      fgets (tmp, LINESIZE, options->parmfile);
    }
  else
    {
      sgets (tmp, LINESIZE, buffer);
    }
  get_next_word(&tmp,",{} \n",&tempstr);
  options->mg[0] = atof(tempstr);
  get_next_word(&tmp,",{} \n",&tempstr);
  options->mg[1] = atof(tempstr);
  myfree(otmp);
  myfree(otempstr);
}


///
/// read the values for start theta from the list in the parmfile
/// called through read_options_master() and also equivalent to read worker
void read_theta (option_fmt * options, char *parmvar, char *varvalue, char **buffer)
{
  (void) parmvar;
   long start = -1;
   long stop = -1;
   char *tempstr;
   char *tmp;
   char *otempstr;
   char *otmp;
   char *v = varvalue;
   tempstr = (char *) mycalloc(LINESIZE,sizeof(char));
   tmp = (char *) mycalloc(LINESIZE,sizeof(char));
   otempstr = tempstr;
   otmp = tmp;
   switch (uppercase (varvalue[0]))
     {
     case 'R':
       options->startguess[0][0] = RANDOMPRIOR;
       break;
     case 'P':
       options->startguess[0][0] = PRIOR;
       get_next_word(&v,":",&tmp);	   
       get_next_word(&v,",{}\n",&tmp);
       while(strchr("{", (tmp)[0]))
	 {
	   get_next_word(&v,",{} \n",&tmp);
	 }
       options->startguess[0][1] = atol(tmp);
       break;
     case 'O':
     case '0':
       options->startguess[0][0] = OWN;
       start = count_char(varvalue,'{');
       stop =  count_char(varvalue,'}');
       if (start == 0 && stop == 0)
	 {
	   get_next_word(&v,":",&tmp);// extract the word OWN: , into tmp
	   add_startparam(options,THETAPRIOR, (float) atof(v));
	 }
       else /* not on one line or using {} notation*/
	 {
	   if (start == 1 && stop == 1)  //all on one line
	     {
	       get_next_word(&v,"{",&tmp);// extract the word OWN: {, into tmp
	       while (v!= NULL)
		 {
		   get_next_word(&v," ,}",&tmp);
		   if(tmp != NULL)
		     add_startparam(options,THETAPRIOR, (float) atof(tmp));
		 }
	     }
	   else // multiple lines
	     {
	       if (start == 1 && stop == 0 )
		 {
		   get_next_word(&v,"{",&tmp);// extract the word OWN: {, into tmp
		   while (v!= NULL)
		     {
		       get_next_word(&v," ,}",&tmp);
		       if (tmp != NULL)
			 add_startparam(options,THETAPRIOR, (float) atof(tmp));
		     }
		   // read the next line 
		   while(stop==0)
		     {
		       if(myID==MASTER)
			 FGETS (varvalue, LINESIZE, options->parmfile);
		       else
			 sgets (varvalue, LINESIZE, buffer);
		     }
		   while (varvalue[0]=='#')
		     {
		       if(myID==MASTER)
			 FGETS (varvalue, LINESIZE, options->parmfile);
		       else
			 sgets (varvalue, LINESIZE, buffer);
		     }
		   //stop =  
		     count_char(varvalue,'}');
		   v = varvalue;
		   while (v!= NULL)
		     {
		       get_next_word(&v," ,}",&tmp);
		       if (tmp != NULL)
			 add_startparam(options,THETAPRIOR, (float) atof(tmp));
		     }		   
		 }
	       else
		 {
	       usererror("theta=OWN starvale not correctly formatted");
		 }
	     }
	 }
       options->numpop = options->startparam.numpop;
       break;
     default:
#ifdef MPI
       fprintf(stderr,"%i>",myID);
#endif
       fprintf(stderr, "Theta= start method is incorrect, should be one of these:");
       fprintf(stderr, "theta=RANDOM or theta=Own:x.x\n");
       fprintf(stderr, "or theta=Own:{x.x, x.x , x.x, .....}\n");
       fprintf(stderr, "or theta=PRIOR:percentage\n");
       fprintf(stderr, "\nI used theta=PRIOR:1 instead\n\n");
       options->startguess[0][0] = PRIOR;
       options->startguess[0][1] = 1;
       break;

     }
#ifdef DEBUG
   show_startparam(stdout,options,THETAPRIOR,VERBOSE);
#endif
   myfree(otmp);
   myfree(otempstr);
}


///
/// read migration rates from file (master) or buffer (worker in MPI code)
void read_mig (option_fmt * options, char *parmvar, char * varvalue, char **buffer)
{
  (void) parmvar;
  long start;
  long stop;
  long test = 0;
  char *tempstr;
  char *tmp;
  char *otempstr;
  char *otmp;
  //int choint;
  //char choice;
  char *v = varvalue;

  tempstr = (char *) mycalloc(LINESIZE,sizeof(char));
  tmp = (char *) mycalloc(LINESIZE,sizeof(char));
  otempstr = tempstr;
  otmp = tmp;
  
  /* 1st example:  1.0 (n-island model)
     2nd example: {1.0} (migration matrix model, all the same start values  
     3rd example: the dashes on the diagonal are NECESSARY, {} are facultativ
     -  1.0 0.1
     1.0  -  2.0
     0.9 1.2  -
     to specify real 0.0 you need to use the custom-migration settings.
     0.0 in the table will be change to SMALLES_MIGRATION
  */
  // we can safely ignore parmvar because this has to be "migration"
  // varvalue must contain something like this
  // OWN:5.0 or OWN:{ - 5.0 2.0 - } or FST or URANDOM:{ a, b}  or NRANDOM:{a,b}
  // if there is a { we search for the last } and then assume after that a new line contains a new option item
   switch (uppercase (varvalue[0]))
     {
     case 'R':
       options->startguess[1][0] = RANDOMPRIOR;
       break;
     case 'P':
       options->startguess[1][0] = PRIOR;
       get_next_word(&v,":",&tmp);// extract the word OWN: , into tmp	   
       get_next_word(&v,",{}\n",&tmp);
       while(strchr("{", (tmp)[0]))
	 {
	   get_next_word(&v,",{} \n",&tmp);
	 }
       options->startguess[1][1] = atol(tmp);
       break;
     case 'O':
     case '0':
       options->startguess[1][0] = OWN;
       options->migration_model = MATRIX;
       start = count_char(varvalue,'{');
       stop =  count_char(varvalue,'}');
       test = count_char(varvalue,'-');
       if (start == 0 && stop == 0)
	 {
	   get_next_word(&v,":",&tmp);// extract the word OWN: , into tmp
	   add_startparam(options,MIGPRIOR, (float) atof(v));
	 }
       else /* not on one line or using {} notation*/
	 {
	   if (start == 1 && stop == 1)  //all on one line
	     {
	       options->migration_model = ISLAND;
	       get_next_word(&v,"{",&tmp);// extract the word OWN: {, into tmp
	       while (v!= NULL)
		 {
		   get_next_word(&v," ,}",&tmp);
		   if ((tmp != NULL) && (tmp[0] != '-'))
		     add_startparam(options,MIGPRIOR, (float) atof(tmp));
		 }
	     }
	   else // multiple lines
	     {
	       if (start == 1 && stop == 0 )
		 {
		   get_next_word(&v,"{",&tmp);// extract the word OWN: {, into tmp
		   while (v!= NULL)
		     {
		       get_next_word(&v," ,}",&tmp);
		       if ((tmp != NULL) && (tmp[0] != '-'))
			 add_startparam(options,MIGPRIOR, (float) atof(tmp));
		     }
		   // read the next line 
		   while(stop==0)
		     {
		       if (myID==MASTER)
			 FGETS (varvalue, LINESIZE, options->parmfile);
		       else
			 sgets (varvalue, LINESIZE, buffer);
		       while (varvalue[0]=='#')
			 {
			   if(myID==MASTER)
			     FGETS (varvalue, LINESIZE, options->parmfile);
			   else
			     sgets (varvalue, LINESIZE, buffer);
			 }
		       stop =  count_char(varvalue,'}');
		       test += count_char(varvalue,'-');
		       v = varvalue;
		       while (v!= NULL)
			 {
			   get_next_word(&v," ,}",&tmp);
			   if ((tmp != NULL) && (tmp[0] != '-'))			   		       
			     add_startparam(options,MIGPRIOR, (float) atof(tmp));
			 }		   
		     }
		 }
	       else
		 {
		   usererror("migration=OWN startvalue not correctly formatted");
		 }
	     }
	 }
      break;
    default:
       options->startguess[1][0] = PRIOR;
       options->startguess[1][1] = 1;
       fprintf(stderr, "migration= start method is incorrect, should be one of these:");
       fprintf(stderr, "migration=RANDOM or migration=Own:x.x\n");
       fprintf(stderr, "or migration=Own:{ -- x.x x.x  x.x -- .....}\n");
       fprintf(stderr, "or migration=PRIOR:percentage\n");
       fprintf(stderr, "\nI used migration=PRIOR:1 instead\n\n");
    }
  myfree(otmp);
  myfree(otempstr);
}


///
/// read the values for start parameter from the list in the parmfile
/// called through read_options_master() and also equivalent to read worker
void read_startparameter (long parampos, option_fmt * options, char *parmvar, char *varvalue, char **buffer)
{
  (void) parmvar;
   long start = -1;
   long stop = -1;
   short prior = -1;
   //long test = 0;
   char *tempstr;
   char *tmp;
   char *otempstr;
   char *otmp;
   char *v = varvalue;
   char errormessage[LINESIZE];
   char ppp[5][10]={"theta","migration","rate","split","splitstd"};
   tempstr = (char *) mycalloc(LINESIZE,sizeof(char));
   tmp = (char *) mycalloc(LINESIZE,sizeof(char));
   otempstr = tempstr;
   otmp = tmp;
   switch(parampos)
     {
     case 0: prior = THETAPRIOR; 
       strcpy(errormessage, "theta=OWN startvalue not correctly formatted");
       break;
     case 1: prior = MIGPRIOR; 
       strcpy(errormessage, "migration=OWN startvalue not correctly formatted");
       break;
     case 2: prior = RATEPRIOR;        
       strcpy(errormessage, "rate=OWN startvalue not correctly formatted");
       break;
     case 3: prior = SPLITPRIOR; 
       strcpy(errormessage, "split=OWN startvalue not correctly formatted");
       break;
     case 4: prior = SPLITSTDPRIOR; 
       strcpy(errormessage, "splitstd=OWN startvalue not correctly formatted");
       break;
     }
   switch (uppercase (varvalue[0]))
     {
     case 'R':
       options->startguess[parampos][0] = RANDOMPRIOR;
       break;
     case 'P':
       options->startguess[parampos][0] = PRIOR;
       get_next_word(&v,":",&tmp);// extract the word OWN: , into tmp	   
       get_next_word(&v,",{}\n",&tmp);
       while(strchr("{", (tmp)[0]))
	 {
	   get_next_word(&v,",{} \n",&tmp);
	 }
       options->startguess[parampos][1] = atol(tmp);
       break;
     case 'O':
     case '0':
       options->startguess[parampos][0] = OWN;
       if(parampos==1) //migration
	 {
	   options->migration_model = MATRIX;
	   //test = count_char(varvalue,'-');
	 }
       start = count_char(varvalue,'{');
       stop =  count_char(varvalue,'}');
       if (start == 0 && stop == 0)
	 {
	   get_next_word(&v,":",&tmp);// extract the word OWN: , into tmp
	   add_startparam(options,prior, (float) atof(v));
	 }
       else /* not on one line or using {} notation*/
	 {
	   if (start == 1 && stop == 1)  //all on one line
	     {
	       if(parampos==1) //migration
		 {
		   options->migration_model = ISLAND;
		 }
	       get_next_word(&v,"{",&tmp);// extract the word OWN: {, into tmp
	       while (v!= NULL)
		 {
		   get_next_word(&v," ,}",&tmp);
		   if((tmp != NULL)  && (tmp[0] != '-'))//THE '-' is for migration
		     add_startparam(options,prior,(float) atof(tmp));
		 }
	     }
	   else // multiple lines
	     {
	       if (start == 1 && stop == 0 )
		 {
		   get_next_word(&v,"{",&tmp);// extract the word OWN: {, into tmp
		   while (v!= NULL)
		     {
		       get_next_word(&v," ,}",&tmp);
		       if ((tmp != NULL)  && (tmp[0] != '-'))
			 add_startparam(options,prior, (float) atof(tmp));
		     }
		   // read the next line 
		   while(stop==0)
		     {
		       if(myID==MASTER)
			 FGETS (varvalue, LINESIZE, options->parmfile);
		       else
			 sgets (varvalue, LINESIZE, buffer);
		     }
		   while (varvalue[0]=='#')
		     {
		       if(myID==MASTER)
			 FGETS (varvalue, LINESIZE, options->parmfile);
		       else
			 sgets (varvalue, LINESIZE, buffer);
		     }
		   //count_char(varvalue,'}');
		    //if (parampos==2)
		    // test += count_char(varvalue,'-');
		   v = varvalue;
		   while (v!= NULL)
		     {
		       get_next_word(&v," ,}",&tmp);
		       if ((tmp != NULL)  && (tmp[0] != '-'))
			 add_startparam(options,prior, (float) atof(tmp));
		     }		   
		 }
	       else
		 {
	       usererror(errormessage);
		 }
	     }
	 }
       //options->numpop = options->startparam.numpop;
       break;
     default:
#ifdef MPI
       fprintf(stderr,"%i>",myID);
#endif
       fprintf(stderr, "%s= start method is incorrect, should be one of these:",ppp[parampos]);
       fprintf(stderr, "%s=RANDOM or theta=Own:x.x\n",ppp[parampos]);
       fprintf(stderr, "or %s=Own:{x.x, x.x , x.x, .....}\n",ppp[parampos]);
       fprintf(stderr, "or %s=PRIOR:percentage\n",ppp[parampos]);
       fprintf(stderr, "\nI used %s=PRIOR:1 instead\n\n",ppp[parampos]);
       options->startguess[2][0] = PRIOR;
       options->startguess[2][1] = 1;
       break;
     }
#ifdef DEBUG
   show_startparam(stdout,options,prior,VERBOSE);
#endif
   myfree(otmp);
   myfree(otempstr);
}



char
skip_space (option_fmt * options)
{
  char ch = (char) getc (options->parmfile);
    while (isspace ((int) ch) || ch == ',')
    {
      ch = (char) getc (options->parmfile);
    }
    if (isalpha (ch))
    {
        ungetc (ch, options->parmfile);
        ch = '\0';
    }
    return ch;
}

#ifdef MPI
char
skip_sspace (char **buffer)
{
    char ch = sgetc (buffer);
    while (isspace ((int) ch) || ch == ',')
    {
        ch = sgetc (buffer);
        if (isalpha (ch))
            return ch;
    }
    return ch;
}
#endif

void
set_profile_options (option_fmt * options)
{
    switch (options->profile)
    {
    case myNONE:
        options->printprofsummary = options->printprofile = FALSE;
        break;
    case ALL:
        options->printprofsummary = options->printprofile = TRUE;
        break;
    case TABLES:
        options->printprofsummary = FALSE;
        options->printprofile = TRUE;
        break;
    case SUMMARY:
        options->printprofsummary = TRUE;
        options->printprofile = FALSE;
        break;
    }
    if (options->profilemethod == 'd')
        options->printprofsummary = FALSE;
}


void add_custom_mig(option_fmt *options, long *z, char ch)
{
  switch (ch)
    {
    case '{': //opens the migration matrix
    case '}': //closes the migration matrix, 
    case ' ': 
    case '\t':
    case '\0':
    case '\n':
      break;
    default:
      options->custm  = (char *) myrealloc (options->custm, sizeof (char) * (size_t) (*z + 2));
      options->custm2 = (char *) myrealloc (options->custm2, sizeof (char) * (size_t) (*z + 2));
      switch (ch)
	{
	case 'C':
	case 'T':
	case 'D':
	case 'M': 
	case 'S':
	  options->custm[(*z)++] = ch;
	  break;
	case 'x':
	case 'X':
	case '*':
	  options->custm[(*z)++] = '*';
	  break;
	case 'c':
	case 'd':
	case 't':
	  options->custm[(*z)++] = ch;
	  break;
	default:
	  // here we need code to adjust for the different keys to improve the 
	  // migration matrix further
	  options->custm[(*z)++] = (char) tolower (ch);
	  break;
	}
    }
}



/* custom-migration:<{> migration matrix and theta on
   diagonal:
   0 means not estimated,
   x means estimated, s means symmetrically
   estimated, m means all are the same, 
   c means remains constant at start value <}>
   example: 
   {* * s
    * c *
    s 0 *}
New version that will supercede the stiff old one if there is a custom-key:
like this:
custom-key={a=m,b=m} #standard schemes do not need to be named
custom-migration={a 0 0 d * a * d b b * 0 b * * *}
standards are : * x d D t T m M s S c 
*/
void
read_custom_migration (FILE * file, option_fmt * options, char *varvalue,
                       long customnumpop,long pos)
{
  (void) customnumpop;
  long i,j,ii;
  long z=pos;
  long zz=0;
  char ch;
  //char *tmp;
  char *v = varvalue;
  //tmp = (char *) mycalloc(LINESIZE,sizeof(char));
  long start = count_char(varvalue,'{');
  long stop =  count_char(varvalue,'}');
  //long test = count_char(varvalue,'-');
  if (start == 0 && stop == 0)
    {
      ch = ' ';
      while(ch!='\0')
	{
	  ch = v[zz++];
	  if(strchr("{ \t",ch))
	    continue;
	  else
	    options->custm[z++] = ch;
	}
    }
  else /* not on one line or using {} notation*/
    {
      if (start == 1 && stop == 1)  //all on one line
	{
	  ch = ' ';
	  while(ch!='\0')
	    {
	      ch = v[zz++];
	      if(strchr("{ ,\t}",ch))
		continue;
	      else
		add_custom_mig(options, &z, ch);
	    }
	}
      else // multiple lines
	{
	  if (start == 1 && stop == 0 )
	    {
	      ch = ' ';
	      while(ch!='\0')
		{
		  ch = v[zz++];
		  if(strchr("{ ,\t}",ch))
		    continue;
		  else
		    add_custom_mig(options, &z, ch);
		}
	      // read the next line 
	      while(stop==0)
		{
		  zz=0;
		  if (myID==MASTER)
		    FGETS (varvalue, LINESIZE, file);
		  else
		    usererror("the parallele migrate version should not read here");
		  while (varvalue[0]=='#')
		    {
		      FGETS (varvalue, LINESIZE, file);
		    }
		  stop =  count_char(varvalue,'}');
		  v = varvalue;
		  ch = ' ';
		  while(ch!='\0')
		    {
		      ch = v[zz++];
		      if(strchr("{ ,\t}",ch))
			continue;
		      else
			add_custom_mig(options, &z, ch);
		    }
		}
	    }
	  else
	    {
	      usererror("custom-migration option not correctly formatted");
	    }
	}
    }
  zz = z;
  options->custm[zz] = '\0';
  long lc = (long) strlen (options->custm);
  long numpop = (long) (sqrt ((MYREAL) lc));
  z = numpop;
  for (i = 0; i < numpop; i++)
    {
      for (j = 0; j < numpop; j++)
        {
	  ii = i * numpop + j;
	  if (i == j)
	    options->custm2[i] = options->custm[ii];
	  else
	    options->custm2[z++] = options->custm[ii];
        }
    }
  options->custm2[z] = '\0';
  specify_migration_type (options);
  if (file != stdin)
    {
      fgetpos(options->parmfile, &thePos);
    }
  while (file != stdin && !(strstr (varvalue, "end") || strchr (varvalue, '=')))
    {
      fgetpos(options->parmfile, &thePos);
      FGETS (varvalue, LINESIZE, file);
    }
  if (file != stdin)
    {
      reset_oneline (options, &thePos);
    }
}

#ifdef MPI
void
read_custom_migration_worker (char **buffer, option_fmt * options,
                              char *value, long customnumpop)
{
  long cc = customnumpop;
    long zz = 0, z = 0;
    char ch = '\0';
    long lc, numpop, i, j, ii;
    char *position;

    if (cc == 0)
        cc = 1000000;
    else
        cc *= cc;

    z = 0;
    zz = 0;
    while (ch != '}')
    {
        ch = value[z];
        switch (ch)
        {
        case '}':
        case '{':
        case ' ':
        case '\t':
            z++;
            break;
        case '\0':
        case '\n':
            z = 0;
            sgets (value, LINESIZE, buffer);
            break;
        default:
            options->custm =
                (char *) myrealloc (options->custm, sizeof (char) * (zz + 2));
            options->custm2 =
                (char *) myrealloc (options->custm2, sizeof (char) * (zz + 2));
            switch (ch)
            {
	    case 'T':
	    case 'D':
            case 'S':
            case 'M':  //do we have code for this?
                options->custm[zz++] = ch;
                break;
            case 'x':
            case 'X':
                options->custm[zz++] = '*';
                break;
	    case '0':
	      options->custm[zz++] = '0';
	      break;
            default:
                options->custm[zz++] = tolower (ch);
                break;
            }
            z++;
        }
    }
    options->custm[zz] = '\0';
    lc = (long) strlen (options->custm);
    numpop = (long) sqrt ((MYREAL) lc);
    z = numpop;
    options->numpop = numpop;
    for (i = 0; i < numpop; i++)
    {
        for (j = 0; j < numpop; j++)
        {
            ii = i * numpop + j;
            if (i == j)
                options->custm2[i] = options->custm[ii];
            else
                options->custm2[z++] = options->custm[ii];
        }
    }
    options->custm2[z] = '\0';
    specify_migration_type (options);
    position = *buffer;
    while (!(strstr (value, "end") || strchr (value, '=')))
    {
        position = *buffer;
        sgets (value, LINESIZE, buffer);
    }
    *buffer = position;
}
#endif /*MPI*/
void
specify_migration_type (option_fmt * options)
{
    long len = (long) strlen (options->custm);
    long ms = 0, xs = 0, ns = 0, ss = 0, len2, i;
    char *p;
    p = options->custm;
    while (*p != '\0')
    {
        switch (*p)
        {
        case 'm':
            ms++;
            break;
	case 'd':
	case 't':
	  options->has_speciation=TRUE;
	  ns++;
	  break;
	case 'D':
	case 'T':
	  options->has_speciation=TRUE;
	  xs++;
	  break;
        case 'x':
        case '*':
            xs++;
            break;
        case '0':
            ns++;
            break;
        case 'S':
        case 's':
            ss++;
            break;
        case 'c':
            break;
        }
        p++;
    }
    if (ms>0 || xs > 0 || ss > 0)
      options->has_migration = TRUE;
    if (ms >= len)
    {
        options->migration_model = ISLAND;
        return;
    }
    if (xs >= len)
    {
        options->migration_model = MATRIX;
        return;
    }
    if (ns >= len && options->has_speciation == FALSE)
    {
        usererror ("Custom migration matrix was completely set to zero?!\n");
        //return;
    }
    len2 = (long) (sqrt ((MYREAL) len));
    if (ms == len2 && xs == len - len2)
    {
        for (i = 0; i < len2; i++)
        {
            if (options->custm[i * len2 + i] != 'm')
            {
                options->migration_model = MATRIX;
                return;
            }
        }
        options->migration_model = MATRIX_SAMETHETA;
        return;
    }
    if (xs == len2 && ms == len - len2)
    {
        for (i = 0; i < len2; i++)
        {
            if (options->custm[i * len2 + i] != '*')
            {
                options->migration_model = MATRIX;
                return;
            }
        }
        options->migration_model = ISLAND_VARTHETA;
        return;
    }
    options->migration_model = MATRIX_ARBITRARY;
}



void
fillup_custm (long len, world_fmt * world, option_fmt * options)
{
    long i, j, ii, z;
    char *tmp;
    len = (long) strlen (options->custm);
    tmp = (char *) mycalloc (1, sizeof (char) * (size_t) (world->numpop2 + 2));
    options->custm =
      (char *) myrealloc (options->custm, sizeof (char) * (size_t) (world->numpop2 + 2));

    options->custm2 =
      (char *) myrealloc (options->custm2, sizeof (char) * (size_t) (world->numpop2 + 2));
    strncpy (tmp, options->custm, (size_t) (world->numpop2 + options->gamma));
    z = world->numpop;
    for (i = 0; i < world->numpop; i++)
    {
        for (j = 0; j < world->numpop; j++)
        {
            ii = i * world->numpop + j;
            if (ii < len)
                options->custm[ii] = tmp[ii];
            else
                options->custm[ii] = '*';
            if (i == j)
                options->custm2[i] = options->custm[ii];
            else
                options->custm2[z++] = options->custm[ii];
        }
    }
    if (options->gamma)
    {
        if (len <= world->numpop2)
        {
            options->custm[world->numpop2] = '*';
            options->custm2[world->numpop2] = '*';
        }
        else
        {
            options->custm[world->numpop2] = tmp[world->numpop2];
            options->custm2[world->numpop2] = tmp[world->numpop2];
        }
        options->custm[world->numpop2 + 1] = '\0';
        options->custm2[world->numpop2 + 1] = '\0';
    }
    else
    {
        options->custm[world->numpop2] = '\0';
        options->custm2[world->numpop2] = '\0';
    }
    strcpy (world->options->custm, options->custm);
    strcpy (world->options->custm2, options->custm2);
    myfree(tmp);
}

void
print_arbitrary_migration_table (FILE * file, world_fmt * world, option_fmt *options,
                                 data_fmt * data)
{
  long i,j,z;
  //char mytext[LONGLINESIZE];
  
  fprintf (file, "Connection matrix:\n");				  
  fprintf (file, "m = average (average over a group of Thetas or M,\n");	  
  fprintf (file, "s = symmetric migration M, S = symmetric 4Nm,\n");	  
  fprintf (file, "0 = zero, and not estimated,\n");			  
  fprintf (file, "* = migration free to vary, Thetas are on diagonal\n");	  
  fprintf (file, "d = row population split off column population\n");	  
  fprintf (file, "D = split and then migration\n");                        
  for (i = 0; i < data->numpop; i++)
    {
      fprintf (file, "%4li %-10.10s     ", options->newpops[i], data->popnames[i]);
      for (j = 0; j < data->numpop; j++)
	{
	  z = world->numpop * (options->newpops[i]-1) + (options->newpops[j]-1);
	  fprintf (file, "%c ", world->options->custm[z]);
	}
      fprintf (file, "\n");
    }
  fprintf (file, "\n\n");
  // print growth table if there is growth:
  if (world->has_growth) // && !(world->options->growpops[0] == 0 && world->options->growpops_numalloc==1))
    {
      for (i = 0; i < data->numpop; i++)
	{
	  fprintf (file, "%4li %-10.10s     ", options->newpops[i], data->popnames[i]);
	  if (world->options->growpops[options->newpops[i]-1]>0)
	    fprintf (file, "%s [%li]\n", "growing/shrinking", world->options->growpops[options->newpops[i]-1]);
	  else
	    fprintf (file, "%s\n", "constant size    ");
	}
      fprintf (file, "\n\n");
    }
}

void
print_distance_table (FILE * file, world_fmt * world, option_fmt * options,
                      data_fmt * data)
{
    long i, j;

    if (options->geo)
    {
        fprintf (file,
                 "   Geographic distance matrix between locations\n      ");
        for (i = 0; i < world->numpop; i++)
        {
            fprintf (file, "%-10.10s     ", data->popnames[i]);
            for (j = 0; j < world->numpop; j++)
            {
                if (i == j)
                    fprintf (file, "   -   ");
                else
                    fprintf (file, "%6.3f ", data->ogeo[i][j]);
            }
            fprintf (file, "\n      ");
        }
        fprintf (file, "\n");
    }
}

void reorder_populations(world_fmt *world, option_fmt *options, data_fmt *data)
{
  //  long newnumpop=0;
  long pop;
  long lim = data->numpop;
  long i;
  for(pop=0;pop<lim;pop++)
    {
      for(i=0; i<world->sumtips;i++)
	{
	  if(world->nodep[i]->actualpop == pop)
	    {
	      world->nodep[i]->actualpop = world->nodep[i]->pop = options->newpops[pop]-1; 
	    }
	}
    }
}

void set_localities(char **value, char **tmp, option_fmt *options)
{
      long z=0;
      options->newpops[0] = 1;
      get_next_word(value,"{}:,; ",tmp);
      while(*tmp != NULL)
	{
	  if(z >= options->newpops_numalloc)
	    {
	      options->newpops_numalloc = z+1;
	      options->newpops = (long*) myrealloc(options->newpops, options->newpops_numalloc * sizeof(long));
	    }
	  options->newpops[z] = atol(*tmp);
	  //printf("%i> population relabel [%li]=%li (input: |%s|)\n",myID,z, options->newpops[z],*tmp);
	  z++;
	  get_next_word(value,"{}:,; ",tmp);
	}
}

void set_growth(char **value, char **tmp, option_fmt *options)
{
      long z=0;
      options->growpops[0] = 0;
      get_next_word(value,"{}:,; ",tmp);
      while(*tmp != NULL)
	{
	  if(z >= options->growpops_numalloc)
	    {
	      options->growpops_numalloc = z+1;
	      options->growpops = (long*) myrealloc(options->growpops, options->growpops_numalloc * sizeof(long));
	    }
	  options->growpops[z] = atol(*tmp);
	  //printf("%i> population-growth relabel [%li]=%li (input: |%s|)\n",myID,z, options->growpops[z],*tmp);
	  z++;
	  get_next_word(value,"{}:,; ",tmp);
	}
}

///
/// delivers sequence submodel from string
/// allowed strings are JC69, K80, TN, HKY, F84, F81
int sequence_submodeltype(char *mutmodel)
{
  int type = uppercase(mutmodel[0]);
  switch(type)
    {
    case 'J':
      return JC69;
    case 'K':
      return K2P;
    case 'F':
      if(mutmodel[2]=='1')
	return F81;
      else 
	return F84;
    case 'H':
      return HKY;
    case 'T':
      return TN;
    default:
      error ("Parmfile contains incompatible mutation model for sequence must be one of  JC69, K80, HKY, F81, F84, TN");
    }
  return TN;
}

void set_updating_choices(double *choices, option_fmt * options, int flag)
{
  //fprintf(stdout,"%i> @@@@ options->has_unassigned=%li [%f]\n",myID, (long) options->has_unassigned, options->unassigned_updatefreq);
  // order in choice vector
  // {treeupdatefreq,bayesupdatefreq,haplotypeupdatefreq, timeparamupdatefreq, assignupdatefreq, seqerrorupdatefreq ml-alpha}
  if(options->tree_updatefreq>0.0)
    choices[TREEUPDATE]=options->tree_updatefreq;
  else
    choices[TREEUPDATE]=0.0;
  if(options->parameter_updatefreq>0.0)
    choices[PARAMETERUPDATE]=options->parameter_updatefreq;
  else
    choices[PARAMETERUPDATE]=0.0;
  if(options->haplotyping && options->haplotype_updatefreq>0.0)
    choices[HAPLOTYPEUPDATE]=options->haplotype_updatefreq;
  else
    choices[HAPLOTYPEUPDATE] = 0.0;
  if(options->skyline_param)
    choices[SKYLINETIMEUPDATE] = options->timeparam_updatefreq;
  else
    choices[SKYLINETIMEUPDATE]=0;
  if(options->has_unassigned)
    choices[ASSIGNMENTUPDATE] = options->unassigned_updatefreq;
  else
    choices[ASSIGNMENTUPDATE]=0;
  if(options->has_estimateseqerror)
    choices[SEQUENCEERRORUPDATE] = options->seqerror_updatefreq;
  else
    choices[SEQUENCEERRORUPDATE]=0;

  if(options->tri_mlalpha==ESTIMATE)
    choices[MITTAGLEFFLERUPDATE] = options->mlalpha_updatefreq;
  else
    choices[MITTAGLEFFLERUPDATE]=0;

  switch(flag)
    {
    case NOPARAMETER:
      choices[1] = 0.0;
      break;
    case NOTREE:
      choices[0] = 0.0;
      break;
    case STANDARD:
    default:
      break;
    }

  double denom = choices[TREEUPDATE] + choices[PARAMETERUPDATE] + choices[HAPLOTYPEUPDATE] +
    choices[SKYLINETIMEUPDATE] + choices[ASSIGNMENTUPDATE] + choices[SEQUENCEERRORUPDATE] + choices[MITTAGLEFFLERUPDATE];

  if(denom>0.0)
    {
      choices[0] = choices[0] / denom;
      choices[1] = choices[0] + choices[1] / denom;
      choices[2] = choices[1] + choices[2] / denom;
      choices[3] = choices[2] + choices[3] / denom;
      choices[4] = choices[3] + choices[4] / denom;
      choices[5] = choices[4] + choices[5] / denom;
      choices[6] = choices[5] + choices[6] / denom;
    }
  else
    {
      usererror("The update probability for trees, parameters, haplotypes, ... is set incorrectly");
    }
}
