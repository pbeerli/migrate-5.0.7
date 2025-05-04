/*! \file main.c */
/*! \mainpage M I G R A T E
*  \section intro Introduction
*  Migrate is a program to estimate population-genetic parameters from genetic data such as
*  electrophoretic marker data, microsatellite data, DNA or RNA sequence data, and single 
*  nucleotide polymorphisms. Migrate uses the concepts of Maximum likelihood and Bayesian
*  inference to infer parameters based on coalescence theory. For more information related
*  to the use of the program or the theory behind check out 
*  http://popgen.scs.fsu.edu/migrate.html
*  
*  \section install Installation
*  to install the program simply do
*  \subsection step1 configure
*  \subsection step2 make
*  \subsection step3 sudo make install
*
*  \section maintainer Maintainer and Copyrights 
*   Peter Beerli,
*   beerli@fsu.edu
*
*   Copyright 1997-2002 Peter Beerli and Joseph Felsenstein, Seattle WA\n
*   Copyright 2003-2016 Peter Beerli, Tallahassee FL\n
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
*  of the Software, and to permit persons to whom the Software is furnished to do
*  so, subject to the following conditions:
* 
*  The above copyright notice and this permission notice shall be included in all copies
*  or substantial portions of the Software.
* 
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
*  INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
*  PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
*  HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
*  CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
*  OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
* 
*/
/* $Id: main.c 2169 2013-08-24 19:02:04Z beerli $*/
#ifndef WINDOWS
#include <sys/resource.h>
#endif
#include "migration.h"
//#include "aic.h"
#include "assignment.h"
#include "autotune.h"
#include "bayes.h"
#include "data.h"
#include "growth.h"
#include "haplotype.h"
#include "heating.h"
#include "tables.h"
#include "histogram.h"
//#include "lrt.h"

#include "marginallike.h"
#include "mcmc.h"
#include "menu.h"
#include "migevents.h"
#include "migrate_mpi.h"
#include "options.h"
#include "../lib/SFMT-src-1.5.1/SFMT.h"
#include "random.h"
#include "reporter.h"
#include "sequence.h"
#include "seqerror.h"
#include "sighandler.h"
#include "skyline.h"
#include "tools.h"
#include "tree.h"
#include "world.h"
#include "priors.h"
#include "pretty.h"
#include "speciate.h"
#include <stdio.h>
#ifdef BEAGLE
#include "calculator.h"
#endif
#include "mutationmodel.h"
#ifdef UEP
#include "uep.h"
#endif
#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif

#ifdef GRANDCENTRAL
#include <dispatch/dispatch.h>
extern dispatch_semaphore_t semaphore;
dispatch_semaphore_t semaphore;
#endif

//use this if you have nice() and need to be nice with others
#ifdef HAVE_NICE
#include <unistd.h>
#endif

#ifdef WIN32
#include "pretty-win32.h"
#endif


/* Definitions for MCMCMC -- heated chains
 * definitions of the first four heated chains, they are ordered from cold to hot */
#define EARTH  universe[0]	/*!< cold chain has always a temperature=1            */
#define VENUS  universe[1]	/*!< second coldest chain                             */
#define MERKUR universe[2]	/*!< thrid coldest chain                              */
#define SUN    universe[3]	/*!< fourth coldest chain                             */

/* GLOBAL VARIABLES +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
extern time_t startseconds; 
extern int simulator;
extern int myID;	  /*!< myID=0 for single-cpu program and master in MPI program, worker nodes myID > 0 */
extern int myRepID;	  /*!< myID=0 for single-cpu program and master in MPI program, worker nodes myID > 0 */
extern int color;	  /*!< either 1 or MPI_UNDEFINED, in anticipation of more complicated MPI schemes     */
extern int numcpu;	  /*!< used if MPI is running otherwise == 1                                          */
extern int locidone;	  /*!< used in MPI workers                                                            */
extern long unique_id_global;
extern int use_beagle_gpu;
extern int use_beagle_dynamicscale;
extern int use_beagle_autoscale;
extern int use_beagle_manualscale;

time_t startseconds;
int simulator;
int myID;	  /*!< myID=0 for single-cpu program and master in MPI program, worker nodes myID > 0 */
int myRepID;	  /*!< myID=0 for single-cpu program and master in MPI program, worker nodes myID > 0 */
int color;	  /*!< either 1 or MPI_UNDEFINED, in anticipation of more complicated MPI schemes     */
int numcpu;	  /*!< used if MPI is running otherwise == 1                                          */
int locidone;	  /*!< used in MPI workers                                                            */
long unique_id_global;
int use_beagle_gpu;
int use_beagle_dynamicscale;
int use_beagle_autoscale;
int use_beagle_manualscale;
#ifdef MPI
#ifdef PARALIO
boolean my_file_open_error;
#endif
#ifdef  USE_MYREAL_FLOAT
#ifdef WINDOWS
#define mpisizeof MPI_FLOAT;
#else
const MPI_Datatype mpisizeof = MPI_FLOAT;
#endif
#else
#ifdef WINDOWS
//#define mpisizeof MPI_DOUBLE;
 MPI_Datatype mpisizeof;
#else
const  MPI_Datatype mpisizeof = MPI_DOUBLE;
#endif
#endif
#endif

#ifdef CAUTIOUS
boolean cautious; /*!< forces a check whether a file exists and will ask before overwriting           */
#endif
#ifdef MEMDEBUG
#include <sys/time.h>
FILE *memfile;
struct timeval memt_start, memt_finish;
double memelapsed;
long totalsize;
#endif
#ifdef MPI
filedb_fmt filedb[30];
long filenum;

MPI_Comm comm_world;		/*!< parallel environment that contains knowledge about all workers and master */
//MPI_Comm comm_workers;		/*!< parallel environment that contains knowledge about all workers */
MPI_Group worker_group;
MPI_Group world_group;
#endif
#ifdef SLOWNET
int profiledone;		/*!< used in MPI workers in profiles on slow networks    */
#endif
#ifdef PTHREADS
tpool_t heating_pool;		/*!< when compiled with PTHREADS then holds all threads */
#else
#ifndef tpool_t
#define tpool_t char
#endif
extern tpool_t heating_pool;
tpool_t heating_pool;
#endif

// random generator related global variables
long *seed;	 /*!< contains the seed of the random number */
long *newseed;	 /*!< contains the new random seed */
char *generator; /*!< string that shows what random number generator is used */
extern sfmt_t sfmt;
extern sfmt_t *sfmtp;
extern sfmt_t ** sfmtH;
sfmt_t sfmt;
sfmt_t *sfmtp;
sfmt_t ** sfmtH;

// for bayesian output counter, initialized in bayes_init()
extern long *mdimfilecount;


int
setup_locus (long locus, world_fmt * world, option_fmt * options,
             data_fmt * data);

void
condense_time (world_fmt * world, long *step, long *j,
               MYREAL * accepted, long *G, long *steps, long oldsteps,
               long rep);


long set_repkind (option_fmt * options);

void
heating_init (world_fmt ** universe, int usize, data_fmt * data,
              option_fmt * options);

void
heating_prepare (world_fmt ** universe, int usize,
                 option_fmt * options, data_fmt * data, long rep);

void heating_prepare2 (world_fmt ** universe, int usize);

long
replicate_number (option_fmt * options, long chain, char type,
                  long replicate);

void set_penalizer (long chain, long chains, char type,
                    world_fmt ** universe);

void
run_sampler (option_fmt * options, data_fmt * data,
             world_fmt ** universe, int usize, long *outfilepos, long *Gmax);

void run_replicate (long locus,
                    long replicate,
                    world_fmt ** universe,
                    option_fmt * options,
                    data_fmt * data,
                    tpool_t * localheating_pool,
                    int usize, long *treefilepos, long *Gmax);


void
run_locus (world_fmt ** universe, int usize, option_fmt * options,
           data_fmt * data, tpool_t * localheating_pool, long maxreplicate,
           long locus, long *treefilepos, long *Gmax);

void
run_loci (world_fmt ** universe, int usize, option_fmt * options,
          data_fmt * data, tpool_t * localheating_pool, long maxreplicate,
          long *treefilepos, long *Gmax);

void run_one_update (world_fmt * world);

void run_updates (world_fmt ** universe,
                  int usize, option_fmt * options,
                  tpool_t * localheating_pool, long inc, long increment,
                  long step, long steps);

void heated_swap (world_fmt ** universe, worldoption_fmt * options);

void
change_chaintype (long locus, char *type, world_fmt * world,
                  long *increment, long *oldsteps, long *chains,
                  option_fmt * options);

void
prepare_next_chain (world_fmt ** universe, worldoption_fmt * options,
                    char type, long chain, long *chains,
                    long *pluschain, long locus, long replicate);



#ifdef GRANDCENTRAL
extern dispatch_queue_t queue;
dispatch_queue_t queue;
#endif

void
print_heating_progress (world_fmt ** universe,
                        worldoption_fmt * options, long stepinc);

boolean analyze_oldbayesdata(world_fmt **universe, option_fmt *options, data_fmt *data, long *outfilepos);
void print_theta0(FILE *file, world_fmt *world, long maxreplicate);
//void profile_tables (option_fmt * options, world_fmt * world, long *gmaxptr);

void finish_mac (option_fmt * options, data_fmt * data);

int
setup_locus (long locus, world_fmt * world, option_fmt * options,
             data_fmt * data);

long set_repkind (option_fmt * options);

void heating_prepare2 (world_fmt ** universe, int usize);

long replicate_number (option_fmt * options, long chain, char type,
                       long replicate);

void print_heating_progress2 (FILE * file, worldoption_fmt * options,
                              world_fmt ** universe);


void get_bayeshist (world_fmt * world, option_fmt * options);
void get_treedata (world_fmt * world, option_fmt * options);
void get_mighistdata (world_fmt * world, option_fmt * options);

void change_longsum_times (world_fmt * world);

boolean  check_parmfile(long argcount, char **arguments, char *parmfilename);
boolean set_usemenu(boolean usemenu, boolean fromparmfile);
void check_bayes_options(option_fmt *options);
void reset_bayesmdimfile(world_fmt *world, option_fmt *options);
long  calculate_newpop_numpop(option_fmt *options, data_fmt *data);

void alloc_sticksize(option_fmt *options, data_fmt *data, size_t size);
void set_skiploci(option_fmt *options, data_fmt *data, world_fmt **universe);
void recalc_skyline_values(world_fmt *world, option_fmt * options, long maxreplicate);
void reset_heated_accept(world_fmt **universe, long unum);
void run_increments (world_fmt ** universe,
		     long usize,
		     option_fmt * options,
		     tpool_t * localheating_pool,
		     long increment, long step, long steps);
void run_steps (world_fmt ** universe,
		long usize,
		option_fmt * options,
		tpool_t * localheating_pool, long increment, long steps, long rep);
void run_chains (world_fmt ** universe,
		 long usize,
		 option_fmt * options,
		 tpool_t * localheating_pool,
		 long replicate,
		 long chains,
		 char type,
		 long increment, long oldsteps, long *treefilepos, long *Gmax);

void print_bf_values(world_fmt * world);
void  print_marginal_order(char *buf, long *bufsize, world_fmt *world);
//##
///
/// specifies the classes of override the menu as a parameter to the program
enum override_enum 
  {
    OVERRIDE_NO, OVERRIDE_MENU, OVERRIDE_NOMENU
  };

#ifdef COREDUMP
#include <sys/resource.h>

static boolean EnableCoreDumps(void)
{
  struct rlimit   limit;

  limit.rlim_cur = RLIM_INFINITY;
  limit.rlim_max = RLIM_INFINITY;
  return setrlimit(RLIMIT_CORE, &limit) == 0;
}
#endif




///
/// the program migrate calculates migration rates and population sizes
/// from genetic data, it allows for a wide variety of data types and many options
//  \callgraph
int
main (int argc, char **argv)
{
  //char      type       = 's';
  int i;
  int       usize      = 1;
  int       usemenu    = OVERRIDE_NO;
  long      Gmax       = 0;
  long      outfilepos = 0;
  data_fmt  *data;
  boolean   restarted_bayes_bool=FALSE;
  option_fmt *options;
  world_fmt **universe;
#ifdef MPI
  int       server;
#ifdef WINDOWS
 mpisizeof = MPI_DOUBLE;
#endif
#endif

 //#ifdef DEBUG
 //proptimes = fopen("proptimes",'w');
 //#endif
 
#ifdef COREDUMP
 EnableCoreDumps();
#endif

  use_beagle_gpu = 0;
  use_beagle_manualscale = 0;
  use_beagle_autoscale = 0;
  use_beagle_dynamicscale = 0;
#ifdef HAVE_NICE
  nice (10); //nice value arbitrarily set to 10, 5 is still nice, 0 is standard
#endif
  
  //---------------------------------------------------------------------------------------------
  // windows specific code for pretty printing
#ifndef MPI
#ifdef WINDOWS
  set_haru_handler();
  printf("no exceptionhandler for the PDF library set\n");
#endif
#endif

    //---------------------------------------------------------------------------------------------
    // set stacksize
#ifndef WINDOWS
    const rlim_t kStackSize = MIN_STACKSIZE_63MB; //64MB fails
    struct rlimit rl;
    int result;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
      {
        if (rl.rlim_cur < kStackSize)
          {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0)
              {
                fprintf(stderr, "Stack not set: setrlimit returned result=%d\n", result);
                result = getrlimit(RLIMIT_STACK, &rl);
                fprintf(stderr, "The stack is still set to result=%d [%d]\n", result, RLIMIT_STACK);
              }
          }
      }
#endif

    //---------------------------------------------------------------------------------------------
    // MPI initialisation
#ifdef MPI
    // parallel version
    filenum = 0; //setting up the filename database so that worker node can write back to master

    MPI_Init (&argc, &argv);
    comm_world = MPI_COMM_WORLD;
    MPI_Comm_size (comm_world, &numcpu);
    MPI_Comm_rank (comm_world, &myID);
    MPI_Comm_group (comm_world, &world_group);
    server = MASTER;		//server ID
    MPI_Group_excl (world_group, 1, &server, &worker_group);
    locidone = 0;
#ifdef FIRSTDEBUG
    if (myID==MASTER)
      {
	printf("Begin of main(), waiting for user to start the master, the worker are not stopped here!\n");
	getc(stdin);
      }
#endif
#ifdef SLOWNET
    // slow network parallel version -- this will turn into the standard version
    //profiledone = 0;
    //which_calc_like (SINGLELOCUS);
#endif

#else /*MPI*/
    //scalar version of migrate
    myID = MASTER;
    numcpu = 1;
#endif /*MPI*/
    simulator = FALSE;
    // debug code for memory problems
#ifdef MEMDEBUG
    totalsize = 0;
    memfile = fopen("memoryfile","w+");
    gettimeofday(&memt_start, NULL);
#endif

    //---------------------------------------------------------------------------------------------
    // initialization of main container and random number parts
    seed = (long *) mymalloc (sizeof (long) * 3);
    newseed = (long *) mymalloc (sizeof (long) * 3);
    generator = (char *) mycalloc (1,sizeof(char) * 80);
    universe = (world_fmt **) mycalloc (HEATED_CHAIN_NUM, sizeof (world_fmt *));
#ifdef GRANDCENTRAL
        queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
#endif
    // set flag to decide whether we use caution when writing files
#ifdef CAUTIOUS
    cautious = FALSE;
#endif
    // try to catch and beautify some error messages
    signalhandling (ON);
    unique_id_global=0;
    // create main data structures
    create_data (&data);
    create_world (&(EARTH), 1L);
    create_options (&options);
    
    // parmfile and menu
    init_options (options);
#ifdef BEAGLE
    //    print_beagle_available_resources(EARTH);
#endif
    //---------------------------------------------------------------------------------------------
    // master initializations, this works for both parallel and single-cpu version
    if (myID == MASTER)
      {
	usemenu = check_parmfile((long) argc,argv,options->parmfilename);
#ifdef DEBUG
	printf("%i> argc=%i\n",myID, argc);
	for(i=0; i < argc; i++)
	  printf("%i> argv[%i]=%s\n",myID, i, argv[i]);
	printf("%i> parmfilename=%s\n",myID, options->parmfilename);
#endif
	printf("before option master\n");
	read_options_master (options);
	printf("after option master\n");
	options->menu = set_usemenu(usemenu, options->menu);
	print_menu_title (stdout, options, EARTH);
	get_menu (options, EARTH, data);
	//check_bayes_options(options);

#ifndef NEWVERSION
	options->haplotyping=FALSE;
	options->haplotyping_report=FALSE;
#endif
	
#ifdef MPI        
	//	MYMPIBARRIER (comm_world);
	//check_bayes_priors(options, data, EARTH);
        broadcast_options_master (options, data);
#endif
      }
    else
      {
#ifdef MPI
	// MYMPIBARRIER (comm_world);
	broadcast_options_worker (options);
#endif
	options->menu = FALSE;
      }
    
    //---------------------------------------------------------------------------------------------
    // data initialization and data reading phase
    usize = (options->heating ? ( MAX (1, (int) options->heated_chains)) : 1 );
    universe = (world_fmt **) myrealloc (universe, usize * sizeof (world_fmt *));
    // opens files on all nodes    
#ifdef DEBUG
    char nowstr[LINESIZE];
    get_time (nowstr, "%c");
    printf("%i> before opening files for master and 'pointers' for workers  at %s\n",myID,nowstr);
#endif
    init_files (EARTH, data, options);
#ifdef DEBUG
    get_time (nowstr, "%c");
    printf("%i> files are now open and ready at %s\n", myID, nowstr);
#endif

    if (options->writelog && myID == MASTER)
      {
        print_menu_title (options->logfile, options, EARTH);
      }
    
    EARTH->repkind = SINGLECHAIN;
    //fprintf(stderr,"myID=%i myRepID=%i\n",myID, myRepID);    
    //---------------------------------------------------------------------------------------------
    // sampling phase
    if (!options->readsum && !options->checkpointing) // all go here except when reading old runs
      {
        run_sampler (options, data, universe, usize, &outfilepos, &Gmax);
      }
    else
      {
	if(options->checkpointing)
	  {
	    // what was the last set of loci worked on (MPI?)
	    long **unfinished = NULL;
#ifdef MPI 
	    long repmax = number_replicates2(options);
#endif
	    boolean recover_needed = FALSE;
	    if (myID==MASTER)
	      {
		recover_needed = checking_bayesallfile(universe[0], options, data, &unfinished);
		options->unfinished = unfinished;
#ifdef MPI
		long mysize = universe[0]->loci * repmax;
		long buffer[3];
		buffer[0]= (long) recover_needed;
		buffer[1]= (long) EARTH->loci;
		buffer[2]= (long) repmax;
		MYMPIBCAST (&buffer[0], 3, MPI_LONG, MASTER, comm_world);
		MYMPIBCAST (&(unfinished[0][0]), mysize, MPI_LONG, MASTER, comm_world);
#endif
	      }
#ifdef MPI
	    else
	      {
		long buffer[3];
		MYMPIBCAST (&buffer[0], 3, MPI_LONG, MASTER, comm_world);
		recover_needed = (boolean) buffer[0];
		EARTH->loci = buffer[1];
		repmax = buffer[2];
		intvec2d(&unfinished, EARTH->loci,repmax);
		long mysize = universe[0]->loci * repmax;
		MYMPIBCAST (&(unfinished[0][0]), mysize, MPI_LONG, MASTER, comm_world);
		options->unfinished = unfinished;
#ifdef DEBUG
		  printf("@Checking bayesallfile, to find last loci and replicates worked on\n");
		  long maxsample1 = options->lsteps - 1;
		  long locus;
		  for (locus=0; locus < EARTH->loci; locus++)
		    {
		      long replicate;
		      for (replicate=0; replicate < repmax; replicate++)
			{
			  long xx = unfinished[locus][replicate];
			  printf("%5li ", xx);
			} 
		      printf("\n");
		    }
		  printf("[%li,%li]\n",maxsample1,mysize);		  
#endif
	      }
#endif /*end MPI*/
	    // set the last parameters like the options so that it can be restarted
	    // this will be tricky for MPI, the master needs to seed the loci correctly
	    //if (TRUE) //recover_needed)
	    //  {
		run_sampler (options, data, universe, usize, &outfilepos, &Gmax);
		//  }
		//else
		//{
	    	//options->readsum = TRUE;
	    	//options->checkpointing = FALSE;
	    	//restarted_bayes_bool = analyze_oldbayesdata(universe, options, data, &outfilepos);
		//}
	  }
	  else
	    {
	      restarted_bayes_bool = analyze_oldbayesdata(universe, options, data, &outfilepos);
	    }
      }
    
    //---------------------------------------------------------------------------------------------
    // combining phase for multiple loci and replicates
#ifdef MPI
    //with MPI:the workers stay here much longer than the Master
    //combine_loci_mpi (type, options, EARTH, &Gmax);
    if(myID != MASTER)
      {
	mpi_maximize_worker (EARTH, options, MULTILOCUS, options->replicatenum);
      }
#else
    //combine_loci_standard (type, options, EARTH, &Gmax);
#endif
    // if bayes intermediate data recording is ON then reset the file
    // for reading for printing and combining, else use the material still in RAM
    if(!restarted_bayes_bool)
      reset_bayesmdimfile(EARTH, options);

    //---------------------------------------------------------------------------------------------
    // printing main results
    if (myID == MASTER)
      {
	bayes_stat (EARTH,data);	   
	get_seqerror(EARTH,options);
#ifdef MPI
	get_assignments (EARTH, options);
#endif
	get_haplotypes(EARTH,options);
	print_haplotype_stat(EARTH,data);
	fflush (EARTH->outfile);
	//---------------------------------------------------------------------------------------------
	// printing additional results
#ifdef MPI
	get_mighistdata (EARTH, options);
#endif
	if(options->skyline)
	  {
	    print_expected_values(EARTH, options);       
	  }
	if(options->mighist || options->skyline)
	  {
	    print_event_values(EARTH);
	    pdf_histogram_legend();
	  }              
#ifdef UEP
	// print UEP probabilities
	if (options->uep)
	  analyze_uep (EARTH);
#endif
      } // end of printing main tables and histograms
    
    #ifdef SLOWNET
    //which_calc_like (PROFILE);
    if (myID == MASTER)
     {
	// release the workers from the mpi_maximize_worker() function.
	// the workers will stop working on locus-likelihoods and advance to profiles
    	mpi_send_stop (EARTH);
     }
    #endif /* SLOWNET */
    
#ifdef MPI
    if (myID == MASTER)
      {
#endif
	if(EARTH->has_estimateseqerror)
	  seqerror_report(EARTH,"seqerror");
	
	print_bayesfactor(universe,options);
	print_mcmc_run_character(EARTH);

	print_finaltree(EARTH, options);

	print_heatingreport(universe, options);
	// print warnings into the PDF and the outfile
	//
	print_stored_warnings(EARTH);
	pdf_print_stored_warnings(EARTH);
	report_unassigned(stdout, EARTH);
	report_unassigned(EARTH->outfile, EARTH);
        pdf_report_unassigned(EARTH);
	print_finish (EARTH, outfilepos);
	pdf_print_end_time(&page_height);
	// write out to PDF file
	pdf_write_file(options);
	// closing all files 
#ifdef MAC
	finish_mac (options, data);
#endif
	
	exit_files (EARTH, data, options);
#ifdef MPI
#    ifndef SLOWNET
	mpi_send_stop (EARTH); //stop workers
	// printf("%i> at barrier in slownet-main before finalize\n", myID);        
#    endif
      }
    MYMPIBARRIER (comm_world);
    MPI_Finalize ();
#endif /*MPI*/
    myfree(seed);
    myfree(newseed);
    myfree(generator);
#ifndef MPIDATAONDEMAND
    destroy_data(data);
#endif
    free_universe(universe, usize, options);
    myfree(data);
    destroy_options(options);
    myfree(options);
    // for debugging with MallocDebug on macosx
#ifdef MACMALLOCDEBUG
    while(TRUE)
      {
	sleep(1);
      } //end for debugging
#endif
    return 0;
}				/* main end */



///
/// calculate the number populations to analyze and resize the newpop array if necessary
long  calculate_newpop_numpop(option_fmt *options, data_fmt *data)
{
  long pop;
  long i;
  long newnumpop=0;
  long numpop = data->numpop;
  //  char *strsep(char **, const char *);
#ifdef WIN32
  // (old?) windows compilers have difficult to assign memory 
  // based on variables
  long temp[5000];//this is __way__ more than migrate ever will be 
                  // able to handle    
#else
  long temp[5000];
  // clang with -Weverything dislikes the variable length initialization:
  // long temp[data->numpop < options->newpops_numalloc ? options->newpops_numalloc : data->numpop];
#endif
  //use order of populations {1,2,3,4,5,....} , do not start with zero!!
  //reset to options->newpops={1,1,2,1,2,3,4,...}
  if (options->newpops_numalloc != numpop)
    {
      // newpops_numalloc > numpop ===> warn and resize
      // newpops_numalloc < numpop ===> resize and fill
      if(options->newpops_numalloc > numpop)
	{
	  warning("Population number in datafile is smaller than in the relabel option!");
	  warning("Relabel option will be adjust to number of populations in datafile");
	  options->newpops_numalloc = numpop;
	}
      options->newpops = (long *) myrealloc(options->newpops,sizeof(long)* (size_t) numpop);
      for(i=options->newpops_numalloc;i<numpop;i++)
	{
	  options->newpops[i] = i+1;
	}
      options->newpops_numalloc = numpop;
    }
  memcpy(temp,options->newpops,sizeof(long)* (size_t) numpop);
  qsort((void *) temp, (size_t) numpop, sizeof(long), longcmp);
  pop = temp[0];
  newnumpop = 1;
  for(i=1;i<numpop;i++)
    {
      if(temp[i] != pop)
	{
	  newnumpop++;
	  pop = temp[i];
	}
    }
  options->newpops_numpop = newnumpop;
  return newnumpop;
}


///
/// the slice sampling stick size is allocated here for the options, this is 
/// done in three different places for: standard, genealogy, MPI
void alloc_sticksize(option_fmt *options, data_fmt *data, size_t size)
{
  (void) data;
  if( options->slice_sticksizes==NULL)
    {
      options->slice_sticksizes = (MYREAL *) mycalloc((size), sizeof(MYREAL));
    }
  else
    {
      options->slice_sticksizes = (MYREAL *) myrealloc(options->slice_sticksizes, (size)* sizeof(MYREAL));
    }
}

void set_skiploci(option_fmt *options, data_fmt *data, world_fmt **universe)
{ 
  long locus;
  long i;
  for(locus=0;locus<data->loci;locus++)
    {
      for (i = 0; i < options->heated_chains; i++)
	universe[i]->data->skiploci[locus] = data->skiploci[locus];
    }
}

///
/// runs the MCMC sampling 
void
run_sampler (option_fmt * options, data_fmt * data, world_fmt ** universe,
             int usize, long *outfilepos, long *Gmax)
{
  MYREAL var;
  long i;
  long maxreplicate;
  long treefilepos;
#ifdef MPI
  MPI_Request *irequests;	// contains requests generated in MPI_Isend
  MPI_Status *istatus;		// conatins stata from MPI_Wait_any
  long minnodes;
  long *twolongs;		// message sent to workers contains locus and replicate number
  
  if (myID < data->loci + 1)
    color = 1;
  else
    color = MPI_UNDEFINED;
  
  twolongs = (long *) mycalloc (TWO, sizeof (long));
#endif
  getseed(options);
#ifdef MPI
  if (myID == MASTER)
    {
#endif

#ifdef DEBUG
      char nowstr[LINESIZE];
      get_time (nowstr, "%c");
      printf("\n\nMaster start reading data at %s\n",nowstr);
#endif
      get_new_data (data->infile, data, options, EARTH);
#ifdef NEWSTUFF
      if (data->has_fastq)
	{
	  get_new_fastq(data->fastqinfile,data,options, EARTH);
	}
#endif
      //check_bayes_priors(options, data, EARTH);
      long locus;
      for(locus=0;locus < data->loci; locus++)
	{
	  //create_mixed_data (data, options);
	  init_sequences_aliases (EARTH, options, data, locus);
	  set_subloci_frequencies_alleles(EARTH, options, data, locus);	
	  finish_mutationmodel(EARTH,data,options,locus);
	}

#ifdef MPI
#ifdef DEBUG
      get_time (nowstr, "%c");
      printf("Master finished reading data at %s\n",nowstr);
#endif
      if (numcpu == 1)
        {
	  error
            ("This program was compiled to use a parallel computer\n and you tried to run it on only a single node.\nThis will not work because it uses a \n\"single_master-many_worker\" architecture \nand needs at least TWO nodes\n");
        }
      //#ifndef MPIDATAONDEMAND      
      broadcast_data_master (data, options);
      //#endif
    } /*end of myID==MASTER loop for MPI*/
  //#ifndef MPIDATAONDEMAND
  else
    {
      broadcast_data_worker (data, options, EARTH);
#ifndef MPIDATAONDEMAND
      long locus;
      for(locus=0;locus < data->loci; locus++)
	{
	  set_subloci_frequencies_alleles(EARTH, options, data, locus);	
	  finish_mutationmodel(EARTH,data,options,locus);
	}
#endif
    }
  //#endif  
#endif
  //set_plot (options);
  EARTH->cold = TRUE; // this is the cold chain when we use heating, the hotter chains have FALSE
  // sticksizes are filled in init_world but initialized here so make sure that 
  // the heated chains do not produce a memory leak by mutiply initializing the sticksizes
  alloc_sticksize(options,data, (size_t) data->numpop);
  // filling of options->slice_sticksizes
  // deferred after world->bayes initialization in init_world() 
  calculate_newpop_numpop(options,data);

  init_world (EARTH, data, options);
  //  set_meanmu(EARTH,options);
  if(EARTH->cold)
    {
      single_chain_var (EARTH, 0L, &var, NULL, NULL);
    }
#ifdef NEXUSTREE
  if (EARTH->options->treeprint != myNONE && myID == MASTER) //was = LASTCHAIN)
    {
      //FPRINTF(EARTH->treefile,"#NEXUS\n\nbegin trees;");
      nexus_treesheader(EARTH, options, data);
    }
#endif 
#ifdef PTHREADS
   tpool_init (&heating_pool, usize, usize, 0);
#endif
  if (options->heating)
    {
      //first in universe is the cold chain[EARTH]
      heating_init (universe, usize, data, options);
    }
  if (myID == MASTER)
    {
  /* report to screen */
      print_menu_options (EARTH, options, data);
      if (options->progress)
	print_data_summary (stdout, EARTH, options, data);
      /* print to outfile */
      pdf_master_init(EARTH, options, data);
      *outfilepos = print_title (EARTH, options);
      print_options (EARTH->outfile, EARTH, options, data);
      print_data_summary (EARTH->outfile, EARTH, options, data);
      print_data (EARTH, options, data);
      print_spectra (EARTH, options, data);
      if(options->murates_fromdata)
	{
	  if (options->progress)
	    print_mutationrate_weights(stdout, options->mu_rates, options->segregs,options->wattersons, data->loci);
	  print_mutationrate_weights(EARTH->outfile, options->mu_rates, options->segregs,options->wattersons, data->loci);
	}

      if (options->printfst)
	print_fst (EARTH, options, data, EARTH->fstparam);

      set_skiploci(options, data, universe);
      if(options->bayes_infer && myID == MASTER)
	{
	  if(options->has_bayesmdimfile)
	    print_bayes_mdimfileheader(EARTH->bayesmdimfile,
				       options->bayesmdiminterval, EARTH, data);
	}
    }
  else
    {
      //print_spectra (EARTH, options, data);// worker only calculate but do not print
    }

   if (options->lchains < 2 && options->replicatenum == 0)
	  {
	  options->replicate = 0;
	}
  maxreplicate = universe[0]->maxreplicate; 
    //(options->replicate
    //		  && options->replicatenum > 0) ? options->replicatenum : 1;
#ifdef MPI
    minnodes = MAX (0, numcpu - data->loci - 1);
    irequests = (MPI_Request *) mycalloc (minnodes + 1, sizeof (MPI_Request));
    istatus = (MPI_Status *) mycalloc (minnodes + 1, sizeof (MPI_Status));
    if (myID == MASTER)
      {
    //--------------------------------------------------------------------------------------------
    // copies the main haplotype structure from data to world->haplotypes for the MASTER
    // the worker node do that in the setup_locus function
	//printf("%i> before MAIN: moving individuals to haplotypes for all %li loci\n",myID, EARTH->loci);
	if(EARTH->data->haplotyping_report)
	  {
	    long locus;
	    for(locus=0;locus<EARTH->loci;locus++)
	      {
		//printf("%i> MAIN: moving individuals to haplotypes for locus %li\n",myID, locus);
		copy_individuals_from_data(data, EARTH->data, locus);
		save_haplotypes(EARTH,locus);
		reset_haplotypes(EARTH,locus);
	      }
	  }
	mpi_runloci_master (data->loci, EARTH->who, EARTH, options, data, options->readsum, options->menu);
      }
    else
      {
        mpi_runloci_worker (universe, usize, options, data,
                            &heating_pool, maxreplicate, &treefilepos, Gmax);
      }
    myfree(istatus);
    myfree(irequests);
#else /*MPI*/
    run_loci (universe, usize, options, data,
              &heating_pool, maxreplicate, &treefilepos, Gmax);    
#endif /*MPI*/

#ifdef BEAGLE
    beagle_stop(universe, usize);
#endif

#ifdef PTHREADS
    tpool_destroy (heating_pool, 1);
#endif
#ifdef MPI
    
    if (myID != MASTER)
    {
#endif
      if (options->heating)
	{
	  for (i = 0; i < options->heated_chains; i++)
	    {
	      //free_tree(universe[i]->root, universe[i]);
	      free_timevector (universe[i]->treetimes);
	    }
	}
      else
	{
	  //free_tree(EARTH->root, EARTH);
	  free_timevector (EARTH->treetimes);
	}
#ifdef MPI
    }
    myfree(twolongs);
#endif
}

/// generate samples of genealogies for all loci
void
run_loci (world_fmt ** universe, int usize, option_fmt * options,
          data_fmt * data, tpool_t * localheating_pool, long maxreplicate,
          long *treefilepos, long *Gmax)
{
    long locus;
    for (locus = 0; locus < data->loci; locus++)
    {
      if(!data->skiploci[locus])
	{
	  run_locus (universe, usize, options, data,
		     localheating_pool, maxreplicate, locus, treefilepos, Gmax);
	}
      free_datapart(data,universe[0],locus);
    }
}


///
/// save all genealogy summaries
void
get_bayeshist (world_fmt * world, option_fmt * options)
{
#ifndef MPI
  (void) world;
  (void) options;
  return;
#else
  (void) options;
  long maxreplicate = world->maxreplicate;    
  if (myID == MASTER)
    {
      mpi_results_master (MIGMPI_BAYESHIST, world, maxreplicate,
                            unpack_bayes_buffer);
    }
  return;
#endif
}

void 	recalc_skyline_values(world_fmt *world, option_fmt * options, long maxreplicate)
{
  (void) options;
  long i, j, locus;
  mighistloci_fmt *aa;
  long * eventnum;
  long npall = world->numpop2 + world->bayes->mu + world->species_model_size * 2 + world->grownum;
  const float invmax = (float) (1./maxreplicate);
  if(world->options->mighist && world->options->skyline)
    {
      for(locus=0; locus < world->loci; locus++)
	{
	  aa = &world->mighistloci[locus];
	  eventnum = world->mighistloci[locus].eventbinnum;
	  for(j=0; j< npall; j++)
	    {
	      for (i = 0; i < eventnum[j]; i++)
		{
		  aa->eventbins[j][i][0] *= invmax;
		  aa->eventbins[j][i][1] *= invmax;
		  aa->eventbins[j][i][2] *= invmax;
		  aa->eventbins[j][i][3] *= invmax;
		  aa->eventbins[j][i][4] *= invmax;
		  aa->eventbins[j][i][5] *= invmax;
		}
	    }
	}
    }
}


/// get migration event time data
/// and skyline data when present
void
get_mighistdata (world_fmt * world, option_fmt * options)
{
#ifndef MPI
  (void) world;
  (void) options;
#else
  long maxreplicate = world->maxreplicate;
  if (myID == MASTER)
    {
        if (options->mighist)
	  {
            //get all mighist data from the workers using unpack_mighist_buffer()
            mpi_results_master (MIGMPI_MIGHIST, world, maxreplicate,
                                unpack_mighist_buffer);
	    if(options->skyline)
	      {
		//get all mighist data from the workers using unpack_mighist_buffer()
		fprintf(stdout,"%i> before unpack skyline\n",myID);
		mpi_results_master (MIGMPI_SKYLINE, world, maxreplicate,
				    unpack_skyline_buffer);
	      }
	  }
	if(/*world->options->treeprint ==BEST &&*/ world->options->treeinmemory==TRUE)
	  {
            mpi_results_master (MIGMPI_TREESPACE, world, maxreplicate,
                                unpack_treespace_buffer);
	  }
	recalc_skyline_values(world, options, maxreplicate);
    }
#endif
    
}

/// swap the tree pointer between chains with different temperatures
void
heated_swap (world_fmt ** universe, worldoption_fmt * options)
{
    long sm = 0, mv = 0, vw = 0;
    long r;
    //debug for heating traverse issues
    if(options->heatedswap_off || (options->heating_count++ % options->heating_interval) != 0)
      return;
    switch (r = RANDINT (0L, options->heated_chains - 2L))
    {
        case 2:
            sm = chance_swap_tree (SUN, MERKUR);
            MERKUR->swapped += sm;
            break;
        case 1:
            mv = chance_swap_tree (MERKUR, VENUS);
            VENUS->swapped += mv;
            break;
        case 0:
            vw = chance_swap_tree (VENUS, EARTH);
            EARTH->swapped += vw;
            break;
        default:
            universe[r]->swapped += chance_swap_tree (universe[r + 1], universe[r]);
    }
}


/// change the type of the chain form short to long and reset the treelist
void
change_chaintype (long locus, char *type, world_fmt * world, long *increment,
                  long *oldsteps, long *chains, option_fmt * options)
{
  (void) locus;
    if (*type == 's')
    {
        *type = 'l';
        create_treetimelist (world, &(world->treetimes));
        set_bounds (increment, oldsteps, chains, options, *type);
        world->increment = *increment;
    }
}

/// prepare the next chain
void
prepare_next_chain (world_fmt ** universe, worldoption_fmt * options,
                    char type, long chain, long *chains, long *pluschain,
                    long locus, long replicate)
{
    long i;
    EARTH->likelihood[0] = EARTH->likelihood[EARTH->G];
    if (options->heating)
    {
        for (i = 1; i < options->heated_chains; i++)
            universe[i]->likelihood[0] = universe[i]->likelihood[universe[i]->G];
    }
    EARTH->treetimes[0].copies = 0;
    if (type == 'l')
    {
        if (options->replicate && options->replicatenum > 0)
            EARTH->chainlikes[locus][replicate] = EARTH->param_like;
        else
            EARTH->chainlikes[locus][chain] = EARTH->param_like;
    }
    EARTH->start = FALSE;
    if (type == 'l')
    {
        if (chain < *chains + *pluschain)
        {
            if (((EARTH->param_like > options->lcepsilon)
		 //|| (options->gelman && EARTH->convergence->gelmanmaxRall > GELMAN_MYSTIC_VALUE)
		 ) 
		&& !(options->replicate && options->replicatenum == 0))
            {
                (*chains)++;
                (*pluschain)--;
            }
        }
    }
}

///
/// resetting of accept and accept_freq for heated chains
void reset_heated_accept(world_fmt **universe, long unum)
{
  long i;
  for(i=0; i < unum; i++)
    {
      universe[i]->swapped = 0;
      universe[i]->accept = 0L;
      universe[i]->accept_freq = 0.;
    }
}


/// print progress about heated chains
void
print_heating_progress (world_fmt ** universe, worldoption_fmt * options,
                        long stepinc)
{
  MYREAL fstepinc = (MYREAL) stepinc;
    if (options->heating)
    {
#ifdef DEBUG_MPI
      printf("%i> stepinc=%li =========\n",myID, stepinc);
#endif
        universe[0]->accept_freq /= fstepinc;
        universe[1]->accept_freq /= fstepinc;
        universe[2]->accept_freq /= fstepinc;
        universe[3]->accept_freq /= fstepinc;
        if (options->progress)
            print_heating_progress2 (stdout, options, universe);
        if (options->writelog)
            print_heating_progress2 (options->logfile, options, universe);
    }
}

///
/// sub-function of print_heating, in MPI mode, the string printing will minimize data transfer.
void
print_heating_progress2 (FILE * file, worldoption_fmt * options,
                         world_fmt ** universe)
{
  long i;
  char *plog;
  long plogsize = 0L;
  char nowstr[STRSIZE];
  char dots[STRSIZE];
  world_fmt *world = universe[0];
  const  long npp = world->numpop2 + (long) world->bayes->mu + 2 * world->species_model_size + world->grownum;
  get_time (nowstr, "%H:%M:%S");
  plog = (char *) mycalloc(LONGLINESIZE,sizeof(char));
#ifdef MPI
  plogsize = sprintf(plog, "[%3i] %s   Sampling Temp[%li]:", myID, nowstr, options->heated_chains);
#else
  plogsize = sprintf(plog, "%s   Sampling Temp[%li]: ", nowstr, options->heated_chains);
#endif
  if (options->heated_chains > 4)
    strcpy(dots,",...");
  else
    dots[0]='\0';
  plogsize += sprintf (plog + plogsize, "(%.4g,%.4g,%.4g,%.4g%s) ",universe[0]->averageheat,
		       universe[1]->averageheat,universe[2]->averageheat,universe[3]->averageheat,dots);
  plogsize += sprintf (plog + plogsize, "Acc(%.2f,%.2f,%.2f,%.2f%s) ",universe[0]->accept_freq,
		       universe[1]->accept_freq,universe[2]->accept_freq,universe[3]->accept_freq,dots);
  plogsize += sprintf (plog + plogsize, "Swap(%li,%li,%li%s)\n                           Param = {",universe[0]->swapped,
			   universe[1]->swapped,universe[2]->swapped /*,universe[3]->swapped*/,dots);
  for (i=0;i<npp;i++)
    {
      plogsize += sprintf (plog + plogsize, " %.4g",universe[0]->param0[i]);
    }
  sprintf (plog + plogsize, "}\n");

  FPRINTF(file,"%s",plog);
  myfree(plog);
}


boolean analyze_oldbayesdata(world_fmt **universe, option_fmt *options, data_fmt *data, long *outfilepos)
{
  (void) outfilepos;
  world_fmt *world = universe[0];
#ifndef MPI
  long pop;
#endif
  long locus;
  char **files = NULL;
  long numfiles=1;
  charvec2d(&files,numfiles,LINESIZE);
  // insert code for number of files and allocation and naming
  strcpy(files[0],options->bayesmdimfilename);
  // read data from bayesfile
  world->cold=TRUE;
  // reads some minimal information form the header of the bayesallfile
  // this only works with files written with migrate 2.5+
#ifdef MPI
  if(myID==MASTER)
    {
      warning("does not work correctly yet");
      read_from_bayesmdim_minimal_info(world->bayesmdimfile, world, options, data);
      read_geofile (data, options, world->numpop);
      alloc_sticksize(options,data, (world->numpop*world->numpop+2L*world->numpop));
      options->newpops_numalloc = world->numpop;
      options->newpops = (long*) mycalloc(options->newpops_numalloc, sizeof(long));
      calculate_newpop_numpop(options,data);
      data->locusweight =
	(MYREAL *) mycalloc (1, sizeof (MYREAL) * (size_t) (data->loci + 1));
      for (locus=0;locus < data->loci; locus++)
	{
	  data->locusweight[locus]=1.0;
	}
      /*options->newpops_numalloc = world->numpop;
      options->newpops = (long*) mycalloc(options->newpops_numalloc, sizeof(long));
      for(pop=0;pop<world->numpop;pop++)
      options->newpops[pop]=pop+1;*/
      init_world (world, data, options);
	read_bayes_fromfile(world->bayesmdimfile, world, options,files, numfiles);
    }
#else
  read_from_bayesmdim_minimal_info(world->bayesmdimfile, world, options, data);
  read_geofile (data, options, world->numpop);
  alloc_sticksize(options,data, (size_t) (data->numpop*data->numpop+2*data->numpop));
  options->newpops_numalloc = world->numpop;
  options->newpops = (long*) mycalloc(options->newpops_numalloc, sizeof(long));
  for(pop=0;pop<world->numpop;pop++)
    options->newpops[pop]=pop+1;
  data->locusweight =
    (MYREAL *) mycalloc (1, sizeof (MYREAL) * (size_t) (data->loci + 1));
  for (locus=0;locus < data->loci; locus++)
    {
      data->locusweight[locus]=1.0;
    }
  init_world (world, data, options);
    read_bayes_fromfile(world->bayesmdimfile, world, options,files, numfiles);
#endif
  //set_meanmu(world,options);
  pdf_master_init(world, options, data);
  for(locus=0;locus<world->loci;locus++)
    {
      if(world->data->skiploci[locus])
	continue;
      if(world->bayes->histogram[locus].results == NULL)
	{
	  world->data->skiploci[locus] = TRUE;
	  continue;
	}
      calc_hpd_credibility(world, locus, world->numpop2, world->numpop2 + world->bayes->mu+2*world->species_model_size + world->grownum);
    }
  myfree(files);
  return TRUE;
}

/// \brief setup a locus
/// 
/// Set up a the structures for a locus. Run-Parameters are copied into the 
/// major structure WORLD from options and data. a first genealogy is build
/// adjusted to the starting parameters (changing times of nodes) and a first
/// conditional likelihood is calculated for all nodes, also the tree-parallele
/// time structure is generated.
/// Several run controllers are checked: replicate, datatype
/// \returns 0 { failed to setup structure for locus because there is no data }
/// \returns 1 { succeeded to setup locus-structure } 
int
setup_locus (long locus, world_fmt * world, option_fmt * options,
             data_fmt * data)
{
  world->locus = locus;
  set_param (world, data, options, locus);
  // replaced with more informative reporting in burnin_bayes(). 
  //if (world->options->progress)
  //  print_menu_locus (stdout, world, locus);
  //if (world->options->writelog)
  //  print_menu_locus (world->options->logfile, world, locus);
  world->start = TRUE;
  reset_growth(world);

  
  if(options->has_unassigned && locus>0)
    {
      reset_all_assigned_nodes(world);
    }
  if(options->haplotyping)
    {
      world->data->haplotyping = TRUE;
      copy_individuals_from_data(data, world->data, locus);
      world->data->numindividuals = data->numindividuals[locus];
      reset_ID_nodelist(locus, world);
      buildtree (world, options, data, locus);
      check_individual_nodes(locus,world);
      set_individuals_request_haplotyping(world, data, locus);
      cleanup_individual_nodes(locus,world);
    }
  else
    {
      buildtree (world, options, data, locus);
    }
#ifdef BEAGLE
      finish_mutationmodel(world, data, options,locus);
  if(world->beagle->instance_handle == NULL)
    {
      init_beagle(world,locus);
    }
  else
    {
      reinit_beagle(world,locus);
    }
  set_beagle_instances(world, locus);
  fill_beagle_instances(world, locus);
#endif
  if (data->skiploci[locus])	/* skip loci with no tips */
    return 0;
  create_treetimelist (world, &world->treetimes);
  fix_times (world, options);
#ifdef DEBUG
  //  printf(" fixed times \n");
#endif
  first_smooth (world, locus);
  world->likelihood[0] = treelikelihood (world);
  world->allikemax = world->likelihood[0]; //-MYREAL_MAX;	/* for best tree option */
  world->treelen = 0.0;
  calc_treelength (world->root->next->back, &world->treelen);

#ifdef UEP
  if (options->uep)
    {
        world->treelen = 0.0;
        if (options->ueprate > 0.0)
	  calc_treelength (world->root->next->back, &world->treelen);
        //world->ueplikelihood = ueplikelihood(world);
        //world->likelihood[0] += world->ueplikelihood;
    }
#endif
  world->has_proposal_details=FALSE;
  return 1;
}

/// condense a single genealogy into a sufficient statistic for the maximization phase
void
condense_time (world_fmt * world, long *step, long *j, MYREAL * accepted,
               long *G, long *steps, long oldsteps, long rep)
{
  static long n = 0;
  world->rep = rep;
  if(world->options->bayes_infer)
    {
      if (world->in_last_chain)
	{
	  if (*step == 0)
	    return;

	  // syncs treeupdate and paramupdate
	  if(world->bayesaccept == -1)
	    {
	      world->bayes->oldval = probg_treetimes(world);
	    }
	  bayes_save (world, *step * world->options->lincr);
	  store_events(world, world->treetimes, world->numpop, rep);
	  *accepted += *j;
	  return;
	}
      else
	{
	  warning("Bayesian analysis: do not run multiple long chains, run 1 long chain! [nothing will be wrong but it is a waste of time]");
	  return;
	}
    }
  if (*step == 0)
    {
      copy_time (world, world->treetimes, FIRSTSTEP, 0L, world->numpop, rep);
      *accepted += *j;
      //for INTEGRATEDLIKE
      n = 0;
    }
  else
    {
      copy_time (world, world->treetimes, *G, *G + (long) (*j > 0),
		 world->numpop, rep);
      if (*step < *steps)
        {
	  //handled in advance_world * G += (long) (*j > 0);
	  *accepted += *j;
        }
    }
  if (*step >= oldsteps - 1 && world->options->movingsteps)
    {
      if (((MYREAL) (*G + 1)) < world->options->acceptfreq * oldsteps)
        {
	  (*steps)++;
        }
    }
}

void print_theta0(FILE *file, world_fmt *world, long maxreplicate)
{
  long i;
  long z;
  long locus;
  long rep;
  long nn = world->numpop2;
  fprintf(file,"Locus Replicate Parameters\n");
  fprintf(file,"-----------------------------------------------------------\n");
  for(locus=0; locus < world->loci; locus++)
    {
      for(rep = 0; rep < maxreplicate; rep++)
	{
	  fprintf(file,"%5li  %5li ",locus, rep+1);
	  for(i=0; i < nn; i++)
	    {
	      if(world->bayes->map[i][1] == INVALID)
		continue;
	      else
		{
		  z  = world->bayes->map[i][1];
		}
	      fprintf(file,"%f ", world->atl[rep][world->locus].param0[z]);
	    }
	  fprintf(file,"\n");
	}
    }
  fprintf(file,"\n\n");
}

/// set type of replication scheme
long
set_repkind (option_fmt * options)
{
    if (options->replicate)
    {
        if (options->replicatenum == 0)
            return MULTIPLECHAIN;
        else
           return MULTIPLERUN;
    }
    return SINGLECHAIN;
}

/// intialize heating scheme
void
heating_init (world_fmt ** universe, int usize, data_fmt * data,
              option_fmt * options)
{
    long chain;
    //char tmp[20];
    for (chain = 1; chain < usize; ++chain)
    {
      //MYSNPRINTF (tmp, 20L, "%10.5f", options->heat[chain]);
        create_world (&universe[chain], data->loci);
	universe[chain]->cold = FALSE;
	universe[chain]->heatid = chain;
        init_world (universe[chain], data, options);

    }
}

/// prepare for heating scheme: Step I
void heating_prepare (world_fmt ** universe, int usize,
                 option_fmt * options, data_fmt * data, long rep)
{
    (void) rep;
    long chain;
    long locus = EARTH->locus;
    EARTH->heat = options->heat[0];
    for (chain = 1; chain < usize; ++chain)
    {
      universe[chain]->heatid = chain;
      klone_mutationmodel(universe[chain],EARTH, data, locus);
      klone (EARTH,universe[chain], options, data, options->heat[chain]);
      if(options->has_unassigned && locus>0)
        {
	  reset_all_assigned_nodes(universe[chain]);
        }
      if(options->haplotyping)
        {
	  universe[chain]->data->haplotyping = TRUE;
	  copy_individuals_from_data(data, universe[chain]->data, locus);
	  buildtree (universe[chain], options, data, locus);
	  check_individual_nodes(locus,universe[chain]);
	  set_individuals_request_haplotyping(universe[chain], data, locus);
	  cleanup_individual_nodes(locus,universe[chain]);
        }
      else
        {
	  buildtree (universe[chain], options, data, locus);
        }

#ifdef BEAGLE
        if(universe[chain]->beagle->instance_handle == NULL)
        {
            init_beagle(universe[chain],locus);
        }
        else
        {
            reinit_beagle(universe[chain],locus);
        }
        set_beagle_instances(universe[chain],locus);
        fill_beagle_instances(universe[chain],locus);
#endif
        //old loc        klone (EARTH,universe[chain], options, data, options->heat[chain]);
	klone_tree_setup(universe[chain], options);
        
    }
}


/// prepare for heating scheme: Step II
void
heating_prepare2 (world_fmt ** universe, int usize)
{
    long chain;
    universe[0]->averageheat = 1.0;
    for (chain = 1; chain < usize; ++chain)
    {
        universe[chain]->G = 0;
        universe[chain]->averageheat = universe[chain]->options->adaptiveheat!=NOTADAPTIVE ?
            0.0 : 1. / universe[chain]->heat;
    }
    for (chain = 1; chain < usize; ++chain)
    {
        clone_polish (EARTH, universe[chain]);
    }
}

/// \brief generates replicate number
///
/// generates replicate number
/// \returns 0 no replicate
/// \returns rep number of replicates
long
replicate_number (option_fmt * options, long chain, char type, long replicate)
{
    long rep = 0;
    if (options->replicate && options->replicatenum > 0)
        rep = replicate;
    else
    {
        if (type == 'l' && options->replicate)
            rep = chain;
        else
            rep = 0;
    }
    return rep;
}

/// generate genealogies for a single locus
//  \callgraph
void
run_locus (world_fmt ** universe, int usize, option_fmt * options,
           data_fmt * data, tpool_t * localheating_pool, long maxreplicate,
           long locus, long *treefilepos, long *Gmax)
{
  long i;
  long convergence_len = universe[0]->numpop2 + 3 * universe[0]->numpop;
  long replicate;
  // DEBUG Hapmap stuff
  if(EARTH->data->skiploci[locus])
    return;
  if(options->bayes_infer)
    convergence_len = universe[0]->numpop2 + 1;
  memset(EARTH->convergence->chain_means,0, (size_t) (maxreplicate * convergence_len) * sizeof(MYREAL));
  memset(EARTH->convergence->chain_s,0, (size_t) (maxreplicate * convergence_len) * sizeof(MYREAL));
  memset(EARTH->convergence->gelmanmeanmaxR,0,(size_t) (maxreplicate * maxreplicate) * sizeof(MYREAL));
  for (replicate = 0; replicate < maxreplicate; replicate++)
    {
      run_replicate (locus, replicate, universe, options, data,
		     localheating_pool, usize, treefilepos, Gmax);
      if(EARTH->cold && EARTH->options->replicatenum > 0)
	{
	  if(options->gelman)
	    {
	      chain_means(&EARTH->convergence->chain_means[replicate * convergence_len], EARTH);
	      calc_chain_s(EARTH->convergence->chain_s, EARTH->convergence->chain_means, EARTH, 
			   replicate);
	      if (replicate > 0)
		{
		  convergence_check_bayes(EARTH, maxreplicate);
		  convergence_progress(stdout, EARTH);
		}
	    }
	}
    }    
#ifdef UEP
  if (options->uep)
    show_uep_store (EARTH);
#endif
  if (options->replicate && options->replicatenum > 0)
    {
      EARTH->repkind = MULTIPLERUN;
      if (options->bayes_infer)
        {
	  if(!options->has_bayesmdimfile)
	    calculate_credibility_interval (EARTH, locus);
	  //	  return;
        }
#ifdef LONGSUM
      change_longsum_times (EARTH);
#endif /*LONGSUM*/
    }
  else
    {
      if(!options->has_bayesmdimfile)
	{
	  //adjust_bayes_bins(EARTH,locus);
	  calculate_credibility_interval (EARTH, locus);
	}
    }
  // here we store the haplotypes into a save container
  //
  if(options->haplotyping_report)
    {
#ifdef MPI
      if(myID!=MASTER)
	{
	  //printf("%i> moving individuals to haplotypes for locus %li",myID, locus);
	  save_haplotypes(EARTH, locus);
	  reset_haplotypes(EARTH, locus);
	}
#else
      save_haplotypes(EARTH, locus);
      reset_haplotypes(EARTH, locus);
#endif
    }
  if (options->heating)
    {
      for (i = 0; i < options->heated_chains; i++)
        {
	  if(options->haplotyping_report)
	    reset_haplotypes(universe[i], locus);
	  free_tree (universe[i]->root, universe[i]);
	  //   free_timevector(universe[i]->treetimes);
	  myfree(universe[i]->nodep);
        }
    }
  else
    {
      free_tree (EARTH->root, EARTH);
      myfree(EARTH->nodep);
      //  free_timevector(EARTH->treetimes);
    }
  if(options->bayes_infer)
    bayes_reset (universe[0]); 
}

/// \brief updates the tree and records acceptance
///
/// updates the tree and records acceptance
/// when in bayesian mode change between changing the tree and the parameters
void
run_one_update (world_fmt * world)
{
  //long count = 0;
  //long np = world->numpop2;
  //fprintf(stderr,"%i> locus=%li heat=%f\n",myID, world->locus, world->heat);
  //traverse_check(crawlback (world->root->next));
  // propose tree, parameter, and haplotype moves
  updating(world);//uses world->options->choices, and world-G
}

/// \brief updates all trees, controls updates and heating
///
/// controls updates and temperatures (threaded or unthreaded)
/// updates the tree and/or parameters
//  \callgraph
void
run_updates (world_fmt ** universe,
             int usize,
             option_fmt * options,
             tpool_t * localheating_pool,
             long inc, long increment, long step, long steps)
{
#ifndef PTHREADS
  (void) usize;
  (void) localheating_pool;
#endif
#ifndef  PTHREADS
#ifndef GRANDCENTRAL
    long ii;
#endif
#endif
    if (options->heating)
    {
#ifdef PTHREADS			/*using threads and running on an SMP machine */
        fill_tpool (*localheating_pool, universe, usize);
        tpool_synchronize (*localheating_pool, usize);
#else /*heating but not using threads */
#ifdef GRANDCENTRAL
        //dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
	//semaphore = dispatch_semaphore_create(1);
        dispatch_apply((size_t) options->heated_chains, queue,
                       ^(unsigned long ii) {
			 sfmtp = sfmtH[ii];
                         run_one_update (universe[ii]);
                       }
                       );
	//dispatch_release(semaphore);
#else
	//count++;
	run_one_update (universe[0]);
        for (ii = 1; ii < options->heated_chains; ii++)
        {
	  run_one_update (universe[ii]);
        }
#endif /*GRANDCENTRAL*/
#endif /* end of not-using threads */
	if(options->bayes_infer && RANDUM() > options->updateratio)
	  {
#ifdef UEP
	    if (options->uep && EARTH->in_last_chain)
	      update_uepanc (EARTH);
#endif
	    return;
	  }
	else
	  {
	    heated_swap (universe, EARTH->options);
	    switch (options->adaptiveheat)
	      {
	      case STANDARD:
		adjust_temperatures (universe, options->heated_chains,
				     inc /*rement */  + step * increment,
				     steps * increment);
		break;
	      case BOUNDED:
		adjust_temperatures_bounded (universe, options->heated_chains,
				     inc /*rement */  + step * increment,
				     steps * increment);
		break;
	      case NOTADAPTIVE:
	      default:
		break;
	      }
	  }
    }
    else
      {		/* no heating */
	run_one_update (EARTH);
      }
#ifdef UEP
    if (options->uep && EARTH->in_last_chain)
      update_uepanc (EARTH);
#endif
}

/// generates trees over the interval increments
void run_increments (world_fmt ** universe,
                long usize,
                option_fmt * options,
                tpool_t * localheating_pool,
                long increment, long step, long steps)
{
  static long treefilepos; /* write position in the treefile */
    long i;
    for (i = 0; i < increment; i++)
    {
        EARTH->actualinc = i;
        run_updates (universe, (int) usize, options,
                     localheating_pool, i, increment, step, steps);
        
        if (EARTH->likelihood[EARTH->G] > EARTH->maxdatallike)
            EARTH->maxdatallike = EARTH->likelihood[EARTH->G];
    }
    if (EARTH->options->treeprint != myNONE)
      print_tree (EARTH, EARTH->G, &treefilepos);
    if(EARTH->has_speciation && options->recorddivtime)
      record_parameters(EARTH);
}

/// \brief run a chain
/// run a chain
void run_steps (world_fmt ** universe,
           long usize,
           option_fmt * options,
           tpool_t * localheating_pool, long increment, long steps, long rep)
{
    long step=0;
    long oldsteps = steps;
    long ii;
    MYREAL var=0.0;
    if(EARTH->cold && EARTH->options->bayes_infer)
      single_chain_var (NULL,  step, &var, NULL, NULL);
#ifdef MPI
    // does not work? check_memory_limit();
#endif
    long start = 0;
    if (options->checkpointing)
      start = options->unfinished[EARTH->locus][rep];
#ifndef MPI
    const boolean progress = EARTH->options->progress && EARTH->cold;
    long delta = ((steps > 100000) ? (steps / 10) : (steps/3));
    boolean reportdone=FALSE;
    long starttime = expected_end_burnin(EARTH,0.0, 0, " ");
    if (delta==0)
      delta=1;
#endif
    for (step = start; step < steps; step++)
    {
#ifdef MPI
      // check_memory_limit() does not work;
#else
      // too much writing for MPI therefore we do not report progress on burnin
      if (!reportdone && progress && (step+1) % delta == 0)
	{
	  reportdone=TRUE;
	  expected_end_burnin(EARTH, 100.0*(step-1.0)/steps, starttime,"Sampling complete ");
	}
#endif
        EARTH->increment = increment;
        run_increments (universe, usize, options, localheating_pool,
                        increment, step, steps);
	
#ifdef UEP
        store_uep (EARTH);
#endif
	//	printf("step=%li inc=%li locus=%li\n",step,increment,EARTH->locus);

        condense_time (EARTH, &step, &EARTH->accept,
                       &EARTH->accept_freq, &EARTH->G, &steps, oldsteps, rep);
	if(EARTH->options->bayes_infer)
	  {
	    calculate_BF(universe,options);
	    if(EARTH->cold)
	      {
		single_chain_var (EARTH, step, &var, EARTH->autocorrelation, EARTH->effective_sample);
	      }
	  }
        if (step > 0)
        {
            advance_clone_like (EARTH, EARTH->accept, &EARTH->G);
            EARTH->accept = 0L;
        }
        // if we do heating
        if (EARTH->options->heating)
        {
            for (ii = 1; ii < EARTH->options->heated_chains; ++ii)
            {
                universe[ii]->accept_freq += universe[ii]->accept;
                advance_clone_like (universe[ii],
                                    universe[ii]->accept, &(universe[ii])->G);
                universe[ii]->accept = 0L;
            }
        }
    }
}

/// run all chains for short and then for long
void run_chains (world_fmt ** universe,
            long usize,
            option_fmt * options,
            tpool_t * localheating_pool,
            long replicate,
            long chains,
            char type,
            long increment, long oldsteps, long *treefilepos, long *Gmax)
{
  (void) chains;
  (void) treefilepos;
  //long i;
  long steps;
  long rep;
  char nowstr[STRSIZE]; // will hold time of day
  //const long locus = EARTH->locus;
  //long kind;
  //long pluschain = 0;
  const MYREAL treeupdateratio = (options->bayes_infer ? options->updateratio : 1.0);

  //for (chain = 0;
  //     chain < chains || (type == 'l' && chain >= chains
  //			  && EARTH->param_like > options->lcepsilon); chain++)
  //  {
  EARTH->in_last_chain = TRUE;
  type = 'l';
  if (options->heating)
    MERKUR->swapped = VENUS->swapped = EARTH->swapped = 0;
  rep = replicate_number (options, 0, type, replicate);
  EARTH->rep = rep;
  precalc_world (EARTH);
  polish_world (EARTH);
  if (options->heating)
    {
      heating_prepare2 (universe, (int) usize);
#ifdef GRANDCENTRAL
      dispatch_apply((size_t) options->heated_chains, queue,
		     ^(unsigned long ii) {
		       sfmtp = sfmtH[ii];		      
		       burnin_chain (universe[ii]);
		     }
		     );
#else      
      int i;
      for(i=0; i<EARTH->options->heated_chains; i++)
        {
          burnin_chain (universe[i]);
        }
#endif
      EARTH->G = 0;
      EARTH->accept_freq = 0.;
      EARTH->accept = 0;
      steps = oldsteps;
      EARTH->maxdatallike = EARTH->likelihood[0];
      reset_heated_accept(universe,EARTH->options->heated_chains);
    }
  else
    {
      burnin_chain (EARTH);
      EARTH->G = 0;
      EARTH->accept_freq = 0.;
      EARTH->accept = 0;
      steps = oldsteps;
      EARTH->maxdatallike = EARTH->likelihood[0];
    }
  get_time (nowstr, "%H:%M:%S");
#ifdef MPI
  FPRINTF(stdout,"[%3i] %8.8s   Sampling of %li steps (Locus: %li/%li, Replicate: %li/%li) \n",myID, nowstr, steps*increment, 1+EARTH->locus, EARTH->loci, 1+EARTH->rep, EARTH->maxreplicate);
#else
  FPRINTF(stdout,"%8.8s   Sampling of %li steps (Locus: %li/%li, Replicate: %li/%li) \n", nowstr, steps*increment, 1+EARTH->locus, EARTH->loci, 1+EARTH->rep, EARTH->maxreplicate);
#endif      
  run_steps (universe, usize, options, localheating_pool,
	     increment, steps, rep);
  
  //decide_plot (EARTH->options, chain, chains, type);
  // prepare for parameter estimation, precompute and copy
  memcpy (EARTH->param00, EARTH->param0,
	  sizeof (MYREAL) * (size_t) EARTH->numpop2);
  EARTH->repkind = SINGLECHAIN;
  //kind = SINGLELOCUS;
  *Gmax = (EARTH->G > *Gmax) ? EARTH->G + 1 : (*Gmax);
#ifdef LONGSUM
  change_longsum_times (EARTH);
#endif /*LONGSUM*/
  if (EARTH->options->heating)
    {
      print_heating_progress (universe, EARTH->options, (long) (increment * steps * treeupdateratio));
      reset_heated_accept(universe,EARTH->options->heated_chains);
    }
  else
    {
      print_progress (EARTH->options, EARTH, rep,
		      (long) (increment * steps * treeupdateratio), (long) EARTH->accept_freq);
      EARTH->accept_freq = 0.;
      EARTH->accept = 0L;
    }
  // cleanup and prepare next, should not be necessary anymore
  // prepare_next_chain (universe, EARTH->options, type, chain,
  //		      &chains, &pluschain, locus, replicate);
  //} //end chains
}



/// \brief run a replicate
///
/// run a replicate
//  \callgraph
void
run_replicate (long locus,
               long replicate,
               world_fmt ** universe,
               option_fmt * options,
               data_fmt * data,
               tpool_t * localheating_pool,
               int usize, long *treefilepos, long *Gmax)
{
    long chains = 0;
    char type = 's';
    long runs;
    long increment = 1;
    long oldsteps = 0;
    /* loop over independent replicates----------------------- */
    /* but make sure that we get single chain estimates, too   */
    EARTH->locus = locus;
    EARTH->repkind = SINGLECHAIN;
    EARTH->replicate = replicate;

    if (options->checkpointing)
      {
	long start = options->unfinished[locus][replicate];
	if (start >= options->lsteps - 1)
	  return;
      }
    if (setup_locus (locus, EARTH, options, data) == 0)
        return;
    if (options->heating)
        heating_prepare (universe, usize, options, data, replicate);
    type = 's';
    runs = 0;
    if(options->bayes_infer)
      {
	type='l';
      }
    /* short and long chains ----------------------------- */
    set_bounds (&increment, &oldsteps, &chains, options, type);
	//trying to find problem with parallel programming and bayesallfile
    //fprintf(stdout,"%i> PROBLEM heat=%f custm=%s custm2=%s\n",myID, EARTH->heat, EARTH->options->custm, EARTH->options->custm2);
    //exit(-1);


    EARTH->increment = increment;
    while (runs-- >= 0)
    {
        if (myID == MASTER)
        {
            print_menu_chain (type, FIRSTCHAIN, oldsteps, EARTH,
                              options, replicate);
            if (EARTH->options->treeprint == ALL)
                print_tree (EARTH, 0, treefilepos);
        }
        EARTH->chains = chains;
        /* loop over chains--------------------------- */
        run_chains (universe, usize, options, localheating_pool, replicate, chains, type, increment, oldsteps, treefilepos, Gmax);	//PB 020203
        change_chaintype (locus, &type, EARTH, &increment, &oldsteps,
                          &chains, options);

        /* evaluate multiple long chain estimates */
        if (runs < 0 && options->replicate && options->replicatenum == 0)
        {
            EARTH->repkind = MULTIPLECHAIN;
#ifdef LONGSUM
            change_longsum_times (EARTH);
#endif	   /*LONGSUM*/
        }
    }				
    // end runs
    // collect correlation information of a locus and stick it into the archive,
    // the autocorrelation is averaged but the ess is summed up assuming that ll loci are independent
    if(EARTH->cold && EARTH->options->bayes_infer)
      {
	  collect_ess_values(EARTH);
	  collect_acceptance(EARTH);
      }
    // cleanup locus: the free_tree() command in buildtree() in setup_locus()
    // fails when a locus i has more subloci than locus i-1 because it tries to 
    // free node data for conditional likelihoods that never was allocated,
    // this should problably more clearly formalized into a cleanup_locus() function.
    free_tree(EARTH->root, EARTH);
    if (options->heating)
      {
	int chain;
	const long hc = options->heated_chains;
	for(chain=0; chain<hc;chain++)
	  {
	    free_tree(universe[chain]->root, universe[chain]);
	  }
      }
}




void print_bf_values(world_fmt * world)
{
  long locus;
  long t;
  const long hc = world->options->heated_chains;
  for(locus=0;locus<world->loci;locus++)
    {
      if(world->data->skiploci[locus])
	continue;
      printf("=====> locus=%li ",locus );
      for(t=0; t < hc; t++)
	{
	  printf("%f ",world->bf[locus*hc+t]);
	}
      printf("\n");
    }
}



  


#ifdef LONGSUM
/// put timepoints into the total tree to allow for bottlenecks and growth 
long
find_chaincutoffs (MYREAL * treetime, world_fmt * world, long locus, long r)
{
    
    long j;
    long T;
    MYREAL totaltime = 0.;
    long copies = 0;
    tarchive_fmt *atl;
    atl = world->atl[r][locus].tl;
    for (j = 0; j < world->atl[r][locus].T; j++)
    {
        T = atl[j].longsumlen - 1;
        copies += atl[j].copies;
        totaltime = atl[j].longsum[T].eventtime / 3.;
        treetime[0] += totaltime;
        treetime[1] += totaltime + totaltime;
        treetime[2] += totaltime + totaltime + totaltime;
    }
    return copies;
}

/// find the cutoffs for the locus, averaging over multiple replicates
void
find_loci_cutoffs (MYREAL * treetime, world_fmt * world)
{
    long r;
    long locus;
    long count = 0;
    long repstart = 0;
    long repstop = 1;
    memset (treetime, 0, sizeof (MYREAL) * 3);
    set_replicates (world, world->repkind, world->rep, &repstart, &repstop);
    
    for (locus = 0; locus < world->loci; locus++)
    {
        for (r = repstart; r < repstop; r++)
        {
            count += find_chaincutoffs (treetime, world, locus, r);
        }
    }
    treetime[0] /= (MYREAL) count;
    treetime[1] /= (MYREAL) count;
    treetime[2] /= (MYREAL) count;
}


/// find the cutoffs avaergin over all replicates 
void
find_replicate_cutoffs (MYREAL * treetime, world_fmt * world)
{
    long r;
    long count = 0;
    long repstart = 0;
    long repstop = 1;
    memset (treetime, 0, sizeof (MYREAL) * 3);
    set_replicates (world, world->repkind, world->rep, &repstart, &repstop);
    
    for (r = repstart; r < repstop; r++)
    {
        count += find_chaincutoffs (treetime, world, world->locus, r);
    }
    treetime[0] /= (MYREAL) count;
    treetime[1] /= (MYREAL) count;
    treetime[2] /= (MYREAL) count;
}


/// change the cutoff times
void
change_longsum_times (world_fmt * world)
{
    long i;
    long numpop3 = 3 * world->numpop;
    long numpop6 = 2 * numpop3;
    MYREAL *treetime;
    treetime = (MYREAL *) mycalloc (3, sizeof (MYREAL));	// we have 3 classes
    if (world->locus < world->loci)
        find_replicate_cutoffs (treetime, world);
    else
        find_loci_cutoffs (treetime, world);
    for (i = numpop3; i < numpop6; i += 3)
    {
        world->flucrates[i] = treetime[0];
        world->flucrates[i + 1] = treetime[1];
        world->flucrates[i + 2] = treetime[2];
    }
    myfree(treetime);
}

#endif /*LONGSUM*/

///
/// check_parmfile checks what is called and also whether we should honor the menu=YES in the file
/// user uses another parmfile than the default when there is a parameter associated with the program call
/// this function checks the first two line of the file to find out whether this is a parmfile or potentially
/// a datafile/ 
/// The addition of an option to override the menu command in the parmfile
/// migrate-n -menu ==> use menu independent of the parmfile setting
/// migrate-n -nomenu ==> does use menu even if parmfile says so
/// migrate-n -version ==> prints the version and quits
/// migrate-n -help   ==> prints a quick help of commandline options
boolean  check_parmfile(long argcount, char **arguments, char *parmfilename)
{
  int argument=0;
  int len;
  boolean  usemenu = OVERRIDE_NO;
  if (argcount > 1)			
    {
      argument = 1;
      while(argument < argcount)
	{
	  //fprintf(stderr,"<|%s|>\n",argv[argument]);
	  //fflush(stderr);
	  if(arguments[argument][0]!='-')
	    {
	      len = (int) (strlen (arguments[argument]) + 1);
	      sprintf (parmfilename, "%-.*s", len, arguments[argument]);		       
	    }
	  else
	    {

	      switch(arguments[argument][1])
		{
		case 's':
		  simulator = TRUE;
		  break;
		case 'm':
		  usemenu = OVERRIDE_MENU;
		  break;
		case 'n':
		  usemenu = OVERRIDE_NOMENU;
		  break;
		case 'v':
		  printf("Migrate-n version: %s, subversion: %s\n", MIGRATEVERSION, MIGRATESUBVERSION);
		  printf("Operating System:  %s\n", MYSYSTEM ); 
#ifdef MPI
#ifdef OMPI_MAJOR_VERSION
		  printf("OPENMPI version: %i.%i\n", OMPI_MAJOR_VERSION, OMPI_MINOR_VERSION);
#endif
#endif
		  exit(-1);
		  //break;
		case 'h':
		  printf("Migrate-n version: %s, subversion: %s\n", MIGRATEVERSION, MIGRATESUBVERSION);
		  printf("(beerli@fsu.edu)\n\n");
		  printf("commandline options:\n");
		  printf("-version  # shows the current version and exits\n");
		  printf("-help     # show this menu and exit\n");
		  printf("\n");
		  printf("-nomenu   # does not display menu, use this for batch jobs\n");
		  printf("-menu     # forces the display of the menu\n");
		  printf("\n");
		  exit(-1);
		  //break;
		}
#ifdef BEAGLE
	      if(strcmp(arguments[argument],"-gpu")==0)
		{
		  use_beagle_gpu = 1;
		}
	      if(strcmp(arguments[argument],"-manualscale")==0)
		{
		  use_beagle_manualscale = 1;
		}
	      if(strcmp(arguments[argument],"-dynamicscale")==0)
		{
		  use_beagle_dynamicscale = 1;
		}
	      if(strcmp(arguments[argument],"-autoscale")==0)
		{
		  use_beagle_autoscale = 1;
		}

#endif
	    }
	  argument++;
	}
    }
  return usemenu;
}

/// set the menu overrider 
boolean set_usemenu(boolean usemenu, boolean fromparmfile)
{
  switch(usemenu)
    {
    case OVERRIDE_MENU:
      return TRUE;
      //break;
    case OVERRIDE_NOMENU:
      return FALSE;
      //break;
    default:
      break;
    }
  return fromparmfile;
} 


/// checks the settings of the number of long an short chain for bayes options and resets useless settings
void check_bayes_options(option_fmt *options)
{
  if (options->bayes_infer)
    {
      if(options->schains>0)
	{
	  //		fprintf(stdout,"NOTICE: with Bayesian inference no short chains are allowed, setting them to ZERO\n");
	  options->schains = 0;
	}
      if(options->lchains>1 || options->lchains < 1)
	{
	  fprintf(stdout,"NOTICE: with Bayesian inference, most efficient runs are with 1 long chain,\nsetting the long chain to 1\n");
	  options->lchains = 1;
	};
    }
}

void reset_bayesmdimfile(world_fmt *world, option_fmt *options)
{
  long locus=0;
  char ** files = NULL;
  long numfiles=0;
  if(world->options->has_bayesmdimfile)
  {
      charvec2d(&files,numcpu,STRSIZE);
      get_bayeshist (world, options); // this is transferring some material but
      // the parameter file is not filled in anymore [-> small data block]
      // save parameters in file and reload
      if(myID==MASTER)
      {
          // first file is alwasy present
          strcpy(files[numfiles++],options->bayesmdimfilename);
#ifdef MPI
#ifdef PARALIO
          fprintf(stdout,"%i> bayesmdimfile opening for analysis\n", myID);
          //this will need to open a whole array of files
          long numfiles=0;
          charvec2d(&files,numcpu,STRSIZE);
          long worker;
          char tmp[STRSIZE];
          MPI_Status status;
          strcpy(files[numfiles++],options->bayesmdimfilename);
          long numelem  = world->numpop2 + (world->options->gamma ? 1 : 0);
          long numelem2 = numelem * 2;
          MYREAL * temp = (MYREAL *) mycalloc(numelem2+2,sizeof(MYREAL));
          temp[0]=MIGMPI_PARALIO;
          for(worker=1;worker<numcpu;worker++)
          {
              MYMPISEND (temp, numelem2 + 2, mpisizeof, worker, worker, comm_world);
          }
          worker=1;
          while(worker<numcpu)
          {
              MYMPIRECV (tmp,STRSIZE, MPI_CHAR, MPI_ANY_SOURCE,
                         MPI_ANY_TAG, comm_world, &status);
              strcpy(files[numfiles++],tmp);
              worker++;
          }
          myfree(temp);
          //	  world->bayesmdimfile = NULL;//fopen(options->bayesmdimfilename,"r+");
#endif /*paralio*/
#endif /*mpi*/
      } /*for the master to get the filenames for bayesallfile */
#ifdef ZNZ
#ifdef PARALIO
      // closing the MPI worker's files and the master
      MPI_File_close(&world->mpi_bayesmdimfile);
#else
      //closing the one and only file (single cpu mode or MPI master alone)
      if(myID==MASTER)
	znzclose(world->bayesmdimfile);
#endif /*paralio*/
      // open the bayesallfile assuming no paralio setting and using the zipped stuff ZNZ
      if(myID==MASTER)
	world->bayesmdimfile = znzopen(options->bayesmdimfilename,"r", options->use_compressed);
#else
      // not ZNZ
#ifdef PARALIO
      // closing the MPI worker's files and the master
      MPI_File_close(&world->mpi_bayesmdimfile);
#else
      //closing the one and only file (single cpu mode or MPI master alone)
      if(myID==MASTER)
	fclose(world->bayesmdimfile);
#endif /*paralio*/
      world->bayesmdimfile = NULL;
      //world->bayesmdimfile = fopen(options->bayesmdimfilename,"r");
      //if(world->bayesmdimfile==NULL)
      //  {
      //    perror("Failed to open file: ");
      //    printf("errno = %d.\n", errno);
      //    exit(1);
      //  }
#endif /*znz*/
      if (myID==MASTER)
      {
	  read_bayes_fromfile(world->bayesmdimfile, world, options,files, numfiles);
          for(locus=0;locus<world->loci;locus++)
          {
              if(world->data->skiploci[locus])
                  continue;
              if(world->bayes->histogram[locus].results == NULL)
              {
                  world->data->skiploci[locus] = TRUE;
                  continue;
              }
              calc_hpd_credibility(world, locus, world->numpop2, world->numpop2 + world->bayes->mu+world->species_model_size*2 + world->grownum);
          }
      }
  }
  else
    {
      // keep everything in RAM
      get_bayeshist (world, options);
    }
}
