/* ! \file menu.c */
/*------------------------------------------------------
  Maximum likelihood estimation
  of migration rate  and effectice population size
  using a Metropolis-Hastings Monte Carlo algorithm
  -------------------------------------------------------
  M E N U   R O U T I N E S

  presents the menu and its submenus.
  Peter Beerli 1996, Seattle
  beerli@fsu.edu

  Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
  Copyright 2003-2007 Peter Beerli, Tallahassee FL

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


  $Id: menu.c 2144 2013-03-12 03:14:13Z beerli $

  -------------------------------------------------------*/
#include "migration.h"
#include "sighandler.h"
#include "tools.h"
#include "sequence.h"

#include "options.h"
#include "priors.h"
#include "migrate_mpi.h"
#include "menu.h"
#include "pretty.h"
#include <string.h>
#include <unistd.h>
#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif
/* prototypes ------------------------------------------- */
//void            print_menu_title(FILE * file, option_fmt * options, world_fmt *world);
//void            print_menu_accratio(long a, long b, world_fmt * world);
//long            print_title(world_fmt * world, option_fmt * options);
//void            get_menu(option_fmt * options);
/* private functions */
void            setup_datatype(char *datatype, option_fmt * options);
void            menuData(option_fmt * options, char datatype[]);
void            menuInput(option_fmt * options);
void            menuParameters(option_fmt * options);
void            menuStrategy(option_fmt * options);
void            menuHeat(option_fmt * options, char *input);
void            display_ml_mcmc(option_fmt * options);
void            display_bayes_mcmc(option_fmt * options);
//void          menuSequences(option_fmt * options);
void            read_custom_menu_migration(option_fmt * options);
void            read_custom_menu_lratio(option_fmt * options);
char           *custom_migration_type(long type);
void            read_heatvalues(option_fmt * options);
void            menuRandom(MYREAL * param, char type);
//void            start_tree_method(option_fmt * options);
void            start_data_method(option_fmt * options);
void menu_msat_submodel(option_fmt *options);
void            how_many_pop(long *numpop);
void            set_menu_localities(option_fmt *options);
void            set_localities_string(char *loc, option_fmt *options);
char            dialog_lrt(long numpop, lratio_fmt * lratio);
char            menuread_lrt_paramvalue(char *value, long *counter, long numpop2);
void            get_plotmenu(option_fmt * options);
boolean         menuStrategy_ml(option_fmt * options);
boolean         menuStrategy_bayes(option_fmt * options);
long            get_prior(char *input);
void            set_prior(char *output, int * prior, boolean without_rate);
void set_proposal(char *output, boolean *proposal, boolean without_rate);


void menu_sequence_submodel(option_fmt *options);
void menu_haplotyping(option_fmt *options);
char * menu_sequence_submodeltype(int type);
char * msat_submodeltype(int type);


void menu_get_filename(char message[], char thedefault[], char *filename);
boolean         setup_categs(option_fmt * options);
  boolean         setup_rcategs(option_fmt * options);
void change_tipdate(char *input, option_fmt * options);
void change_mutationrate(char *input, option_fmt * options);
void change_generationtime(char *input, option_fmt * options);
void change_randomsubset(char *input, option_fmt * options);
void change_assignment(char *input, option_fmt * options);
void change_inheritance(char *input, option_fmt * options);
void current_datatype_text(char *text, option_fmt *options);
void display_inheritance_option(option_fmt *options, const int inheritance);
void display_random_subset_option(option_fmt *options, const int randomsubset);
void display_assignment_option(option_fmt *options, const int assignsubset);
void display_sampledate_option(char * text, option_fmt * options, 
			       const int tipdate, 
			       const int mutationrate,
			       const int generationtime);
void setup_starttree(char *starttree, option_fmt * options);
void set_menu_ownstartvalue(option_fmt *options, long numpop);
void display_allele_options(char *text, char *starttree, option_fmt *options);
void display_brownian_options(char *text, char *starttree, option_fmt *options);
void display_stepwise_options(char *text, char *starttree, option_fmt *options);
void display_seq_mutationmodel(char *text, char *starttree, option_fmt *options);
void set_menu_priorpercent(option_fmt *options, float value);
void set_menu_randomprior(option_fmt *options);
void            print_bottom_menu_part(void);
void set_mult_prior(int paramgroupm, prior_fmt *prior);
void set_exp_prior(int paramgroupm, prior_fmt *prior);
void set_gamma_prior(int paramgroupm, prior_fmt *prior);
void set_wexp_prior(int paramgroupm, prior_fmt *prior);
void set_uni_prior(int paramgroupm, prior_fmt *prior);
void set_binning(prior_fmt *prior);
char * is_prior(int priorkind, prior_fmt *p, int priorset);
boolean set_proposal_menu(int paramgroup, option_fmt *options);
prior_fmt * set_theta_priormenu(option_fmt *options);
prior_fmt * set_mig_priormenu(option_fmt *options);
prior_fmt * set_split_priormenu(option_fmt *options, int type);
boolean set_prior_menu(int paramgroup, option_fmt *options);
void set_autotuning(option_fmt *options);
void menuProposal(option_fmt * options);
void menuPrior(option_fmt * options);
void display_JC69(option_fmt * options);
void  display_K2P(option_fmt * options);
void  display_F81(option_fmt * options);
void  display_F84(option_fmt * options);
void  display_HKY(option_fmt * options);
void   display_TN(option_fmt * options);
void sequence_modelparameters(int type, char * text, option_fmt *options);

long  is_priortype(prior_fmt *p, long pnum, int priortype);

extern time_t startseconds;
//##

#ifdef POPMODEL
#include "popmodel.h"
#endif
#define MI_INFILE 1
#define MI_RAND 2
#define MI_TITLE  3
#define MI_SUMREAD 4
#define MI_PROGRESS 5
#define MI_PRINTDATA 6
#define MI_OUTFILE 7
#define MI_TERSEPDF 8
#define MI_TREES 12
//#define MI_SUMFILE 14
#define MI_LOGFILE 15
#ifdef UEP
#define MI_UEPFILE 16
#define MI_UEPRATE 17
#define MI_UEPFREQ 18
#endif
#define MI_MIGHISTOGRAM 19
#define MI_SKYLINE 20

///
///Enumeration of possible population models
enum popmodel_menu 
{
  WRIGHTFISHER, CANNING, MORAN
};

///
///Enumeration of possible choices for ML menu
enum ml_menu 
{
  MLSTRATEGY, MLSHORTCHAINS, MLSHORTSKIP, MLSHORTSAMPLES, MLLONGCHAINS, MLLONGSKIP, MLLONGSAMPLES,
  MLBURNIN, MLREPLICATE, MLHEAT, MLMOVINGSTEPS, MLEPSILON, MLGELMAN
};

///
///Enumeration of possible choices for Datatype menu
enum dtseq_menu 
{
  DTSEQZERO, DTSEQTYPE, DTSEQHAP, DTSEQTRATIO, DTSEQFREQ, DTSEQSITECATEGS, DTSEQRATES, DTSEQCORR, DTSEQWEIGHT, DTSEQINTERLEAVED, DTSEQERROR, DTSEQFAST, DTSEQTREE, DTSEQINHERITANCE, DTSEQRANDOMSUBSET, DTSEQASSIGN, DTSEQTIPDATE, DTSEQMUTRATE, DTSEQGENERATION 
};

enum dtmsat_menu {
  DTMSATZERO, DTMSATTYPE, DTMSATXXX, DTMSATMISS, DTMSATTRESH, DTMSATTREE, DTMSATINHERITANCE, DTMSATRANDOMSUBSET, DTMSATASSIGN, DTMSATTIPDATE, DTMSATMUTRATE, DTMSATGENERATION
};

enum dtbrown_menu 
{
  DTBROWNZERO, DTBROWNTYPE, DTBROWNXXX, DTBROWNMISS, DTBROWNTREE, DTBROWNINHERITANCE, DTBROWNRANDOMSUBSET, DTBROWNASSIGN, DTBROWNTIPDATE, DTBROWNMUTRATE, DTBROWNGENERATION
};
enum dtep_menu {
  DTEPZERO, DTEPTYPE, DTXXXX, DTEPMISS, DTEPTREE, DTEPINHERITANCE, DTEPRANDOMSUBSET, DTEPASSIGN, DTEPTIPDATE, DTEPMUTRATE, DTEPGENERATION
};

///
///Enumeration of possible choices for Bayes menu
enum bayesian_menu {
  BAYESSTRATEGY, BAYESOUT, BAYESMDIMOUT, BAYESBINNING, BAYESPRETTY, BAYESFREQ, BAYESPROPOSAL, BAYESPRIOR, BAYESLCHAINS, BAYESSKIP, BAYESSAMPLES, BAYESBURNIN,
  BAYESREPLICATE, BAYESHEAT, BAYESMOVINGSTEPS, BAYESGELMAN, BAYESPRIORALONE 
};



///
///get user supplied filenames or fill the filenames with default values
void menu_get_filename(char message[], char thedefault[], char *filename)
{
  char            input[LINESIZE];
  printf("%s\n[Default: %s]\n===> ", message, thedefault);
  fflush(stdout); FGETS(input, LINESIZE, stdin);
  if (input[0] == '\0')
    strcpy(filename, thedefault);
  else
    strcpy(filename, input);
}

long print_menu_title(FILE * file, option_fmt * options, world_fmt *world)
{
  long filepos = -1;
  char nowstr[LINESIZE];
  if (!(options->menu || options->progress))
    return -1;

  //if(options->menu && file!=world->outfile)
  //  clearscreen();
  FPRINTF(file, " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  FPRINTF(file, " +                                                                +\n");
  FPRINTF(file, " +   POPULATION SIZE, MIGRATION, DIVERGENCE, ASSIGNMENT, HISTORY  +\n");
  FPRINTF(file, " +   Bayesian inference using the structured coalescent           +\n");
  FPRINTF(file, " +                                                                +\n");
  FPRINTF(file, " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
#ifdef MPI
  FPRINTF(file, "  Compiled for a PARALLEL COMPUTER ARCHITECTURE\n");
  FPRINTF(file, "  One master and %i compute nodes are available.\n", numcpu-1);
#endif
#ifdef OSTYPE
  FPRINTF(file, "  Operating System: %s\n", OSTYPE);
#endif
#ifdef AVX
  FPRINTF(file, "  Using Intel AVX (Advanced Vector Extensions)\n");
#endif
#ifdef THREAD
  FPRINTF(file, "  Compiled for a SYMMETRIC multiprocessors\n");
#endif
#ifdef GRANDCENTRAL
  FPRINTF(file, "  Compiled for a SYMMETRIC multiprocessors (GrandCentral)\n");
#endif
#ifdef FAST_EXP
  FPRINTF(file, "  Fast approximation to Exp() and Log() used\n");
#endif
#ifdef PRETTY
  FPRINTF(file, "  PDF output enabled [%s]\n",
#ifdef A4PAPER
	  "A4-size"
#else
	  "Letter-size"
#endif
	  );
#endif  
  FPRINTF(file, "  Version %s ", MIGRATEVERSION);
  if(strlen(MIGRATESUBVERSION)>0)
    { 
      FPRINTF(file, "  [%s]\n", MIGRATESUBVERSION);
    }
  else
    FPRINTF(file, "\n");
  get_time(nowstr, "  %c");
  startseconds = time(0);
  if (nowstr[0] != '\0')
    FPRINTF(file, "  Program started at %s\n", nowstr);
  if(file == world->outfile)
    {
      filepos = ftell(world->outfile);
      /*
       * DO NOT REMOVE the blanks in the following print
       * statement, they are overwritten at the program end
       * with the timestamp and the runtime
       */
      FPRINTF(world->outfile,
	      "                                                                        \n");
    }
  FPRINTF(file, "\n\n");
  return filepos;
}

void
print_menu_accratio(long a, long b, world_fmt * world)
{
  char            buffer[STRSIZE];
  boolean         writelog = world->options->writelog;
  boolean         progress = world->options->progress;

  if (writelog || progress) 
    {
      sprintf(buffer, "           Acceptance-ratio = %li/%li (%f)\n\n", a, b,
	      (MYREAL) a / (MYREAL) b);
      if (progress)
	FPRINTF(stdout, "%s", buffer);
      if (writelog)
	FPRINTF(world->options->logfile, "%s", buffer);
    }
}

long print_title(world_fmt * world, option_fmt * options)
{
  long            len = 45, i, filepos = -1;
  if (!world->options->simulation) 
    {
      FPRINTF(world->outfile, " ");
      if (options->title[0] != '\0') 
	{
	  len = (long) MAX(strlen(options->title), 66);
	  for (i = 0; i < len; i++)
	    FPRINTF(world->outfile, "+");
	  FPRINTF(world->outfile, "\n%-*s\n", (int) len,options->title);
	}
      filepos = print_menu_title(world->outfile,options,world);
    }
  return filepos;
}

///
/// main menu display, if counter is bigger than 100 the program will exit, to 
/// make sure that batch jobs are not filling up the available harddisk space
void
get_menu(option_fmt * options, world_fmt *world, data_fmt *data)
{
#ifdef MPI
  	long bufsize = 9;
	char * menubuffer = (char *) mycalloc (bufsize, sizeof (char));
#endif
	
  static boolean first=TRUE;
  long           counter=0L;
  char            input[LINESIZE];
  char           *datatype;
  datatype = (char *) mycalloc(LINESIZE, sizeof(char));
  if (options->menu) {
    setup_datatype(datatype, options);
    do {
      if (first==TRUE)
	{
	  first = FALSE;
	}
      else
	{
	  if(options->menu)
	    clearscreen();
	}
      printf("  =============================================================\n");
      printf("  MAIN MENU\n");
      printf("  =============================================================\n\n");
      printf("  D       Data type currently set to: %-30.30s\n", datatype);
      printf("  I       Input/Output formats and Event reporting\n");
      printf("  P       Parameters  [start, migration model]\n");
      printf("  S       Search strategy\n");
      printf("  W       Write a parmfile\n");
      printf("  Q       Quit the program\n");
      printf("\n\n");
      printf("  To change the settings type the letter for the menu to change\n");
      printf("  Start the program with typing Yes or Y\n===> ");fflush(stdout);
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      printf("\n");
      switch (uppercase(input[0])) {
      case 'D':
	if(options->menu)
	  clearscreen();
	menuData(options, datatype);
	break;
      case 'I':
	if(options->menu)
	  clearscreen();
	menuInput(options);
	break;
      case 'P':
	if(options->menu)
	  clearscreen();
	menuParameters(options);
	break;
      case 'S':
	if(options->menu)
	  clearscreen();
	menuStrategy(options);
	break;
      case 'W':
	save_parmfile(options, world, data);
	break;
      case 'Q':
#ifdef MPI
        //MPI_Abort(comm_world,0);
	sprintf(menubuffer,"MENUSTOP");
	MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
	MYMPIBCAST (menubuffer, bufsize, MPI_CHAR, MASTER, comm_world);
	myfree(menubuffer);
	MPI_Finalize();
#endif /*MPI*/
	exit(0);
      default:
	break;
      }
      if(counter++ > 100)
	{
	  warning("Batch job uses menu ON option, and gets non-sensical input, set menu=NO in parmfile\nProgram continues but check!!!!!");
	  input[0]='Y';
	  break;
	}
    }
    while (uppercase(input[0]) != 'Y');
    if (options->usertree && !strchr(SEQUENCETYPES, options->datatype))
      options->usertree = FALSE;
    if (options->dist && !strchr(SEQUENCETYPES, options->datatype))
      options->dist = FALSE;
    //prevents incompatability
    // of tree and distance designation, only for
    //sequence data the datalines and tree - tips can match
    // msats or EP data is ambiguos concerning the tree.
      }
  if (options->datatype == 'g')
    //prevents overwriting sumfile
    options->writesum = FALSE;

#ifdef UEP

  if (options->uep)
    options->randomtree = TRUE;
#endif
  myfree(datatype);
#ifdef MPI
  myfree(menubuffer);
#endif
}
/*
 * private
 * functions---------------------------------------------------
 */

boolean         setup_categs(option_fmt * options) {
  boolean         didchangecat = FALSE;
  if              (options->categs > 1)
    didchangecat = TRUE;
  else
    didchangecat = FALSE;
  return didchangecat;
}

boolean         setup_rcategs(option_fmt * options) {
  boolean         didchangercat = FALSE;
  if              (options->rcategs > 1)
    didchangercat = TRUE;
  else
    didchangercat = FALSE;

  return didchangercat;
}

void            setup_datatype(char *datatype, option_fmt * options) 
{
  switch (options->datatype) 
    {
    case 'a':
      strcpy(datatype, "Allele model");
      break;
    case 'm':
      if(options->msat_option == SINGLESTEP)
	strcpy(datatype, "Microsatellite data");
      else
	strcpy(datatype, "Microsatellite data");
      break;
    case 'b':
      strcpy(datatype, "Microsatellite data");
      break;
    case 's':
      strcpy(datatype, "DNA sequence model");
      break;
    case 'n':
      strcpy(datatype, "SNP model");
      break;
    case 'h':
      strcpy(datatype, "SNP model (Hapmap data)");
      break;
      //case 'u':
      //strcpy(datatype, "Panel-SNP (panel is population 1)");
      //break;
      //case 'f':
      //strcpy(datatype, "ancestral state reconstruction method");
      //break;
      //case 'g':
      //strcpy(datatype, "Genealogy summary");
      //options->readsum = TRUE;
      //break;
    default:
      //if (options->readsum) {
      //options->datatype = 'g';
      //strcpy(datatype, "Genealogy summary");
      //} else {
      options->datatype = 's';
      strcpy(datatype, "DNA sequence model");
	//}
      break;
    }
}

///
/// set up start tree
void setup_starttree(char *starttree, option_fmt * options) 
{
  if (options->usertree && strchr(SEQUENCETYPES, options->datatype))
    sprintf(starttree, "is supplied in %s", options->utreefilename);
  else 
    {
      if (options->dist && strchr(SEQUENCETYPES, options->datatype))
	sprintf(starttree, "generates using %s", options->distfilename);
      else 
	{
	  if (options->randomtree)
	    strcpy(starttree, "random start genealogy");
	  else
	    strcpy(starttree, "UPGMA based start genealogy");
	}
    }
}


///
/// checks whether the data type can use a start tree or not
/// and also allows to change the datatype
//void            check_type_tree(char *input, option_fmt * options) {
// switch (options->datatype) 
//  {
//  case 'a':
//    switch (atol(input)) {
//    case DTEPTREE:
//start_tree_method(options);
//break;    
//    case DTEPTYPE:
//      start_data_method(options);
//      break;
//    default:
//break;
//    }
//    break;
//  case 'm':
//    switch (atol(input)) 
//{
//case DTMSATTREE:
//  start_tree_method(options);
//  break;
////case DTMSATTYPE:
//  start_data_method(options);
//  break;        
//default:
//  break;
//}
//    break;
//  case 'b':
//    switch (atol(input)) 
//{
//case DTBROWNTREE:
//  start_tree_method(options);
//  break;
//case DTBROWNTYPE:
//  start_data_method(options);
//  break;
//default:
//  break;
//}
///    break;
// case 's':
//  case 'n':
//  case 'h':
//  case '@':
//    switch (atol(input)) 
//{
//case DTSEQTREE:
//  start_tree_method(options);
//  break;
//case DTSEQTYPE:
//  start_data_method(options);
//  break;
///default:
//  break;
//}
//    break;
//    //case 'u':
//    //case 'f':
//  case 'g':
//    switch (atol(input)) {
//    case 1:
//start_data_method(options);
//break;
////    default:
//break;
//    }
//    break;
////default:
//break;
//}
//}


///
/// changing tipdate treatment and filename
void change_tipdate(char *input, option_fmt * options)
{
	printf(" Samples (Tips) where sampled at different dates? [No]\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	if (uppercase(input[0]) == 'Y') 
	  {
	    options->has_datefile = TRUE;
	    menu_get_filename(" What is the filename that contains the sample dates?", TIPDATEFILE, options->datefilename);
	  } 
	else 
	  {
	    options->has_datefile = FALSE;
	  }
	input[0] = 'X';
}

///
/// changing mutation rate associated with tipdate
void change_mutationrate(char *input, option_fmt * options)
{
  char *in;
  char *tmp;
  long count=0;
  printf(" Give the number of loci or\n when all loci have the same mutation rate, the rate\nor a C to cancel\n===> ");
  fflush(stdout); FGETS(input, LINESIZE, stdin);
  if (uppercase(input[0]) == 'C')
    return;
  if (input[0] == '0') 
    {
      options->mutationrate_year_numalloc = 1;
      options->mutationrate_year[0] = atof(input);
    } 
  else 
    {
      options->mutationrate_year_numalloc = atol(input);
      options->mutationrate_year = (MYREAL *) myrealloc(options->mutationrate_year,
							sizeof(MYREAL)* (size_t) options->mutationrate_year_numalloc);
      while(count < options->mutationrate_year_numalloc)
	{
	  printf(" Give the mutation rate per year for each locus\n[either in groups or singly use a space to separate]\n===> ");
	  fflush(stdout); FGETS(input, LINESIZE, stdin);
	  in = input;
	  tmp = strsep(&in,";, ");
	  while(tmp != NULL)
	    {
	      options->mutationrate_year[count++] = atof(tmp);
	      if(in==NULL)
		break;
	      tmp = strsep(&in,";, ");
	    }
	}
    }
  input[0] = 'X';
}

///
/// changing generation time associated with tipdates
void change_generationtime(char *input, option_fmt * options)
{
  printf(" Give the number of generation per year\nor a C to cancel\n===> ");
  fflush(stdout); FGETS(input, LINESIZE, stdin);
  if (uppercase(input[0]) == 'C')
    return;
  options->generation_year = atof(input);
  input[0] = 'X';
}

///
/// changing hwo many individuals are used for the run, subsetting the data randomly
void change_randomsubset(char *input, option_fmt * options)
{
  printf(" How may individuals should be used\nGive a number or A for all or a C to cancel\n===> ");
  fflush(stdout); FGETS(input, LINESIZE, stdin);
  if (uppercase(input[0]) == 'C')
    return;
  if (uppercase(input[0]) == 'A')
    {
      options->randomsubset = 0;
      options->randomsubsetseed = 0;
    }
  else
    {
      options->randomsubset = (long) atol(input);
      printf(" If you want to generate a subset that be regenerated\nin another run of migrate you need to give\na random number seed (a number between 1 and 2147483648),\nif you give the same numbers in other runs migrate will pick the same set\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      options->randomsubsetseed= (unsigned long) atol(input);
    }
  input[0] = 'X';
}

///
/// turn on and off assignment, the data must contain a ? as the first
/// character in the individual name
void change_assignment(char *input, option_fmt * options)
{
  if (options->has_unassigned)
    options->has_unassigned = FALSE;
  else
    {
      options->has_unassigned = TRUE;
      
    }
  printf(" Assignment of individuals to populations: Use {F)requency or (U)niform\n===> ");
  fflush(stdout); FGETS(input, LINESIZE, stdin);
  if (uppercase(input[0]) == 'F')
    options->has_unassignedfreq=TRUE;
  else
    options->has_unassignedfreq=FALSE;
  input[0] = 'X';
}

///
/// changing inheritance scalar for each locus
void change_inheritance(char *input, option_fmt * options)
{
  char *in;
  char *tmp;
  long count=0;
  printf(" Inheritance scalars: Give the number of loci or\nor a C to cancel\n===> ");
  fflush(stdout); FGETS(input, LINESIZE, stdin);
  if (uppercase(input[0]) == 'C')
    return;
  
  options->inheritance_scalars_numalloc = atol(input);
  options->inheritance_scalars = (MYREAL *) myrealloc(options->inheritance_scalars,
						      sizeof(MYREAL) * (size_t) options->inheritance_scalars_numalloc);
  while(count < options->inheritance_scalars_numalloc)
    {
      printf(" Give the inheritance scaler for each locus\n[either in groups or singly use a space to separate]\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      in = input;
      tmp = strsep(&in,";, ");
      while(tmp != NULL)
	{
	  options->inheritance_scalars[count++] = atof(tmp);
	  if(in==NULL)
	    break;
	  tmp = strsep(&in,";, ");
	}
    }
  input[0] = 'X';
}

///
/// print current datatype into menu
void current_datatype_text(char *text, option_fmt *options)
{
  text[0]='\0';
  switch(options->datatype)
    {
    case 'a':
      sprintf(text,"'Infinite' allele model");
      break;
    case 'm':
      sprintf(text, "%22.22s", msat_submodeltype(options->msat_option));
      break;
    case 'b':
      sprintf(text, "%22.22s", "Brownian motion model");
      break;
    case 'u':
    case 'n':
    case 'h':
    case 's':
    case 'f':
      sprintf(text, "%22.22s", menu_sequence_submodeltype(options->sequence_model));
      break;
    case 'g':
    default:
      sprintf(text,  "Model not specified"); 
      break;
    }
}

///
/// displays the inheritance menu
void display_inheritance_option(option_fmt *options, const int inheritance)
{
      printf(" %2i   Inheritance scalar set                                      %3.3s\n",inheritance, 
	     options->inheritance_scalars_numalloc>1 ? "YES" : " NO");
}

///
/// displays the random subset menu
void display_random_subset_option(option_fmt *options, const int randomsubset)
{
      if(options->randomsubset > 0)
	printf(" %2i   Pick random subset per population of individuals             %4li\n", randomsubset,  options->randomsubset);
      else
	printf(" %2i   Pick random subset per population of individuals             NO\n", randomsubset);
}
///
/// displays the random subset menu
void display_assignment_option(option_fmt *options, const int assignsubset)
{
      if(options->has_unassigned)
	printf(" %2i   Assign individuals [?name in data] to population            YES\n",assignsubset);
      else
	printf(" %2i   Assign individuals [?name in data] to population             NO\n",assignsubset);
}

///
/// shows all the option for the dampling dates,
/// user needs to give a datefile and also generation time and year
void display_sampledate_option(char * text, option_fmt * options, 
			       const int tipdate, 
			       const int mutationrate,
			       const int generationtime)
{
  if(options->has_datefile)
    {
      printf(" %2i   Tip date file %49.49s\n", tipdate, options->datefilename);
      if(options->mutationrate_year_numalloc > 1)
	sprintf(text,"multiple rates");
      else
	sprintf(text,"%.12f",options->mutationrate_year[0]);
      printf(" %2i   Mutation rate per locus and year %30.30s\n", mutationrate, text);
      sprintf(text,"%10.4f",options->generation_year);
      printf(" %2i   How many generations per year  %30.30s\n", generationtime, text);
    }
  else
    {
      printf(" %2i   Tip date file                     None, all tips a contemporary\n", tipdate);
    }
}

///
/// displays the Electrophoretic marker data menu
void display_allele_options(char *text, char *starttree, option_fmt *options)
{
  (void) starttree;
  printf(" %2i   Discard missing data      %37.37s\n", DTEPMISS, options->include_unknown ? " NO" : "YES");
  //printf(" %2i   Start genealogy            %37.37s\n", DTEPTREE, starttree);
  display_inheritance_option(options, DTEPINHERITANCE);
  display_random_subset_option(options,DTEPRANDOMSUBSET);
  display_assignment_option(options,DTEPASSIGN);
  display_sampledate_option(text, options, DTEPTIPDATE, DTEPMUTRATE, DTEPGENERATION);
}

///
/// display the options for the brownnian motion model
void display_brownian_options(char *text, char *starttree, option_fmt *options)
{
  (void) starttree;
  printf(" %2i   Discard missing data:     %37.37s\n", DTBROWNMISS, options->include_unknown ? " NO" : "YES");
  //  printf(" %2i   Start genealogy            %37.37s\n", DTBROWNTREE, starttree);
  display_inheritance_option(options, DTBROWNINHERITANCE);
  display_random_subset_option(options,DTBROWNRANDOMSUBSET);
  display_assignment_option(options,DTBROWNASSIGN);
  display_sampledate_option(text, options, DTBROWNTIPDATE, DTBROWNMUTRATE, DTBROWNGENERATION);
}  

///
/// display the options for the stepwise mutation model
void display_stepwise_options(char *text, char *starttree, option_fmt *options)
{
  (void) starttree;
  printf(" %2i   Discard missing data:     %37.37s\n", DTMSATMISS, options->include_unknown ? " NO" : "YES");
  //  printf(" %2i   Start genealogy            %37.37s\n", DTMSATTREE, starttree);
  display_inheritance_option(options, DTMSATINHERITANCE);
  display_random_subset_option(options,DTMSATRANDOMSUBSET);
  display_assignment_option(options,DTMSATASSIGN);
  display_sampledate_option(text, options, DTMSATTIPDATE, DTMSATMUTRATE, DTMSATGENERATION);
}  

///
/// displayes the sequence options
void display_seq_mutationmodel(char *text, char *starttree, option_fmt *options)
{
  (void) starttree;
  //long z = 0;
  //  long numchar = 0;      
  printf(" %2i   One category of sites?           %30.30s\n", DTSEQSITECATEGS, options->categs == ONECATEG ?
	 "One category" :
	 "More than one category of sites");
  printf(" %2i   One region of substitution rates?", DTSEQRATES);
  if (options->rcategs == 1)
    printf("                           YES\n");
  else {
    printf("  %4li categories of regions\n", options->rcategs);
    printf(" %2i   Rates at adjacent sites correlated?", DTSEQCORR);
    if (!options->autocorr)
      printf("    NO, they are independent\n");
    else
      printf("YES, mean block length =%6.1f\n",
	     1.0 / options->lambda);
  }
  printf(" %2i   Sites weighted?            %36.36s\n", DTSEQWEIGHT,
	 (options->weights ? "YES" : "NO"));
  //printf(" %2i   Input sequences interleaved? %34.34s\n", DTSEQINTERLEAVED,
  //	 (options->interleaved ? "YES" : "NO, sequential"));
  if(options->has_estimateseqerror)
    printf(" %2i   Sequencing error rate?                             ESTIMATE\n", DTSEQERROR);
  else
    printf(" %2i   Sequencing error rate?                [%4.3f %4.3f %4.3f %4.3f]\n", DTSEQERROR,
	   options->seqerror[0], options->seqerror[1],options->seqerror[2],options->seqerror[3]);
  //  printf(" %2i   Slow but safer Data likelihood calculation %20.20s\n", DTSEQFAST,
  //	 options->fastlike ? "NO" : "YES");
  //  printf(" %2i   Start genealogy:          %37.37s\n", DTSEQTREE, starttree);
  
  
  //      while (options->ttratio[z] > 0.0)
  //numchar += sprintf(text + numchar, "%8.4f ", options->ttratio[z++]);
  //
  //
  //printf(" %2i   Transition/transversion ratio:   %31s\n", DTSEQTRATIO, text);
  //printf(" %2i   Use empirical base frequencies?  %30s\n", DTSEQFREQ, (options->freqsfrom ? "YES" : "NO"));
  //printf(" %2i   Start genealogy           %37.37s\n", DTSEQTREE, starttree);
  
  display_inheritance_option(options, DTSEQINHERITANCE);
  display_random_subset_option(options,DTSEQRANDOMSUBSET);
  display_assignment_option(options,DTSEQASSIGN);
  display_sampledate_option(text, options, DTSEQTIPDATE, DTSEQMUTRATE, DTSEQGENERATION); 
}

///
///let the user set all data related options
void
menuData(option_fmt * options, char datatype[]) {
  static boolean  didchangecat, didchangercat;

  long            i;
  //long            z = 0;
  //long            numchar = 0;

  char            input[LINESIZE];
  char            starttree[LINESIZE];

  char           *text = (char *) mycalloc(LINESIZE, sizeof(char));
  char           *text2 = (char *) mycalloc(LINESIZE, sizeof(char));

  didchangecat = setup_categs(options);
  didchangercat = setup_rcategs(options);

  do {
    setup_datatype(datatype, options);
    setup_starttree(starttree, options);
    current_datatype_text(text, options);
    printf("  ===================================================================\n");
    printf("  DATATYPE AND DATA SPECIFIC OPTIONS\n");
    printf("  ===================================================================\n\n");
    printf("  D   change Datatype, currently:%36.36s\n", datatype);
    printf(" %2i   change Mutation model, currently:%30.30s\n", DTSEQTYPE, text);
  if(options->haplotyping)
    sprintf(text, "YES:%s",options->haplotyping_report ? "reporting haplotypes" : "no reporting of haplotype" );
  else
    sprintf(text, "NO");

    printf(" %2i   Haplotyping is turned on:     %33.33s\n", DTSEQHAP,  text);
    switch (options->datatype) 
      {
      case 'a':
	display_allele_options(text, starttree, options);  // allelic data ----------------------------------------
	break;
      case 'b':
	display_brownian_options(text, starttree, options);// brownian motion data --------------------------------
	break;
      case 'm':
	display_stepwise_options(text, starttree, options);// stepwise mutation model  data -----------------------
	break;
      case 'u':
      case 'n':
      case 'h':
      case 's':
      case 'f':
	display_seq_mutationmodel(text, starttree, options);
	break;
	
      case 'g':
	printf("       [Reanalyze an old run]\n");
	break;
      }
    printf("\n\n");
    printf("  Are the settings correct?\n");
    printf("  (Type Y or the number of the entry to change)\n===> ");
    input[0] = '\0';
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if (strlen(input) == 0)
      continue;
    if (uppercase(input[0])=='D')
      {
	start_data_method(options);      
	continue;
      }
    //check_type_tree(input, options);
    //start tree option for all and switch datatype
    if (strchr("ab", options->datatype)) 
      // this only works if brownian and allele datatype have the same number of options
      {
	switch (atol(input)) 
	  {
	  case DTSEQTYPE:
	    menu_msat_submodel(options);
	    break;
	  case DTSEQHAP:
	    menu_haplotyping(options);
	    break;
	  case DTEPMISS:	/* discard unknowns */
	    printf("  Discard missing Alleles (YES|NO)?\n");
	    printf("  [Default is YES, this is the best choice for most situations]\n===> ");
	    fflush(stdout); FGETS(input, LINESIZE, stdin);
	    if (strchr("Yy", input[0]))
	      options->include_unknown = FALSE;
	    else
	      options->include_unknown = TRUE;
	    input[0] = '\0';
	    break;
	  case DTEPINHERITANCE:
	    change_inheritance(input,options);
	    break;
	  case DTEPRANDOMSUBSET:
	    change_randomsubset(input,options);
	    break;
	  case DTEPASSIGN:
	    change_assignment(input,options);
	    break;
	  case DTEPTIPDATE: 
	    change_tipdate(input,options);
	    break;
	  case DTEPMUTRATE:
	    change_mutationrate(input,options);
	    break;
	  case DTEPGENERATION:
	    change_generationtime(input,options);
	    break;
	  default:
	    break;
	  }
	continue;
      }
    if (options->datatype == 'm') 
      {
	switch (atol(input)) {
	case DTSEQTYPE:
	  menu_msat_submodel(options);
	  break;
	case DTMSATMISS:	/* discard unknowns */
	  printf("  Discard missing Alleles?\n");
	  printf("  [Default is YES, this is the best choice for most situations]\n===> ");
	  fflush(stdout); FGETS(input, LINESIZE, stdin);
	  if (strchr("Yy", input[0]))
	    options->include_unknown = TRUE;
	  else
	    options->include_unknown = FALSE;
	  break;
	case DTMSATTRESH:	/* micro-threshold */
	  printf("  What is the threshold value? [needs to be even!]\n");
	  printf
	    ("  E.g. if your allele is 24 and the threshold is 10\n");
	  printf("  there is some probability that the allele 24 can\n");
	  printf
	    ("  change to allele 14 (or 38), but there is a probability\n");
	  printf("  of 0.0 (ZERO) to go to 13 (39),\n");
	  printf
	    ("  if you choose this too small, than the program will fail\n");
	  printf
	    ("  The default is set 100, if the biggest difference in the data is smaller the value\n");
	  printf("  will be adjust to that maximal difference\n");
	  printf
	    ("  [the bigger the longer the run and the more\n accurate the estimate]\n===> ");
	  fflush(stdout); FGETS(input, LINESIZE, stdin);
	  options->micro_threshold = atol(input);
	  if(options->micro_threshold % 2 != 0)
	    options->micro_threshold += 1;
	  input[0] = '\0';
	  break;
	case DTMSATTIPDATE: 
	  change_tipdate(input,options);
	  break;
	case DTMSATMUTRATE:
	  change_mutationrate(input,options);
	  break;
	case DTMSATGENERATION:
	  change_generationtime(input,options);
	  break;
	  case DTMSATINHERITANCE:
	    change_inheritance(input,options);
	    break;
	  case DTMSATRANDOMSUBSET:
	    change_randomsubset(input,options);
	    break;
	  case DTMSATASSIGN:
	    change_assignment(input,options);
	    break;
	default:
	  break;
	}
	continue;
      }
    if (strchr(SEQUENCETYPES, options->datatype)) {
      switch (atol(input)) {
      case DTSEQTYPE:
	menu_sequence_submodel(options);
	break;
      case DTSEQHAP:
	menu_haplotyping(options);
	break;
      case DTSEQTRATIO:
	initratio(options);
	break;
      case DTSEQFREQ:
	options->freqsfrom = !options->freqsfrom;
	if (!options->freqsfrom) {
	  initfreqs(&options->freqa, &options->freqc,
		    &options->freqg, &options->freqt);
	}
	break;
      case DTSEQSITECATEGS:
	if (options->categs == ONECATEG) {
	  options->categs = MANYCATEGS;
	  didchangecat = TRUE;
	} else {
	  options->categs = ONECATEG;
	  didchangecat = FALSE;
	}
	break;
      case DTSEQRATES:
	if (didchangercat) {
	  options->autocorr = FALSE;
	  didchangercat = FALSE;
	  options->rcategs = ONECATEG;
	} else {
	  printf("\n  Regional rates:\n");
	  initcatn(&options->rcategs);
	  options->probcat =
	    (MYREAL *) myrealloc(options->probcat,
				 options->rcategs * sizeof(MYREAL));
	  options->rrate =
	    (MYREAL *) myrealloc(options->rrate,
				 options->rcategs * sizeof(MYREAL));
	  didchangercat = TRUE;
	  options->gammarates =
	    initcategs(options->rcategs, options->rrate,
		       options->probcat);
	  if (!options->gammarates)
	    {
	      initprobcat(options->rcategs, &options->probsum,
			  options->probcat);
	      constrain_rates(options->rcategs, options->rrate,
			      options->probcat);
	    }
	    FPRINTF(stdout,
		    "\n\nRegion type     Rate of change    Probability\n");
	  FPRINTF(stdout,
		  "---------------------------------------------\n");
	  for (i = 0; i < options->rcategs; i++)
	    FPRINTF(stdout, "%9ld%16.3f%17.3f\n",
		    i + 1, options->rrate[i], options->probcat[i]);
	  FPRINTF(stdout, "\n");
	}
	break;
      case DTSEQCORR:
	options->autocorr = !options->autocorr;
	if (options->autocorr)
	  initlambda(options);
	break;
      case DTSEQWEIGHT:
	options->weights = !options->weights;
	break;
      case DTSEQINTERLEAVED:
	options->interleaved = !options->interleaved;
	break;
      case DTSEQERROR:
	FPRINTF(stdout,
		"Enter the sequencing error rates for sites with A, C, G, T or E for ESTIMATE\n[Good values are 0.0 (=no error), or 0.01 (1%% error):\nEnter either one value (all nulceotides have the same error)\n or four values separated by spaces\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	if (input[0]=='E' || input[0]=='e')
	  {
	    options->has_estimateseqerror=TRUE;
	  }
	else
	  {
	    options->has_estimateseqerror=FALSE;
	    long vals = sscanf(input,"%lf%lf%lf%lf",
			       &options->seqerror[0],&options->seqerror[1],&options->seqerror[2],&options->seqerror[3]);
	    if (vals==1)
	      {
		options->seqerror[1] =options->seqerror[0];
		options->seqerror[2] =options->seqerror[0];
		options->seqerror[3] =options->seqerror[0];
	      }
	    else
	      {
		if (vals!=4)
		  warning("Sequencing error input incorrect");
	      }
	  }
	input[0] = '\0';
	break;
	//case DTSEQFAST:
	//options->fastlike = !options->fastlike;
	//break;
      case DTSEQTIPDATE: 
	change_tipdate(input,options);
	break;
      case DTSEQMUTRATE:
	change_mutationrate(input,options);
	break;
      case DTSEQGENERATION:
	change_generationtime(input,options);
	break;
      case DTSEQINHERITANCE:
	change_inheritance(input,options);
	break;
      case DTSEQRANDOMSUBSET:
	change_randomsubset(input,options);
	break;
      case DTSEQASSIGN:
	change_assignment(input,options);
	break;	
      default:
	break;
      }
      if (!didchangercat) {
	options->probcat =
	  (MYREAL *) myrealloc(options->probcat, sizeof(MYREAL) * 2);
	options->rrate =
	  (MYREAL *) myrealloc(options->rrate, sizeof(MYREAL) * 2);
	options->rrate[0] = 1.0;
	options->probcat[0] = 1.0;
      }
      if (!didchangecat) {
	options->rate =
	  (MYREAL *) myrealloc(options->rate, sizeof(MYREAL) * 2);
	options->rate[0] = 1.0;
      }
    }
    if (options->datatype != 'g')
      options->readsum = FALSE;
  }
  while (uppercase(input[0]) != 'Y');

  myfree(text);
  myfree(text2);

}

void
menuInput(option_fmt * options) {
  long            test = 0;
  //long            len = (long) strlen(options->title);
  unsigned long   timeseed;
  char            input[LINESIZE];
  char            treeinc[LINESIZE];
  char            outputstring[LINESIZE];
  char           *stringstep, *string2, *string3, *string4;
  char            treestr[4][20] = {"None", "All", "Best", "Last chain"};
  char            progress[3][8] = {"Verbose", "YES", "NO"};
  char           *extension;
  //  int             retval;
  stringstep = (char *) mymalloc(sizeof(char) * 128);
  string2 = (char *) mymalloc(sizeof(char) * 128);
  string3 = (char *) mymalloc(sizeof(char) * 128);
  string4 = (char *) mymalloc(sizeof(char) * 128);

  do {
    printf("  ====================================================================\n");
    printf("  INPUT/OUTPUT FORMATS (for %s)\n",options->bayes_infer ? "Bayesian inference" : "DONOTUSEMaximum likelihood approach");
    printf("  ====================================================================\n\n");
    printf("  Input:\n");
    printf("  %2i   Datafile name is %46.46s\n", MI_INFILE,
	   options->infilename);

    sprintf(treeinc,"[every %li]",options->treeinc);

    switch (options->autoseed) {
    case AUTO:
      sprintf(stringstep, "YES");
      break;
    case NOAUTO:
      sprintf(stringstep, "NO, use seedfile");
      break;
    case NOAUTOSELF:
      sprintf(stringstep, "NO, seed=%li ", options->inseed);
      break;
    default:
      options->autoseed = AUTO;
      sprintf(stringstep, "YES");
      break;
    }

    printf("  %2i   Use automatic seed for randomisation?  %24.24s\n",
	   MI_RAND, stringstep);
    //    printf("  %2i   Title of the analysis is", MI_TITLE);
    //
    //if (len == 0 || (len == 5 && strstr(options->title, "Migration analysis")))
    //  printf("%39.39s\n", "<no title given>");
    //else
    //  printf("%39.39s\n", options->title);

    if (options->readsum && !options->bayes_infer) {
      printf("  %2i    Summary of genealogies are read from %s\n",
	     MI_SUMREAD, options->sumfilename);
    }
    else
      {
	if (options->readsum && options->bayes_infer) {
	  printf("  %2i    Bayesian output data are read from %s\n",
		 MI_SUMREAD, options->bayesmdimfilename);
	}
      }
    printf("\n  OUTPUT:\n");

    printf("  %2i   Print indications of progress of run? %25.25s\n",
	   MI_PROGRESS, options->progress ? (options->
					     verbose ? progress[0] :
					     progress[1]) : progress[2]);
    printf("  %2i   Print the data?%48.48s\n",
	   MI_PRINTDATA, options->printdata ? "YES" : "NO");

    printf("  %2i   Outputfile name is  %43.43s\n", MI_OUTFILE,
	   options->outfilename);
#ifdef PRETTY
    printf("       PDF outfile name is %43.43s\n", options->pdfoutfilename);
    printf("  %2i   PDF outputfile is terse  %38.38s\n", MI_TERSEPDF,
	   options->tersepdf ? "YES" : "NO");
#endif

    switch (options->treeprint) 
      {
      case ALL:
	printf("  %2i   Print genealogies?     %40.40s\n",
	       MI_TREES, treestr[1]);
	break;
      case BEST:
	printf("  %2i   Print genealogies?     %40.40s\n",
	       MI_TREES, treestr[2]);
	break;
      case LASTCHAIN:
	printf("  %2i   Print genealogies?  %19.19s %20.20s\n",
	       MI_TREES, treeinc ,treestr[3]);
	break;
      case myNONE:
      default:
	printf("  %2i   Print genealogies?     %40.40s\n",
	       MI_TREES, treestr[0]);
	break;
      }

    if (options->writelog)
      printf("  %2i   Save logging information into %33.33s\n", MI_LOGFILE,
	     options->logfilename);
    else
      printf("  %2i   Save logging information?  %36.36s\n", MI_LOGFILE, "NO");
#ifdef UEP

    if (options->uep) {
      printf("  %2i   Read unique polymorphism? From file %30.30s\n",
	     MI_UEPFILE, options->uepfilename);
      printf("  %2i   Mutation rate for UEP is %10.5f %s\n",
	     MI_UEPRATE, options->ueprate, " x mutation rate");
      printf("  %2i   Base frequency for UEP is \"0\"=%f, \"1\"=%f\n",
	     MI_UEPFREQ, options->uepfreq0, options->uepfreq1);
    } else
      printf("  %2i   Read unique polymorphism? %20.20s\n", MI_UEPFILE, "NO");
#endif
    if (options->mighist) 
      {
	sprintf(outputstring,"%s %s", options->mighistfilename,
		options->mighist_all ? "(all events)" : "(migration events)");
	printf("  %2i   Show event statistics%42.42s\n", MI_MIGHISTOGRAM, outputstring);
	if(options->mighist_increment > 1)
	  sprintf(outputstring,"every %li %s",  options->mighist_increment,"sample steps");
	else
	  sprintf(outputstring,"every %s", "sample step");
	printf("       Events are recorded every     %32.32s\n",outputstring);
	sprintf(outputstring,"%f",(double) options->eventbinsize);
	printf("       Histogram bin width            %32.32s\n",outputstring);

      } 
    else 
      {
	sprintf(outputstring,"NO");
	printf("  %2i   Show event statistics          %32.32s\n", MI_MIGHISTOGRAM, outputstring);
    
      }
    if (options->skyline) 
      {
	remove_trailing_blanks(&options->skylinefilename);
	sprintf(outputstring,"%s", options->skylinefilename);
	printf("  %2i   Record parameter change through time?%26.26s\n", MI_SKYLINE, outputstring);
	sprintf(outputstring,"%f",(double) options->eventbinsize);
	printf("       Histogram bin width            %32.32s\n",outputstring);
      } 
    else 
      {
	sprintf(outputstring,"NO");
	printf("  %2i   Record parameter change through time?%26.26s\n", MI_SKYLINE, outputstring);
      }

    printf("\n\n  Are the settings correct?\n");
    printf
      ("  (type Y to go back to the main menu or the letter for the entry to change)\n===> ");
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    test = atoi(input);
    switch (test) 
      {
      case MI_INFILE:
	printf("  What is the datafile name?\n[Default: %s]\n===> ", options->infilename);
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	if (input[0] != '\0')
	  {
	    strcpy(options->infilename, input);
	  }
	break;
      case MI_RAND:
	do 
	  {
	    printf("  (A)utomatic or (S)eedfile or (O)wn\n");
	    printf("  Start value for Random-generator seed\n===> ");
	    fflush(stdout); FGETS(input, LINESIZE, stdin);
	    switch (uppercase(input[0])) 
	      {
	      case 'A':
		options->autoseed = AUTO;
#ifndef MAC
		timeseed = (unsigned long) time(NULL) / 4;
#else
		timeseed = (unsigned long) clock() / 4;
#endif
		options->inseed = (unsigned long) timeseed + 1;
		break;
	      case 'S':
		menu_get_filename(" What is the filename that contains the random number seed?", SEEDFILE, options->seedfilename);
		openfile(&options->seedfile, options->seedfilename, "r", NULL);
		if (options->seedfile) {
		  options->autoseed = NOAUTO;
		  fscanf(options->seedfile, "%ld%*[^\n]",
			 &options->inseed);
		  fclose(options->seedfile);
		} else
		  printf("\n\n  There is no seedfile present\n");
		break;
	      case 'O':
		options->autoseed = NOAUTOSELF;
		printf("  Random number seed (best values are x/4 +1)?\n");
		scanf("%ld%*[^\n]", &options->inseed);
		break;
	      }
	  }
	while (options->autoseed < AUTO || options->autoseed > NOAUTOSELF);
	break;
      case MI_TITLE:
	printf("  Enter a title? [max 80 Characters are used]\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	if (input[0] == '\0')
	  options->title[0] = '\0';
	else
	  sprintf(options->title,"%80.80s", input);
	break;
      case MI_SUMREAD:
	if(!options->bayes_infer)
	  {
	    printf
	      (" What is the filename for the summary of genealogies\n[Default: %s]\n===> ",
	       SUMFILE);
	    fflush(stdout); FGETS(input, LINESIZE, stdin);
	    if (input[0] == '\0')
	      strcpy(options->sumfilename, SUMFILE);
	    else {
	      strcpy(options->sumfilename, input);
	    }
	  }
	else
	  {
	    printf
	      (" What is the filename of the recorded Bayes data\n[Default: %s]\n===> ",
	       BAYESMDIMFILE);
	    fflush(stdout); FGETS(input, LINESIZE, stdin);
	    if (input[0] == '\0')
	      strcpy(options->bayesmdimfilename, BAYESMDIMFILE);
	    else {
	      strcpy(options->bayesmdimfilename, input);
	    }
	    unpad(options->bayesmdimfilename, " ");
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
      case MI_PROGRESS:
	printf("  Progress report during the run? <YES | NO>\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	switch (tolower(input[0])) {
	case 'n':
	  options->progress = FALSE;
	  options->verbose = FALSE;
	  break;
	case 'v':
	  options->progress = TRUE;
	  options->verbose = TRUE;
	  break;
	case '\0':
	case 'y':
	default:
	  options->progress = TRUE;
	  options->verbose = FALSE;
	  break;
	}
	input[0] = 'X';
	break;
      case MI_PRINTDATA:
	options->printdata = !options->printdata;
	break;
      case MI_OUTFILE:
	menu_get_filename("  What is the output filename?", OUTFILE, options->outfilename);
#ifdef PRETTY
	strcpy(options->pdfoutfilename,options->outfilename);
	strcat(options->pdfoutfilename,".pdf");
#endif
	break;
      case MI_TERSEPDF:
	options->tersepdf = !options->tersepdf;
	break;
      case MI_TREES:
	do 
	  {
	    printf("  Print genealogies:\n");
	    printf("  (N)one, (A)all [!], (B)est\n===> ");
	    fflush(stdout); FGETS(input, LINESIZE, stdin);
	  }
	while (strchr("NAB", (int) uppercase(input[0])) == NULL);
	switch (uppercase(input[0])) 
	  {
	  case 'N':
	    options->treeprint = myNONE;
	    break;
	  case 'B':
	    options->treeprint = BEST;
	    menu_get_filename(" What is the filename for storing recorded genealogies?", TREEFILE, options->treefilename);
	    break;
	  case 'A':
	    options->treeprint = ALL;
	    menu_get_filename(" What is the filename for storing recorded genealogies?", TREEFILE, options->treefilename);
	    printf("Give the increment to record genealogies\n===>");
	    fflush(stdout); FGETS(input, LINESIZE, stdin);
	    sscanf(input,"%li",&options->treeinc);
	    break;
	  default:
	    options->treeprint = myNONE;
	    break;
	  }
	break;
	//case MI_SUMFILE:
	//printf(" Save genealogy summaries? [NO]\n===> ");
	//fflush(stdout); FGETS(input, LINESIZE, stdin);
	//if (uppercase(input[0]) == 'Y') 
	//  {
	//    options->writesum = TRUE;
	//    menu_get_filename(" What is the filename for the summary of genealogies?", SUMFILE, options->sumfilename);
	//  } 
	//else 
	//  {
	//    options->writesum = FALSE;
	//  }
	//input[0] = 'X';
	//break;
      case MI_LOGFILE:
	printf(" Save logging information? [NO]\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	if (uppercase(input[0]) == 'Y') 
	  {
	    options->writelog = TRUE;
	    menu_get_filename(" What is the filename for logging?", LOGFILE, options->logfilename);
	  }
	else 
	  {
	    options->writelog = FALSE;
	  }
	input[0] = 'X';
	break;
#ifdef UEP
      case MI_UEPFILE:
	printf(" UEP estimation? [NO]\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	if (uppercase(input[0]) == 'Y') 
	  {
	    options->uep = TRUE;
	    menu_get_filename(" What is the filename of the UEP data?", UEPFILE, options->uepfilename);
	  } 
	else 
	  {
	    options->uep = FALSE;
	  }
      input[0] = 'X';
      break;
      case MI_UEPRATE:
      do 
	{
	  FPRINTF(stdout, "What is mutation rate for the UEP locus\n");
	  FPRINTF(stdout,
		  "[expressed as a ratio of the point mutation rate\n");
	  FPRINTF(stdout, "Reasonable values are 0.0 =< x <=1.0]\n===> ");
	  fflush(stdout); FGETS(input, LINESIZE, stdin);
	  options->ueprate = atof(input);
	}
      while (options->ueprate < 0.0);
      break;
      case MI_UEPFREQ:
      do 
	{
	  FPRINTF(stdout, "What is the base frequency for the \"1\" allele\n");
	  FPRINTF(stdout, "Reasonable values are 0.0 < x <1.0]\n===> ");
	  fflush(stdout); FGETS(input, LINESIZE, stdin);
	  options->uepfreq1 = atof(input);
	  options->uepfreq0 = 1. - options->uepfreq1;
	}
      while (options->uepfreq1 < 0.0);
      break;
#endif
      case MI_MIGHISTOGRAM:
	printf(" Display event statistic and save events into file?\n [options are: NO, ALL events, MIGration events]\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	if (uppercase(input[0]) == 'N') 
	  {
	    options->mighist = FALSE;
	  } 
	else 
	  {
	    options->mighist = TRUE;
	    if (uppercase(input[0]) == 'A')
	      options->mighist_all = TRUE;
	    else
	      options->mighist_all = FALSE;

	    /*	    printf(" Thin the raw data and sample only every x sample steps,\n Enter increment [default is 1]:\n===> ");
	    fflush(stdout); FGETS(input, LINESIZE, stdin);
	    options->mighist_increment = atol(input);
	    if(options->mighist_increment < 1)
	      options->mighist_increment = 1;
	    */
	    menu_get_filename(" Specify a filename for the raw event data!", MIGHISTFILE, options->mighistfilename);
	    
	    do
	      {
		printf("Events are recorded per time intervals.\nWhat is the time interval size?\n");
		printf("[Current value: %f, try values that are about 10x smaller than the largest Theta]\n===> ",
		       (double) options->eventbinsize);
		fflush(stdout); FGETS(input, LINESIZE, stdin);
		if(input[0]=='\0' && options->eventbinsize > 0.0f)
		  break;
		options->eventbinsize = (float) atof(input);
	      } 
	    while (options->eventbinsize <= 0.0f);

	  }
	input[0] = 'X';
	break;
      case MI_SKYLINE:
	printf(" Display time dependent parameters and save them into file? [Choices are NO or  YES]\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	if (uppercase(input[0]) == 'N') 
	  {
	    options->skyline = FALSE;
	  } 
	else 
	  {
	    options->skyline = TRUE;
	    do
	      {
		printf("What is the time interval size?\n");
		printf("[Current value: %f, try values that are about 10x smaller than the largest Theta]\n===> ",
		       (double) options->eventbinsize);
		fflush(stdout); FGETS(input, LINESIZE, stdin);
		if(input[0]=='\0' && options->eventbinsize > 0.0f)
		  break;
		options->eventbinsize = (float) atof(input);
	      } 
	    while (options->eventbinsize <= 0.0f);
	    menu_get_filename(" Specify a filename for the raw skyline output!", SKYLINEFILE, options->skylinefilename);
	    if(!options->mighist)
	      {
		printf("The skyline option needs the coalescence and migration event recording\n");
		menu_get_filename(" Specify a filename for the raw event data!", MIGHISTFILE, options->mighistfilename);		
		options->mighist_all = TRUE;
		options->mighist_increment = 0;
	      }

	  }
	input[0] = 'X';
	break;
      default:
	break;
      }
  }
  while (uppercase(input[0]) != 'Y');
  myfree(stringstep);
  myfree(string2);
  myfree(string3);
  myfree(string4);
}


///
/// menu parameter related functions
void set_menu_priorpercent(option_fmt *options, float value)
{
  int i;
  long v = (long) value;
  if (v < 0 ||  v > 100)
    v = 1;
  for (i=0;i<STARTGUESSNUM;i++)
    {
      options->startguess[i][0] = PRIOR;
      options->startguess[i][1] = v;
    }
}

void set_menu_randomprior(option_fmt *options)
{
  int i;
  for (i=0;i<STARTGUESSNUM;i++)
    {
      options->startguess[i][0] = RANDOMPRIOR;
    }
}

void set_menu_ownstartvalue(option_fmt *options, long numpop)
{
  int i;
  int j;
  int z=0;
  int check;
  float value;
  char *input;
  input = (char *) mycalloc(LINESIZE,sizeof(char));
  if (numpop == 0) {
    printf("  I do not know yet how many populations the data set contains? Specify the number of populations > ");
    fflush(stdout);
    scanf("%ld", &numpop);
    scanf("%*[^\n]");
  }
  realloc_startparam(options,numpop);
  printf("THETA: Mutation scaled population size\n  Give the initial Theta estimates:\n");
  for (i = 0; i < numpop; i++) 
    {
      printf("Population %3i> ", i+1);fflush(stdout);
      scanf("%f", &value);
      fill_startparam(options,THETAPRIOR,i,value);
    }
  //options->startparam.numtheta = numpop;
  options->startguess[0][0] = OWN; 
  printf("MIGRATION: Mutation scaled migration rate\n");
  printf("Give initial migration rate estimate\n[default parameters are mutation-scaled migration rates\n, if you estimate Nm instead of M the values entered are in terms of xNm where x depends on the data (see manual)]");
  z = 0;
  for (i = 0; i < numpop; i++) 
    {
      for (j = 0; j < numpop; j++) 
	{
	  if (j == i) 
	    {
	      printf("Population %3i              > --\n",i + 1);
	      continue;
	    }
	  printf("From population %-3i to %-3i> ", j + 1, i + 1);fflush(stdout);
	  check = scanf("%f", &value);
	  scanf("%*[^\n]");
	  if (check == 0)
	    break;
	  fill_startparam(options,MIGPRIOR,z++,value);
	}
    }
  options->startguess[1][0] = OWN; 
  fflush(stdout); 
  //FGETS(input, LINESIZE, stdin);
  printf
    ("\n     You typed in the following migration matrix\n\n     ");
  show_migownparam(stdout,options, &input);
  printf("\n     [Press <Return> to continue]\n");
  printf("\n\n");
  getchar();
  myfree(input);
}


///
///allows to set all parameter related issues:start parameters, mutation type, migration model
void
menuParameters(option_fmt * options) 
{
  float           value;
  char            input[LINESIZE];
  char            custmexplain[LINESIZE];
  char            localities[LINESIZE];
  char           *temp;
  char           *outputstring;
  long            count = 0;
  long            numpop = 0;
  long            i;
  MYREAL          sum = 0.;
  //double          tmpdouble;
  //  int             retval;
  outputstring = (char *) mycalloc(LINESIZE, sizeof(char));


  do {
    printf("  ===================================================================\n");    
    printf("  PARAMETERS\n");
    printf("  ===================================================================\n\n");
    printf
      ("  Start parameters:\n");
    //switch_startparam(outputstring, options);
    strcpy(outputstring,"combined startparameter report not finished");
    printf("  1   First parameter values are?\n      %64s\n", outputstring);
    printf
      ("\n Gene flow parameter and Mutation rate variation among loci:\n");
    printf("  2   Use M for the gene flow parameter   %28.28s\n",
	   options->usem ? "YES [M=m/mu]" : "NO [Theta * M]");
    if(options->murates)
      {
	if(options->murates_fromdata)
	  {
	    printf("  3   Mutation rate is  %46.46s\n","Varying (from data)");
	  }
	else
	  {
	    printf("  3   Mutation rate is  %46.46s\n","Varying (from user)");
	  }
      }
    else
      {
	if(options->bayesmurates)
	  {
	    printf("  3   Mutation rate is  %46.46s\n","Estimated (see prior menu: Rate)");
	  }
	else
	  {
	    printf("  3   Mutation rate is  %46.46s\n","Constant");
	  }
      }
    strcpy(custmexplain, custom_migration_type(options->migration_model));
    set_localities_string(&localities[0],options);
    printf("\n  Structured coalescent model and combination of localities:\n");
    printf("  4   Sampling localities %43.43s\n", localities);
    printf("  5   Model is set to %48.48s\n", custmexplain);
    if(options->geo)
      printf("  6   Geographic distance matrix: YES:%34.34s\n", options->geofilename);
    else
      printf("  6   Geographic distance matrix: %36.36s\n", "NO");
    printf("\n\n  Are the settings correct?\n");
    printf
      ("  (Type Y to go back to the main menu or the letter for an entry to change)\n===> ");
    fflush(stdout); 
    FGETS(input, LINESIZE, stdin);
    switch (input[0]) 
      {
      case '1': // change start parameters
	printf("  Which method?]\n");
	printf("  P    Value at x%% of prior\n");
	printf("  R    Random from prior distribution\n");
	printf("  O    You enter your OWN start values\n");
	printf("  \n===> ");
	fflush(stdout); 
	FGETS(input, LINESIZE, stdin);
	switch (uppercase(input[0])) {
	case 'P':
	  do
	    {
	      printf("  Enter the %% point [0-100] of the prior distribution\n");
	      printf("  Good values are closer to the minimum, for example 10%% \n");
	      printf("  \n===> ");
	      fflush(stdout); 
	      FGETS(input, LINESIZE, stdin);
	      value = atol(input);
	    } while (value < 0 && value > 100);
	  set_menu_priorpercent(options, value);
	  break;
	case 'R':
	  set_menu_randomprior(options);
	  break;
	case 'O':
	  set_menu_ownstartvalue(options, numpop);	  
	  break;
      default:
	set_menu_priorpercent(options, 1.0);
	break;
      }
      break;
    case '2':
      options->usem = !options->usem;
      break;
    case '3':
      printf("Mutation rate among loci\n");
      printf("(C)onstant    All loci have the same mutation rate [default]\n");
      printf("(E)stimate    Mutation rate \n");
      printf("(V)arying     Mutation rates are different among loci [user input]\n");
      printf("(R)elative    Mutation rates estimated from data\n===> ");
      fflush(stdout); 
      FGETS(input, LINESIZE, stdin);
      switch (uppercase(input[0])) {
      case 'E':
	options->bayesmurates = TRUE;
	options->murates = FALSE;
	options->murates_fromdata=FALSE;
	break;
      case 'V':
	printf
	  ("Enter the number of loci and a rate for each of them\n[the rates are normalized to average to 1.0]\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	temp = strtok(input, ";,;\n");
	if (temp == NULL) {
	  options->murates = FALSE;
	  options->gamma = FALSE;
	  options->murates_fromdata=FALSE;
	  options->bayesmurates = FALSE;
	} else {
	  options->gamma = FALSE;
	  options->murates = TRUE;
	  options->murates_fromdata=FALSE;
	  options->muloci = atol(temp);
	  options->bayesmurates = FALSE;
	  options->mu_rates = (MYREAL *) myrealloc(options->mu_rates, sizeof(MYREAL) * (size_t) options->muloci);
	  //  printf("%i> menu.c: opt realloc murate size %li\n",myID,  options->muloci * sizeof (MYREAL));			
	  sum = 0.;
	  count = 0;
	  while (temp != NULL) {
	    temp = strtok(NULL, ";,;\n");
	    if (temp == NULL)
	      break;
	    options->mu_rates[count] = atof(temp);
	    sum += options->mu_rates[count];
	    count++;
	  }
	  for (i = count; i < options->muloci; i++) {
	    options->mu_rates[i] = options->mu_rates[count - 1];
	    sum += options->mu_rates[count];
	  }
	  sum /= options->muloci;
	  for (i = 0; i < options->muloci; i++) {
	    options->mu_rates[i] /= sum;
	  }
	}
	break;
      case 'R':
	options->gamma = FALSE;
	options->murates = TRUE;
	options->murates_fromdata=TRUE;
	options->bayesmurates = FALSE;
	break;
      default:
	options->gamma = FALSE;
	options->murates = FALSE;
	options->murates_fromdata=FALSE;
	options->bayesmurates = FALSE;
	break;
      }
      break;
    case '4':
      set_menu_localities(options);
      break;
    case '5':	/* fill in the custom migration
		 * matrix */
      read_custom_menu_migration(options);
      break;
    case '6':
      options->geo = !options->geo;
      if(options->geo)
	{
	  menu_get_filename("  What is the filename for the distance matrix file?",
		     GEOFILE, options->geofilename);
	}
      break;
    }
  }
  while (uppercase(input[0]) != 'Y');
  myfree(outputstring);
}

void set_menu_localities(option_fmt *options)
{
  char *input;
  char *tmp;
  int i;
  input = (char *) mycalloc(LINESIZE,sizeof(char));
  tmp = (char *) mycalloc(LINESIZE,sizeof(char));
  printf("Associate sampling locations with populations\n");
  printf("---------------------------------------------\n");
  printf("This menu allows to combine sample locations into populations\n");
  printf("For example there are 4 locations: 1, 2, 3, 4\n");
  printf("They can be combined into 2 populations\n");
  printf("by mapping the 4 positions 1, 2, 3, 4 to 1, 1, 1, 2\n");
  printf("Migrate will now combine the first three locations\n\n");
  printf("Give (1) the number of populations <return> then (2) the mapping\n\n");
  printf("How many localities are in the data set?\n");
  printf("[Default: every sampling location is a population]\n");
  printf("> ");fflush(stdout);
  FGETS(input,LINESIZE,stdin);
  if(input[0]!='\0')
    {
      options->newpops_numalloc = atoi(input);
      options->newpops = (long*) myrealloc(options->newpops, options->newpops_numalloc * sizeof(long));
      printf("Enter now the remappings (little checking is done with this, enter exactly the number of locations)\n");
      for(i = 1; i <= options->newpops_numalloc; i++)
	{
	  switch((int) (log10((double) options->newpops_numalloc)))
	    {
	    case 0: 
	    case 1:
	      printf("%2i",i); break;
	    case 2:
	      printf("%3i",i); break;
	    default:  
	      printf("%4i",i); break;
	    }
	}
      printf("\n>");fflush(stdout);
      FGETS(input,LINESIZE,stdin);
      set_localities(&input,&tmp,options);
    }
  else
    {
      printf("All locations are treated as populations.");
      options->newpops = (long*) myrealloc(options->newpops, sizeof(long));
      options->newpops[0] = 1;
      options->newpops_numalloc  = 1;
    }
  myfree(input);
  myfree(tmp);
}

///print the standard question "are the settings correct? ....." at the end of the menu
void            print_bottom_menu_part(void) 
{
  printf("\n\n  Are the settings correct?\n");
  printf("  (Type Y to go back to the main menu or the number for a menu to change)\n===>");
}

///
///displays and makes changes to the options that manage the analysis strategy and run condition
///
void            menuStrategy(option_fmt * options) {
  boolean         done = FALSE;
  do {
    printf("  ===================================================================\n");    
    printf("  SEARCH STRATEGY\n");
    printf("  ===================================================================\n\n");
    //if (!options->bayes_infer) {
    //  display_ml_mcmc(options);
    //  print_bottom_menu_part();
    //  done = menuStrategy_ml(options);
    //} else {
      display_bayes_mcmc(options);
      print_bottom_menu_part();
      done = menuStrategy_bayes(options);
      //}
  } while (done == FALSE);
}

//-------------------------------------------------------------------------------------
// prior menu functions

/// \brief Menu for multiplier proposal
/// Menu for multiplier proposal
void set_mult_prior(int paramgroupm, prior_fmt *prior)
{
  (void) paramgroupm;
  char input[LINESIZE];
  long count = 0;
  do {
    printf("Specify the MINIMUM, MAXIMUM and maximal MULTIPLIER\n[Current setting: {%f %f %f}\n===> ",
	   prior->min, prior->max, prior->delta);
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if (input[0] == '\0')
      break;
#ifdef USE_MYREAL_FLOAT
    count = sscanf(input, "%f%f%f", &prior->min, &prior->max, &prior->delta);
#else
    count = sscanf(input, "%lf%lf%lf", &prior->min, &prior->max, &prior->delta);
#endif
  } while (count != 3 && count > 0);
}

/// \brief Menu for exponential proposal
/// Menu for exponential proposal
void set_exp_prior(int paramgroupm, prior_fmt *prior)
{
  (void) paramgroupm;
  char input[LINESIZE];
  long count = 0;
  do {
    printf("Specify the MINIMUM, MEAN, MAXIMUM\n[Current setting: {%f %f %f}\n===> ",
	   prior->min, prior->mean, prior->max);
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if (input[0] == '\0')
      break;
#ifdef USE_MYREAL_FLOAT
    count = sscanf(input, "%f%f%f", &prior->min, &prior->mean, &prior->max);
#else
    count = sscanf(input, "%lf%lf%lf", &prior->min, &prior->mean, &prior->max);
#endif
  } while (count != 3 && count > 0);
}

/// \brief Menu for gamma proposal
/// Menu for gamma proposal, mean = alpha * beta = 1 -> 
void set_gamma_prior(int paramgroupm, prior_fmt *prior)
{
  (void) paramgroupm;
  char input[LINESIZE];
  long count = 0;
  do {
    printf("Specify the MINIMUM, MEAN, MAXIMUM, ALPHA\n[Current setting: {%f %f %f %f}\n===> ",
	   prior->min, prior->mean, prior->max, prior->alpha);
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if(input[0]=='\0')
      return;
#ifdef USE_MYREAL_FLOAT
    count = sscanf(input, "%f%f%f%f", &prior->min, &prior->mean, &prior->max, &prior->alpha);
#else
    count = sscanf(input, "%lf%lf%lf%lf", &prior->min, &prior->mean, &prior->max, &prior->alpha);
#endif
  } while (count != 4 && count > 0);
}

/// \brief Menu for exponential with window proposal
/// Menu for exponential with window proposal
void set_wexp_prior(int paramgroupm, prior_fmt *prior)
{
  (void) paramgroupm;
  char input[LINESIZE];
  long count = 0;
  do {
    printf("Specify the MINIMUM, MEAN, MAXIMUM and DELTA\nUse for DELTA about 1/10 of the MIN-MAX range\n[Current setting: {%f %f %f %f}\n===> ",
	   prior->min, prior->mean, prior->max, prior->delta);
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if(input[0]=='\0')
      return;
#ifdef USE_MYREAL_FLOAT
    count = sscanf(input, "%f%f%f%f", &prior->min, &prior->mean, &prior->max, &prior->delta);
#else
    count = sscanf(input, "%lf%lf%lf%lf", &prior->min, &prior->mean, &prior->max, &prior->delta);
#endif
  } while (count != 4 && count > 0);
}

/// \brief Menu for uniform proposal
/// Menu for uniform proposal
void set_uni_prior(int paramgroupm, prior_fmt *prior)
{
  (void) paramgroupm;
  char input[LINESIZE];
  long count = 0;
  do {
    printf("Specify the MINIMUM, MAXIMUM, DELTA\nUse for DELTA about 1/10 of the MIN-MAX range\n[Current setting: {%f %f %f}\n===> ",
	   prior->min, prior->max, prior->delta);
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if(input[0]=='\0')
      return;
#ifdef USE_MYREAL_FLOAT
    count = sscanf(input, "%f%f%f", &prior->min, &prior->max, &prior->delta);
#else
    count = sscanf(input, "%lf%lf%lf", &prior->min, &prior->max, &prior->delta);
#endif
  } while (count != 3 && count > 0);
  prior->mean = (prior->max + prior->min) / 2.;
  //prior->delta = (prior->max - prior->min) / 10.;
}


/// \brief Menu for binning
/// Menu for uniform proposal
void set_binning(prior_fmt *prior)
{
  char input[LINESIZE];
  long count = 0;
  do {
    printf("Specify the maximal number of bins used\n[Current setting: %li]\n===> ",
	   prior->bins);
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if(input[0]=='\0')
      return;
    count = sscanf(input, "%li", &prior->bins);
  } while (count != 1);
}

/// \brief return "set" "-" depending on prior set
/// returns a short string that shows whether this prior distribution is set
char * is_prior(int priorkind, prior_fmt *p, int priorset)
{
  long count = 0;
  if (p->kind == priorkind && p->type == priorset)
    count++;
  if (count > 0)
    {
      return "Set";
    }
  return  "-";
}

/// \brief returns priortype sting
/// returns a string that shows what prior distribution is set
char * is_proposaltype(boolean proposalset)
{
  switch(proposalset)
    {
    case TRUE: return  "Slice sampling"; 
    default: return "Metropolis sampling";
    }
}


/// \brief Menu to set all prior distributions
/// menu to set proposal distribution, returns TRUE when standard course of action
/// FALSE when a problem occured.
boolean set_proposal_menu(int paramgroup, option_fmt *options)
{
  char            input[LINESIZE];
  boolean proposal;

  switch(paramgroup)
    {
    case THETAPRIOR:
      printf("Proposal setting for all population sizes [Theta]:\n");
      proposal = options->slice_sampling[THETAPRIOR];
      break;
    case MIGPRIOR:
      printf("Proposal setting for all migration rates [%s]:\n",options->usem ? "M" : "Theta*M");
      proposal = options->slice_sampling[MIGPRIOR];
      break;
    case RATEPRIOR:
      printf("Proposal setting for mutation rate modifier [r]:\n");
      proposal = options->slice_sampling[RATEPRIOR];
      break;
    case SPLITPRIOR:
      printf("Proposal setting for divergence mean [Delta]:\n");
      proposal = options->slice_sampling[SPECIESTIMEPRIOR];
      break;
    case SPLITSTDPRIOR:
      printf("Proposal setting for the divergence spread  [Std]:\n");
      proposal = options->slice_sampling[SPECIESSTDPRIOR];
      break;
    default:
      return FALSE;
    }
    input[0] = '\0';
    printf("Which proposal distribution? [current: %11s] (choices: Slice or Metropolis)\n===> ", 
    	   is_proposaltype(proposal));
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    switch (uppercase(input[0])) 
      {
      case 'M':
	options->slice_sampling[paramgroup] = FALSE;
	break;
      case 'S':
      default:
	options->slice_sampling[paramgroup] = TRUE;
	break;
      }
    input[0]='\0';
    return TRUE;
  }


/// \brief returns priortype string
/// returns a string that shows what prior distribution is set
long  is_priortype(prior_fmt *p, long pnum, int priortype)
{
  long c=7;
  long i;
  long *count;
  count = (long *) calloc(8,sizeof(long));
  for (i=0;i<pnum;i++)
    {
      
      if (p[i].type == priortype)
	{
	  count[0]++;
	  switch(p[i].kind)
	    {
	    case MULTPRIOR: count[1]++; break;
	    case EXPPRIOR: count[2]++; break;
	    case WEXPPRIOR: count[3]++; break;
	    case GAMMAPRIOR: count[4]++; break;
	    case UNIFORMPRIOR: count[5]++; break;
	    case NORMALPRIOR:  count[6]++; break;
	    default:  count[7]++; break;
	    }
	}
    }
  if (count[0]==0)
    {
      free(count);
      return c;
    }
  for (i=1;i<8;i++)
    {
      if (count[i]>count[c])
	c=i;
    }
  free(count);
  return c;
}

prior_fmt * set_theta_priormenu(option_fmt *options)
{
  char input[LINESIZE];
  long pop;
  double mini;
  double mean;
  double maxi;
  prior_fmt *prior=NULL;
  mini = SMALLEST_THETA;
  mean = strchr(SEQUENCETYPES,options->datatype) ? DNA_GUESS_THETA : ALLELE_GUESS_THETA; 
  maxi = mean * 10.0;
  while(prior==NULL)
    {
      printf("N change number of populations (currently set to %li)\n",options->numpop);
      printf("1 Set the same prior for all Thetas\n");
      printf("2 Set prior for Theta of population i: first population in data is 1, seconde is 2 etc.\n");
      printf("Enter either N, 1, or 2, or <enter> to return from to the parent menu\n>>> ");
      fflush(stdout);
      FGETS(input, LINESIZE, stdin);
      switch(input[0])
	{
	case 'n':
	case 'N':
	  printf("\n  Give the number of populations\n===> ");
	  fflush(stdout); 
	  FGETS(input, 1024L, stdin);
	  options->numpop = atol(input);
	  printf("\n  NEW NUMBER OF POPULATIONS = %li\n",options->numpop);
	  break;
	case '1':
	  //pop = 0;
	  prior =  find_prior_menu(-1, -1, THETAPRIOR, options);
	  break;
	case '2':
	  printf("Enter the population label number (1,2,3,...)\n");
	  FGETS(input, LINESIZE, stdin);
	  if (input[0]=='\0')
	    break;
	  pop = atol(input) - 1;
	  if (pop >= options->numpop)
	    {
	      printf("\nLabel is larger than number of populations! Is the number of populations set?\n\n");
	      break;
	    }
	  else
	    {
	      prior =  find_prior_menu(pop,pop, THETAPRIOR, options);
	    }
	}
    }
  return prior;
}

prior_fmt * set_mig_priormenu(option_fmt *options)
{
  char input[LINESIZE];
  char * fromstr = (char *) calloc(2 * LINESIZE, sizeof(char));
  char * tostr = fromstr + LINESIZE;
  long pop;
  double mini;
  double mean;
  double maxi;
  prior_fmt *prior=NULL;
  mini = SMALLEST_MIGRATION;
  mean = strchr(SEQUENCETYPES,options->datatype) ? DNA_GUESS_MIG : ALLELE_GUESS_MIG; 
  maxi = mean * 10.0;
  while(prior == NULL)
    {
      printf("N change number of populations (currently set to %li)\n",options->numpop);
      printf("1 Set the same prior for all Migration rates\n");
      printf("2 Set prior for migration from pop j to pop i\n");
      FGETS(input, LINESIZE, stdin);
      switch(input[0])
	{
	case 'n':
	case 'N':
	  printf("\n  Give the number of populations\n===> ");
	  fflush(stdout); 
	  FGETS(input, 1024L, stdin);
	  options->numpop = atol(input);
	  printf("\n  NEW NUMBER OF POPULATIONS = %li\n",options->numpop);
	  break;
	case '1':
	  //pop = 0;
	  prior =  find_prior_menu(-1,-1, MIGPRIOR, options);
	  break;
	case '2':  
	  printf("Set prior for %s [specify the population label for the from- and to- population\n Enter two numbers with a space]\n",options->usem ? "M" : "Theta*M");
	  FGETS(input, LINESIZE, stdin);
	  if (input[0]=='\0')
	    break;
	  sscanf(input,"%s%s", fromstr, tostr);
	  long from = atol(fromstr) -1;
	  long to   = atol(tostr) -1;
	  pop = mm2m(from, to, options->numpop);
	  if (pop >= options->numpop*options->numpop)
	    {
	      printf("\nLabel is larger than number of populations**2! Is the number of populations set?\n\n");
	      break;
	    }
	  else
	    {
	      prior =  find_prior_menu(from,to, MIGPRIOR, options);	      
	    }
	}
    }
  free(fromstr);
  return prior;
}

prior_fmt * set_split_priormenu(option_fmt *options, int type)
{
  char input[LINESIZE];
  char * fromstr = (char *) calloc(2 * LINESIZE, sizeof(char));
  char * tostr = fromstr + LINESIZE;
  long pop;
  double mini;
  double mean;
  double maxi;
  prior_fmt  *prior=NULL;
  char thistype[14] = "splitting";
  int thisdefault = 3;
  long splitstart;
  long mu = (options->bayesmurates == TRUE ? 1 : 0);
  if (type == SPLITSTDPRIOR)
    {
      strcpy(thistype,"splitting std"); 
      thisdefault = 4;
    }
  splitstart = options->numpop*options->numpop + mu;
  mini = SMALLEST_SPLIT;
  mean = strchr(SEQUENCETYPES,options->datatype) ? DNA_GUESS_SPLIT : ALLELE_GUESS_SPLIT; 
  maxi = mean * 10.0;
  while(prior == NULL)
    {
      printf("N change number of populations (currently set to %li)\n",options->numpop);
      printf("1 Set the same prior for all split times (%s)\n",thistype);
      printf("2 Set prior for %s from pop j to pop i\n",thistype);
      FGETS(input, LINESIZE, stdin);
      switch(input[0])
	{
	case 'n':
	case 'N':
	  printf("\n  Give the number of populations\n===> ");
	  fflush(stdout); 
	  FGETS(input, 1024L, stdin);
	  options->numpop = atol(input);
	  printf("\n  NEW NUMBER OF POPULATIONS = %li\n",options->numpop);
	  break;
	case '1':
	  prior =  find_prior_menu(-1,-1, SPLITPRIOR, options);	      
	  break;
	case '2':  
	  printf("Set prior for %s [specify the population label for the ancestor- and descendent-population\n Enter two numbers with a space]\n",thistype);
	  FGETS(input, LINESIZE, stdin);
	  if (input[0]=='\0')
	    break;
	  sscanf(input,"%s%s", fromstr, tostr);
	  long from = atol(fromstr) -1;
	  long to   = atol(tostr) -1;
	  pop = mm2m(from, to, options->numpop);
	  if (pop >= options->numpop*options->numpop)
	    {
	      printf("\nLabel is larger than number of populations**2! Is the number of populations set?\n\n");
	      break;
	    }
	  else
	    {
	      prior =  find_prior_menu(from,to, SPLITPRIOR, options);	      
	    }
	}
    }
  free(fromstr);
  return prior;
}


/// \brief Menu to set all prior distributions
/// menu to set all prior distributions, returns TRUE when standard course of action
/// FALSE when a problem occured.
/// This menu will do this:
///     - ask the number of populations (enter means default -- the value show in the question)
///     - previous menu shows the parameter groups, now after knowing the number of populations
///     - menu shows:
///                - all parameter in group have same prior: YES/NO
///                - if no, specify parameter ID with population numbers (1 ,, n)
///                - if no, for migration specify pair of parameters from, to
///                - if no, for splits and std specify ancestor and descendent
      
boolean set_prior_menu(int paramgroup, option_fmt *options)
{
  int        tmp,i;
  char       input[LINESIZE];
  prior_fmt *p=NULL, *prior;
  switch(paramgroup)
    {
    case THETAPRIOR:
      printf("Prior setting for all population sizes [Theta]:\n");
      prior = set_theta_priormenu(options);
      break;
    case MIGPRIOR:
      printf("Prior setting for all migration rates [%s]:\n",options->usem ? "M" : "Theta*M");
      prior = set_mig_priormenu(options);
      break;
    case SPLITPRIOR:
      printf("Prior setting for all splitting events:\n");
      prior = set_split_priormenu(options,SPLITPRIOR);
      break;
    case SPLITSTDPRIOR:
      printf("Prior setting for all stand. dev. of the splittings:\n");
      prior = set_split_priormenu(options,SPLITSTDPRIOR);
      break;
    case RATEPRIOR:
      i=0;
      printf("Prior setting for mutation rate modifier [r]:\n");
      p = &options->bayes_priors[0];
      while (p->type != RATEPRIOR && i< options->bayes_priors_num)
	p = &options->bayes_priors[i++];
      prior = p;
      break;
    default:
      return FALSE;
    }
  if (prior==NULL)
    return FALSE;
  do{
    input[0] = '\0';
    //    printf("\nSet the prior distribution and its parameters\n\n");
    //printf("  %i   Set Slice sampler?                             %11s\n", SLICE, 
    //	   is_prior(SLICE,options->bayesprior[paramgroup]));
    //printf("  %i   Set Multiplication prior distribution?         %11s\n", MULTPRIOR, 
    //	   is_prior(MULTPRIOR, prior , paramgroup));
    printf("  %i   Set Exponential prior distribution?            %11s\n", EXPPRIOR, 
	   is_prior(EXPPRIOR,prior,paramgroup));
    printf("  %i   Set Exponential prior with window distribution?%11s\n", WEXPPRIOR,
	   is_prior(WEXPPRIOR,prior,paramgroup));
    printf("  %i   Set Gamma prior distribution?                  %11s\n", GAMMAPRIOR, 
    	   is_prior(GAMMAPRIOR,prior,paramgroup));
    printf("  %i   Set Uniform prior distribution?                %11s\n", UNIFORMPRIOR,
	   is_prior(UNIFORMPRIOR,prior,paramgroup));
    printf("  B   Set binning of prior and posterior              %11li\n\n",
	   prior->bins);
    printf("  (Type Y to go back to the main menu or the letter for an entry to change)\n===> ");
    fflush(stdout); 
    FGETS(input, LINESIZE, stdin);
    if(uppercase(input[0]) == 'Y')
      return TRUE;
    input[1] = '\0';
    if (uppercase(input[0])=='B')
      tmp = -99;
    else
      tmp = atoi(input);
    switch (tmp) 
      {
      case UNIFORMPRIOR:
	prior->kind = UNIFORMPRIOR;
	set_uni_prior(paramgroup, prior);
	break;
      case EXPPRIOR:
	prior->kind = EXPPRIOR;
	set_exp_prior(paramgroup, prior);
	break;
      case WEXPPRIOR:
	prior->kind = WEXPPRIOR;
	set_wexp_prior(paramgroup, prior);
	break;
	//case MULTPRIOR:
	//set_mult_prior(paramgroup, prior);
	//break;
      case GAMMAPRIOR:
	prior->kind = GAMMAPRIOR;
	set_gamma_prior(paramgroup, prior);
	break;
      case -99:
	set_binning(prior);
	break;
      default:
	prior->kind = UNIFORMPRIOR;
	set_uni_prior(paramgroup, prior);
	break;
      }
  }    while (uppercase(input[0]) != 'Y');
  input[0]='\0';
  return TRUE;
}

void set_autotuning(option_fmt *options)
{

  char input[LINESIZE];
  printf("Autotuning for Metropolis-Hastings sampling (YES or NO)\n");
  printf("===> ");
  fflush(stdout);
  FGETS(input, LINESIZE, stdin);
  if (input[0] != '\0')  
    {
      if(input[0]=='Y' || input[0]=='y')
	options->has_autotune = TRUE;
      else
	options->has_autotune=FALSE;
    }
  if (options->has_autotune)
    {
      printf("Give a target acceptance ratio [suggested 0.44]\n");
      printf("===> ");
      fflush(stdout);
      FGETS(input, LINESIZE, stdin);
      if (input[0] != '\0')  
	options->autotune = atof(input);
      else
	{
	  if (options->autotune>1.0 || options->autotune<0.01)
	    options->autotune=AUTOTUNEDEFAULT;
	}
    }
}

///
/// sets proposal distribution
void menuProposal(option_fmt * options) {
  char            input[LINESIZE];
  do 
    {
      input[0] = '\0';
      printf ("  Proposal distribution setting [You still need to set the PRIOR distribution!]:\n");
      printf("  1   Set proposal distribution for Theta?          %11s\n",is_proposaltype(options->slice_sampling[THETAPRIOR]));
      printf("  2   Set proposal distribution for Migration?      %11s\n",is_proposaltype(options->slice_sampling[MIGPRIOR]));
      printf("  3   Set proposal distribution for Divergence?     %11s\n",is_proposaltype(options->slice_sampling[SPECIESTIMEPRIOR]));
      printf("  4   Set proposal distribution for Div. Spread?    %11s\n",is_proposaltype(options->slice_sampling[SPECIESSTDPRIOR]));
      if(options->bayesmurates)
	{      
	  printf("  5   Set mutation rate modifier prior distribution?%11s\n",
	     is_proposaltype(options->slice_sampling[RATEPRIOR]));
	  if(options->has_autotune)
	      printf("  6   Autotuning of the acceptance ratio        %7s %.2f\n","YES",options->autotune);
	  else
	      printf("  6   Autotuning of the acceptance ratio        %11s\n","NO");
	}
      else
	{
	  if(options->has_autotune)
	      printf("  5   Autotuning of the acceptance ratio        %7s %.2f\n","YES",options->autotune);
	  else
	      printf("  5   Autotuning of the acceptance ratio        %11s\n","NO");
	}
      printf("\n  (Type Y to go back to the main menu or the letter for a menu to change)\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      switch (input[0]) 
	{
	case '1':
	  set_proposal_menu(THETAPRIOR, options);
	  break;
	case '2':
	  set_proposal_menu(MIGPRIOR, options);
	  break;
	case '3':
	  set_proposal_menu(SPECIESTIMEPRIOR, options);
	  break;
	case '4':
	  set_proposal_menu(SPECIESSTDPRIOR, options);
	  break;
	case '5':
	  if(options->bayesmurates)
	    {      
	      set_proposal_menu(RATEPRIOR, options);
	    }
	  else
	    {
	      set_autotuning(options);
	    }
	  break;
	case '6':
	  set_autotuning(options);
	  break;
	default:
	  break;
	}
    } while (uppercase(input[0]) != 'Y');;
}
///
/// sets prior distribution parameters
void menuPrior(option_fmt * options) {
  char            input[LINESIZE];
  const char text[8][20]={"All", "Multiplier", "Exponential", "Exp window", "Gamma", "Uniform", "Truncated Normal", "-"};
  prior_fmt *p = options->bayes_priors;
  long pnum = options->bayes_priors_num;
  do 
    {
      p = options->bayes_priors;
      pnum = options->bayes_priors_num;
      input[0] = '\0';
      printf ("  Prior distribution setting:\n");
      printf("  N   Change the number of populations, currently   %11li\n",options->numpop);
      printf("  1   Set Theta prior distribution?                 %11s\n", 
	     text[is_priortype(p,pnum, THETAPRIOR)]);
      printf("  2   Set Migration prior distribution?             %11s\n",
	     text[is_priortype(p, pnum, MIGPRIOR)]);
      printf("  3   Set mutation rate modifier prior distribution?%11s\n",
	     text[is_priortype(p,pnum, RATEPRIOR)]);
      printf("  4   Set splitting parameter prior distribution?   %11s\n",
	     text[is_priortype(p, pnum, SPECIESTIMEPRIOR)]);
      printf("  5   Set splitting std dev. prior distribution?    %11s\n",
	     text[is_priortype(p, pnum, SPECIESSTDPRIOR)]);
      printf("\n  (Type Y to go back to the main menu or the letter for a menu to change)\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      switch (uppercase(input[0])) 
	{
	  case 'N':
	    printf("\n  Give the number of populations\n===> ");
	    fflush(stdout); 
	    FGETS(input, 1024L, stdin);
	    options->numpop = atol(input);
	    printf("\n  NEW NUMBER OF POPULATIONS = %li\n",options->numpop);
	    break;
	case '1':
	  set_prior_menu(THETAPRIOR, options);
	  break;
	case '2':
	  set_prior_menu(MIGPRIOR, options);
	  break;
	case '3':
	      set_prior_menu(RATEPRIOR, options);
	case '4':
	      set_prior_menu(SPECIESTIMEPRIOR, options);
	  break;
	case '5':
	      set_prior_menu(SPECIESSTDPRIOR, options);
	  break;
	default:
	  break;
	}
    } while (uppercase(input[0]) != 'Y');;
}

///
///reads bayes menu choices and sets options accordingly
boolean menuStrategy_bayes(option_fmt * options) 
{
  //long            pop;
  char            input[LINESIZE];
  char           *priorstring;
  //long            count = 0;
#ifdef ZNZ
    char           *extension=NULL;
#endif
  //read user input
  fflush(stdout); FGETS(input, LINESIZE, stdin);
  //if user types yes or YES or YES exit and tell menuStrategy that we are done
    if (strchr("Yy", input[0]))
      return TRUE;

  priorstring = (char *) mycalloc(LINESIZE, sizeof(char));
  //make changes to options
  switch (atoi(input)) {
  case BAYESSTRATEGY:
    //strategy change not allowed
    options->bayes_infer = TRUE;
    options->lchains=1;
    input[0] = '\0';
    break;
  case BAYESOUT:
    printf(" Save posterior distribution (frequency histogram values)? [YES]\n===> ");
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if (uppercase(input[0]) != 'N') 
      {
	options->has_bayesfile = TRUE;
	menu_get_filename("  What is the filename for the posterior distribution (frequency histogram)?",
		     BAYESFILE, options->bayesfilename);
      }
    else
      {
	options->has_bayesfile = FALSE;
      }
    break;
  case BAYESMDIMOUT:
    printf(" Save posterior distribution (all raw parameter values)? [NO]\n===> ");
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if (uppercase(input[0]) == 'Y') 
      {
	options->has_bayesmdimfile = TRUE;
#ifdef ZNZ
	menu_get_filename("What is the filename for the complete posterior distribution (raw parameter values)?\n[if the extension of the file is .gz then the will be compressed]",
		     BAYESMDIMFILE, options->bayesmdimfilename);
	unpad(options->bayesmdimfilename, " ");
	extension = strrchr(options->bayesmdimfilename,'.');
	if(extension!=NULL && !strncmp(extension,".gz",3))
	  {
	    options->use_compressed = 1;
	  }
	else
	  {
	    options->use_compressed = 0;
	  }  	
#else
	menu_get_filename("What is the filename for the complete posterior distribution (raw parameter values)?",
		     BAYESMDIMFILE, options->bayesmdimfilename);
	options->use_compressed = 0;
#endif
	do
	  {
	    printf("Sampling interval for the raw parameter values\n");
	    printf("[Small values can result in a HUGE file!!]\n");
	    printf("[Default is the same as the sampling increment]\n");
	    printf("[Examples: 1 --> saving all parameters every sampling increment]\n");
	    printf("[          2 --> saving every second sampling increment]\n");
	    printf("[        100 --> saving only very hundredth sample]\n");
	    printf("[                if long-inc=50 and here we have here 100 then only]\n");
	    printf("[                every 5,000th step is saved to file]\n");
	    printf("[                with many loci and large run this is still a large file]\n");
	    printf("Current setting: %li\n===> ",options->bayesmdiminterval);
	    fflush(stdout); FGETS(input, LINESIZE, stdin);
	    if(input[0] == '\0' && options->bayesmdiminterval > 0)
	      break;
	    else
	      options->bayesmdiminterval = atol(input);
	  }
	while(options->bayesmdiminterval<1);
      }
    else
      {
	options->has_bayesmdimfile = FALSE;
      }
    break;
  case BAYESBINNING:
    fprintf(stdout,"\n\n\nBinning of posterior has to be set in the prior specification\n\n\n");
#ifdef WINDOWS
    fprintf(stdout,"press return to continue");
    FGETS(input, LINESIZE, stdin);
    input[0]='\0';
#else
    sleep(5);
#endif
    break;
  case BAYESPRETTY:
    printf("Specify the method for plotting the posterior distribution!\n");
    printf("Valid options are:\n");
    printf("ALL     use the range of the prior distribution\n");
    printf("P99     do not plot the values higher than the 99%% percentile for each parameter\n");
    printf("MAXP99  plot up to the maximum 99%% percentile of all parameters.\n");
    printf("TOTAL   plot up to the 100%% percentile of all parameters.\n");
    printf("Current setting: %s\n",
	   (options->bayespretty == PRETTY_P99 ? "P99" :
	   (options->bayespretty == PRETTY_MAX ? "ALL" : 
	   (options->bayespretty == PRETTY_P100 ? "TOTAL" : "MAXP99"))));
    printf("\n\n==> ");
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if (input[0] == '\0')
      break;
    switch(uppercase(input[0]))
      {
      case 'A':
	options->bayespretty = PRETTY_MAX;
	break;
      case 'P':
	options->bayespretty = PRETTY_P99;
	break;
      case 'T':
	options->bayespretty = PRETTY_P100;
	break;
      case 'M':
      default:
	options->bayespretty = PRETTY_P99MAX;
	break;
      }
    break;
  case BAYESFREQ:
    do {
      printf("What is the frequency of the tree updates vs the parameter updates?\n");
      printf("[A frequncy of 1.0 updates only trees, frequency of 0.0 updates only parameters]\n");
      printf("[Suggestion for small problems: 0.5, actual value: %f]\n===> ", options->updateratio);
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      if (input[0] != '\0')  
	  options->updateratio = atof(input);
    } while (options->updateratio > 1);
    break;
  case BAYESPROPOSAL:
    menuProposal(options);
    break;
  case BAYESPRIOR:
    menuPrior(options);
    break;
  case BAYESLCHAINS:
    do {
      printf("  How many Long Chains?\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      options->lchains = atoi(input);
      if (options->lchains < 0)
	printf("  Must be non-negative\n");
    }
    while (options->lchains < 0);
    break;
	   case BAYESSKIP:
	   do {
	     printf("  How many steps (tree changes, parameter changes) to skip?\n===> ");
	     fflush(stdout); FGETS(input, LINESIZE, stdin);
	     options->lincrement = atoi(input);
	     if (options->lincrement <= 0)
	       printf("  Must be positive\n");
	   }
	   while (options->lincrement <= 0);
	   break;
	   case BAYESSAMPLES:
	   do {
	     printf("  How many trees and parameter steps to record?\n===> ");
	     fflush(stdout); FGETS(input, LINESIZE, stdin);
	     options->lsteps = atoi(input);
	     if (options->lsteps <= 0)
	       printf("  Must be a positive integer\n");
	   }
	   while (options->lsteps <= 0);
	   break;
	   case BAYESBURNIN:
	   do {
	     printf("  How many steps to discard? (burn-in)\n");
	     //printf("  For standard fixed burn-in simply give a number, for example 50000\n");
	     //printf("  With Metropolis-Hastings sampling (ignored with Slice-sampling)\n"); 
	     //printf("  append a space and a character (a, e, or t) to invoke autostopping\n");
	     //printf("  for example like this: 100000 a\n");
	     //printf("  a  using variance ratios every 1000 steps to stop: abs(1-var(n)/var(n-1))<0.01\n"); 
	     //printf("  e  using the effective samples size to stop: min(ess(param_i)>%li\n",(long) ESSMINIMUM);
	     //printf("  t  using the average acceptance ratio: %f<avg(acceptance(param_i)<%f\n",(MYREAL) options->autotune-0.05, (MYREAL) options->autotune+0.05);
	     printf("===> ");
	     fflush(stdout); 
	     FGETS(input, LINESIZE, stdin);
	     // a autostop using variance ratios every 1000 steps
	     // e autostop using effective sample size of run > ESSMINIMUM (of worst parameter)
	     // t autostop using when target acceptance ratio is reached (of the average)
	     if(strchr(input,'a'))
	       {
		 options->burnin_autostop = 'a';
	       }
	     else
	       {
		 if(strchr(input,'t'))
		   {
		     options->burnin_autostop = 't';
		     if(!options->has_autotune)
		       {
			 options->has_autotune=TRUE;
			 options->autotune=AUTOTUNEDEFAULT;
		       }
		   }
		 else
		   {
		     if(strchr(input,'e'))
		       {
			 options->burnin_autostop = 'e';
		       }
		     else
		       {
			 options->burnin_autostop = ' ';
		       }
		   }
	       }
	     options->burn_in = atoi(input);
	     if (options->burn_in <= 0)
	       printf("  Must be a positive integer or zero (0)\n");
	   }
	   while (options->burn_in < 0);
	   break;
	   case BAYESREPLICATE:
	   do {
	     printf("  Summarize over multiple chains? [YES | NO] \n===> ");
	     fflush(stdout); FGETS(input, LINESIZE, stdin);
	     if (uppercase(input[0]) == 'Y') 
	       {
		 options->replicate = TRUE;
		 printf("  How many independent runs\n===> ");
		 fflush(stdout); FGETS(input, LINESIZE, stdin);
		 options->replicatenum = ATOL(input);
		 if (options->replicatenum < 1)
		   printf("  Enter a number >= 1\n");
	       } 
	     else 
	       {
		 options->replicate = FALSE;
		 options->replicatenum = 0;
	       }
	   }
	   while (options->replicatenum < 0);
	   break;
	   case BAYESHEAT:
	   printf
	   ("  Heating scheme? < NO | YES | STATIC | BOUNDED_ADAPTIVE >\n===> ");
	   fflush(stdout); FGETS(input, LINESIZE, stdin);
	   switch (tolower(input[0])) {
	   case 'b':
	     options->heating = 1;
	     options->adaptiveheat = BOUNDED;
	     options->heating_interval = 1;
	     menuHeat(options, input);
	     break;
	   case 'y':
	   case 's':
	     options->heating = 1;
	     options->adaptiveheat = NOTADAPTIVE;
	     options->heating_interval = 1;
	     menuHeat(options, input);
	     break;
	   case '\0':
	   case 'n':
	   default:
	     options->heating = 0;
	     options->adaptiveheat = NOTADAPTIVE;
	     options->heated_chains = 1;
	     break;
	   }
	   break;
 /* do we need this
  case BAYESMOVINGSTEPS:
    options->movingsteps = !options->movingsteps;
    if (options->movingsteps) {
      do {
	printf
	  ("  How big should the fraction of new genealogies\n");
	printf
	  ("  of the originally proposed number of samples be?\n[Use this option only when acceptance ratio is small (<0.01)\nRuntime can increase tremendeously]\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	options->acceptfreq = atof(input);
	if (options->acceptfreq < 0 || options->acceptfreq > 1)
	  printf("  Range should be between 0 - 1, and not %f\n",
		 options->acceptfreq);
      }
      while (options->acceptfreq < 0 || options->acceptfreq > 1);
    }
    break;
  case BAYESGELMAN:
    printf("  Convergence statistic options: Pairs , Sum, NO\n");
    printf("  [Currently set to %15s]\n===> ",
	   options->gelman ? (options->gelmanpairs ? "YES:Pairs" : "YES:Summary" ) : " NO");
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if (input[0] == '\0')
      break;
    switch(uppercase(input[0]))
      {
      case 'N':
	options->gelman = FALSE;
	options->gelmanpairs = FALSE;
	break;
      case 'P':
	options->gelman = TRUE;
	options->gelmanpairs = TRUE;
	break;
      default:
	options->gelman = TRUE;
	options->gelmanpairs = FALSE;
	break;
      }
    break;
 */
  case BAYESPRIORALONE:
    printf
      ("  Show only the prior distributions? < NO | YES >\nwith YES no data is analyzed!!!\n===> ");
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    switch (tolower(input[0])) {
    case 'y':
      options->prioralone = TRUE;
      break;
    case 'n':
    default:
      options->prioralone = FALSE;
    }
    break;
  default:
    break;
  }
  myfree(priorstring);
  return FALSE;
}


boolean         menuStrategy_ml(option_fmt * options) 
{
  char            input[LINESIZE];
  //read user input
  fflush(stdout); FGETS(input, LINESIZE, stdin);
  //if user types yes or Yes or YES exit and tell menuStrategy that we are done
  if              (strchr("Yy", input[0]))
    return TRUE;
  
  //make changes to options
  switch          (atoi(input)) {
  case MLSTRATEGY://strategy change not allowed
    options->bayes_infer = TRUE;
    input[0] = '\0';
    break;
  case MLSHORTCHAINS:
    do {
      printf("  How many Short Chains?\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      options->schains = atoi(input);
      if (options->schains < 0)
	printf("  Must be non-negative\n");
    }
    while           (options->schains < 0);
    break;
  case MLSHORTSKIP:
    do {
      printf("  How many trees to skip?\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      options->sincrement = atoi(input);
      if (options->sincrement <= 0)
	printf("  Must be positive\n");
    }
    while (options->sincrement <= 0);
    break;
  case MLSHORTSAMPLES:
    do {
      printf("  How many trees to record?\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      options->ssteps = atoi(input);
      if (options->ssteps <= 0)
	printf("  Must be a positive integer\n");
    }
    while (options->ssteps <= 0);
    break;
  case MLLONGCHAINS:
    do {
      printf("  How many Long Chains?\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      options->lchains = atoi(input);
      if (options->lchains < 0)
	printf("  Must be non-negative\n");
    }
    while (options->lchains < 0);
    break;
  case MLLONGSKIP:
    do {
      printf("  How many trees to skip?\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      options->lincrement = atoi(input);
      if (options->lincrement <= 0)
	printf("  Must be positive\n");
    }
    while (options->lincrement <= 0);
    break;
  case MLLONGSAMPLES:
    do {
      printf("  How many trees to record?\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      options->lsteps = atoi(input);
      if (options->lsteps <= 0)
	printf("  Must be a positive integer\n");
    }
    while (options->lsteps <= 0);
    break;
  case MLBURNIN:
    do {
      printf("  How many genealogies to discard?\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      options->burn_in = atoi(input);
      if (options->burn_in <= 0)
	printf("  Must be a positive integer or zero (0)\n");
    }
    while (options->burn_in < 0);
    break;
  case MLREPLICATE:
    options->replicate = !options->replicate;
    if (options->replicate) {
      do {
	printf("  Summarize over (L)ong chains?\n");
	printf("  or (M)ultiple runs ?\n");
	printf("  [Default is (L)]\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	if (uppercase(input[0]) == 'M') {
	  printf("  How many independent runs\n===> ");
	  fflush(stdout); FGETS(input, LINESIZE, stdin);
	  options->replicatenum = ATOL(input);
	  if (options->lcepsilon < 0)
	    printf("  Enter a number >= 1\n");
	} else
	  options->replicatenum = 0;
      }
      while (options->replicatenum < 0);
    } else {
      options->replicatenum = 0;
    }
    break;
  case MLHEAT:
    printf
      ("  Heating scheme? < NO | YES | ADAPTIVE | BOUNDED_ADAPTIVE>\n===> ");
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    switch (tolower(input[0])) {
    case 'a':
      options->heating = 1;
      options->adaptiveheat = STANDARD;
      options->heating_interval = 1;
      menuHeat(options, input);
      break;
    case 'b':
      options->heating = 1;
      options->adaptiveheat = BOUNDED;
      options->heating_interval = 1;
      menuHeat(options, input);
      break;
    case 'y':
      options->heating = 1;
      options->adaptiveheat = NOTADAPTIVE;
      options->heating_interval = 1;
      menuHeat(options, input);
      break;
    case '\0':
    case 'n':
    default:
      options->heating = 0;
      options->adaptiveheat = NOTADAPTIVE;
      break;
    }
    break;
  case MLMOVINGSTEPS:
    options->movingsteps = !options->movingsteps;
    if (options->movingsteps) {
      do {
	printf
	  ("  How big should the fraction of new genealogies\n");
	printf
	  ("  of the originally proposed number of samples be?\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	options->acceptfreq = atof(input);
	if (options->acceptfreq < 0 || options->acceptfreq > 1)
	  printf("  Range should be between 0 - 1, and not %f\n",
		 options->acceptfreq);
      }
      while (options->acceptfreq < 0 || options->acceptfreq > 1);
    }
    break;
  case MLEPSILON:
    do {
      printf("  Parameter likelihood epsilon?\n[INF is Default]\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      if (uppercase(input[0]) == 'I')
	options->lcepsilon = LONGCHAINEPSILON;
      else
	options->lcepsilon = atof(input);
      if (options->lcepsilon <= 0)
	printf
	  ("  Must be a positive value, be warned: too small values will run the program forever\n");
    }
    while (options->lcepsilon <= 0);
    break;
  case MLGELMAN:
      printf("  Convergence statistic options: Pairs , Sum, NO\n");
      printf("  [Currently set to %10s]\n===> ",
	 options->gelman ? (options->gelmanpairs ? "YES:Pairs" : "YES:Summary" ) : " NO");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      if (input[0] == '\0')
	break;
      switch(uppercase(input[0]))
	{
	case 'N':
	  options->gelman = FALSE;
	  options->gelmanpairs = FALSE;
	  break;
	case 'P':
	  options->gelman = TRUE;
	  options->gelmanpairs = TRUE;
	  break;
	default:
	  options->gelman = TRUE;
	  options->gelmanpairs = FALSE;
	  break;
	}
      break;
  }
  return FALSE;
}

void            menuHeat(option_fmt * options, char *input) 
{
  if (options->heating > 0) {
    printf("Enter the number of different \"heated\" chains.\nMinimum is 4\n===> ");
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    options->heated_chains = atol(input);
    if (options->heated_chains < 4)
      options->heated_chains = HEATED_CHAIN_NUM;
    printf("Enter the interval between swapping trees\nEnter 0 (zero) for NO swapping\n[Current interval is %li]\n===> ",
	    options->heating_interval);
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    //    if (input[0] == '\0')
    //  return;
    if (atol(input) > 0) {
      options->heating_interval = atol(input);
      options->heating = 1;
    } else {
      options->heating = 1;
      options->heating_interval = 0;
      options->heatedswap_off=TRUE;
    }
    read_heatvalues(options);
  }
}

void
display_ml_mcmc(option_fmt * options) {
  char            temp[LINESIZE];

  //  printf("  0   Strategy: %42.42s\n", options->bayes_infer ?
  // "Bayesian Inference" :
  // "Maximum Likelihood");
  printf("  1   Number of short chains to run?                %6ld\n",
	 options->schains);
  if              (options->schains > 0) {
    printf
      ("  2   Short sampling increment?                     %6ld\n",
       options->sincrement);
    printf
      ("  3   Number of recorded genealogies in short chain?%6ld\n",
       options->ssteps);
  }
  printf("  4   Number of long chains to run?                 %6ld\n",
	 options->lchains);
  if (options->lchains > 0) {
    printf
      ("  5   Long sampling increment?                      %6ld\n",
       options->lincrement);
    printf
      ("  6   Number of recorded genealogies in long chain? %6ld\n",
       options->lsteps);
  }
  switch (options->burnin_autostop)
    {
    case 'a':
      sprintf(temp,"(stopping crit: variance) %10li", options->burn_in);
      break;
    case 'e':
      sprintf(temp,"(stopping crit: ESS)      %10li", options->burn_in);
      break;
    case ' ':
    default:
      sprintf(temp,"                          %10li", options->burn_in);
    }
  printf("  7   Burn-in for each chain:%30.30s\n",temp);
  if (!options->replicate)
    printf
      ("  8   Combine chains or runs for estimates?             NO\n");
  else {
    if (options->replicatenum == 0)
      printf
	("  8   Combine chains for estimates?       YES, long chains\n");
    else
      printf
	("  8   Combine chains for estimates?     YES, over %3li runs\n",
	 options->replicatenum);
  }
  if (options->heating == 0)
    printf
      ("  9   Heating:                                          NO\n");
  else {
    if (options->adaptiveheat!=NOTADAPTIVE)
      sprintf(temp, "Adaptive (%s %3li chains, swap interval is %3li)\n",
	      options->adaptiveheat == STANDARD ? "Standard" : "Bounded", 
	      options->heated_chains, options->heating_interval);
    else
      sprintf(temp, "     YES (%3li chains, swap interval is %3li)\n",
	      options->heated_chains, options->heating_interval);
    printf("  9   Heating: %s", temp);
  }
  printf
    ("\n -------------------------------------------------------------\n");
  printf(" Obscure options (consult the documentation on these)\n\n");
  if (options->movingsteps)
    printf
      (" 10   Sample a fraction of %2.2f new genealogies?       YES\n",
       options->acceptfreq);
  else
    printf
      (" 10   Sample at least a fraction of new genealogies?    NO\n");
  if (options->lcepsilon < LONGCHAINEPSILON)
    sprintf(temp, "%17.2f", options->lcepsilon);
  else
    sprintf(temp, "%17.17s", "infinity");
  printf(" 11   Epsilon of parameter likelihood    %s\n", temp);
  printf(" 12   Use Gelman's convergence criterium?    %13s\n",
	 options->gelman ? (options->gelmanpairs ? "YES:Pairs" : "YES:Summary" ) : " NO");
}

void
display_bayes_mcmc(option_fmt * options) 
{
  long count = 0;
  prior_fmt *p;
  char           *outputstring;
  outputstring = (char *) mycalloc(LINESIZE, sizeof(char));

  //printf("  %2i   Strategy:%54.54s\n", BAYESSTRATEGY, options->bayes_infer ?
  // "Bayesian Inference" :
  // "Maximum Likelihood");
  if(options->has_bayesfile)
    sprintf(outputstring,"%s%s",  "YES:", options->bayesfilename);
  else
    sprintf(outputstring,"NO");
  printf(  "  %2i   File for recording posterior distribution?%21s\n", BAYESOUT, outputstring);
  if(options->has_bayesmdimfile)
    sprintf(outputstring,"%s%s",  "YES:", options->bayesmdimfilename);
  else
    sprintf(outputstring,"NO");
  printf(  "  %2i   File for recording all parameter values?  %21s\n", BAYESMDIMOUT, outputstring);
  if(options->has_bayesmdimfile)
    {
      sprintf(outputstring,"[save every %li x sample-increment = %li]", options->bayesmdiminterval, 
	      options->bayesmdiminterval*options->lincrement);
      printf(  "            %58.58s\n", outputstring);
    }
  p = &options->bayes_priors[0];
  count = 0;
  while(p!=NULL)
    {
      count += sprintf(outputstring+count, "%li ",p->bins); 
      p = p->next;
      //      if (count > 20)
      //{
      //  count += sprintf(outputstring+count, "..."); 
      //  break;
      //}
    }
  printf("  %2i   Number of bins of posterior? %35.35s\n", BAYESBINNING, outputstring);
  printf("  %2i   Plotting type of posterior distribution?%23s\n", BAYESPRETTY, 
	 (options->bayespretty == PRETTY_P99 ? "up to 99% percentile" :
	 (options->bayespretty == PRETTY_MAX ? "prior distr. range" :
	 (options->bayespretty == PRETTY_P100 ? "up to ~100% percentile" : "up to maximal 99%"))));
  printf("  %2i   Frequency of tree updates vs. parameter updates?         %6.2f\n",
	 BAYESFREQ, options->updateratio);
  set_proposal(outputstring, options->slice_sampling, !options->bayesmurates);
  printf("  %2i   Proposal distribution?%41.41s\n",
	 BAYESPROPOSAL, outputstring);
  printf("  %2i   Prior distribution submenu\n",
	 BAYESPRIOR);
  printf("  %2i   Number of long chains to run?                            %6ld\n",
	 BAYESLCHAINS, options->lchains);
  printf("  %2i   Sampling increment?                                      %6ld\n",
	 BAYESSKIP, options->lincrement);
  printf("  %2i   Number of recorded steps in chain                  %12ld\n",
	 BAYESSAMPLES, options->lsteps);
  switch (options->burnin_autostop)
    {
    case 'a':
      sprintf(outputstring,"(stopping crit: variance) %10li", options->burn_in);
      break;
    case 'e':
      sprintf(outputstring,"(stopping crit: ESS)      %10li", options->burn_in);
      break;
    case ' ':
    default:
      sprintf(outputstring,"                          %10li", options->burn_in);
    }
  printf("  %2i   Burn-in for each chain:%40.40s\n",
	 BAYESBURNIN,outputstring);

  if(!options->replicate)
    sprintf(outputstring,"NO");
  else
    sprintf(outputstring,"YES (%li independent chains)",options->replicatenum);
  printf("  %2i   Running multiple replicates:  %33.33s\n", BAYESREPLICATE, outputstring);

  if(options->heating == 0)
    printf("  %2i   Heating:            %43.43s\n", BAYESHEAT, "NO");
  else
    printf("  %2i   Heating:                       %10s (%3li parallel chains)\n",
	   BAYESHEAT, (options->adaptiveheat==STANDARD ? "ADAPTIVE" : (options->adaptiveheat==BOUNDED ? "BOUNDED" : "STATIC")) , options->heated_chains);
  /* do we need this see action items further down
  sprintf(outputstring,"%f",options->acceptfreq);
  printf("  %2i   Sampling at least fraction of new genealogies:       %10.10s\n", BAYESMOVINGSTEPS, outputstring);
  sprintf(outputstring,"%s",options->gelman ? (options->gelmanpairs ? "YES:Pairs" : "YES:Summary") : "NO");
  printf("  %2i   Convergence diagnostic for replicates:            %13.13s\n", BAYESGELMAN, outputstring);
  */

  printf("  %2i   Run analysis without data:                        %13.13s\n", BAYESPRIORALONE, options->prioralone ? "YES" : "NO");

  myfree(outputstring);
}


void
how_many_pop(long *numpop) {
  char            input[LINESIZE];
  if              (*numpop == 0) {
    do {
      printf("  How many populations?\n===> ");
      fflush(stdout); FGETS(input, 1024, stdin);
      *numpop = atoi(input);
    }
    while           (*numpop <= 0 && *numpop < 100);
  } else
    printf("  The data set contains %li populations\n", *numpop);
}


void
read_custom_menu_migration(option_fmt * options) 
{
  char            input[LINESIZE];
  long            z = 0, numpop = 0, numpop2;
  printf("  Specify the migration model as an {n x n} matrix\n");
  printf("  Theta values are on the diagonal, migration rates are\n");
  printf("  off-diagonal, spaces (\" \"), \"{\", \"}\", or newlines\n");
  printf("  are allowed, but not necessary.\n");
  printf("\n  Syntax:\n");
  printf("      * = independent parameter\n");
  printf("      0 = (zero) not estimated]\n");
  printf("      c = (constant) not estimated, taken from start-parameter]\n");
  printf("      s = symmetric migration rates (M=m/mu)\n");
  printf("      S = symmetric migration rates (xNm) \n");
  printf("      m = average of each label group [not c, or s]\n"); 
  printf("      d = population split with no immigration\n"); 
  printf("      D = population split with immigration\n"); 
  //  printf("      a,b,c,... = average of each label group [not k, or s]\n"); 
  //  printf("      1,2,3,... = if both migration rates are labeled then combine populations\n"); 
  //printf("      M = average, either ALL thetas and/or ALL migration rates\n");
 
  numpop = options->startparam.numpop;
  printf("\n  NUMBER OF POPULATIONS = %li\n",numpop);
  printf("\n  If the population number is different, type a \"N\" at the prompt!\n");
  printf("\n"); 
  numpop2 = numpop * numpop;
  //printf("\n  You must give %li values or \"N\"\n===> ", numpop2);
  while (z < numpop2) 
    {
      printf("\n  You must give %li more values out of %li\n===> ", numpop2-z, numpop2);
      fflush(stdout); 
      FGETS(input, 1024L, stdin);
      while (input[0]=='N' || input[0]=='n')
	{
	    printf("\n  Give the number of populations\n===> ");
	    fflush(stdout); 
	    FGETS(input, 1024L, stdin);
	    numpop = atol(input);
	    printf("\n  NEW NUMBER OF POPULATIONS = %li\n",numpop);
	    numpop2 = numpop * numpop;
	    printf("\n  You must now give %li more values out of %li\n===> ", numpop2-z, numpop2);
	    fflush(stdout); 
	    FGETS(input, 1024L, stdin);
	}
      read_custom_migration(stdin, options, input, numpop, z);
      z = (long) strlen(options->custm);
    }
}

char * custom_migration_type(long type) 
{
  switch (type) 
    {
    case MATRIX:
      return "Full migration matrix model";
      //break;
    case MATRIX_SYMMETRIC:
      return "Symmetric migration matrix model";
      //break;
    case MATRIX_SAMETHETA:
      return "Full migration matrix model (same Theta)";
      //break;
    case MATRIX_ARBITRARY:
      return "User specified migration matrix model";
      //break;
    case ISLAND:
      return "N-population island model";
      //break;
    case ISLAND_VARTHETA:
      return "N-population island model (variable Theta)";
      //break;
    case STEPSTONE:
      return "Stepping stone model";
      //break;
    case CONTINUUM:
      return "Continuum model";
      //break;
    case NEIGHBOR:
      return "Isolation by distance model";
      //break;
    default:
      return "Illegal migration model";
      //break;
    }
}


void
read_heatvalues(option_fmt * options) {
  long            i, z;
  char           *tmp;
  char            input[LINESIZE];
  MYREAL          diff = 0.;
  MYREAL          x;
  FPRINTF(stdout, " ");
  printf("Enter %li \"temperatures\"\n", options->heated_chains);
  printf("[The coldest temperature, which is the first, has to be 1\n");
  printf(" For example: 1.0 1.5 3.0 1000000.0]\n");
 
  FPRINTF(stdout,
	  "OR give a range of values  [linear increase:      1 - 10]\n");
  FPRINTF(stdout,
	  "                           [exponential increase: 1 @ 10]\n");
  FPRINTF(stdout,
	  "or, most lazily, let me suggest a range [simply type a #]\n");
  FPRINTF(stdout,"@@@@@ For model comparison, the range of temperatures @@@@@\n");
  FPRINTF(stdout,"@@@@@ MUST include a very hot chain (>100000.0)       @@@@@\n");
  printf("===> ");
  fflush(stdout); FGETS(input, LINESIZE, stdin);
  if (strstr(input, "-")) 
    {
      tmp = strstr(input, "-");
      options->heat[options->heated_chains - 1] = fabs(atof(tmp + 1));
      options->heat[0] = 1.0;
      diff = options->heat[options->heated_chains - 1] - options->heat[0];
      diff /= options->heated_chains - 1.;
      for (i = 1; i < options->heated_chains; i++)
	options->heat[i] = options->heat[i - 1] + diff;
    } 
  else 
    {
      if (strstr(input, "#")) 
	{
	  options->heat[0] = 1.0;
	  diff = 1./(options->heated_chains - 1);
	  x = 1 - diff;
	  for (i = 1; i < options->heated_chains-1; i++)
	    {
	      options->heat[i] = 1./ x;
	      x -= diff;
	    }
	  options->heat[i] = 1000000.;
	}
      else
	{
	  if (strstr(input, "@")) 
	    {
	      tmp = strstr(input, "@");
	      options->heat[options->heated_chains - 1] = fabs(atof(tmp + 1));
	      options->heat[0] = 1.0;
	      diff = 2. * (options->heat[options->heated_chains - 1] -
			   options->heat[0]);
	      diff /= pow(2.0, (MYREAL) options->heated_chains) - 2.;
	      for (i = 1; i < options->heated_chains; i++)
		options->heat[i] =
		  options->heat[i - 1] + pow(2.0, i - 1.) * diff;
	    } 
	  else 
	    {
	      z = 0;
	      tmp = strtok(input, " ,");
	      if (tmp != NULL) 
		{
		  options->heat[z++] = atof(tmp);
		}
	      for (;;) 
		{
		  tmp = strtok(NULL, " ,");
		  if (tmp != NULL) 
		    {
		      options->heat[z++] = atof(tmp);
		    } 
		  else 
		    {
		      if (z < options->heated_chains) 
			{
			  printf("%li more values are needed\n===> ",
				  options->heated_chains - z);
			  fflush(stdout); FGETS(input, LINESIZE, stdin);
			  tmp = strtok(input, " ,");
			  if (tmp != NULL) 
			    {
			      options->heat[z++] = atof(tmp);
			    }
			} 
		      else
			break;
		    }
		}
	    }
	}
    }
  printf("Chain         Temperature\n");
  printf("-------------------------\n");
  for (i = options->heated_chains - 1; i >= 0; --i)
    printf("%4li %20.5f\n", i+1, options->heat[i]);
  if(fabs(options->heat[0] - 1.0) > EPSILON)
    {
      printf("The temperature of chain is is NOT 1.0 -- proceed at your own risk\n");
    }
  printf("\n\nIs this correct? [YES or NO]\n===> ");
  fflush(stdout); FGETS(input, LINESIZE, stdin);
  if ('y' != tolower(input[0]))
    read_heatvalues(options);
}

void
menuRandom(MYREAL * param, char type) {
  //long            check;
  //  int retval;
  if(type == 'N')
    {
      printf
	("Specify a MEAN and a STANDARD DEVIATION\nThis will be used to generate a random value from a Normal distribution\n");
    }
  else
    {
      printf
	("Specify a MINIMUM and a MAXIMUM\nThis will be used to generate a random value from a Uniform distribution\n");
    }
#ifdef USE_MYREAL_FLOAT
  scanf("%f%f", &param[0], &param[1]);
#else
 scanf("%lf%lf", &param[0], &param[1]);
#endif
}

//void
//start_tree_method(option_fmt * options) {
//  char            input[LINESIZE];
//  char            compstring[LINESIZE];
//  printf("  Start genealogy is created with\n");
//  printf("      (a)utomatic\n");
//  if              (strchr(SEQUENCETYPES, options->datatype)) {
//    printf("      (u)sertree\n");
//    printf("      (d)istancematrix\n");
//  }
//  printf("      (r)andom\n\n===> ");
//  strcpy(compstring, (strchr(SEQUENCETYPES, options->datatype) ?
//		      "daur" : "dar"));
// do {
//  fflush(stdout); FGETS(input, LINESIZE, stdin);
//}
//while (strchr(compstring, (int) (lowercase(input[0]))) == NULL);
//switch (lowercase(input[0])) {
//case 'u':
//  options->usertree = TRUE;
//  options->dist = FALSE;
//  break;
//case 'd':
//  options->usertree = FALSE;
//  options->dist = TRUE;
//  options->randomtree = FALSE;
//  break;
//case 'a':
//  options->usertree = FALSE;
//  options->dist = FALSE;
//  options->randomtree = FALSE;
//  break;
//case 'r':
//  options->usertree = FALSE;
//  options->dist = FALSE;
//  options->randomtree = TRUE;
//  break;
//}
//}


void display_JC69(option_fmt * options)
{
  options->sequence_model = JC69;
  options->datamodel = JC69;
  options->sequence_model_numparam = 0;
  options->freqsfrom = FALSE; 
}

void  display_K2P(option_fmt * options)
{
  char *input = (char *) mycalloc(LINESIZE,sizeof(char));
  options->sequence_model = K2P;
  options->datamodel = K2P;
  options->sequence_model_numparam = 1;
  fprintf(stdout,"Enter the transition tranversion ratio for the Kimura 2-parameter model\n[Current value = %f]\n===> ", 
	  options->sequence_model_parameters[0]);
  fflush(stdout);
  FGETS(input,LINESIZE,stdin);
  options->sequence_model_parameters[0] = atof(input);
  options->ttratio[0] = options->sequence_model_parameters[0];
  options->freqsfrom = FALSE; 
  myfree(input);
}

void  display_F81(option_fmt * options)
{
  char *input = (char *) mycalloc(LINESIZE,sizeof(char));
  options->sequence_model = F81;
  options->datamodel = F81;
  options->sequence_model_numparam = 0;
    printf("Use empirical base frequencies?  [currently set to: %15s]\n==>", (options->freqsfrom ? "YES" : "NO"));
  fflush(stdout);
  FGETS(input,LINESIZE,stdin);
  //printf("debug: @%s@\n",input);
  if (input[0]=='\0')
    {
      myfree(input);
    }
  else
    {
      if (input[0]=='Y' || input[0]=='y')
	{
	  options->freqsfrom = TRUE;
	}
      else
	{
	  options->freqsfrom = FALSE;
	  initfreqs(&options->freqa, &options->freqc,
		    &options->freqg, &options->freqt);
	}
      myfree(input);
    }
}

void  display_F84(option_fmt * options)
{
  char *input = (char *) mycalloc(LINESIZE,sizeof(char));
  options->sequence_model = F84;
  options->datamodel = F84;
  options->sequence_model_numparam = 3;
  fprintf(stdout,"Enter the transition-tranversion ratio for the Felsenstein 84  model\n[Current value = %f]\n===> ", 
	  options->sequence_model_parameters[0]);
  fflush(stdout);
  FGETS(input,LINESIZE,stdin);
  options->sequence_model_parameters[0] = 1.0 + atof(input); 
  options->sequence_model_parameters[1] =   1.0 + options->sequence_model_parameters[0];
  options->sequence_model_parameters[2] = 1.0; 
  options->ttratio[0] = options->sequence_model_parameters[0];
  options->ttratio[1] = options->sequence_model_parameters[1];
  options->ttratio[2] = options->sequence_model_parameters[2];
  printf("Use empirical base frequencies?  [currently set to: %15s]\n==>\n", (options->freqsfrom ? "YES" : "NO"));
  fflush(stdout);
  FGETS(input,LINESIZE,stdin);
  if (input[0]=='\0')
    {
      myfree(input);
    }
  else
    {
      if (input[0]=='Y' || input[0]=='y')
	{
	  options->freqsfrom = TRUE;
	}
      else
	{
	  options->freqsfrom = FALSE;
	  initfreqs(&options->freqa, &options->freqc,
		    &options->freqg, &options->freqt);
	}
      myfree(input);
    }
}

void  display_HKY(option_fmt * options)
{
  char *input = (char *) mycalloc(LINESIZE,sizeof(char));
  options->sequence_model = HKY;
  options->datamodel = HKY;
  options->sequence_model_numparam = 3;
  fprintf(stdout,"Enter the transition-tranversion ratio for the Hasegawa-Kishino-Yano model\n[Current value = %f]\n===> ", 
	  options->sequence_model_parameters[0]);
  fflush(stdout);
  FGETS(input,LINESIZE,stdin);
  options->sequence_model_parameters[0] = atof(input); 
  options->sequence_model_parameters[1] = options->sequence_model_parameters[0]; 
  options->sequence_model_parameters[2] = 1.0; 
  options->ttratio[0] = options->sequence_model_parameters[0];
  options->ttratio[1] = options->sequence_model_parameters[1];
  options->ttratio[2] = options->sequence_model_parameters[2];
  printf("Use empirical base frequencies?  [currently set to: %15s]\n==>\n", (options->freqsfrom ? "YES" : "NO"));
  fflush(stdout);
  FGETS(input,LINESIZE,stdin);
  if (input[0]=='\0')
    {
      myfree(input);
    }
  else
    {
      if (input[0]=='Y' || input[0]=='y')
	{
	  options->freqsfrom = TRUE;
	}
      else
	{
	  options->freqsfrom = FALSE;
	  initfreqs(&options->freqa, &options->freqc,
		    &options->freqg, &options->freqt);
	}
      myfree(input);
    }
}

void   display_TN(option_fmt * options)
{
  char *input = (char *) mycalloc(LINESIZE,sizeof(char));
  options->sequence_model = TN;
  options->datamodel = TN;
  options->sequence_model_numparam = 3;
  fprintf(stdout,"Enter the transition-tranversion ratio and the ratio of transitions for the Tamura-Nei model\n[Current values: Tv/Tt=%f, Ti=%f]\n===> ", 
	  options->sequence_model_parameters[0],options->sequence_model_parameters[2]);
  fflush(stdout);
  FGETS(input,LINESIZE,stdin);
  sscanf(input,"%lf%lf",  &options->sequence_model_parameters[0],  &options->sequence_model_parameters[2]); 
  options->sequence_model_parameters[1] = 1.0; 
  options->ttratio[0] = options->sequence_model_parameters[0];
  options->ttratio[1] = options->sequence_model_parameters[1];
  options->ttratio[2] = options->sequence_model_parameters[2];
  printf("Use empirical base frequencies?  [currently set to: %15s]\n==>\n", (options->freqsfrom ? "YES" : "NO"));
  fflush(stdout);
  FGETS(input,LINESIZE,stdin);
  if (input[0]=='\0')
    {
      myfree(input);
    }
  else
    {
      if (input[0]=='Y' || input[0]=='y')
	{
	  options->freqsfrom = TRUE;
	}
      else
	{
	  options->freqsfrom = FALSE;
	  initfreqs(&options->freqa, &options->freqc,
		    &options->freqg, &options->freqt);
	}
      myfree(input);
    }
}

///
/// delivers string for particular sequence submodel
char * menu_sequence_submodeltype(int type)
{
  switch(type)
    {
    case JC69:
      return "         Jukes Cantor";
    case K2P:
      return "   Kimura 2-parameter";
    case F81:
      return "       Felsenstein 81";
    case F84:
      return "       Felsenstein 84";
    case HKY:
      return "Hasegawa-Kishino-Yano";
    case TN:
      return "           Tamura-Nei";
    case GTR:
      return "      Not implemented";
    default:
      return "      Not implemented";
    }
}


///
/// delivers string for particular sequence submodel
void sequence_modelparameters(int type, char * text, option_fmt *options)
{
  switch(type)
    {
    case JC69:
      strcpy(text,"None");
      break;
    case K2P:
      sprintf(text,"Tv/Ti=%f",options->sequence_model_parameters[0]);
      break;
    case F84:
      sprintf(text,"Tv/Ti=%f",options->sequence_model_parameters[0]);
      break;
    case F81:
      strcpy(text,"None");
      break;
    case HKY:
      sprintf(text,"Tv/Ti=%f",options->sequence_model_parameters[0]);
      break;
    case TN:
      sprintf(text,"Tv/Ti=%f, Ti-rat=%f",options->sequence_model_parameters[0], options->sequence_model_parameters[1]);
      break;
    case GTR:
    default:
      strcpy(text, "Not implemented");
    }
}

void menu_haplotyping(option_fmt *options)
{
  char text[LINESIZE];
  char input[LINESIZE];
  if(options->haplotyping)
    sprintf(text, "YES:%s",options->haplotyping_report ? "reporting haplotypes" : "no reporting of haplotype" );
  else
    sprintf(text, "NO");
  
  fprintf(stdout,"Turn on haplotyping [default is set to:%s]\n[use NO or YES]===> ",text);
  fflush(stdout); FGETS(input, LINESIZE, stdin);
  if (input[0] == '\0')
    return;      
  if(strchr("Yy", input[0]))
    {
#ifdef NEWVERSION
      options->haplotyping = TRUE;
      fprintf(stdout, "Report haplotypes [Default: %s]\n[use NO or YES]\n",text+4); 
      fprintf(stdout, "   [Remember: data with larger numbers of unphased blocks will\n");
      fprintf(stdout, "   demand huge computer resources I suggest not to report the\n");
      fprintf(stdout, "   haplotype states.]\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      if (input[0] == '\0')
	return;      
      if(strchr("Yy", input[0]))
	{
	  options->haplotyping_report = TRUE;
	}
      else
	{
	  options->haplotyping_report = FALSE;
	}
#else
      fprintf(stdout, "\nHaplotyping is not turned on in this version\n");
#endif
    }
  else
    {
      options->haplotyping = FALSE;
      options->haplotyping_report = FALSE;
    }
}


/// \brief Menu for sequence submodels
/// Menu sequence submodels
void menu_sequence_submodel(option_fmt *options)
{
  char input[LINESIZE];
  char text[LINESIZE];
  char text2[LINESIZE];
  //char val;
  long numval;
  //double tune = 0.;
  //double pinc = 0.5;
  //long count = 0;
  
  do {
    text[0]='\0';
    sprintf(text,"Current sequence model: %s\n", menu_sequence_submodeltype(options->sequence_model));
    sequence_modelparameters(options->sequence_model, text2, options);
    printf("\nChoose a sequence model from the list?\n");
    printf(" 1  Jukes-Cantor model\n");
    printf(" 2  Kimura 2-parameter model\n");
    printf(" 3  Felsenstein 1981 (F81) model\n");
    printf(" 4  Felsenstein 1984 (F84) model\n");
    printf(" 5  Hasegawa-Kishino-Yano model\n");
    printf(" 6  Tamura-Nei model\n");
    printf(" Y  leave this submenu\n\n");
    printf("\n===> ");
    fflush(stdout); 
    FGETS(input, LINESIZE, stdin);
    if (input[0] == '\0')
      break;
    if(uppercase(input[0])=='Y')
      break;
    sscanf(input, "%li", &numval);
    switch(numval-1)
      {
      case JC69:
	display_JC69(options);
	break;
      case K2P:
	display_K2P(options);
	break;
      case F81:
	display_F81(options);
	break;
      case F84:
	display_F84(options);
	break;
      case HKY:
	display_HKY(options);
	break;
      case TN:
	display_TN(options);
	break;
      default:
	printf("Wrong model, try again\n");
      }
  } while (uppercase(input[0])!='Y');
}



///
/// delivers string for particular msat submodel
char * msat_submodeltype(int type)
{
  switch(type)
    {
    case MULTISTEP:
      return "Multistep method";
    case SINGLESTEP:
      return "Singlestep method";
    default:
      return "Singlestep method";
    }
}

/// \brief Menu for msat submodel
/// Menu for microsatellite submodel (using Watkins methods)
void menu_msat_submodel(option_fmt *options)
{
  char input[LINESIZE];
  char val;
  double tune = 0.;
  double pinc = 0.5;
  long count = 0;
  printf("Microsatellite data can be analyzed using\n");
  printf("(a) infinite allele model\n");
  printf("(b) Brownian motion model (choose this one!)\n");
  printf("(s) Single-step mutation model\n");
  printf("(m) Multi-step mutation model (tune=%f, p_increase=%f)\n",options->msat_tuning[0], options->msat_tuning[1]);
  printf("             Single Step model is the standard stepwise mutation model\n");
  printf("             The multistep model is between the infinite model and and the\n");
  printf("             singlestep model using two additional parameters:\n");
  printf("              - chance of increasing repeat length\n");
  printf("              - tuning between single step (tune=0) and infinite allele (tune=1)\n");
  printf("\n===> ");
  fflush(stdout); 
  FGETS(input, LINESIZE, stdin);
  val = uppercase(input[0]);
  switch(val)
    {
    case 'A':
      options->datatype='a';
      break;
    case 'B':
      options->datatype = 'b';
      break;
    case  'S':
      options->datatype = 'm';
      options->msat_option = SINGLESTEP;
      break;
    case 'M':
      options->datatype = 'm';
      options->msat_option = MULTISTEP;
      do
	{
	  printf("MULTISTEP model\nGive two numbers separated by a space:\nFor the tune (range 0 to 1) and\nfor the probability of repeat number increase (0 to 2/3)\n===> ");
	  fflush(stdout); FGETS(input, LINESIZE, stdin);
	  count = sscanf(input, "%lf%lf", &tune, &pinc);
	  options->msat_tuning[0] = (MYREAL) tune;
	  options->msat_tuning[1] = (MYREAL) pinc;
	} while (count != 2 && count > 0);
      //count = 1;
      break;
    default:
      break;
    }
}


void
start_data_method(option_fmt * options) {
  char            input[LINESIZE];
  do {
    //printf("  (a)llele model\n");
    //printf("  (m)icrosatellite model [SLOW! Ladder model: %s]\n", msat_submodeltype(options->msat_option));
    //printf("  (b)rownian microsatellite model [Brownian motion model]\n");
    //printf("  (s)equence model\n");
    //printf("  (n)ucleotide polymorphism (SNP)\n");
    printf("  (a)llele data\n");
    printf("  (m)icrosatellite datamodel\n");
    printf("  (b)icrosatellite datamodel (Brownian motion)\n");
    printf("  (s)equence data (DNA or RNA sequences)\n");
    printf("  (n)ucleotide polymorphism (SNP)\n");
    printf("  (h)apmap data nucleotide polymorphism (SNP)\n");
    //printf("  (u)nlinked nucleotide polymorphism (SNP) using a PANEL\n");
    //printf("  (g)enealogy summaries\n\n===> ");
    printf("\n  Are the settings correct? (current datatype is %c)\n", options->datatype);
    printf("  (Type Y or the number of the entry to change)\n");
    printf("\n===> ");
    fflush(stdout); 
    FGETS(input, LINESIZE, stdin);
  }
  while           (strchr("abmnshy", (int) (lowercase(input[0]))) == NULL);
  //while           (strchr("abmnsgh", (int) (lowercase(input[0]))) == NULL)
    if (lowercase(input[0])!='y')
      {
	options->datatype = input[0];
	switch(options->datatype)
	  {
	  case 'm':
	    menu_msat_submodel(options);
	    break;
	  case 's':
	  case 'n':
	    break;      
	  case 'a':
	  case 'b':
	    break;
	  default:
	    break;
	  }
      }
  if (!strchr(SEQUENCETYPES, options->datatype))
    options->usertree = FALSE;
}

long            get_prior(char *input) {
  switch (toupper(input[0])) {
    //case 'S':
    //return SLICE;
    //break;
  case 'M':
    return MULTPRIOR;
    //break;
  case 'E':
    return EXPPRIOR;
    //break;
  case 'W':
    return WEXPPRIOR;
    //break;
  case 'U':
    return UNIFORMPRIOR;
  default:
    return UNIFORMPRIOR;
    //break;
  }
}

///
/// prints type of proposal in the menu
void set_proposal(char *output, boolean *proposal, boolean without_rate)
{
  int i;
  int l=0;
  const char types[][LINESIZE]={"Theta","Mig","Rate"};
  
  for(i=THETAPRIOR; i <= RATEPRIOR; i++)
    {
      if(without_rate && i==RATEPRIOR)
	continue;
      switch (proposal[i]) {
      case TRUE:
	l+=sprintf(output+l, " %s:Slice", types[i]);
	break;
      case FALSE:
      default:
	l+=sprintf(output+l, " %s:MH", types[i]);
	break;
	//      default:
	//l+=sprintf(output+l, " %s:?", types[i]);
	//break;
      }
    }
}


void set_localities_string(char *loc, option_fmt *options)
{
  long z;
  long i;
  if(options->newpops_numalloc>1)
    {
      z=0;
      for(i=0; i<options->newpops_numalloc; i++)
	{
	  z += sprintf(loc+z,"%li,",options->newpops[i]);
	}
    }
  else
    {
      sprintf(loc,"default");
    }
}
