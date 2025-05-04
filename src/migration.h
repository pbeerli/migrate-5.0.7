/*! \file migration.h*/
#ifndef MIGRATION_HEADER /*migration.h */
#define MIGRATION_HEADER

#include "definitions.h"
/*-----------------------------------------------------------------
  Bayesian inference of population genetic forces: drift, migraiton, divergence
  allowing for the n-coalescent, the f-coalescent, and the BSC-coalescent
 
  Peter Beerli
  Department of Scientific Computing
  Florida State University
  Tallahassee FL 32306-4120
  beerli@fsu.edu
 
  Copyright 2002 Peter Beerli and Joseph Felsenstein, Seattle WA
  Copyright 2003 Peter Beerli, Tallahassee FL
  Copyright 2017 Peter Beerli, Tallahassee FL

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

 $Id: migration.h 2165 2013-08-24 16:33:45Z beerli $
  *----------------------------------------------------------------
*/
/**
 * Definitions used for the program migrate
 * this file holds structures and typedefs
 * and defines
 */

/* Altivec support [NOT WORKING YET]*/
#ifdef ALTIVEC
#include "altivec.h"
#endif /*ALTIVEC*/

#ifdef ZNZ
#include "znzlib.h"
#endif


typedef char allele_type[DEFAULT_ALLELENMLENGTH]; //!< allele type string

typedef long twin_fmt[2];
typedef long quad_fmt[4];
typedef MYREAL pair[2];
typedef long longpair[2];
typedef float tetra[9];
typedef float duo[2];
#ifdef ALTIVEC
typedef FloatVec sitelike;
#else
#ifdef VARMUT
//  best model for condlike? chunks of numcategs * (condlike_site * numsites)
//  accessing a single 1D array?
//
typedef MYREAL * sitelike; //this allows to tread each site as a different datatype with its own mutation model
#else
#ifdef GAP
typedef MYREAL sitelike[5];
#else
typedef MYREAL sitelike[4];
#endif /*GAP*/
#endif /*VARMUT*/
#endif /*ALTIVEC*/

typedef sitelike *ratelike;
typedef ratelike *phenotype;
typedef MYREAL contribarr[MAXCATEGS]; //!< \todo conditional likelihoods should change to 1-D vector
typedef short seqval[MAXCATEGS]; 

typedef char allele_fmt[LINESIZE];
//typedef char allele_fmt[DEFAULT_ALLELENMLENGTH];

/// tipdate format
typedef struct tipdate_fmt 
{
  MYREAL date;
  long   id;
  char * name;
} 
tipdate_fmt;

/// records when the automatic stop of burnin happened
typedef struct _burnin_record_fmt
{
  long locus;
  long replicate;
  long stopstep;
  MYREAL ess;
  MYREAL accept;
  MYREAL variance;
  MYREAL oldvariance;
  long worker;
} burnin_record_fmt;


#ifdef BEAGLE
#include <beagle.h>
typedef struct _beagle_fmt {
  boolean instance;
  int *instance_handle; // indicator to the workorder for the GPU/CPU likelihood calculation
  int numallocbranches;
  int numbranches;
  double *branch_lengths; // branchlengths for (for each edge leading to a child)
  int *branch_indices;    // indicator for transition matices
  int numallocoperations;
  int numoperations;
  BeagleOperation * operations;  // holds the operation for each parent,leftchild,lefttransitionmatrix,rightchild,rightt.. 
  int *scalingfactorsindices;
  int scalingfactorscount;
  boolean scaling;
  boolean ievectrans;
  boolean logscalers;
  boolean eigencomplex;
  boolean dynamicscaling;
  boolean autoscaling;
  boolean requireDoublePrecision;
  boolean requireSSE;
  int scalingIndex1;
  int scalingIndex2;
} beagle_fmt;
#endif /*BEAGLE*/


/// likelihood ratio test data structure
typedef struct lr_data_fmt
{
    short type;   //!< type of test \todo check whether this is obsolete
    long elem;    //!< number of elements in value1 and value2 
    char *value1; //!< values to test using a mix of numbers and m,c, etc
    char *value2; //!< not used, reserved fo pairwise tests 
    char *connect;//!< custom migration matrix for value1 [no user interaction]
}
lr_data_fmt;

/// likelihood ratio test structure
typedef struct lratio_fmt
{
    long alloccounter; //!< number of elements allocated in data
    long counter;      //!< number of elements in data
    lr_data_fmt *data; //!< contains the LRT data
}
lratio_fmt;

/* used in the tree structure */
/// defines a node in the tree structure

/// temporary structure to hold precomputed values for the conditional likelihood calculation
typedef struct valrec
{
    MYREAL rat, ratxi, ratxv, zz, z1, y1, ww1, zz1, ww2, zz2, z1zz, z1yy, xiz1,
    xiy1xv, ww1zz1, vv1zz1, ww2zz2, vv2zz2;
}
valrec;

/// 2-D array to pointers to valrec structures 
typedef valrec ***tbl_fmt;

/// union to hold the conditional likelihood data in tree nodes
typedef union xarray_fmt
{
    MYREAL *a;  //!< used for allelic types
    phenotype s;//!< used for sequence data types 
}
xarray_fmt;

#ifdef UEP
/// union to hold unique event polymorphism conditional likelihood data 
typedef union ueparray_fmt
{
    MYREAL *a; //!< holds probabilities for many state data \todo check the usage of UEP union
    pair *s;   //!< holds probabilities for 0/1 data
}
ueparray_fmt;
#endif
#ifdef VARMUT
typedef struct varmut_fmt
{
  //sequence model
  MYREAL freqa, freqt, freqg, freqc;      //!< base frequencies
  MYREAL freqr, freqy;                    //!< frequency of purins
  MYREAL freqar, freqcy, freqgr, freqty;  //!< frequency of transitions
  MYREAL aa, bb;                          //
  long oldsite;                           //!< number of site patterns (raw, with SNP adjustment
  long endsite;                           //!< number of site patterns 
  MYREAL xi;                              //!< ratio of transitions
  MYREAL xv;                              //!< ratio of transversions
  MYREAL ttratio;                         //!< ration between transitions and transversions
  MYREAL fracchange;                      //
  //msat model

  //allele model

} varmut_fmt;
#endif /*VARMUT*/

typedef struct _unassigned_fmt
{
  char *key; //individual name
  long index;
  struct _unassigned_fmt * next;
  MYREAL *probloc; //loci * numpop
  struct _node **thenodes; //working node in tree of particular locus
  long ploidy;
  long lastpop;
  long lastlocus;
  long lastreplicate;
  long *accept;
  long *trials;
} unassigned_fmt;

typedef struct _hash_fmt 
{
  char *key;
  struct _hash_fmt * next;
  long value;
  size_t numhash;
} hash_fmt;


/* used in the tree structure */
/// defines a node in the tree structure
typedef struct _node
{
    struct _node *next, *back;
    boolean tip;
    char type;
    long number;
    long pop;
    long actualpop;
  long truepop;
    long id;
  long bid;
  boolean visited;
    xarray_fmt *x;
  //for asignment:
  double *freqs;
#ifdef UEP

    int *uep;
    ueparray_fmt ux;
    MYREAL uepscale;
#endif
    MYREAL *s;
  MYREAL *probloc;
    MYREAL **scale;
    char *nayme;
  char *truename;
    boolean top;
    boolean dirty;
  MYREAL v;
  MYREAL tyme;
    MYREAL length;
  char *sequence;
}
node;

/// sequence data structure
typedef struct seqmodel_fmt
{
  // MYREAL freqa, freqt, freqg, freqc;      //!< base frequencies
  // MYREAL freqr, freqy;                    //!< frequency of purins
  // MYREAL freqar, freqcy, freqgr, freqty;  //!< frequency of transitions
  // BASEFREQLENGTH is defined in definitions.h
  MYREAL basefrequencies[BASEFREQLENGTH];             //!< base frequencies (a,c,g,t), gap, purins (r,y), transitions (ar,cy,gr,ty)
  // rate matrix that should be in basefrequencies [g=nucleotide g, G=gap]
  // pa pc pg pt pG
  // rac rag rat raG rcg rct rcG rgt rgG rtG
  // model vector (Huelsenbeck and Alfaro) augmented to allow the inclusion of gaps but using an additional class 0 --> set to zero.
  //          abcdefghij
  // no gaps:
  // minimal: 1110110100
  // JC/F81:  1110110100
  // F84:     1110110100
  // maximal: 1230450600
  // full model:
  // minimal: 1111111111
  // JC/F81:  1110110100
  // maximal: 123456789a
  // 
  MYREAL aa, bb;                          //
  long oldsite;                           //!< number of site patterns (raw, with SNP adjustment
  long endsite;                           //!< number of site patterns 
  MYREAL xi;                              //!< ratio of transitions
  MYREAL xv;                              //!< ratio of transversions
  MYREAL ttratio;                         //!< ration between transitions and transversions
  MYREAL fracchange;                      //
  long *sites;                            //
  long *alias;                            //
  long *ally;                             //
  long *category;                         //
  short *weight;                          //
  long weightsum;                         //
  double *aliasweight;                      //
  long *location;                         //
  long addon;                             //
  boolean *links;                         //!< information about link status of locus
  double *savealiasweight;                 //!< keeps a savecopy of the aliasweights
  boolean done;
}
seqmodel_fmt;


// definition of dists are in definitions.h
typedef struct _species_fmt {
  char type;
  int  dist;
  long id;
  double mu;
  double min;
  double max;
  double sigma;
  double sigmamin;
  double sigmamax;
  long   from;
  long   to;
  long paramindex_mu;
  long paramindex_sigma;
  long size;
  long allocsize;
  double *data;
} species_fmt;

enum siteclass_enum { SITEWORD, SITECHARACTER };

//defines the mutations models used in the extended "Ican do everything mutation model string thing"
typedef struct _mutationmodel {
  boolean from_infile;
  boolean finished;
  char datatype; //specifices the model
  enum siteclass_enum dataclass;
  int model;
  long numpatterns; // number of site patterns
  long numsites;
  long startsite;
  long numstates; // number of states in model: DNA=4, DNA+gap=5, msat>2
  long numsiterates;
  double lambda;
  double *parameters; //model parameters
  double *basefreqs;
  double *siterates;  // site rates and site probs
  double *siteprobs;
  // transition matrix material
  double *eigenvectormatrix;
  double *inverseeigenvectormatrix;
  double *eigenvalues;
  // pattern material
  long *alias;                            //
  long *ally;                             //
  short *weight;                          //
  long weightsum;                         //
  double *aliasweight;                      //
  long *location;                         //
  double *savealiasweight;
  contribarr *contribution;
  long *category;
  valrec ***tbl;
  MYREAL xi;
  MYREAL xv;
  MYREAL ttratio;
  MYREAL fracchange;
  // allelic data
  MYREAL freq;
  MYREAL freqlast;
  long maxalleles;
  // brownian
  MYREAL browniandefault;
  // msats
  long micro_threshold;
  long microrange;
  long microstart;
  MYREAL **steps;
  // snps
  long addon;
  // hard coded categories
  long numcategs;
  MYREAL *rate; // for categs not relative site rate variation -- needs change of name
  MYREAL scoring_error[4];
  boolean estimateseqerror;
  boolean seqerrorcombined;
  boolean scaling;
#ifdef BEAGLE
  long numallocpartials;
  double *partials;
  long numpartials;
  //  long scalingfactorscount;
  //int *scalingfactorsindices;
#endif
} mutationmodel_fmt;


typedef char * site_fmt; 

typedef struct _region_fmt 
{
  long startlocus;
  long endlocus;
} region_fmt;

typedef struct _haplotyping_orderform 
{
  char **pick;
  char *key;
} haplotyping_orderform_fmt;  

typedef struct _haplocount 
{
  char *haplotype;
  long count;
} haplocount_fmt;  

// ready to delete ? May 11 2015
/// define the holding structure haplotype individual names in the data object 
//typedef struct _individualHOLD
//{
//  
//} individualHOLD;

/// define the structure to hold the individual name and haplotype information
typedef struct _individualDB 
{
  long id;
  char name[100];
  long *ind; 
  long nodenum;
  node ** nodep;
  long region;
  int *difference;
  long checksum;
  int *states1;
  int *states2;
  long numcounts;
  hash_fmt *hash;
  size_t numhash;
  float total1;
  float total2;
  long numstates;
  long targeted;
} individualDB_fmt;


typedef struct _fastadata
{
  char * indname;
  char * fasta;
  char * fastq;
} fastadata_fmt;


/// defines the data structure read from infile
typedef struct _data
{
    FILE *infile;           //!< data file
    FILE *utreefile;        //!< user tree input file
    FILE *weightfile;       //!< site weighting file
    FILE *catfile;          //!< site categories file
    FILE *sumfile;          //!< intermediate output summary file
    FILE *distfile;         //!< distance file among individuals, instead of treefile
    FILE *geofile;          //!< geographic distance among sampling locations
    FILE *divfile;          //!< fixed divergence time among populations
    FILE *datefile;         //!< sample date for each individual
#ifdef UEP
    FILE *uepfile;          //!< unique event polymorphism file
    int **uep;              //!< uep data holder
    long uepsites;          //!< number of UEP sites
    MYREAL uepfreq0;        //!< base frequency for 0
    MYREAL uepfreq1;        //!< base frequency for  1
#endif
    site_fmt *****yy;           //!< data holder: each site is a string minimal 1 character.
#ifdef MPIDATAONDEMAND
  site_fmt **datapart;      // used in data on demand  
#endif
    MYREAL *geo;            //!< geographic distance data holder
    MYREAL *lgeo;           //!< log of geo
    MYREAL **ogeo;          //!< original of geo as read from the geofile
    char ***allele;    //
    long *maxalleles;       //!< maximal number of alleles in dataset
    boolean *skiploci;      //!< number of loci to skip -- loci with no data
    char **popnames;        //!< array of the population names
    long **numind;          //!< array of the number of individuals
    long **numalleles;      //
  tipdate_fmt ***sampledates;
  MYREAL maxsampledate;
    char ****indnames;       //!< array of individual names per population
    long numpop;            //!< number of populations
    long loci;              //!< number of loci
  long allsubloci;
  long *position;           //!< position on chromosome [what to do with mutiple chromosomes?]
    seqmodel_fmt **seq;      //!< sequence data model    
    MYREAL freq;            //
    MYREAL freqlast;        //
    char dlm;               //!< delimiter for allelic data
    boolean hasghost;       //!< whether there is a population with no data or not \todo check if hasghost is still needed
  long ***shuffled;         //!< holds the random subset of individuals
  boolean oneliner;         //!< if a () is used in the sites lines it is assumed that all loci are on one line
  char *datatype;
  long numdatatypealloc;
  char *locitypes;          //!< holds all datatypes for each sublocus, options->datatype is the default
  long numlocitypesalloc;   //!< allocation counter for locitypes 
  char **locusname;

  long *subloci;            //!< holds the number of strictly linked "loci" that make up the locus for migrate
  long *totalsites;         //!< holds the number all sites within a locus, these sites can have different meaning 
  boolean oldsyntax;        //!< either single type data or combinations of data (new)
  long numsublocialloc;
#ifdef MPI
  long *sublocistarts;               // pointer to world(cold)->sublocistarts
  mutationmodel_fmt *mutationmodels; // pointer to world(cold)->mutationmodels
#endif
  long *numrepeatnumbers;
  long **repeatnumbers;
  long *repeatlength;
  boolean has_repeats;
  long numregions;
  region_fmt * regions;

  long *numindividuals;
  individualDB_fmt **individuals;
  boolean haplotyping;
  boolean haplotyping_report;
  haplotyping_orderform_fmt *haplotyping_list;
  long numhaplotyping;
  MYREAL *locusweight; //invariant loci treatment
  boolean has_fastq;
  FILE *fastqinfile;
  fastadata_fmt *fastadata;
}
data_fmt;



/// holds histogram and associated statistics for each locus
typedef struct _bayeshistogram 
{
    long binsum; //sum of all elements of all parameters and bins
  boolean *smoothed;
    long *bins;        // number of bins for each parameter
    double *results;   // contains potentially smoothed histogram, size is bins*numparam
    double *results2;   // contains unsmoothed histogram, size is bins*numparam
    char *set95; // holds in/out HPD set for 95% probability
    char *set50; // holds in/out HPD set for 50% probability [this is attached to cred95]
    // on a per parameter basis
    // structure has a data storage vectors and the following are all pointers into it
    long numparam;    // number of parameters
    MYREAL *datastore; // data storage, size is numparam*10
    // pointers into data storage
    MYREAL *minima;    // contains minimal values for each parameter
    MYREAL *maxima;    // contains maximal values for each parameter
    MYREAL *adjmaxima;// holds maxima values used in histogram [are smaller than maxima]
    MYREAL *cred50l;    // holds 50%-credibility margins (<all lower values>, 
    MYREAL *cred50u;   // <all high values>)
    MYREAL *cred95l;    // holds 95%-credibility margins (<all lower values>)
    MYREAL *cred95u;   // <all high values>)
    MYREAL *modes;    // holds 95%-credibility margins (<all lower values>, <all high values>)
    MYREAL *medians;
    MYREAL *means;
    MYREAL *stds;
    MYREAL **covariance; // holds the covariance structure per locus
  long n; //holds the number of samples taken, this seem to be missing for the reading from mdimfile, but available otherwise
} bayeshistogram_fmt;

typedef struct _hyper_fmt
{
  MYREAL alpha;
  MYREAL alphastd;
  long alphan;
  MYREAL mean;
  MYREAL meanstd;
  long meann;
  long count;
} hyper_fmt; 

typedef struct _rules
{
  long **indices;
  long *indicesnum;
  char **custm2rules;
  long rulesnum;
} rules_fmt;


/// holds data for Bayesian approach
typedef struct _bayes
{
  long count;
  long hypercount;
  long hyperinterval;
  boolean hyperprior;
  MYREAL hyperfactormean;
  MYREAL hyperfactoralpha;
  hyper_fmt *hyperp; //the name hyper is ineligible with visual studio
  long numpop2;
  long mapsize;
  longpair *map; //mapping of parameter using custom migration matrix
  MYREAL *datastore; //data store for
  // pointers into datastore
  MYREAL *priormean; // holds the means of priors
  MYREAL *priorstd; // holds the standard deviation of priors [not always filled]
  MYREAL *prioralpha; // holds alpha of prior distribution [not always filled]
  MYREAL *delta; // change scale [how big can this be
  MYREAL *minparam;// minimal allowed parameter value
  MYREAL *maxparam; // maximal allowed parameter value
  MYREAL *meanparam; // mean values for expprior
  MYREAL *alphaparam; // alpha value for gamma
  MYREAL *alphaorigparam; // alpha value for gamma (original [used for Hyperp])
  MYREAL *betaparam; // beta value for gamma
  // record all changes of parameters
  MYREAL *params;  // save for parameter vectors
  long allocparams; //number of allocated parameter vectors
  long numparams;  //number of saved parameter vectors
  long paramnum;  // which param one is working with
  MYREAL starttime; // start time that was used to change the lineage on tree
  MYREAL stoptime; // stop time that was used to change the lineage on the tree
  MYREAL oldval;    //holds prob(G|param)
  long *datastore2; // datastore for long integer variables
  // pointers into datastore2
  long *accept;  //holds numbers of accepted changes for each parameter AND the tree
  long *trials; //holds how many time was tried to change a parameter
  // holds histogram for each locus and summary statistics
  bayeshistogram_fmt *histogram;
  MYREAL *deltahist; // delta for all histograms per parameters [important for multilocus case 
  char * custm2; //pointer to the world->options->custm2;
  boolean mu;
  long mdiminterval ; // interval at which the whole parameter list is printed to file
  int prettyhist;
  MYREAL *scaling_factors;
  MYREAL *histtotal;
  MYREAL maxmaxvala;
  // speed improvement to beat the slow malloc for large string allocation
  // measures size of first string and uses that currently size * 2;
  boolean has_linesize;
  long linesize;
  long progresslinesize;
  long *mdimfilecount;
  MYREAL *priors;
  rules_fmt *rules;
}
bayes_fmt;


/// storage for bayesian prior related parameters

typedef float (randomptr)(float []);
typedef float (cdfptr)(float, float []);

typedef struct prior_fmt
{
  struct prior_fmt * next;
  int kind; // gamma, exp, uni, etc
  int type; // theta, mig, species, rate
  char ptypename[13];
  long number;
  MYREAL min;
  MYREAL mean;
  MYREAL std;
  MYREAL max;
  MYREAL delta;
  MYREAL alpha;
  MYREAL beta;
  MYREAL updatefreq;
  long   bins;
  long from;
  long to;
  float v[4]; // holds lower upper ,mean ,std except for gamma: lower, upper alpha, beta
  randomptr *random;
  cdfptr *cdf;
  //long from;  //these holds the parameter from value e.g. M_i->j i, 
  //long to;    //these holds the parameter to value e.g. M_i->j j,
  // for theta and rates from and to are the same
} prior_fmt;

typedef struct startparam_fmt
{
  boolean allocated;
  long numpop;
  long numtheta;
  float *theta;
  long nummig;
  float *mig;
  long numsplit;
  float *split;
  long numsplitstd;
  float *splitstd;
  long numrate;
  float *rate; 
} startparam_fmt;

/// data storage for all options
typedef struct _option
{
    char **buffer;  //pointer to buffer when MPI-worker otherwise NULL
    FILE *parmfile;
    FILE *seedfile;
    FILE *logfile;
    FILE *mixfile;
    /*general options */
    //
    // name length
    long nmlength;  /* length of individual names */
    long popnmlength;  /* length of population names */
    long allelenmlength;  /* length of allele names */
    //
    // custom migration matrix setup
    char *custm;   /* custom migration matrix: theta are first */
    char *custm2;   /*custom migration matrix: theta on diagonal */
    //long symn;
    //long sym2n;
    //long zeron;
    //long constn;
    //long tmn;
    //long mmn;
   // long *zeroparam;
   // long *constparam;
   // twin_fmt *symparam;
   // quad_fmt *sym2param;
   // long *mmparam;
    
    /*input/output options */
    int menu;
    boolean progress;
    boolean verbose;
    boolean writelog;
    boolean uep;
    MYREAL ueprate;
    MYREAL uepmu;
    MYREAL uepnu;
    MYREAL uepfreq0;
    MYREAL uepfreq1;
    boolean movingsteps;
    MYREAL acceptfreq;
    boolean printdata;
  boolean printcov;
    boolean usertree;
    boolean usertreewithmig;
    boolean randomtree;
    short treeprint;
  long treeinc;
    boolean usem;
  boolean recordedusem;
    short migvar;
    boolean plot;
    boolean plotnow;
    short plotmethod;
    short plotvar;
    short plotscale;
    MYREAL plotrange[4];
    long plotintervals;
    MYREAL *plotxvalues;
    MYREAL *plotyvalues;
    boolean simulation;
    boolean mighist;
  boolean mighist_all;
  long mighist_counter;
  long mighist_increment;
  boolean skyline;
  boolean skyline_param;
  long timeelements;
  float eventbinsize;
  boolean mixplot;
  boolean dist;
  boolean geo;
  boolean div;
  char *parmfilename;
  char *infilename;
  char *outfilename;
  char *pdfoutfilename;
  char *logfilename;
  char *mathfilename;
  char *treefilename;
  char *utreefilename;
  char *catfilename;
  char *weightfilename;
  char *sumfilename;
  char *mighistfilename;
  char *skylinefilename;
  char *distfilename;
  char *geofilename;
  char *divfilename;
  char *bootfilename;
  char *seedfilename;
  char *mixfilename;
#ifdef UEP
  char *uepfilename;
#endif
  char *datefilename;
  char *bayesfilename;
  char *bayesmdimfilename;
  boolean mdimdelete;
  int use_compressed;
  boolean prioralone;
  FILE *divtimefile;
  char *divtimefilename;
  char title[LINESIZE + 1];
  lratio_fmt *lratio;
  boolean aic;
  boolean fast_aic;
  MYREAL aicmod;
  FILE *aicfile;
  char * aicfilename;
  char aictype[3];
  char fsttype;
  boolean printfst;
  short profile;
  char profilemethod;
  long df;
  boolean qdprofile;
  boolean printprofsummary;
  boolean printprofile;
  short profileparamtype;
  
  /* data options */
  char datatype;
  int datamodel;
  boolean include_unknown;
  short migration_model;
  char dlm;
  long micro_stepnum;
  long micro_threshold;
  int msat_option;
  pair msat_tuning;
  MYREAL **steps;
  boolean interleaved;
  MYREAL seqerror[4];
  boolean seqerrorcombined;
  boolean has_estimateseqerror;
  MYREAL *ttratio;
  boolean fastlike;
  boolean freqsfrom;
  long categs;
  MYREAL *rate;
  long rcategs;
  MYREAL *rrate;
  MYREAL *probcat;
  long seqrate_gamma_num; /*number of loci for gamma values*/
  MYREAL *seqrate_gamma; /* for gamma values for many loci*/
  MYREAL probsum;
  boolean gammarates;
  boolean autocorr;
  MYREAL lambda;
  boolean weights;
  MYREAL freqa;
  MYREAL freqc;
  MYREAL freqg;
  MYREAL freqt;
  int sequence_option;
  int sequence_model;
  int sequence_model_numparam;
  MYREAL * sequence_model_parameters;
  /* random number options */
    short autoseed;
    unsigned long inseed;
    unsigned long saveseed;
    /* mcmc options */
    boolean bayes_infer;
  boolean integrated_like;
  /* slice sampling*/
  boolean slice_sampling[PRIOR_SIZE];
  MYREAL *slice_sticksizes;
  MYREAL updateratio; //needs revision
  MYREAL tree_updatefreq;
  MYREAL parameter_updatefreq;
  MYREAL haplotype_updatefreq;
  MYREAL timeparam_updatefreq;
  MYREAL unassigned_updatefreq;
  MYREAL seqerror_updatefreq;
  MYREAL mlalpha_updatefreq;
  //int bayesprior[4];
  int bayespretty;
  long bayes_priors_num;
  prior_fmt *bayes_priors; 
  boolean automatic_bins;
  long bayes_posterior_bins[PRIOR_SIZE];
  bayes_fmt *bayes;
  boolean has_bayesfile;
  boolean has_bayesmdimfile;
  long bayesmdiminterval;
  twin_fmt startguess[STARTGUESSNUM];//theta,migration,split, splitstd
  // PRIOR, RANDOM, PRIOR can take up the second arg for percentage
  startparam_fmt startparam; // holds all start values 
    boolean gamma;
  boolean murates;
  boolean murates_fromdata;
  boolean bayesmurates;
    MYREAL alphavalue;
    MYREAL *mu_rates;
    long muloci;
  long *newpops;
  long newpops_numalloc;
  long newpops_numpop;
  long *growpops;
  long growpops_numalloc;
  long growpops_numpop;
  //  boolean poprelabeled;
  boolean has_inheritance;
  MYREAL *inheritance_scalars;
  long inheritance_scalars_numalloc;
  boolean has_datefile;
#ifdef LONGSUM

    boolean fluctuate;
    MYREAL *flucrates;
#endif
  long burn_in;
  char burnin_autostop;
  short heating;
  MYREAL *thetag;
  long numthetag;
  MYREAL *mg;
  long nummg;
  MYREAL *splitg;
  long numsplitg;
  MYREAL *splitstdg;
  long numsplitstdg;
  long schains;
  long sincrement;
  long ssteps;
  long lchains;
  long lincrement;
  long lsteps;
  MYREAL lcepsilon;
  boolean gelman;
  boolean gelmanpairs;
  long pluschain;  //how many chains are allowed , currently HARDCODED
  boolean replicate;
  long replicatenum;
  long gridpoints;
  long numpop;
  MYREAL heat[50];
  boolean adaptiveheat;
  long heating_interval;
  long heated_chains;
  /* save genealogy summary options */
  boolean readsum;
  boolean checkpointing;
  boolean writesum;
  /* threading over loci */
  int cpu;
  //
  // fatal attraction of zero resistance
  MYREAL minmigsumstat;
  //
  MYREAL *mutationrate_year; // keeps mutation rate per year for each locus 
  long mutationrate_year_numalloc; // 
  MYREAL generation_year;
  long randomsubset;
  unsigned long randomsubsetseed;
#ifdef SEASON
  MYREAL *timings; // vector of time points (options)
  long numtimings; // number of elements in the timings vector
  MYREAL **timemodifiers; // vector of parameter modifier vectors
#endif
  boolean heatedswap_off;
  boolean haplotyping;
  boolean haplotyping_report;
  boolean has_autotune;
  MYREAL autotune;
  boolean has_unassigned;
  boolean has_unassignedfreq;
  boolean has_unassignedpoplist;
  size_t unassignedpoplistnum;
  long *unassignedpoplist;
  boolean oneliner;
  long totalsites; //for snp data
  MYREAL *wattersons;
  long *segregs;
  boolean onlyvariable;  // analyze only variable sequence loci
  MYREAL locusweight;    // use a weight to calculate the invariant locus
  boolean has_variableandone; // use only one invariant locus and reweight
  long firstinvariant; // locus number of first invariant locus to reweight using locusweight
  boolean has_speciation;
  boolean has_migration;
  boolean tersepdf; // prints only summaries and not per locus information
  boolean allposteriors; // plot all loci-posteriors
  long **unfinished;
  boolean hyperprior;
  long hyperinterval;
  MYREAL hyperfactormean;
  MYREAL hyperfactoralpha;
  boolean recorddivtime;
  long species_model_dist;
  int tri_mlalpha;
  double mlalpha;
  double mlinheritance;
  int smoothing[2];
  long smoothing_window[2];
}
option_fmt;


/// defines a node in the the timelist structure with reference to the tree node
typedef struct vtlist
{
    node *eventnode;  /* node with age=tyme */
    MYREAL age;   /* tyme from top nodelet */
    MYREAL interval;  /* interval t[i+1].age - t[i].age */
    long *lineages;
    long from;
    long to;
    /*  long pop; */
    long slice;
  long timeslice;
  boolean visited;
  char type;
}
vtlist;

/// holds the tree, access should be done only through root
typedef struct tree_fmt
{
    node **nodep;
    node *root;
    long pop;
    long tips;
}
tree_fmt;

/// holds timeslist that is used in the MCMC run
typedef struct timelist_fmt
{
    long numpop;
    long copies;
    long allocT;
    long T;
    long oldT;
  vtlist *tl;
  long *lineages;
}
timelist_fmt;

///
/// holds minimal statistic for the longsum version, this allows to calculate changes in size etc through time
/// holds all tree information: ordered: list of lineages, frompop, topop, eventtime, eventtype
#ifdef LONGSUM
typedef struct longsum_fmt
{
    long *lineages;
    long *lineages2;
    long fromto; //indicator into migration matrix mm2m(from,to)
    long to;
    MYREAL eventtime;
    MYREAL interval;
    char eventtype;
}
longsum_fmt;
#endif /*LONGSUM*/

///
/// time archive, holds minimal statistic for a single tree
/// is a compressed version of the longsum_fmt with many things precalculated
typedef struct tarchive_fmt
{
    long copies;
    MYREAL lcopies;
    MYREAL *data;   //hold content
#ifdef ALTIVEC

    FloatVec *vdata;
#endif

    MYREAL *point;  // points into data
    MYREAL *wait;   // points into data
    MYREAL *kt;   //points to wait
    MYREAL *km;   //points to wait
    MYREAL *p;   // points to point
    MYREAL *mindex;  // points to point
#ifdef LONGSUM

    longsum_fmt *longsum; // holds all tree information: ordered: list of lineages, frompop, topop, eventtime, eventtype
    long longsumlen; // holds list of groups in longsum.
#endif /*LONGSUM*/

}
tarchive_fmt;

///
/// time archive holds all statistics for all trees
typedef struct timearchive_fmt
{
    long allocT;
    long T;
    long numpop;
    long sumtips;
    MYREAL param_like;
    MYREAL thb;
    MYREAL alpha;
#ifdef ALTIVEC
    // FloatVec *lcopiesvec;
    // FloatVec *data;
#endif /*ALTIVEC*/

    tarchive_fmt *tl;         // holds each tree (summary)
    MYREAL *parameters;  // holds all the parameter data
    MYREAL *param;  // pointer into parameters
    MYREAL *param0;  //pointer into parameters
    MYREAL *lparam0;  //pointer into parameters
    MYREAL *likelihood;  //pointer into parameters
    long trials;
    MYREAL normd;
}
timearchive_fmt;


/// holds pltting parameters
typedef struct _plotmax
{
    MYREAL x1;
    MYREAL y1;
    MYREAL l1;
    MYREAL x2;
    MYREAL y2;
    MYREAL l2;
}
plotmax_fmt;

/// holds quantile information
typedef struct _quantile_fmt
{
    char *name;
    MYREAL *param;
}
quantile_fmt;

/// format for migration events
/// holds time, from, to, and sumtips
///typedef MYREAL migevent_fmt[4];
typedef struct _migevent_fmt
{
  float age;
  int from;
  int to;
  long sumlines;
} migevent_fmt;

/// holds migration events for migration event time histogram
typedef struct _mighist_fmt
{
  long allocsize;
  long copies;
  long weight;
  long migeventsize;
  migevent_fmt *migevents;
}
mighist_fmt;

/// holds all migration histograms for all loci
typedef struct _mighistloci_fmt
{
  // records migration events into mighist container
  mighist_fmt *mighist;
  long mighistnum;
  long allocsize;

  // records all events  parameters per bin derived from timelist
  // should make communication over MPI much smaller
  duo **migeventbins;
  long *migeventbinnum;
  // records all expected  parameters per bin derived from timelist
  tetra **eventbins; 
  long *eventbinnum;
  MYREAL eventbinsize; // used for eventbins and migeventbins

}
mighistloci_fmt;

/// \todo what does this format do?
typedef struct _histogram_fmt
{
    long count;
    MYREAL *time;
    long *weight;
}
histogram_fmt;

/// reduced set of options that are used to run the MCMC chain
typedef struct _worldoption
{
  boolean kernelsmooth;
  boolean allposteriors; // plot all loci-posteriors
    boolean gamma;
    MYREAL alphavalue;
    boolean murates;
  boolean murates_fromdata;
    long muloci;
    MYREAL *mu_rates;
    MYREAL *lmu_rates;
  MYREAL *meanmu;
  boolean has_inheritance;
  MYREAL *inheritance_scalars;
  boolean prioralone;
  MYREAL *heat;
  MYREAL *averageheat;
#ifdef LONGSUM

    boolean fluctuate;
#endif

    short migration_model;
    char *custm;
    char *custm2;
    MYREAL *thetag;
    MYREAL *mg;
    long zeron;
    long *zeroparam;
    long constn;
    long *constparam;
    long symn;
    long sym2n;
    twin_fmt *symparam;
    quad_fmt *sym2param;
    long tmn;
    long mmn;
    long *mmparam;
    boolean mixplot;
    boolean progress;
    boolean writelog;
    FILE *logfile;
    boolean plotnow;
    boolean verbose;
  boolean tersepdf; // prints only summaries and not per locus information
    boolean replicate;
    boolean gelman;
  boolean gelmanpairs;
    MYREAL lcepsilon;
    boolean simulation;
    char datatype;
    long lsteps;
    long lincr;
    MYREAL loglsteps;
    long treeprint;
  long treeinc;
    long movingsteps;
    MYREAL acceptfreq;
    long rcategs;
    long categs;
    short heating;
    long heated_chains;
    long heating_interval;
  long heating_count;
    boolean adaptiveheat;
    char profilemethod;
    boolean printprofile;
    boolean printprofsummary;
    long profileparamtype;
    long df;
    long lchains;
    long replicatenum;
    long *micro_threshold;
    long micro_stepnum;
  pair msat_tuning;
    MYREAL ***steps;
    MYREAL *rrate;
    MYREAL *rate;
    MYREAL *probcat;
    long pluschain;
    boolean mighist;
  boolean mighist_all;
  long mighist_counter;
  long mighist_increment;
  boolean skyline;
  boolean skyline_param;
  long timeelements;
  float eventbinsize;
    long burn_in;
  char burnin_autostop;
    boolean usem;
    short migvar;
    boolean plot;
    long plotmethod;
    long plotintervals;
    MYREAL plotrange[4];
    short plotscale;
    long plotvar;
    MYREAL *plotxvalues;
    MYREAL *plotyvalues;
    lratio_fmt *lratio;
    boolean aic;
    boolean fast_aic;
    MYREAL aicmod;
    FILE *aicfile;
    char aictype[3];
    MYREAL lambda;
    FILE *mixfile;
#ifdef UEP

    boolean uep;
    MYREAL ueprate;
    MYREAL uepmu;
    MYREAL uepnu;
    MYREAL uepfreq0;
    MYREAL uepfreq1;
#endif

    boolean fastlike;
    boolean bayes_infer;
  boolean slice_sampling[PRIOR_SIZE];
  MYREAL *slice_sticksizes;
  //@@@@@@@  MYREAL updateratio;
  double choices[6];
  boolean has_bayesfile;
  boolean has_bayesmdimfile;
  long bayesmdiminterval;
    MYREAL minmigsumstat;
  boolean has_datefile;
  MYREAL *mutationrate_year; // keeps mutation rate per year for each locus 
  long mutationrate_year_numalloc; // 
  MYREAL generation_year;
  boolean treeinmemory;
  boolean heatedswap_off;
  boolean has_autotune;
  MYREAL autotune;
  boolean onlyvariable;  // analyze only variable sequence loci
  MYREAL locusweight;    // use a weight to calculate the invariant locus
  boolean has_variableandone; // use only one invariant locus and reweight
  long firstinvariant; // locus number of first invariant locus to reweight using locusweight
  long *growpops; // copy of options->growthpop to specify whether populations have different growth rate [growth is directly under world]
  long growpops_numalloc;
  long growpops_numpop;
  int smoothing[2];
  long smoothing_window[2];
}
worldoption_fmt;

/// defines the time range of a unique event polymorphism
typedef struct _ueptime
{
    long size;
    long *populations;  // population the nodes are in
    MYREAL *ueptime;  // first elementis bottom most time , last elem is
    // topmost time,in between are for each migration events
}
ueptime_fmt;

/// reduced set of parameters related to the data used to run MCMC
typedef struct _worlddata
{
  boolean *skiploci;
  MYREAL *geo;
  MYREAL *lgeo;
  long *maxalleles;
  seqmodel_fmt **seq;
  FILE *sumfile;
  MYREAL freq;
  MYREAL freqlast;
  tipdate_fmt ***sampledates;
  MYREAL maxsampledate;
  long **numind;
#ifdef UEP

    long uepsites;
#endif
  long numindividuals;
  individualDB_fmt *individuals;
  boolean haplotyping;
  boolean haplotyping_report;
  long allsubloci;
  //  long *haplotyping_numind;
  //char ***haplotyping_buffer;
  //long **haplotyping_buflen;
  MYREAL *locusweight; //invariant loci treatment
}
worlddata_fmt;

//typedef struct _proposal_fmt proposal_fmt;

typedef struct _convergence
{
  MYREAL *gelmanmeanmaxR;
  MYREAL gelmanmeanRall;
  MYREAL gelmanmaxRall;
  MYREAL *chain_s;
  MYREAL *chain_means;
  long *chain_counts;
} convergence_fmt;

#ifdef MPI
typedef struct _mpirequest_fmt 
{
  char *tempstr;
  int sender;
  int tag;
} mpirequest_fmt;
#endif
/// holds all the parameters etc to run the program

//typedef struct proposal_fmt;
typedef struct _divtime_fmt 
{
  double age;
  long from;
  long to;
} divtime_fmt; 



typedef struct _world
{
  /* generalities */
  char *name;
  char *worldname;
  worldoption_fmt *options;
  worlddata_fmt *data;
  
#ifdef MPI
  int *who;
  int *mpistack;
  long mpistack_numalloc;
  mpirequest_fmt *mpistack_request;
  long mpistack_request_numalloc;
  long mpistack_requestnum;
  long mpistacknum;
#ifdef SLOWNET
  
  int *profilewho;
#endif
#endif /*MPI*/
  long allocbufsize;
  char *buffer;   // buffer for profiles, not needed before profiletables()
  long loci;
  long skipped;   /*loci with no data */
  long locus;   /* the current locus, if single then set to 0 */
  long thislocus;  /* the real current locus, the whole locus scheme
			needs revision */
  long numpop;
  long numpop2;
  long sumtips;
  /* migration parameter array  starts and ends */
  int *mstart;
  int *mend;
  
  /* time archives, contains the data/results for summarizing */
  timearchive_fmt **atl;
  
  /* migration histogram reporter */
  mighistloci_fmt *mighistloci;
  long mighistlocinum;
  
  /*tree material */
  node **nodep;
  node *root;
  MYREAL treelen;
  long unique_id;
  tbl_fmt tbl;
  contribarr *contribution;
  /* parameter */
  MYREAL *param0;
  MYREAL *param00;
  MYREAL **fstparam;
  /*allows for 3 different rate changes at specific times*/
#ifdef LONGSUM
  MYREAL *flucrates;
  MYREAL *lflucrates;
#endif
  /* mcmc related */
  long *lineages;
  timelist_fmt *treetimes;
  MYREAL *mig0list;
  MYREAL **migproblist;
  MYREAL **speciesproblist;
  long *design0list;
  
  /*nr related */
  MYREAL ***apg0;  /* part-loglikelihoods of param0 */
  MYREAL ***apg;  /* part-loglikelihoods of param */
  
  /*heating scheme */
  boolean cold; /* true for cold chain, false for all others*/
  boolean has_proposal_details;
  boolean has_proposal_first;
  boolean in_burnin;
  long actualinc;
  long increment;
  //  MYREAL heatratio;
  MYREAL heat;
  MYREAL averageheat;
  MYREAL varheat;
  long heatid;
  struct _proposal_fmt *proposal;
  MYREAL essminimum;
  MYREAL logprior;
    long treeswapcount;
    /* reporting time */
    time_t starttime;
    long treesdone;
    long treestotal;

  burnin_record_fmt *burnin_stops;
  long burnin_stops_alloc;
  long burnin_counter;
  long burnin_z;
    /* gelman R reporting/checking "gelman-convergence option" */
  convergence_fmt *convergence;

    long chains;
    boolean start;  //I need this in estimateParameter()
    /* reporting */
    MYREAL *likelihood;  /* data likelihoods */
    long alloclike; /* how many likelihoods can we store */
    long numlike;   /* how many likelihoods maximally are stored */
    plotmax_fmt **plotmax;
    quantile_fmt *quantiles;
    boolean **percentile_failed; /*indicator whether the percentile calculation failed*/
    boolean percentile_some_failed; /*some precentiles failed*/
    MYREAL maxdatallike;  /* the maximum log likelihood of a chain */
    MYREAL allikemax;  /* the maximum log likelihood  the best tree */
    boolean in_last_chain; /*for last CHAIN option in print-tree */
    MYREAL param_like;  /* current parameter likelihood [=current chain] */
    MYREAL **chainlikes;  /* parameter likelihood of replicate/last chains all loci */
    long trials;
    MYREAL normd;
    long repkind;
    long rep;   // defines single cahin estimators
    long replicate;  //holdss the replicate we are in
    long repstop;
    long lsteps;
    MYREAL ***cov;
    long migration_counts;
    char ****plane;
#ifdef UEP

    MYREAL ueplikelihood;
    MYREAL **ueplike;  // contains uep like per tree pop x ueps
    MYREAL ****ueplikestore; // constains all ueplike long-steps x pop x uep
    ueptime_fmt *ueptime;  //contains ueptime per tree and mutation
    ueptime_fmt ***ueptimestore; //contains all ueptimes
    long ***ueprootstore;  //stores uep status at root
    long *oldrootuep;
    long **uepanc;
#endif
    long accept;
    MYREAL accept_freq;
    long swapped;
    long G;
  long bayesaccept;
  bayes_fmt *bayes;
  FILE *bayesfile;
#ifdef ZNZ
  znzFile bayesmdimfile;
#else
  FILE *bayesmdimfile;
#endif
  FILE *outfile;
  FILE *pdfoutfile;
  FILE *treefile;
  FILE *mathfile;
  FILE *mighistfile;
  FILE *skylinefile;
  char **treespace;
  long *treespacenum;
  long *treespacealloc;
  MYREAL *besttreelike;
  // Bayes factor
  MYREAL *bf;
  MYREAL *bfscale;
  MYREAL *hmscale;
  MYREAL *amscale;
  MYREAL *hm;
  MYREAL *am;
  // autocorrelation and ESS
  MYREAL *autocorrelation;
  MYREAL *effective_sample;
  MYREAL *auto_archive;
  MYREAL *ess_archive;
  long archive_n;
  node **node_collection;
  long node_collection_count;
  long node_collection_allocated;
#ifdef MPI
#ifdef PARALIO
  MPI_File mpi_bayesmdimfile;
#endif
#endif
  long *numsubloci;//per locus
  long *sublocistarts;//per locus
  long *maxnumpattern;//per locus
  mutationmodel_fmt *mutationmodels;
  //
  float page_width;
  float page_height;
  //
  // haplotyping
  individualDB_fmt **haplotypes;
  long *numhaplotypes;
  //
#ifdef BEAGLE
  beagle_fmt *beagle;
#endif
#ifdef SEASON
  MYREAL *timings; // vector of time points
  long numtimings; // number of elements in the timings vector
  MYREAL **timemodifiers; // vector of parameter modifier vectors
#endif
  char * warning;
  long warningsize;
  long warningallocsize;
#ifdef MPI
  char ****indnames;//otherwise I need the data structure for packing, only used in workers
#endif
  long timeelements;
  MYREAL *times;
  MYREAL *timek;
  long numtimesalloc;
  long numtimes;
  float **recording_times;
  // assignment
  boolean has_unassigned;
  boolean has_unassignedfreq;
  long unassignednum;
  unassigned_fmt **unassigned;
  boolean has_unassignedpoplist;
  size_t unassignedpoplistnum;
  long *unassignedpoplist;
  long assigncount;
  long *seqerrorcount;
  MYREAL **seqerrorrates;
  boolean has_estimateseqerror;
  long *seqerrorratesnum;
  long *seqerrorallocnum;
  long *seqerrorsteps;
  boolean seqerrorcombined;
  long maxreplicate;
  long *accept_archive;
  long *trials_archive;
  boolean has_speciation;
  boolean has_migration;
  species_fmt *species_model;
  long species_model_size;
  long species_model_dist;
  divtime_fmt *divtime;
  long divtime_alloc;
  long divtime_num;
  FILE *divtimefile;
  int tri_mlalpha;
  double mlalpha; //mittag-leffler
  double mlinheritance;
  double *steppingstones;
  double *steppingstone_scalars;
  double *steppingstone_counters;
  boolean has_growth;
  double *growth; // contains growth values: growpops={1,1,1,1} => growth={x},growpop={1,2,1} => growth={x1,x2}
  double *savegrowth;
  long grownum;
}
world_fmt;



/// helper for the ML maximizer, nr originally came from newton-raphson
typedef struct _nr_fmt
{
    long partsize;  /*number of part-variables, fixed per model */
    MYREAL *parts;  /* parts of the first and second derivatives */
    MYREAL *d;   /* first derivates */
    MYREAL *od;
    MYREAL *dv;
    MYREAL *delta;
    MYREAL *gdelta;
    MYREAL *param;  /* changed values of param */
    MYREAL *lparam;  /* changed values of log(param) */
    MYREAL *values;  /* profile values */
    MYREAL *locilikes;

  //  MYREAL *saved;  /* saved first derivates (used in MPI code) */
    /* migration parameter array  starts and ends */
    int *mstart;   //link to world->mstart
    int *mend;   //link to world->mend
    MYREAL normd;
    long repkind;
    long repstart;
    long repstop;
    //  long numg;
    timearchive_fmt **atl;
    world_fmt *world;
    MYREAL **dd;   /* second derivatives */
    MYREAL llike;   /* parameter LOGlikelihood */
    MYREAL lastllike;  /* parameter LOGlikelihood */
    MYREAL ollike;  /* old "   " */
    MYREAL *datalike;  /*P(D|G) */
    MYREAL ***apg0;  /* part-loglikelihoods of param0 */
    MYREAL ***apg;  /* part-loglikelihoods of param */
    MYREAL *apg_max;  /* maxvalue of apg */
  long seqrate_gamma_num; /*number of loci for gamma values*/
  MYREAL *seqrate_gamma; /* for gamma values for many loci*/
    MYREAL *rate;   /* rates for gamma rates */
    MYREAL *probcat;  /* probabilities for gamma rates */
    long categs;   /* #categories for gamma rates */
    MYREAL alpha;
    long numpop;   /* number of populations */
    long numpop2;   /* 2*numpop */
    MYREAL *PGC;   /* uncorrected llike */
    MYREAL *oPGC;   /* uncorrected ollike */
    long copy_nr;   /* number of genealogies */
    boolean *skiploci;
    long profilenum;  /* number of profile parameters */
    long *profiles;  /* which profile parameters */
    long *indeks;   /* which noprofile parameters */
}
nr_fmt;

/// helper format for the maximizer and the profile calculator
typedef struct helper_fmt
{
    long locus;
    nr_fmt *nr;
    timearchive_fmt **atl;
    //MYREAL *param;
    long which;
    MYREAL weight;
    MYREAL ll;
    boolean multilocus;
    boolean boolgamma;
    long analystype;
    MYREAL *dv;
    MYREAL *xv;
    MYREAL *expxv;
    MYREAL sign;
    MYREAL lamda;
}
helper_fmt;


/// migration table 
typedef struct _migr_table_fmt
{
  long from;
  long to;
  MYREAL time;
  char event;
}
migr_table_fmt;

typedef struct _div_fmt 
{
  longpair *divlist;
  long div_allocsize;
  long div_elem;
} div_fmt;


///
/// during the treeupdates we need a scratchpad, proposal_fmt
/// holds all temporary data
typedef  struct _proposal_fmt
{
  world_fmt *world;
  char datatype;
  long sumtips;
  long numpop;
  long endsite;
  MYREAL fracchange;
  MYREAL *param0;
  MYREAL *param0save;
  node *root;
  short migration_model;
  boolean mig_removed;
  MYREAL rr;
  node **nodedata;      // holds abovenodes and bordernodes
  node *origin;            // pointers into nodedata
  node *target;            // .... 
  node *realtarget;
  node *tsister;
  node *realtsister;
  node *osister;
  node *realosister;
  node *ocousin;
  node *realocousin;
  node *oback;
  node *realoback;
  long line_f_allocsize;
  long line_t_allocsize;
  node **line_f;
  node **line_t;
  node *connect;
  MYREAL likelihood;
  MYREAL ueplikelihood;
  MYREAL **ueplike;
  MYREAL time;
  MYREAL v;
  MYREAL vs;
  xarray_fmt *xt;
  xarray_fmt *xf;
  MYREAL **mf;
  MYREAL **mt;
  //could be used to calculate p(gn|go):   MYREAL nu;
#ifdef UEP
  ueparray_fmt ut;
  ueparray_fmt uf;
  node *firstuep;
  MYREAL *umf;
  MYREAL *umt;
#endif
  long aboveorigin_allocsize;
  node **aboveorigin;
  long divlist_allocsize;
  div_fmt *divlist;
  long bordernodes_allocsize;
  node **bordernodes;
  migr_table_fmt *migr_table;
  migr_table_fmt *migr_table2;
  long migr_table_counter;
  long migr_table_counter2;
  long old_migr_table_counter;
  long old_migr_table_counter2;
  long timeslice;
  MYREAL *mig0list;
  long *design0list;
  MYREAL treelen;
  
#ifdef BEAGLE
  long parentid;
  long leftid;
  long rightid;
#endif
} proposal_fmt;


/// splines are not yet used in the profiles
typedef struct spline_fmt
{
    long ntab;
    long nwork;
    MYREAL *param;
    MYREAL *like;
    MYREAL *diff;
    MYREAL *diff2;
    long *constr;
    long *diagn;
    MYREAL *work;
}
spline_fmt;

/// ?????????????????
typedef struct locusdata_fmt
{
    long locus;
    world_fmt **universe;
    option_fmt *options;
    data_fmt *data;
}
locusdata_fmt;


/// holds the minimal statistic for AIC for each model
typedef struct _aic_fmt
{
    MYREAL aic;
    char *pattern;
    long numparam;
    MYREAL mle;
    MYREAL lrt;
    MYREAL prob;
    MYREAL probcorr;
}
aic_fmt;

/// holds all models that are used for the AIC
typedef struct _aic
{
    aic_fmt * aicvec;
    long aicnum;
    MYREAL *param0;
}
aic_struct;

///
/// plot options
typedef struct _plotfield_fmt
{
    boolean print;
    char type;   // 'a' = ascii plot (see mig-histogram.c) , 'p' = PDF plotting (see pretty.c)
    long xsize;   // width  printpositions in chars
    long ysize;   // height printpositions in chars
    char xaxis[255];  //xaxis label
    char yaxis[255];  //yaxis label
    char yfaxis[255];  //frequency yaxis label
    char title[255];
    float *yfreq;
    long *y;
    char **data;   //the plotplane
}
plotfield_fmt;

/// 3 longs used for random number
typedef long longer[3];  /* used in random.c */


//function pointers tupe definitions
typedef MYREAL (*logpriorratioptr)(MYREAL, MYREAL, bayes_fmt *, long);
typedef MYREAL (*logpriorptr)(world_fmt *, long);
typedef MYREAL (*logprior1ptr)(world_fmt *, long, MYREAL);
typedef MYREAL (*profuncptr)(MYREAL, long, world_fmt *, MYREAL *);
typedef MYREAL (*hastratioptr)(MYREAL , MYREAL, MYREAL, MYREAL, bayes_fmt *, long);

typedef void (*nuview_function) (mutationmodel_fmt *, long, long , node *, world_fmt *, long);
typedef double (*prob_micro_function) (MYREAL , long , world_fmt *, mutationmodel_fmt *, pair *);
typedef void (*pseudonuview_function) (mutationmodel_fmt *, proposal_fmt *, xarray_fmt *, MYREAL *, MYREAL , xarray_fmt *, MYREAL *, MYREAL, long);

typedef double (*time_to_speciate_func) (world_fmt *, long, double, char *, long *, long *);
typedef double (*log_prob_wait_speciate_func)(double, double, double, double, species_fmt *);
typedef double (*log_point_prob_speciate_func)(double, double, double, species_fmt *);

extern long *seed;
extern long *newseed;
#ifdef SLOWNET
extern MYREAL (*calc_like) (helper_fmt *, MYREAL *, MYREAL *);
extern void (*calc_gradient) (nr_fmt *, helper_fmt *, MYREAL *);
extern void (*setup_param0) (world_fmt *, nr_fmt *, long, long, long, long,
                                 long, boolean);

#endif
//#ifdef MEMDEBUG
//#include <memcheck.h>
//#endif
#include "sort.h"
//#include "migrate_mpi.h"

#ifdef WINDOWS
typedef struct mycomplex {
  double r;
  double i;
} mycomplex;
#endif

#endif
