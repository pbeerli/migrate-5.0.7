/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 S E Q U E N C E S   R O U T I N E S 
 
 
 Peter Beerli 1996, Seattle
 beerli@fsu.edu
 
 Copyright 2002 Peter Beerli and Joseph Felsenstein
 
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
 
 
 $Id: sequence.c 2170 2013-09-19 12:08:27Z beerli $
 
-------------------------------------------------------*/
/* \file sequence.c

*/
#include "migration.h"
#include "sighandler.h"
#include "data.h"
#include "tools.h"
#include "migrate_mpi.h"
#include "watterson.h"
#include "tree.h"
#include "mutationmodel.h"
#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif
/* prototypes ------------------------------------------- */
//void make_sequences_old (world_fmt * world, option_fmt * options, data_fmt * data,
//                     long locus);
void init_sequences (world_fmt * world, option_fmt * options, data_fmt * data,
                     long locus);
void init_sequences2 (world_fmt * world, mutationmodel_fmt *s);
void init_sequences2_old (world_fmt * world, seqmodel_fmt * seq, long locus);

void initratio (option_fmt * options);
void initfreqs (MYREAL *freqa, MYREAL *freqc, MYREAL *freqg, MYREAL *freqt);
void initcatn (long *categs);
boolean initcategs (long categs, MYREAL *rate, MYREAL *probcat);
void initprobcat (long categs, MYREAL *probsum, MYREAL *probcat);
void init_tbl (world_fmt * world, long locus);
void print_weights (FILE * outfile, world_fmt * world, option_fmt * options,
                    long locus);
void print_tbl (FILE * outfile, world_fmt * world, option_fmt * options,
                long locus);
MYREAL treelike_seq (mutationmodel_fmt *s, long sublocus, world_fmt * world, long locus);
MYREAL treelike_snp (mutationmodel_fmt *s, long sublocus, world_fmt * world, long locus);
void find_rates_fromdata(data_fmt * data, option_fmt * options, world_fmt * world);

/*private functions */
void getbasefreqs (option_fmt * options, seqmodel_fmt * seq, long locus);
void empiricalfreqs (world_fmt * world, option_fmt * options,
                     mutationmodel_fmt * s, long sublocus);
void makeweights (world_fmt * world, data_fmt * data, option_fmt *options, long locus);
void makevalues_seq (world_fmt * world, option_fmt * options, data_fmt * data,
                     long locus);
void make_invarsites (world_fmt * world, data_fmt * data, long locus);
void make_invarsites_unlinked (world_fmt * world, data_fmt * data,
                               long locus);
void sitecombine2_old (world_fmt * world, data_fmt * data, option_fmt *options, long sites,
                   long locus);
void sitesort2_old (world_fmt * world, data_fmt * data, option_fmt *options, long sites, long locus);
void sitescrunch2_old (world_fmt * world, long sites, long i, long j, long locus);

void sitecombine2 (world_fmt * world, data_fmt * data, option_fmt *options, long sites,
                   long locus, long sublocus);
void sitesort2 (world_fmt * world, data_fmt * data, option_fmt *options, long sites, long locus, long sublocus);
void sitescrunch2 (world_fmt * world, long sites, long i, long j, long sublocus);
void inputoptions (world_fmt * world, option_fmt * options, data_fmt * data,
                   long locus);
void inputweights (mutationmodel_fmt * s, data_fmt * data, long chars);
void inputcategs (long a, long b, mutationmodel_fmt * s, option_fmt * options,
                  data_fmt * data);
void initlambda (option_fmt * options);
void printweights (FILE * outfile, world_fmt * world, option_fmt * options,
                   short inc, long chars, short *weight, char *letters);
void print_seqfreqs (FILE * outfile, world_fmt * world, option_fmt * options);

void snp_invariants (mutationmodel_fmt *s, long sublocus, contribarr invariants, world_fmt *world, long locus, phenotype x1, MYREAL *scale);

void makevalues_snp(world_fmt * world, option_fmt * options, data_fmt * data,
                     long locus);

void init_sequences_aliases (world_fmt * world, option_fmt * options,
                             data_fmt * data, long locus);
void constrain_rates(long categs, MYREAL *rate, MYREAL *probcat);
void makeweights_old (world_fmt * world, data_fmt * data, option_fmt *options, long locus);
void set_nucleotide(MYREAL *treedata, const char nucleotide, const MYREAL *seqerr);
MYREAL treelike_snp_unlinked (mutationmodel_fmt *s, long xs, world_fmt * world, long locus);
void copy_seq (world_fmt * original, world_fmt * kopie);
void find_rates_fromdata_alleles(data_fmt * data, option_fmt * options, world_fmt *world, MYREAL mean);
void free_seq(seqmodel_fmt **seq, long seqnum);

//##

void check_basefreq (option_fmt * options);

extern void swap (contribarr *a, contribarr *b);




void
init_sequences_aliases (world_fmt * world, option_fmt * options,
                        data_fmt * data, long locus)
{
    inputoptions (world, options, data, locus);
    makeweights (world, data, options, locus);
}


/* menu material ----------------------------------------- */
void
initratio (option_fmt * options)
{
    long z = 0;
    char *tmp;
    char input[LINESIZE];
    printf
    ("Transition/transversion ratio?\nEnter a value for each locus, spaced by blanks or commas\n[Minium value > 0.5]");
    FGETS (input, LINESIZE, stdin);
    tmp = strtok (input, " ,\n");
    while (tmp != NULL)
    {
        options->ttratio[z++] = atof (tmp);
        tmp = strtok (NULL, " ,;\n");
        options->ttratio =
	  (MYREAL *) myrealloc (options->ttratio, sizeof (MYREAL) * (size_t) (z + 1));
        options->ttratio[z] = 0.0;
    }
}

void
initfreqs (MYREAL *freqa, MYREAL *freqc, MYREAL *freqg, MYREAL *freqt)
{
    char input[LINESIZE];
    int scanned;
    MYREAL summ = 0;

    printf
    ("Base frequencies for A, C, G, T/U\n (use blanks to separate, if all are equal use a sign [=])?\n");
    for (;;)
    {
        FGETS (input, LINESIZE, stdin);
        if (input[0] == '=')
        {
            scanned = 4;
            *freqa = *freqc = *freqg = *freqt = 0.25;
        }
        else
#ifdef USE_MYREAL_FLOAT
            scanned = sscanf (input, "%f%f%f%f%*[^\n]", freqa, freqc, freqg, freqt);
#else
        scanned = sscanf (input, "%lf%lf%lf%lf%*[^\n]", freqa, freqc, freqg, freqt);
#endif
        if (scanned == 4)
            break;
        else
            printf ("Please enter exactly 4 values.\n");
    };
    // adjust frequencies to a total of 1
    summ = *freqa + *freqc + *freqg + *freqt;
    if (summ != 1.0)
        printf ("Frequency values were adjusted to add up to 1.0.\n");
    *freqa /= summ;
    *freqc /= summ;
    *freqg /= summ;
    *freqt /= summ;
    printf ("Nucleotide frequencies: A=%f, C=%f, G=%f, T=%f\n", *freqa, *freqc,
            *freqg, *freqt);
}


void
initcatn (long *categs)
{    /* initialize category number */
  //int retval;
    do
    {
        printf ("Number of categories (1-%d)?\n", MAXCATEGS);
        scanf ("%ld%*[^\n]", categs);
        getchar ();
    }
    while (*categs > MAXCATEGS || *categs < 1);
}


boolean
initcategs (long categs, MYREAL *rate, MYREAL *probcat)
{    /* initialize rate categories */
    long i;
    char input[LINESIZE];
    char rest[LINESIZE];
    int scanned;
    boolean done;

    for (;;)
    {
        printf
	  ("Either enter the Shape parameter alpha for Gamma deviated rates\n*OR* enter the rates for each category (use a space to separate)\n===>");fflush(stdout);
        FGETS (input, LINESIZE, stdin);
        done = TRUE;
        if (count_words (input) == 1)
        {
            gamma_rates (rate, probcat, categs, input);
            return TRUE;
        }
        for (i = 0; i < categs; i++)
        {
#ifdef USE_MYREAL_FLOAT
            scanned = sscanf (input, "%f %[^\n]", &rate[i], rest);
#else
            scanned = sscanf (input, "%lf %[^\n]", &rate[i], rest);
#endif
            if ((scanned < 2 && i < (categs - 1))
                    || (scanned < 1 && i == (categs - 1)))
            {
                printf ("Please enter exactly %ld values.\n", categs);
                done = FALSE;
                break;
            }
            strcpy (input, rest);
        }
        if (done)
            break;
    }
    return FALSE;
}

void
initprobcat (long categs, MYREAL *probsum, MYREAL *probcat)
{
    long i;
    boolean done;
    char input[LINESIZE];
    char rest[LINESIZE];
    int scanned;

    do
    {
        printf ("Probability for each category?");
        printf (" (use a space to separate)\n");
        FGETS (input, LINESIZE, stdin);
        done = TRUE;
        for (i = 0; i < categs; i++)
        {
#ifdef USE_MYREAL_FLOAT
            scanned = sscanf (input, "%f %[^\n]", &probcat[i], rest);
#else
            scanned = sscanf (input, "%lf %[^\n]", &probcat[i], rest);
#endif
            if ((scanned < 2 && i < (categs - 1))
                    || (scanned < 1 && i == (categs - 1)))
            {
                done = FALSE;
                printf ("Please enter exactly %ld values.\n", categs);
                break;
            }
            strcpy (input, rest);
        }
        if (!done)
            continue;
        *probsum = 0.0;
        for (i = 0; i < categs; i++)
            *probsum += probcat[i];

        if (fabs (1.0 - (*probsum)) > 0.001)
        {
            for (i = 0; i < categs; i++)
                probcat[i] /= *probsum;
            printf ("Probabilities were adjusted to add up to one\n");
            for (i = 0; i < categs; i++)
                printf ("  %li> %f\n", i + 1, probcat[i]);
            printf ("\n\n");
        }
    }
    while (!done);
}

///
/// constrains the arbitrary site rate variation to an average of 1.0
void constrain_rates(long categs, MYREAL *rate, MYREAL *probcat)
{
  char input[LINESIZE];
  long i;
  MYREAL mean;
  printf("By default the rates will be constrained to an have an average rate of 1.0\nPress return if that is OK (preferred option) or enter the word RAW\n===>"); fflush(stdout);
  FGETS (input, LINESIZE, stdin);
  if(input[0] == '\0')
    {
      mean = 0.0;
      for (i = 0; i < categs; i++)
	{
	  mean += probcat[i] * rate[i];
	}
      for (i = 0; i < categs; i++)
	{
	  rate[i] /= mean;
	}
    }
}

/*data read material ===================================== */
///
/// read sequence data and linked SNP data 
/*
void
make_sequences_old (world_fmt * world, option_fmt * options, data_fmt * data,
                long locus)
{
    if (world->sumtips==0)
      {
	data->skiploci[locus] = TRUE;
        world->data->skiploci[locus] = TRUE;
        world->skipped += 1;
      }

  makevalues_seq (world, options, data, locus);
  //if (options->freqsfrom)
  //  {
  //    empiricalfreqs (world, options, world->data->seq[0], locus);
  //    getbasefreqs (options, world->data->seq[0], locus);
  // }   
}

*/

/* private functions================================== */
void
getbasefreqs (option_fmt * options, seqmodel_fmt * seq, long locus)
{
  long l;

  register MYREAL freqa, freqc, freqg, freqt, freqr, freqy;
  register MYREAL /*freqar, freqcy,*/ freqgr, freqty;

  MYREAL aa, bb;
  if (locus == 0)
    seq->ttratio = options->ttratio[0];
  else
    {
      for (l = 1; l <= locus; l++)
        {
	  if (options->ttratio[l] == 0.0)
            {
	      seq->ttratio = options->ttratio[l - 1];
	      break;
            }
	  seq->ttratio = options->ttratio[l];
        }
      if (l > locus)
	seq->ttratio = options->ttratio[locus];
    }
  check_basefreq (options);
  
  seq->basefrequencies[NUC_A] = options->freqa;
  seq->basefrequencies[NUC_C] = options->freqc;
  seq->basefrequencies[NUC_G] = options->freqg;
  seq->basefrequencies[NUC_T] = options->freqt;    
  freqa = seq->basefrequencies[NUC_A];
  freqc = seq->basefrequencies[NUC_C];
  freqg = seq->basefrequencies[NUC_G];
  freqt = seq->basefrequencies[NUC_T];
  seq->basefrequencies[NUC_R] = freqa + freqg;
  seq->basefrequencies[NUC_Y] = freqc + freqt;
  freqr = seq->basefrequencies[NUC_R];
  freqy = seq->basefrequencies[NUC_Y];
  seq->basefrequencies[NUC_AR]= freqa / freqr;
  seq->basefrequencies[NUC_CY]= freqc / freqy;
  seq->basefrequencies[NUC_GR]= freqg / freqr;
  seq->basefrequencies[NUC_TY]= freqt / freqy;
  //freqar = seq->basefrequencies[NUC_AR];
  //freqcy = seq->basefrequencies[NUC_CY];
  freqgr = seq->basefrequencies[NUC_GR];
  freqty = seq->basefrequencies[NUC_TY];

  aa =
    seq->ttratio * (freqr) * (freqy) - freqa * freqg -
    freqc * freqt;
  bb = freqa * (freqgr) + freqc * (freqty);
  seq->xi = aa / (aa + bb);
  seq->xv = 1.0 - seq->xi;
  if (seq->xi <= 0.0)
    {
      warning ("This transition/transversion ratio (%f)\n",seq->ttratio);
      warning ("is impossible with these base frequencies (%f, %f, %f, %f)!\n",freqa,freqc,freqg,freqt);
      seq->xi = 0.00001; // do not set this to zero because of the 1/(fracchange=xi*(...))
      seq->xv = 0.99999;
      seq->ttratio =
	(freqa * freqg +
	 freqc * freqt) / ((freqr) * (freqy));
      
      warning (" Transition/transversion parameter reset\n");
      warning ("  so transition/transversion ratio is %10.6f\n\n",
	       (seq->ttratio));
    }
  // use 1/frac as precomputation speed up
  seq->fracchange = 1. / (
			  (seq->xi) * (2. * freqa * (freqgr) +
				       2. * freqc * (freqty)) + (seq->xv) * (1.0 -
										      freqa *
										      freqa -
										      freqc *
										      freqc -
										      freqg *
										      freqg -
										      freqt *
										      freqt));
}

/*===================================================*/

void makeweights_old (world_fmt * world, data_fmt * data, option_fmt *options, long locus)
{
    /* make up weights vector to avoid duplicate computations */
    long i;
    seqmodel_fmt *seq = world->data->seq[0];
    world->data->seq[0]->endsite = 1;
    for (i = 0; i < seq->sites[locus]; i++)
    {
        seq->alias[i] = i + 1;
        seq->ally[i] = 0;
        seq->aliasweight[i] = seq->weight[i];
        seq->location[i] = 0;
    }
    sitesort2_old (world, data, options, seq->sites[locus], locus);
    sitecombine2_old (world, data, options, seq->sites[locus], locus);
    sitescrunch2_old (world, seq->sites[locus], 1, 2, locus);
    for (i = 1; i <= seq->sites[locus]; i++)
    {
        if (seq->aliasweight[i - 1] > 0)
            seq->endsite = i;
    }
    for (i = 1; i <= seq->endsite; i++)
    {
        seq->location[seq->alias[i - 1] - 1] = i;
        seq->ally[seq->alias[i - 1] - 1] =
            seq->alias[i - 1];
    }
    init_sequences2_old (world, seq, locus);
    memcpy(seq->savealiasweight, seq->aliasweight, sizeof(long) * (size_t) seq->endsite);
}    /* makeweights */

void
makeweights (world_fmt * world, data_fmt * data, option_fmt *options, long locus)
{
  long sublocus;
  mutationmodel_fmt *s;
  long numsites;
  long i;
  //for sublocus s
  //
  long sublocistart = world->sublocistarts[locus];
  long sublociend = world->sublocistarts[locus+1];
  for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
    { 
      s = &(world->mutationmodels[sublocus]);
      set_subloci_basedefaults(s, world, options, data, sublocus);
      numsites = s->numsites;
      s->numpatterns = 1;
      for (i = 0; i < numsites; i++)
	{
	  s->alias[i] = i + 1;
	  s->ally[i] = 0;
	  s->aliasweight[i] = s->weight[i];
	  s->location[i] = 0;
	}
      sitesort2 (world, data, options, numsites, locus, sublocus);
      sitecombine2 (world, data, options, numsites, locus, sublocus);
      sitescrunch2 (world, numsites, 1, 2, sublocus);
      for (i = 1; i <= numsites; i++)
	{
	  if (s->aliasweight[i - 1] > 0)
            s->numpatterns = i;
	}
      for (i = 1; i <= s->numpatterns; i++)
	{
	  s->location[s->alias[i - 1] - 1] = i;
	  s->ally[s->alias[i - 1] - 1] =
            s->alias[i - 1];
	}
      init_sequences2 (world, s);
      memcpy(s->savealiasweight, s->aliasweight, sizeof(long) * (size_t) s->numpatterns);
      set_siterates(sublocus,world,options);
    }
}    /* makeweights */

void
init_sequences2_old (world_fmt * world, seqmodel_fmt * seq, long locus)
{
  (void) locus;
    if (world->contribution == NULL)
      {
        world->contribution =
            (contribarr *) mymalloc ((4 + seq->endsite) * sizeof (contribarr));
	//	printf("%i> temp=%f size=%li: world->contribution allocated\n",myID, world->heat,(4 + seq->endsite) * sizeof (contribarr));
      }
    else
      {
        world->contribution =
	  (contribarr *) myrealloc (world->contribution,
                                    (4 + seq->endsite) * sizeof (contribarr));
	//	printf("%i> temp=%f size=%li: world->contribution REallocated\n",myID, world->heat,(4 + seq->endsite) * sizeof (contribarr));
      }
}
void
init_sequences2 (world_fmt * world, mutationmodel_fmt * s)
{
  (void) world;

    if (s->contribution == NULL)
      {
        s->contribution =
            (contribarr *) mymalloc ((4 + s->numpatterns) * sizeof (contribarr));
	//printf("%i> temp=%f size=%li: world->contribution allocated\n",myID, world->heat,(4 + s->numpatterns) * sizeof (contribarr));
      }
    else
      {
        s->contribution =
	  (contribarr *) myrealloc (s->contribution,(4 + s->numpatterns) * sizeof (contribarr));
	//printf("%i> temp=%f size=%li: world->contribution REallocated\n",myID, world->heat,(4 + s->numpatterns) * sizeof (contribarr));
      }
}


// set the conditional tip likelihoods using sequencing error 
// an errorrate for each nucleotide must be given, no error a vector of zeroes is fine.
//
void set_nucleotide(MYREAL *treedata, const char nucleotide, const MYREAL *seqerr)
{
  long b;

  const MYREAL a  = seqerr[NUC_A];
  const MYREAL c  = seqerr[NUC_C];
  const MYREAL g  = seqerr[NUC_G];
  const MYREAL t  = seqerr[NUC_T];
  const MYREAL ac = a + c;
  const MYREAL ag = a + g;
  const MYREAL at = a + t;
  const MYREAL cg = c + g;
  const MYREAL ct = c + t; 
  const MYREAL gt = g + t;
  
  treedata[NUC_A] = a;
  treedata[NUC_C] = c;
  treedata[NUC_G] = g;
  treedata[NUC_T] = t;

  switch (nucleotide)
    {
    case 'A':
      treedata[NUC_A] = 1.0 - cg - t;
      break;
      
    case 'C':
      treedata[NUC_C] = 1.0 - gt - a;
      break;
      
    case 'G':
      treedata[NUC_G] = 1.0 - ac - t;
      break;
      
    case 'T':
    case 'U':
      treedata[NUC_T] = 1.0 - ac - g;
      break;
      
    case 'M':
      treedata[NUC_A] = 1.0 - gt;
      treedata[NUC_C] = 1.0 - gt; 
      break;
      
    case 'R':
      treedata[NUC_A] = 1.0 - ct;
      treedata[NUC_G] = 1.0 - ct;
      break;
      
    case 'W':
      treedata[NUC_A] = 1.0 - cg;
      treedata[NUC_T] = 1.0 - cg;
      break;
      
    case 'S':
      treedata[NUC_C] = 1.0 - at;
      treedata[NUC_G] = 1.0 - at;
      break;
      
    case 'Y':
      treedata[NUC_C] = 1.0 - ag;
      treedata[NUC_T] = 1.0 - ag;
      break;
      
    case 'K':
      treedata[NUC_G] = 1.0 - ac;
      treedata[NUC_T] = 1.0 - ac;
      break;
      
    case 'B':
      //treedata[NUC_A] = a; //seqerr;
      treedata[NUC_C] = 1.0 - a; // oneseqerr3;
      treedata[NUC_G] = 1.0 - a; //oneseqerr3;
      treedata[NUC_T] = 1.0 - a; //oneseqerr3;
      break;
      
    case 'D':
      treedata[NUC_A] = 1.0 - c; //oneseqerr3;
      //treedata[NUC_C] = c;//seqerr;
      treedata[NUC_G] = 1.0 - c; //oneseqerr3;
      treedata[NUC_T] = 1.0 - c; //oneseqerr3;
      break;
      
    case 'H':
      treedata[NUC_A] = 1.0 - g;
      treedata[NUC_C] = 1.0 - g;
	//treedata[NUC_G] = g
      treedata[NUC_T] = 1.0 - g;
      break;
      
    case 'V':
      treedata[NUC_A] = 1.0 - t;
      treedata[NUC_C] = 1.0 - t;
      treedata[NUC_G] = 1.0 - t;
      //treedata[NUC_T] = t;
      break;
      
    case 'N':
      for (b = 0; b < 4; b++)
	treedata[b] = 1.0;
      break;
      
    case 'X':
      for (b = 0; b < 4; b++)
	treedata[b] = 1.0;
      break;
      
    case '?':
      for (b = 0; b < 4; b++)
	treedata[b] = 1.0;
      break;
      
    case 'O':
      for (b = 0; b < 4; b++)
	treedata[b] = 1.0;
      break;
#ifdef GAP                            
    case '-':
      treedata[4] = 1.0;
      break;
#else
    case '-':
      for (b = 0; b < 4; b++)
	treedata[b] = 1.0;
      break;
#endif
    }
}



void
sitesort2_old (world_fmt * world, data_fmt * data, option_fmt *options, long sites, long locus)
{
  long gap, i, i1, j, jj, jg, k, kk, kkk, itemp, pop, z = 0;
    boolean flip, tied, samewt;
    seqmodel_fmt *seq;
    long *tempsum, *temppop;
    long numind;
    MYREAL a1,a2;
    tempsum = (long *) mycalloc (data->numpop, sizeof (long));
    temppop = (long *) mycalloc (data->numpop, sizeof (long));
    for (i = 0; i < data->numpop; i++)
    {
      if(options->randomsubset > 0)
	{
	  numind = (options->randomsubset < data->numind[i][locus] ? options->randomsubset : data->numind[i][locus]);
	}
      else
	{
	  numind = data->numind[i][locus];
	}
      if (numind > 0)
        {
	  temppop[z] = i;
	  if (z == 0)
	    tempsum[z] = numind;
	  else
	    tempsum[z] = tempsum[z - 1] + numind;
            z++;
        }
    }
    seq = world->data->seq[0];
    gap = sites / 2;
    while (gap > 0)
    {
        for (i = gap + 1; i <= sites; i++)
        {
            j = i - gap;
            flip = TRUE;
            while (j > 0 && flip)
            {
                jj = seq->alias[j - 1];
                jg = seq->alias[j + gap - 1];
                samewt = ((seq->weight[jj - 1] != 0)
                          && (seq->weight[jg - 1] != 0))
                         || ((seq->weight[jj - 1] == 0) && (seq->weight[jg - 1] == 0));
                tied = samewt
                       && (seq->category[jj - 1] == seq->category[jg - 1]);
                flip = ((!samewt) && (seq->weight[jj - 1] == 0))
                       || (samewt
                           && (seq->category[jj - 1] > seq->category[jg - 1]));
                k = 0;
                pop = 0;
                kk = -1;
                while (k < world->sumtips && tied)
                {
                    if (k == tempsum[pop])
                    {
                        kk = 0;
                        pop++;
                    }
                    else
                    {
                        kk++;
                    }
		    i1 = temppop[pop];
		    kkk = data->shuffled[i1][locus][kk];
		    a1 = data->yy[i1][kkk][locus][0][jj - 1][0];
		    a2 = data->yy[i1][kkk][locus][0][jg - 1][0];
                    flip =
                        (a1 > a2);
                    tied = (tied
                            && (fabs(a1-a2) <= (double) FLT_EPSILON)); 
                    k++;
                }
                if (!flip)
                    break;
                itemp = seq->alias[j - 1];
                seq->alias[j - 1] = seq->alias[j + gap - 1];
                seq->alias[j + gap - 1] = itemp;
                itemp = (long) seq->aliasweight[j - 1];
                seq->aliasweight[j - 1] = seq->aliasweight[j + gap - 1];
                seq->aliasweight[j + gap - 1] = itemp;
                j -= gap;
            }
        }
        gap /= 2;
    }
    myfree(tempsum);
    myfree(temppop);
}    /* sitesort2 */

void
sitesort2 (world_fmt * world, data_fmt * data, option_fmt *options, long sites, long locus, long sublocus)
{
  long gap, i, i1, j, jj, jg, k, kk, kkk, itemp, pop, z = 0;
  boolean flip, tied, samewt;
  mutationmodel_fmt *s;
  long *tempsum, *temppop;
  long numind;
  MYREAL a1,a2;
  long sumtips = 0;
  tempsum = (long *) mycalloc (data->numpop, sizeof (long));
  temppop = (long *) mycalloc (data->numpop, sizeof (long));
  for (i = 0; i < data->numpop; i++)
    {

      if(options->randomsubset > 0)
	{
	  numind = (options->randomsubset < data->numind[i][locus] ? options->randomsubset : data->numind[i][locus]);
	}
      else
	{
	  numind = data->numind[i][locus];
	}
      //      if (numind > 0)
      //  {
	  temppop[z] = i;
	  if (z == 0)
	    tempsum[z] = numind;
	  else
	    tempsum[z] = tempsum[z - 1] + numind;
          z++;
	    //  }
	  sumtips += numind;
    }
  if (world->sumtips == 0 && sumtips > 0)
    world->sumtips = sumtips;
  
  s = &world->mutationmodels[sublocus];
    //    seq = world->data->seq[0];
  gap = sites / 2;
    while (gap > 0)
    {
        for (i = gap + 1; i <= sites; i++)
        {
            j = i - gap;
            flip = TRUE;
            while (j > 0 && flip)
            {
                jj = s->alias[j - 1];
                jg = s->alias[j + gap - 1];
                samewt = ((s->weight[jj - 1] != 0)
                          && (s->weight[jg - 1] != 0))
                         || ((s->weight[jj - 1] == 0) && (s->weight[jg - 1] == 0));
                tied = samewt
                       && (s->category[jj - 1] == s->category[jg - 1]);
                flip = ((!samewt) && (s->weight[jj - 1] == 0))
                       || (samewt
                           && (s->category[jj - 1] > s->category[jg - 1]));
                k = 0;
                pop = 0;
                kk = -1;
                while (k < sumtips && tied)
                {
                    if (k == tempsum[pop])
                    {
                        kk = 0;
                        pop++;
                    }
                    else
                    {
                        kk++;
                    }
		    i1 = temppop[pop];
		    //debug1
		    if (data->numind[i1][locus]>0)
		      {
			kkk = data->shuffled[i1][locus][kk];
			a1 = data->yy[i1][kkk][sublocus][0][jj - 1][0];
			a2 = data->yy[i1][kkk][sublocus][0][jg - 1][0];
			flip =
			  (a1 > a2);
			tied = (tied
				&& (fabs(a1-a2) <= (double) FLT_EPSILON));
		      }
                    k++;
                }
                if (!flip)
                    break;
                itemp = s->alias[j - 1];
                s->alias[j - 1] = s->alias[j + gap - 1];
                s->alias[j + gap - 1] = itemp;
                itemp = (long) s->aliasweight[j - 1];
                s->aliasweight[j - 1] = s->aliasweight[j + gap - 1];
                s->aliasweight[j + gap - 1] = itemp;
                j -= gap;
            }
        }
        gap /= 2;
    }
    myfree(tempsum);
    myfree(temppop);
}    /* sitesort2 */


void
sitecombine2_old (world_fmt * world, data_fmt * data, option_fmt *options, long sites, long locus)
{
  long i, i1, j, k, kk, pop, z = 0;
    boolean tied, samewt;
    seqmodel_fmt *seq;
    long *tempsum, *temppop;
    long numind;
    long kkk;
    tempsum = (long *) mycalloc (data->numpop, sizeof (long));
    temppop = (long *) mycalloc (data->numpop, sizeof (long));
    if(options->randomsubset > 0)
      {
	tempsum[0] = (options->randomsubset < data->numind[0][locus] ? options->randomsubset : data->numind[0][locus]);
      }
    else
      {
	tempsum[0] = data->numind[0][locus];
      }
    
    for (i = 0; i < data->numpop; i++)
    {
      numind = ((options->randomsubset > 0) && (options->randomsubset < data->numind[i][locus])) ? options->randomsubset : data->numind[i][locus];
      //debug2
      //      if (numind > 0)
      //  {
	  temppop[z] = i;
	  if (z == 0)
	    tempsum[z] = numind;
	  else
                tempsum[z] = tempsum[z - 1] + numind;
            z++;
	    //  }
    }

    seq = world->data->seq[0];
    i = 1;
    while (i < sites)
    {
        j = i + 1;
        tied = TRUE;
        while (j <= sites && tied)
        {
	  samewt = ((seq->aliasweight[i - 1] > (double) FLT_EPSILON)
                      && (seq->aliasweight[j - 1] > (double) FLT_EPSILON))
                     || ((seq->aliasweight[i - 1] < (double) FLT_EPSILON)
                         && (seq->aliasweight[j - 1] < (double) FLT_EPSILON));
            tied = samewt
	      && (seq->category[seq->alias[i - 1] - 1] ==
			seq->category[seq->alias[j - 1] - 1]);
            k = 0;
            pop = 0;
            kk = -1;
            while (k < world->sumtips && tied)
            {
                if (k == tempsum[pop])
                {
                    kk = 0;
                    pop++;
                }
                else
                {
                    kk++;
                }
		i1 = temppop[pop];
		if (data->numind[i1][locus]>0)
		  {
		    kkk = data->shuffled[i1][locus][kk];
		    tied = (tied
			    && data->yy[i1][kkk][locus][0][seq->
							   alias[i - 1] -
							   1][0] ==
			    data->yy[i1][kkk][locus][0][seq->alias[j - 1] -
							1][0]);
		  }
		k++;
	    }
	    if (!tied)
	      break;
	    seq->aliasweight[i - 1] += seq->aliasweight[j - 1];
	    seq->aliasweight[j - 1] = 0;
	    seq->ally[seq->alias[j - 1] - 1] = seq->alias[i - 1];
	    j++;
	}
	i = j;
    }
    myfree(temppop);
    myfree(tempsum);
}    /* sitecombine2 */


void
sitecombine2 (world_fmt * world, data_fmt * data, option_fmt *options, long sites, long locus, long sublocus)
{
    long i, j, k, kk, pop, z = 0;
    boolean tied, samewt;
    //seqmodel_fmt *seq;
    long *tempsum, *temppop;
    long numind;
    long kkk;
    mutationmodel_fmt *s;
    long sumtips = 0;
    tempsum = (long *) mycalloc (1, sizeof (long) * (size_t) data->numpop);
    temppop = (long *) mycalloc (1, sizeof (long) * (size_t) data->numpop);
    if(options->randomsubset > 0)
      {
	tempsum[0] = (options->randomsubset < data->numind[0][locus] ? options->randomsubset : data->numind[0][locus]);
      }
    else
      {
	tempsum[0] = data->numind[0][locus];
      }
    
    for (i = 0; i < data->numpop; i++)
    {
      numind = ((options->randomsubset > 0) && (options->randomsubset < data->numind[i][locus])) ? options->randomsubset : data->numind[i][locus];
      if (numind > 0)
        {
	  temppop[z] = i;
	  if (z == 0)
	    tempsum[z] = numind;
	  else
                tempsum[z] = tempsum[z - 1] + numind;
            z++;
        }
      sumtips += numind;
    }
    s = &world->mutationmodels[sublocus];
    //seq = world->data->seq[0];
    i = 1;
    while (i < sites)
    {
        j = i + 1;
        tied = TRUE;
        while (j <= sites && tied)
        {
            samewt = ((s->aliasweight[i - 1] > (double) FLT_EPSILON)
                      && (s->aliasweight[j - 1] > (double) FLT_EPSILON))
                     || ((s->aliasweight[i - 1] < (double) FLT_EPSILON)
                         && (s->aliasweight[j - 1] < (double) FLT_EPSILON));
            tied = samewt
	      && (s->category[s->alias[i - 1] - 1] ==
		  s->category[s->alias[j - 1] - 1]);
            k = 0;
            pop = 0;
            kk = -1;
            while (k < sumtips && tied)
            {
                if (k == tempsum[pop])
                {
                    kk = 0;
                    pop++;
                }
                else
                {
                    kk++;
                }
		kkk = data->shuffled[pop][locus][kk];
                tied = (tied
                        && data->yy[temppop[pop]][kkk][sublocus][0][s->
                                                                alias[i - 1] -
                                                                1][0] ==
                        data->yy[temppop[pop]][kkk][sublocus][0][s->alias[j - 1] -
                                                             1][0]);
                k++;
            }
            if (!tied)
                break;
            s->aliasweight[i - 1] += s->aliasweight[j - 1];
            s->aliasweight[j - 1] = 0;
            s->ally[s->alias[j - 1] - 1] = s->alias[i - 1];
            j++;
        }
        i = j;
    }
    myfree(temppop);
    myfree(tempsum);
}    /* sitecombine2 */


void
sitescrunch2_old (world_fmt * world, long sites, long i, long j, long locus)
{
  (void) locus;
    /* move so positively weighted sites come first */
    long itemp;
    boolean done, found;
    seqmodel_fmt *seq;
    seq = world->data->seq[0];
    done = FALSE;
    while (!done)
    {
        //found = FALSE;
        if (seq->aliasweight[i - 1] > 0)
            i++;
        else
        {
            if (j <= i)
                j = i + 1;
            if (j <= sites)
            {
                //found = FALSE;
                do
                {
                    found = (seq->aliasweight[j - 1] > 0);
                    j++;
                }
                while (!(found || j > sites));
                if (found)
                {
                    j--;
                    itemp = seq->alias[i - 1];
                    seq->alias[i - 1] = seq->alias[j - 1];
                    seq->alias[j - 1] = itemp;
                    itemp = (long) seq->aliasweight[i - 1];
                    seq->aliasweight[i - 1] = seq->aliasweight[j - 1];
                    seq->aliasweight[j - 1] = itemp;
                }
                else
                    done = TRUE;
            }
            else
                done = TRUE;
        }
        done = (done || i >= sites);
    }
}    /* sitescrunch2 */

void
sitescrunch2 (world_fmt * world, long sites, long i, long j, long sublocus)
{
    /* move so positively weighted sites come first */
    long itemp;
    boolean done, found;
    //seqmodel_fmt *seq;
    mutationmodel_fmt *s;
    //    seq = world->data->seq[0];
    s = &world->mutationmodels[sublocus];
    done = FALSE;
    while (!done)
    {
      //found = FALSE;
        if (s->aliasweight[i - 1] > 0)
            i++;
        else
        {
            if (j <= i)
                j = i + 1;
            if (j <= sites)
            {
                //found = FALSE;
                do
                {
                    found = (s->aliasweight[j - 1] > 0);
                    j++;
                }
                while (!(found || j > sites));
                if (found)
                {
                    j--;
                    itemp = s->alias[i - 1];
                    s->alias[i - 1] = s->alias[j - 1];
                    s->alias[j - 1] = itemp;
                    itemp = (long) s->aliasweight[i - 1];
                    s->aliasweight[i - 1] = s->aliasweight[j - 1];
                    s->aliasweight[j - 1] = itemp;
                }
                else
                    done = TRUE;
            }
            else
                done = TRUE;
        }
        done = (done || i >= sites);
    }
}    /* sitescrunch2 */

void
inputoptions (world_fmt * world, option_fmt * options, data_fmt * data,
              long locus)
{
    long i;
    mutationmodel_fmt * s;
    long numsites;
    long sublocusstart;
    long sublocusend;
    long sublocus;
    sublocusstart = world->sublocistarts[locus];
    sublocusend   = world->sublocistarts[locus+1];

    for(sublocus=sublocusstart; sublocus < sublocusend; sublocus++)
      {
	s = &world->mutationmodels[sublocus];
	numsites = s->numsites;
	for (i = 0; i < numsites; i++)
	  {
	    s->category[i] = 1;
	    s->weight[i] = 1;
	  }
	if (options->weights)
	  inputweights (s, data, numsites);
	s->weightsum = 0;
	for (i = 0; i < numsites; i++)
	  s->weightsum += s->weight[i];
	if (options->categs > 1)
	  {
	    inputcategs (0, numsites, s, options, data);
	  }
      }
}    /* inputoptions */

void
inputweights (mutationmodel_fmt * s, data_fmt * data, long chars)
{
    /* input the character weights, 0-9 and A-Z for weights 0 - 35 */
    char ch;
    long i;
	char input[1024];

	ch = (char) getc (data->weightfile);
    while (ch == '#')
    {
        FGETS(input, LINESIZE,data->weightfile);
        ch = (char) getc (data->weightfile);
    }
    ungetc (ch, data->weightfile);
    for(i=0;i<chars;i++)
      {
	s->weight[i] = 1;
        ch = (char) getc (data->weightfile);
	
	if (isdigit ((int) ch))
	  s->weight[i] = ch - '0';
	else 
	  {
	    if (isalpha ((int) ch))
	      {
		ch = uppercase (ch);
		s->weight[i] = (short) (ch - 'A' + 10);
	      }
	    else
	      {
		if(isspace((int) ch))
		  {
		    i--;
		    continue;
		  }
		else
		  printf ("ERROR: Bad weight character: %c\n", ch);
		exit (EXIT_FAILURE);
	      }
	  }
      }
}    /* inputweights */

void inputcategs (long a, long b, mutationmodel_fmt * s, option_fmt * options, data_fmt * data)
{
    /* input the categories, 1-9 */
    char ch;
    long i;
    char input[LINESIZE];
    //int retval;
    ch = (char) getc (data->catfile);
    while (ch == '#')
    {
        FGETS(input, LINESIZE,data->catfile);
        ch = (char) getc (data->catfile);
    }
    ungetc (ch, data->catfile);
    fscanf (data->catfile, "%ld", &options->categs);
    s->numcategs = options->categs;
    s->rate =
      (MYREAL *) myrealloc (s->rate, sizeof (MYREAL) * (size_t) s->numcategs);
    for (i = 0; i < options->categs; i++)
    {
#ifdef USE_MYREAL_FLOAT
        fscanf (data->catfile, "%f", &s->rate[i]);
#else
        fscanf (data->catfile, "%lf", &s->rate[i]);
#endif
    }

    for (i = a; i < b; i++)
    {
        if ((ch >= '1') && (ch <= ('0' + options->categs)))
            s->category[i] = ch - '0';
        else
        {
			if(isspace((int) ch))
			  {
				i--;
				continue;
			  }
			else
			  {
				printf
            ("BAD CATEGORY CHARACTER: %c -- CATEGORIES ARE CURRENTLY 1-%ld\n",
             ch, options->categs);
            exit (EXIT_FAILURE);
			  }
        }
    }
}    /* inputcategs */




void
empiricalfreqs (world_fmt * world, option_fmt * options, mutationmodel_fmt * s,
                long sublocus)
{
  (void) options;
  
    /* Get empirical base frequencies from the data */
    long i, j;
    MYREAL summ, suma, sumc, sumg, sumt, w;
    MYREAL freqa, freqc, freqg, freqt;
    long snps = s->datatype == 'u' ? 5 : 1;
    xarray_fmt xx;
    freqa = 0.25;
    freqc = 0.25;
    freqg = 0.25;
    freqt = 0.25;
    suma = 0.0;
    sumc = 0.0;
    sumg = 0.0;
    sumt = 0.0;
    for (i = 0; i < world->sumtips; i++)
      {
	xx = world->nodep[i]->x[sublocus];
	for (j = 0; j < s->numpatterns * snps; j += snps)
	  {
	    w = (MYREAL) s->aliasweight[j / snps];
	    summ = (freqa) * xx.s[j][0][NUC_A] + (freqc) * xx.s[j][0][NUC_C]
	      + (freqg) * xx.s[j][0][NUC_G]
	      + (freqt) * xx.s[j][0][NUC_T];
	    suma += w * (freqa) * xx.s[j][0][NUC_A] / summ;
	    sumc += w * (freqc) * xx.s[j][0][NUC_C] / summ;
	    sumg += w * (freqg) * xx.s[j][0][NUC_G] / summ;
	    sumt += w * (freqt) * xx.s[j][0][NUC_T] / summ;
	    
	  }
      }
    if(suma<EPSILON)
      {
	suma += EPSILON; 
      }
    if(sumc<EPSILON)
      {
	sumc += EPSILON; 
      }
    if(sumg<EPSILON)
      {
	sumg += EPSILON; 
      }
    if(sumt<EPSILON)
      {
	sumt += EPSILON; 
      }
    summ = suma + sumc + sumg + sumt;
    freqa = suma / summ;
    freqc = sumc / summ;
    freqg = sumg / summ;
    freqt = sumt / summ;
    s->basefreqs[NUC_A] = freqa;
    s->basefreqs[NUC_C] = freqc;
    s->basefreqs[NUC_G] = freqg;
    s->basefreqs[NUC_T] = freqt;
}    /* empiricalfreqs */


void
initlambda (option_fmt * options)
{
  //int retval;
    while (options->lambda <= 1.0)
    {
        printf
        ("Mean block length of sites having the same rate\n (needs to be greater than 1)?\n");
#ifdef USE_MYREAL_FLOAT
        scanf ("%f%*[^\n]", &options->lambda);
#else
        scanf ("%lf%*[^\n]", &options->lambda);
#endif
        getchar ();
    }
    options->lambda = 1.0 / options->lambda;
}



void init_tbl (world_fmt * world, long locus)
{
    /* Define a lookup table. Precompute values and print them out in tables */
    long i, j;
    MYREAL sumrates;
    long categs = world->options->categs;
    long rcategs;// = world->options->rcategs;
    mutationmodel_fmt *s;
    long sublocusstart;
    long sublocusend;
    long sublocus;
    sublocusstart = world->sublocistarts[locus];
    sublocusend   = world->sublocistarts[locus+1];

    for(sublocus=sublocusstart; sublocus < sublocusend; sublocus++)
      {
	s = &world->mutationmodels[sublocus];
	if(strchr(SEQUENCETYPES,s->datatype))
	  {
	    rcategs = s->numsiterates;
	    s->tbl = (valrec ***) mymalloc (rcategs * sizeof (valrec **));
	    for (i = 0; i < rcategs; i++)
	      {
		s->tbl[i] = (valrec **) mymalloc (categs * sizeof (valrec *));
		for (j = 0; j < categs; j++)
		  s->tbl[i][j] = (valrec *) mymalloc (sizeof (valrec));
	      }
	    for (i = 0; i < rcategs; i++)
	      {
		for (j = 0; j < categs; j++)
		  {
		    s->tbl[i][j]->rat = s->siterates[i] * world->options->rate[j];
		    s->tbl[i][j]->ratxi = s->tbl[i][j]->rat * s->xi;
		    s->tbl[i][j]->ratxv = s->tbl[i][j]->rat * s->xv;
		  }
	      }
	    sumrates = 0.0;
	    for (i = 0; i < s->numpatterns; i++)
	      {
		long ccc = s->numcategs > 0 ? (s->category[s->alias[i] - 1] - 1) : 0;
		for (j = 0; j < rcategs; j++)
		  sumrates +=
		    s->aliasweight[i] * s->siteprobs[j] *
		    s->tbl[j][ccc]->rat;
	      }
	    sumrates /= (MYREAL) s->numsites;
	    for (i = 0; i < rcategs; i++)
	      for (j = 0; j < categs; j++)
		{
		  s->tbl[i][j]->rat /= sumrates;
		  s->tbl[i][j]->ratxi /= sumrates;
		  s->tbl[i][j]->ratxv /= sumrates;
		}
	  }
      }
}    /* inittable */

void
print_weights (FILE * outfile, world_fmt * world, option_fmt * options,
               long locus)
{
    if (options->weights)
    {
        if ((options->printdata) || (options->progress && outfile == stdout))
        {
            printweights (outfile, world, options, 0,
                          world->data->seq[0]->sites[locus],
                          world->data->seq[0]->weight, "Sites");
        }
    }
}



void
printweights (FILE * outfile, world_fmt * world, option_fmt * options,
              short inc, long chars, short *weight, char *letters)
{
  (void) world;
    /* print out the weights of sites */
    long i, j;
    FPRINTF (outfile, "\n    %s are weighted as follows:\n", letters);
    for (i = 0; i < chars; i++)
    {
        if (i % 60 == 0)
        {
            FPRINTF (outfile, "\n");
            for (j = 1; j <= options->nmlength + 3; j++)
               FPRINTF (outfile, " ");
        }
        FPRINTF (outfile, "%hd", weight[i + inc]);
        if ((i + 1) % 10 == 0 && (i + 1) % 60 != 0)
            FPRINTF (outfile, " ");
    }
    FPRINTF (outfile, "\n\n");
}    /* printweights */



void
print_seqfreqs (FILE * outfile, world_fmt * world, option_fmt * options)
{
  mutationmodel_fmt *s;
  if(outfile==NULL)
    return;
  if (options->freqsfrom)
    FPRINTF (outfile, "\nEmpirical ");
  FPRINTF (outfile, "Base Frequencies\n");
  FPRINTF (outfile,
	   "------------------------------------------------------------\n");
  FPRINTF (outfile,
	   "Locus     Sublocus  Nucleotide                        Model parameters/\n");
  FPRINTF (outfile,
	   "                    ------------------------------ \n");
  FPRINTF (outfile, "                    A       C       G       T(U)\n");
  FPRINTF (outfile,
	   "----------------------------------------------------------------------\n");
  long locus;
  long sublocus;
  long startlocus;
  long endlocus;
  if(outfile==stdout || outfile==world->options->logfile)
    {
      startlocus=world->locus;
      endlocus=world->locus+1;
    }
  else
    {
      startlocus=0;
      endlocus=world->loci;
    }
  for(locus=startlocus; locus < endlocus; locus++)
    {
      long sublocistart = world->sublocistarts[locus];
      long sublociend   = world->sublocistarts[locus+1];
      for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
	{
	  s = &world->mutationmodels[sublocus];

	  if(strchr(SEQUENCETYPES, s->datatype))
	    {
	      FPRINTF (outfile, "%4li      %4li      %6.4f  %6.4f  %6.4f  %6.4f\n",
		       locus + 1, sublocus + 1, s->basefreqs[NUC_A],s->basefreqs[NUC_C],s->basefreqs[NUC_G],s->basefreqs[NUC_T]
		       );
	    }
	}
    }
  if (outfile == stdout)
    FPRINTF (outfile, "\n");
}

MYREAL
treelike_seq (mutationmodel_fmt *s, long sublocus, world_fmt * world, long locus)
{
  (void) locus;
    const MYREAL freqa = s->basefreqs[NUC_A];
    const MYREAL freqc = s->basefreqs[NUC_C];
    const MYREAL freqg = s->basefreqs[NUC_G];
    const MYREAL freqt = s->basefreqs[NUC_T];
    const long numpatterns= s->numpatterns;
    const long rcategs = s->numsiterates;
    contribarr tterm;
    contribarr like;
    contribarr nulike;
    contribarr clai;

    MYREAL summ    = 0.0;
    MYREAL sum2    = 0.0;
    MYREAL sumc    = 0.0;
    MYREAL sumterm = 0.0; 
    MYREAL lterm;
    long i, j, k, lai;
    MYREAL scale;
    node *p;
    sitelike *x1;
    p = crawlback (world->root->next);
    summ = 0.0;

    if (s->numsiterates == 1)
    {
        for (i = 0; i < numpatterns; i++)
        {
            x1 = &(p->x[sublocus].s[i][0]);
            scale = p->scale[sublocus][i];
            tterm[0] =
                freqa * (*x1)[0] + freqc * (*x1)[1] +
                freqg * (*x1)[2] + freqt * (*x1)[3];
            summ += s->aliasweight[i] * (LOG (tterm[0]) + scale);
        }
    }
    else
    {
        for (i = 0; i < numpatterns; i++)
        {
            sumterm = 0.0;
            scale = p->scale[sublocus][i];
            for (j = 0; j < rcategs; j++)
            {
                x1 = &(p->x[sublocus].s[i][j]);
                tterm[j] =
                    freqa * (*x1)[0] + freqc * (*x1)[1] +
                    freqg * (*x1)[2] + freqt * (*x1)[3];
                sumterm += s->siteprobs[j] * tterm[j];
	    }
            lterm = LOG (sumterm) + scale;
            for (j = 0; j < s->numsiterates; j++)               
	      clai[j] = tterm[j] / sumterm;
            swap (&clai, &s->contribution[i]);
            summ += s->aliasweight[i] * lterm;
            if(MYISNAN((float)summ))
	      {
		// error("summ is not a number, should not happen\n");
		summ = (double) -HUGE;
	      }
        }
        for (j = 0; j < s->numsiterates; j++)
            like[j] = 1.0;
        for (i = 0; i < s->numsites; i++)
        {
            sumc = 0.0;
            for (k = 0; k < s->numsiterates; k++)
                sumc += s->siteprobs[k] * like[k];
            sumc *= s->lambda;
            if ((s->ally[i] > 0) && (s->location[s->ally[i] - 1] > 0))
            {
                lai = s->location[s->ally[i] - 1];
                swap (&s->contribution[lai - 1], &clai);
                //memcpy (clai, world->contribution[lai - 1], size);
                for (j = 0; j < s->numsiterates; j++)
                    nulike[j] = ((1.0 - s->lambda) * like[j] + sumc) * clai[j];
            }
            else
            {
                for (j = 0; j < s->numsiterates; j++)
                    nulike[j] = ((1.0 - s->lambda) * like[j] + sumc);
            }
            swap (&nulike, &like);
            //memcpy (like, nulike, size);
        }
        sum2 = 0.0;
        for (i = 0; i < s->numsiterates; i++)
            sum2 += s->siteprobs[i] * like[i];
        summ += LOG (sum2);
    }
    //printf("liketerm(summ)=%f\n",summ);
    return summ;
}    /* treelikelihood */


void
snp_invariants (mutationmodel_fmt *s, long xs, contribarr invariants, world_fmt *world, long locus, phenotype x1, MYREAL *scale)
{
  (void) xs;
  (void) world;
  (void) locus;

    //worldoption_fmt *opt;
    const long rcategs = s->numsiterates;
    //MYREAL summ;
    const long numpatterns = s->numpatterns;
    const long addon =  s->numpatterns + s->addon;
    sitelike *val;
    long i, j;
    register MYREAL freqa, freqc, freqg, freqt;
    freqa = s->basefreqs[NUC_A];
    freqc = s->basefreqs[NUC_C];
    freqg = s->basefreqs[NUC_G];
    freqt = s->basefreqs[NUC_T];

    //opt = world->options;
    //summ = 0.0;
    /* snp invariants  set to 4 because s->addon may be zero */
    MYREAL invariant_vec[FOUR]={0.0};
    for (j = 0; j < rcategs; j++)
      {
	invariants[j] = 0.0;
	MYREAL max_invariant = (double) -HUGE;
	long z;
	for (i = numpatterns; i < addon; i++)
	  {
	    z = i-numpatterns;
	    val = &(x1[i][j]);
	    invariant_vec[z] = log(freqa * (*val)[0] + freqc * (*val)[1] +
					       freqg * (*val)[2] + freqt * (*val)[3]) + scale[i];
	    if(max_invariant < invariant_vec[z])
	      {
		max_invariant = invariant_vec[z];
	      }
	  }
	for(i = 0; i < s->addon; i++)
	  {
	    invariants[j] += exp(invariant_vec[i] - max_invariant);
	  }
	invariants[j] = 1.0 - (invariants[j] * exp(max_invariant));
      }
    //    printf("invariant0=%f\n",invariants[0]);
}    /* snp_invariants*/


MYREAL
treelike_snp (mutationmodel_fmt *s, long xn, world_fmt * world, long locus)
{
  const long endsite = s->numpatterns;
  const long rcategs = s->numsiterates;
  worldoption_fmt *opt;
  MYREAL scale;
  contribarr tterm;
  contribarr invariants;
  contribarr like;
  contribarr nulike;
  MYREAL summ, sum2, sumc, sumterm, lterm;
  long i, j, k, lai;
  
  node *p;
  sitelike *x1;
  register MYREAL freqa, freqc, freqg, freqt;
  
  freqa = s->basefreqs[NUC_A];
  freqc = s->basefreqs[NUC_C];
  freqg = s->basefreqs[NUC_G];
  freqt = s->basefreqs[NUC_T];
  
  opt = world->options;
  p = crawlback (world->root->next);
  summ = 0.0;
  /* snp invariants */
  snp_invariants (s, xn, invariants,world, locus, p->x[xn].s, p->scale[xn]);
  for (j = 0; j < rcategs; j++)
    {
      if(invariants[j] <= 0.0)
	{
#ifdef DEBUG
	  printf("invariants are funny %f\n", invariants[j]);
#endif
	  return -MYREAL_MAX;
	}
    }
  if(rcategs==1)
    {
      for (i = 0; i < endsite; i++)
	{
	  scale = p->scale[xn][i];
	  x1 = &(p->x[xn].s[i][0]);
	  tterm[0] =
	    (freqa * (*x1)[0] + freqc * (*x1)[1] +
	     freqg * (*x1)[2] + freqt * (*x1)[3]) /invariants[0];
	  // fixing positive logL
	  double tt = (LOG (tterm[0]) + scale);
	  if (tt>0.0)
	    tt=0.0;
	  summ += s->aliasweight[i] * tt;
	}
    }
  else
    {
      for (i = 0; i < endsite; i++)
	{
	  sumterm = 0.0;
	  scale = p->scale[xn][i];
	  for (j = 0; j < rcategs; j++)
	    {
	      x1 = &(p->x[xn].s[i][j]);
	      tterm[j] =
		(freqa * (*x1)[0] + freqc * (*x1)[1] +
		 freqg * (*x1)[2] + freqt * (*x1)[3]) /invariants[j];
	      //fixing positive logL
	      if (tterm[j] > 1.0)
		tterm[j] = 1.0;
	      sumterm += s->siteprobs[j] * tterm[j];
	    }
	  lterm = LOG (sumterm) + scale;
	  for (j = 0; j < rcategs; j++)
	    {
	      //fprintf(stderr,"snp+: %f %f=log(%f)+%f\n",tterm[j],lterm,sumterm, scale);
	      s->contribution[i][j] = tterm[j] / sumterm;
	    }
	  summ += s->aliasweight[i] * lterm;
	}    /* over endsite without the  4[snp-invariants] */
      for (j = 0; j < rcategs; j++)
	like[j] = 1.0;
      for (i = 0; i < endsite; i++)
	{
	  sumc = 0.0;
	  for (k = 0; k < rcategs; k++)
	    sumc += s->siteprobs[k] * like[k];
	  sumc *= opt->lambda;
	  if ((s->ally[i] > 0) && (s->location[s->ally[i] - 1] > 0))
	    {
	      lai = s->location[s->ally[i] - 1];
	      for (j = 0; j < rcategs; j++)
		nulike[j] =
		  ((1.0 - s->lambda) * like[j] +
		   sumc) * s->contribution[lai - 1][j];
	    }
	  else
	    {
	      for (j = 0; j < rcategs; j++)
		nulike[j] = ((1.0 - s->lambda) * like[j] + sumc);
	    }
	  swap (&nulike, &like);
	  //memcpy (like, nulike, size);
	}
      sum2 = 0.0;
      for (i = 0; i < rcategs; i++)
	sum2 += s->siteprobs[i] * like[i];
      summ += LOG (sum2);
    }
  return summ;
}    /* treelike_snp */


MYREAL treelike_snp_unlinked (mutationmodel_fmt *s, long xs, world_fmt * world, long locus)
{
    //worldoption_fmt *opt;

    MYREAL scale;
    contribarr tterm, invariants;
    MYREAL summ, datasum = 0, lterm, result = 0;
    long i, ii;

    node *p;
    sitelike *x1;

    register MYREAL freqa, freqc, freqg, freqt;

    freqa = s->basefreqs[NUC_A];
    freqc = s->basefreqs[NUC_C];
    freqg = s->basefreqs[NUC_G];
    freqt = s->basefreqs[NUC_T];

    //    opt = world->options;
    p = crawlback (world->root->next);
    summ = 0.0;
    /* snp invariants */
    snp_invariants (s, xs,invariants, world, locus, p->x[xs].s, p->scale[xs]);
    /* no rate categories used */
    for (i = 0; i < s->numpatterns; i++)
    {
        ii = i / 5;
        scale = p->scale[xs][i];
        x1 = &(p->x[xs].s[i][0]);
        tterm[0] =
            (freqa * (*x1)[0] + freqc * (*x1)[1] +
             freqg * (*x1)[2] + freqt * (*x1)[3]);
        if (i % 5 == 0)
        {
            lterm = LOG (tterm[0]) + scale;
            summ = 0;
            datasum = s->aliasweight[ii] * lterm;
        }
        else
            summ += pow (tterm[0], (MYREAL) s->aliasweight[ii]);
        if (((i + 1) % 5) == 0 && i != 0)
	  {
	    if(summ > 0.)
	      {
		result +=
		  datasum + LOG ((1 - EXP (LOG (summ) - datasum)) / invariants[0]);
	      }
	  }
    }
    return result;
}    /* treelike_snp_unlinked */

void
check_basefreq (option_fmt * options)
{

    if (options->freqa == 0. || options->freqc == 0. || options->freqt == 0.
            || options->freqg == 0.)
    {
        options->freqa = 0.25;
        options->freqc = 0.25;
        options->freqg = 0.25;
        options->freqt = 0.25;
    }
}


void copy_seq (world_fmt * original, world_fmt * kopie)
{
    size_t sites;
    seqmodel_fmt *kseq;
    seqmodel_fmt *oseq;
    kseq = kopie->data->seq[0];
    oseq = original->data->seq[0];
    sites = (size_t) oseq->sites[original->locus];
    memcpy(kseq->basefrequencies,oseq->basefrequencies,sizeof(MYREAL)*BASEFREQLENGTH);
    kseq->aa = oseq->aa;
    kseq->bb = oseq->bb;
    kseq->endsite = oseq->endsite;
    kseq->xi = oseq->xi;
    kseq->xv = oseq->xv;
    kseq->ttratio = oseq->ttratio;
    kseq->fracchange = oseq->fracchange;
    memcpy (kseq->sites, oseq->sites, sizeof (long) * (size_t) original->loci);
    memcpy (kseq->alias, oseq->alias, sizeof (long) * sites);
    memcpy (kseq->ally, oseq->ally, sizeof (long) * sites);
    memcpy (kseq->category, oseq->category, sizeof (long) * sites);
    memcpy (kseq->weight, oseq->weight, sizeof (short) * sites);
    kseq->weightsum = oseq->weightsum;
    memcpy (kseq->aliasweight, oseq->aliasweight, sizeof (long) * sites);
    memcpy (kseq->location, oseq->location, sizeof (long) * sites);
    kseq->addon = oseq->addon;
}

void find_rates_fromdata(data_fmt * data, option_fmt * options, world_fmt * world)
{
  mutationmodel_fmt *s;
  long sublocus;
  long pop;
  long ind;
  long site;
  long i;
  long n;
  long ss;
  long **v;
  long maxsites = 1;
  long *numind;
  site_fmt *reference;
  site_fmt *indseq;
  long subloci = data->allsubloci;
  long numpop = data->numpop;
  long segreg;
  MYREAL mean=0.;
  MYREAL delta = 0.;
  long locus = world->locus;
  //long newsite;
  numind = (long *) mycalloc (subloci, sizeof (long));
  options->segregs = (long *) mycalloc (subloci, sizeof (long));
  options->wattersons = (MYREAL *) mycalloc (subloci, sizeof (MYREAL));
  for (sublocus=0; sublocus < subloci; sublocus++)
    {
      s = &world->mutationmodels[sublocus];
      maxsites = (s->numsites > maxsites ? s->numsites : maxsites);
      for(pop=0;pop < data->numpop; pop++)
	{
	  numind[sublocus] += data->numind[pop][locus];
	}
    }

  v = (long **) mycalloc (subloci, sizeof (long *));
  v[0] = (long *) mycalloc (subloci * (size_t) maxsites, sizeof (long));
  for (i = 1; i < subloci; i++)
    {
      v[i] = v[0] + maxsites * i;
    }
  
  for (sublocus=0; sublocus< subloci; sublocus++)
    {
      s = &world->mutationmodels[sublocus];
      reference = data->yy[0][0][sublocus][0];
      for(pop=0;pop < numpop; pop++)
	{
	  for(ind=0; ind < data->numind[pop][locus];ind++)
	    {
	      indseq=data->yy[pop][ind][sublocus][0];
	      // if(strcmp(indseq,reference)!=0)
	      //{
		  for(site=0; site < s->numsites; site++)
		    {
		      //ss = 0;
		      if(v[sublocus][site] == 0)
			{
			  ss = (long)(indseq[site][0] != reference[site][0]);
			  v[sublocus][site] += ss;
			}
		    }
		  //	}
	    }
	}
    }
  //n = 0;
  //mean = 0.;
  if(options->mu_rates == NULL)
    {
      options->mu_rates = (MYREAL *) mycalloc( subloci, sizeof(MYREAL));
      //      printf("%i> opt murate size %li\n",myID,subloci * sizeof (MYREAL));
    }
  n = 0;
  mean = 0.;
  options->muloci = subloci;
  for (sublocus=0; sublocus < subloci; sublocus++)
    {
      s = &world->mutationmodels[sublocus];
      segreg = 0;
      for(site=0; site < s->numsites; site++)
	{
	  segreg += (long) (v[sublocus][site] > 0); 
	}
      options->segregs[sublocus] = segreg;
      fprintf(stdout,"sequence.c: segreg=%li\n",segreg);
      if(segreg==0)
	{
	  if (options->onlyvariable)
	    {
	      data->skiploci[locus] = TRUE;
	      data->locusweight[locus] = 0.0;
	    }
	  else
	    {
	      if (options->has_variableandone)
		{
		  if (options->firstinvariant == -1)
		    {
		      options->firstinvariant=locus;
		      data->locusweight[locus] = 1.0;
		    }
		  else
		    {
		      data->skiploci[locus] = TRUE;
		      data->locusweight[locus] = 0.0;
		      data->locusweight[options->firstinvariant] +=1;
		    }
		}
	    }
	}
      // originally I used just the plain watterson estimates but
      // this leads to strange results with highly unequal sequences
      // now the watterson's theta is adjusted per site
      // and that is used to generate the relative rate from the data
      // PB April 3 2011
      options->wattersons[sublocus] = watterson(segreg,numind[locus]);
      options->mu_rates[sublocus] = (options->wattersons[sublocus] + 0.0000001)/s->numsites;
      if (data->seq[0]->sites[locus]>0)
	options->wattersons[sublocus] /= s->numsites;
      n = n + 1;
      delta = options->mu_rates[sublocus] - mean;
      mean += delta/n;
    }
  if(mean > 0.0)
    {
      for (sublocus=0; sublocus < subloci; sublocus++)
	{
	  options->mu_rates[sublocus] /= mean;
	  //s->mu_rate = options->mu_rates[sublocus]; 
	}
    }
  else
    {
      options->murates = FALSE;
      options->murates_fromdata = FALSE;
      warning("Relative mutation rates estimation from data was requested but attempt failed --> reset to all-equal rates");
    }
  myfree(v[0]);
  myfree(v);
  myfree(numind);
}

void find_rates_fromdata_alleles(data_fmt * data, option_fmt * options, world_fmt *world, MYREAL mean)
{
  (void) world;
  long locus;
  options->segregs = (long *) mycalloc(data->loci, sizeof(long));
  for (locus=0;locus < data->loci; locus++)
    {
      options->segregs[locus] = (long) (options->mu_rates[locus] *  mean + 0.5);
    }

}


void free_seq(seqmodel_fmt **seq, long seqnum)
{
  long i;
  for(i=0;i<seqnum;i++)
    {
      if(seq[i]->links!=NULL)
	{
	  myfree(seq[i]->links);
	}
      myfree(seq[i]->sites);
    }
  myfree(seq[0]);
  myfree(seq);
}
