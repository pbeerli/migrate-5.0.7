/*------------------------------------------------------
 Bayesian Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 D A T A   R O U T I N E S 
    
 creates data structures,
 read data (Electrophoretic loci, sequences, microsats),
 prints data,
 destroys data.

Copyright 1997-2002 Peter Beerli and Joseph Felsenstein 
Copyright 2003-2008 Peter Beerli
Copyright 2009-2012 Peter Beerli and Michal Palczewski

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
 
$Id: data.c 2169 2013-08-24 19:02:04Z beerli $

-----------------------------------
add this to the microsatellite data, this will allow to use
a known repeat lengths to read number of sites instead of 
number of repeats, the datafile will need an additional line
perhaps use something like #@ 2 2 2 3  .... on the second line
so that other non microsatellite data models 
can ignore it
y = Map[If[(a = Mod[#, rlen]) != 0, 
      If[Random[] < 0.5, # - (rlen + a), # + rlen - a], #] &, x]/
   rlen  - Min[y] + 1 
-------------------------------------------------------*/
/*! \file data.c

Data manipulation routines

*/


#include <string.h>
#include <wchar.h>
#ifndef WIN32
#include <locale.h>
#endif
#include "migration.h"
#include "sighandler.h"
#include "tools.h"
#include "tree.h"
#include "migrate_mpi.h"
#include "hash.h"
#include "data.h"
#include "sequence.h"
#include "random.h"
#include "haplotype.h"
#include "mutationmodel.h"
#include "pretty.h"

extern long number_genomes (int type);
extern void jumble (long *s, long n);
extern void jumble_ownseed (long *s, long n);

#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif
/* prototypes ----------------------------------------- */
//void create_data (data_fmt ** data);
void shuffle(long **shuffled, long n, long shuffle_ON);
void shuffle_data(data_fmt *data, option_fmt *options);
void decide_full_variable_dataset(data_fmt * data, option_fmt *options, long locus);
//void get_data (FILE * infile, data_fmt * data, option_fmt * options);
//void print_data (world_fmt * world, option_fmt * options, data_fmt * data);
//long find_missing(data_fmt *data, long pop, long locus);
//void print_data_summary (FILE * file, world_fmt * world, option_fmt * options,
//                         data_fmt * data);
//short findAllele (data_fmt * data, char s[], long locus);
//void free_datapart (data_fmt * data, option_fmt * options, long locus);
/*private functions */
//void init_data_structure1 (data_fmt ** data, option_fmt * options);
void read_header (FILE * infile, data_fmt * data, option_fmt * options);
void read_sites (data_fmt * data, world_fmt *world, option_fmt *options);
//void init_data_structure2 (data_fmt ** data, option_fmt * options, long pop);
//void init_data_structure3 (data_fmt * data, world_fmt *world);
CLANG_ANALYZER_NORETURN
void print_utf_warning(void);
void check_ascii(FILE *infile);
void read_popheader (FILE * infile, data_fmt * data, world_fmt *world, long pop, long genomes);
void read_indname (FILE * file, data_fmt * data, long pop, long ind,
                   long locus, long nmlength);
void read_popdata (FILE * file, data_fmt * data, long pop,
                   option_fmt * options, world_fmt *world);
void read_microalleles (FILE * infile, data_fmt * data, option_fmt *options, long pop, long ind);
void read_alleles (FILE * infile, data_fmt * data, long pop, long ind);
long read_ind_seq (FILE * infile, data_fmt * data, option_fmt * options,
                   long locus, long pop, long ind, long baseread);
long read_ind_seq_oneliner (FILE * infile, data_fmt * data, option_fmt * options, world_fmt *world,
			    long pop, long ind, long baseread, long startlocus, long endlocus);
void read_data_line (FILE *infile, data_fmt *data, option_fmt * options, world_fmt * world, long pop, long ind, long region);
long check_list(long la1, long nla1, long locus, data_fmt *data);
void len2repeat(char *a1, long rlen);
void read_hapmap (FILE * infile, data_fmt * data, option_fmt * options, long pop);
void read_hapmap_genotypes (FILE * infile, data_fmt * data, option_fmt * options, long pop);
void read_distance_fromfile (FILE * dfile, long tips, long nmlength, MYREAL **m);
void finish_read_seq (FILE * infile, data_fmt * data, option_fmt * options,
                      long pop, long baseread);
void print_alleledata (world_fmt * world, data_fmt * data,
                       option_fmt * options);
void print_microdata (world_fmt * world, data_fmt * data,
                       option_fmt * options);
void print_seqdata (world_fmt * world, option_fmt * options, data_fmt * data);

void print_header (FILE * outfile, long pop, world_fmt * world,
                   option_fmt * options, data_fmt * data);
MYREAL findleastsquare(MYREAL *rawdata, long total, long repeatlength, long shift, MYREAL *startvalue);
void find_allele_repeatlength(data_fmt *data, option_fmt *options, long locus);
MYREAL create_alleles (world_fmt *world,data_fmt * data, option_fmt *options);
void addAllele (data_fmt * data, char s[], long locus, long *z);
void set_numind (data_fmt * data);
void print_seq_pop (long locus, long pop, world_fmt * world,
                    option_fmt * options, data_fmt * data);
void print_seq_ind (long locus, long pop, long ind, world_fmt * world,
                    option_fmt * options, data_fmt * data);
void print_locus_head (long locus, world_fmt * world, option_fmt * options,
                       data_fmt * data);
void find_delimiter (char *title, char *dlm);
void set_datatype(char input, long locus, char ** locitypes, long *alloc);
void setup_all_nonDNA_subloci(data_fmt * data, world_fmt *world);
void setup_all_oldDNA_subloci(data_fmt * data, world_fmt *world);
void handle_bracket(char *c, char *input, char **newinput, long *allocbufsize, long *z, long *z2, long *pos);
void read_geofile (data_fmt * data, option_fmt * options, long numpop);
MYREAL read_date_fromfile (FILE * datefile, data_fmt *data, option_fmt *options, long nmlength);
void read_datefile (data_fmt * data, option_fmt * options, long numpop);
void read_uep_fromfile (FILE * uepfile, long tips, long nmlength, int **uep,
                        long *uepsites, long datatype);
void read_uepfile (data_fmt * data, option_fmt * options, long numpop);

long read_sublocus(char *input, long *z, long locus, data_fmt * data, world_fmt *world);
void read_sites_new (data_fmt * data, world_fmt *world, option_fmt *options);
void read_site(FILE *infile, char *site, long locus, long sublocus, long pop, long ind, world_fmt *world, data_fmt *data, option_fmt *options);
void read_comments(option_fmt *options, data_fmt *data, world_fmt *world);
MYREAL create_mixed_data (data_fmt * data, option_fmt *options);
void print_random_subset(FILE * file, data_fmt * data, option_fmt *options);

/*=====================================================*/
/// creates the data structure
void
create_data (data_fmt ** data)
{
    (*data) = (data_fmt *) mycalloc (1, sizeof (data_fmt));
}


///
/// free the data module 
void
destroy_data (data_fmt * data)
{
  long ind;
  long indalloc;
  long locus;
  long pop;
  long loci = data->loci;
  //long allsubloci = data->allsubloci;
  long numpop = data->numpop;

  // free data from init_data_structure3
  for (locus = 0; locus < loci; locus++)
    {
      myfree(data->allele[locus]);
    }
  myfree(data->allele);
  myfree(data->subloci);
  myfree(data->maxalleles);
  myfree(data->skiploci);
  myfree(data->locusweight);

  // free data from init_data_structure2
  for(pop=0; pop < numpop; pop++)
    {
      indalloc = -1;
      for(locus=0; locus < loci; locus++)
	{
	  if(indalloc < data->numind[pop][locus])
	    indalloc = data->numind[pop][locus];
	}
      for (ind = 0; ind < indalloc; ind++)
	{
	  for(locus=0; locus < loci; locus++)
	    {
	      myfree(data->indnames[pop][ind][locus]);
	    }
	  myfree(data->indnames[pop][ind]);
	  myfree(data->yy[pop][ind]);
	}
      myfree(data->indnames[pop]);
      myfree(data->yy[pop]);
    }
  myfree(data->indnames);
  myfree(data->yy);

  // data->yy were already freed in free_datapart()

  // free data from init_structure_1
  myfree(data->locitypes);
  myfree(data->locusname);
  myfree(data->popnames[0]);
  myfree(data->numind[0]);
  myfree(data->numalleles[0]);
  myfree(data->popnames);
  myfree(data->numind);
  myfree(data->numalleles);
  myfree(data->seq[0]->sites);
  myfree(data->seq[0]->links);
  if(data->position!=NULL)
    myfree(data->position);
  myfree(data->datatype);
  myfree(data->geo);
  myfree(data->lgeo);
  if(data->ogeo != NULL)
    {
      myfree(data->ogeo[0]);
      myfree(data->ogeo);
    }
  myfree(data->totalsites);
  myfree(data->numrepeatnumbers);
  myfree(data->repeatlength);
  // free sampledates
  for(locus=0;locus<data->loci;locus++)
    {
      myfree(data->repeatnumbers[locus]);
      for(pop=0;pop < data->numpop; pop++)
	{
	  myfree(data->sampledates[pop][locus]);
	}
    }
  myfree(data->repeatnumbers);
  myfree(data->sampledates[0]);
  myfree(data->sampledates);
  //myfree(data);
}
 
///
/// shuffles (jumbles) the individuals so that we can take the first y to subsample the dataset
void shuffle(long **shuffled, long n, long shuffle_ON)
{
  long i;
  for(i=0;i<n;i++)
    {
      (*shuffled)[i] = i;
    }
  switch(shuffle_ON)
    {
    case 0:
      break;
    case -1:
      jumble(*shuffled, n);
      break;
    default:
      jumble_ownseed(*shuffled, n);
    }
}

void shuffle_data(data_fmt *data, option_fmt *options)
{
  long pop;
  long locus;
  data->shuffled = (long ***) mycalloc(data->numpop,sizeof(long **));
  for(pop=0; pop < data->numpop; pop++)
    {
      data->shuffled[pop] = (long **) mycalloc(data->loci, sizeof(long *));
      for(locus=0; locus < data->loci; locus++)
	{
	  data->shuffled[pop][locus] = (long *) mycalloc(data->numind[pop][locus], sizeof(long));
	  if(options->randomsubset>0 && options->randomsubset <  data->numind[pop][locus])
	    {	
	      shuffle(&data->shuffled[pop][locus],data->numind[pop][locus], (long) options->randomsubsetseed);
	    }
	  else
	    {
	      shuffle(&data->shuffled[pop][locus],data->numind[pop][locus],0);
	    }
	}
    }
}


long max_shuffled_individuals(option_fmt *options, data_fmt *data, long pop, long locus)
{
  long top;
  if(options->randomsubset>0 && options->randomsubset <  data->numind[pop][locus])
    {	
      top = options->randomsubset;
    }
  else
    {
      top = data->numind[pop][locus];
    }
  return top;
}


void decide_full_variable_dataset(data_fmt * data, option_fmt *options, long locus)
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


//---------------------------------------------------------------
// ..number_pop.number_loci.<delimiter>.<title [no numbers to start]
// locu_specification {sNumber_of_sites}{n{number of sites}{m1}{a1}
// [parentheses () mark link loci
// example:
// (n3 s120 m1)(m1 m1 m1 a1 s100)(n1)(n2)(m1)
// contains 5 linked loci, the first two are compounds
// the system is to have for each individual haplotype a set of 
// loci so that it may look like this
// ind11......ACA ATTAGA 13 15 5 3 2 A ATACG A AT 13
// ind12......ACA ATTCGA 12 15 6 3 2 A ATACG T CT 5
void
get_new_data (FILE * infile, data_fmt * data, option_fmt * options, world_fmt *world)
{
  //long locus;
  MYREAL mean=1.0;
  long pop;
  long locus;
  long genomes=1;
  boolean saved_option_murate=FALSE;
  data->hasghost = FALSE;
  // read how many populations and loci and delimiter and title
  // if we have the new haplotype input then there should be no delimiter
  // needed
  read_header (infile, data, options);
  data->haplotyping = options->haplotyping;
  data->haplotyping_report = options->haplotyping_report;
  init_data_structure1 (&data);
  init_mutationmodel_first(world, data);
  read_sites_new(data, world, options);
  if(data->oldsyntax)
    {   
      genomes = number_genomes(options->datatype);
    }
  init_mutationmodel_second(world,data);
  if (options->progress)
    fprintf (stdout, "\n\n");
  if (options->writelog)
    fprintf (options->logfile, "\n\n");

  read_comments(options, data, world);// read comments or hapltypes or mutation model

  for (pop = 0; pop < data->numpop; pop++)
    {
      read_popheader (infile, data, world, pop, genomes);
      if (options->progress)
	fprintf (stdout, "Reading (%li) %s ...\n", pop+1, data->popnames[pop]);
      if (options->writelog)
	fprintf (options->logfile, "Reading (%li) %s ...\n", pop+1, data->popnames[pop]);
      init_data_structure2 (&data, options, world, pop);
      read_popdata (infile, data, pop, options, world);
    }
  read_geofile (data, options, data->numpop);
#ifdef UEP

    read_uepfile (data, options, data->numpop);
#endif
    read_datefile(data, options, data->numpop);
    if (options->progress)
        fprintf (stdout, "\n\n");
    init_data_structure3 (data, options, world);
    if ((options->onlyvariable || options->has_variableandone) && !options->murates_fromdata)
      {
	saved_option_murate = TRUE;
	options->murates_fromdata = TRUE;
      }
    switch (options->datatype)
    {
    case 'a':
        mean = create_alleles (world, data, options);
        break;
    case 'b':
        mean = create_alleles (world, data, options);
        for (pop = 0; pop < data->loci; pop++)
            data->maxalleles[pop] = XBROWN_SIZE;
        break;
    case 'm':
        mean = create_alleles (world,data, options);
        for (pop = 0; pop < data->loci; pop++)
            data->maxalleles[pop] = options->micro_stepnum;
        break;
    case '@':
      // mixed input
      //printf ("mixed input");
      //mean = create_mixed_data(data, options);
      break;
    case 'h':
    case 'n':
    case 'u':
      mean = create_alleles (world, data, options);
      break;
    default: /*DNA types*/
      for (pop = 0; pop < data->loci; pop++)
            data->maxalleles[pop] = 4;
    }

    shuffle_data(data, options);
    if (options->murates_fromdata)
      {
	if (strchr (SEQUENCETYPES, options->datatype) || options->datatype=='@')
	  {
	    find_rates_fromdata(data, options, world);
	  }
	else
	  {
	    find_rates_fromdata_alleles(data, options, world, mean);
	  }
      }
    if(saved_option_murate)
      {
	for (locus=0; locus < data->loci; locus++)
	  {
	    options->mu_rates[locus] = 1.0;
	  }
	options->murates_fromdata=FALSE;
	//saved_option_murate = FALSE;
      }
}

void
get_data (FILE * infile, data_fmt * data, option_fmt * options, world_fmt *world)
{
  long locus;
  long pop;
  long genomes=1;
  MYREAL mean = -1;
  boolean saved_option_murate=FALSE;
    data->hasghost = FALSE;
    read_header (infile, data, options);
    genomes =  number_genomes (options->datatype);
    init_data_structure1 (&data);
    data->datatype[0] = options->datatype;
    switch (options->datatype)
    {
    case 's':
    case 'f':// read standard sequence data
      read_sites (data,world, options);
        break;
    case 'n': // read snp data
      read_sites (data,world, options);
        data->seq[0]->addon = 4;
        break;
    case 'h': //read data from hapmap allele frequency files
      // there is no sites line in these files!
      for(locus=0;locus<data->loci;locus++)
	data->seq[0]->sites[locus] = 1;
      data->seq[0]->addon = 4;
        break;
    case 'u': // read unlinked snp
      read_sites (data,world, options);
        data->seq[0]->addon = 4;
        break;
    default:
      read_sites(data,world, options); // checks whether there is a line for microsat repeat numbers
      data->seq[0]->fracchange = 1.0;
      break;
    }
    if (options->progress)
        fprintf (stdout, "\n\n");
    if (options->writelog)
        fprintf (options->logfile, "\n\n");
    for (pop = 0; pop < data->numpop; pop++)
    {
      read_popheader (infile, data, world, pop, genomes);
        if (options->progress)
            fprintf (stdout, "Reading (%li) %s ...\n", pop+1, data->popnames[pop]);
        if (options->writelog)
            fprintf (options->logfile, "Reading (%li) %s ...\n", pop+1, data->popnames[pop]);
        init_data_structure2 (&data, options, world, pop);
        read_popdata (infile, data, pop, options,world);
    }
    read_geofile (data, options, data->numpop);
#ifdef UEP

    read_uepfile (data, options, data->numpop);
#endif
    read_datefile(data, options, data->numpop);
    if (options->progress)
        fprintf (stdout, "\n\n");
    init_data_structure3 (data, options, world);
    if ((options->onlyvariable || options->has_variableandone) && !options->murates_fromdata)
      {
	saved_option_murate = TRUE;
	options->murates_fromdata = TRUE;
      }
    switch (options->datatype)
    {
    case 'a':
        mean=create_alleles (world, data, options);
        break;
    case 'b':
        mean=create_alleles (world, data, options);
        for (pop = 0; pop < data->loci; pop++)
            data->maxalleles[pop] = XBROWN_SIZE;
        break;
    case 'm':
        mean=create_alleles (world, data, options);
        for (pop = 0; pop < data->loci; pop++)
            data->maxalleles[pop] = options->micro_stepnum;
        break;
    case 'h':
    case 'n':
    case 'u':
      mean=create_alleles (world, data, options);
      break;
    default: /*DNA types*/
      for (pop = 0; pop < data->loci; pop++)
            data->maxalleles[pop] = 4;
    }

    shuffle_data(data, options);

    if(options->murates_fromdata)
      {
	if (strchr (SEQUENCETYPES, options->datatype))
	  {
	    find_rates_fromdata(data, options, world);
	  }
	else
	  {
	    find_rates_fromdata_alleles(data, options, world, mean);
	  }
      }
    if(saved_option_murate)
      {
	for (locus=0; locus < data->loci; locus++)
	  {
	    options->mu_rates[locus] = 1.0;
	  }
	options->murates_fromdata=FALSE;
	//saved_option_murate = FALSE;
      }
}

/* private functions ========================================== */

void
init_data_structure1 (data_fmt ** data)
{
    long pop;
    long locus;
    long numpop = (*data)->numpop;
    long loci   = (*data)->allsubloci;
    
    (*data)->ogeo = NULL;
    (*data)->geo = NULL;

    // subloci initialization first time!
    if ((*data)->yy == NULL)
    {
      (*data)->yy = (site_fmt *****) mymalloc (sizeof (site_fmt ****) * (size_t) numpop);
        (*data)->seq = (seqmodel_fmt **) mycalloc (1, sizeof (seqmodel_fmt *));
        (*data)->seq[0] = (seqmodel_fmt *) mycalloc (1, sizeof (seqmodel_fmt));
        (*data)->popnames =(char **) mymalloc (sizeof (char *) * (size_t) numpop);
        (*data)->popnames[0] =(char *) mycalloc (numpop * STRSIZE,sizeof(char));
        (*data)->indnames = (char ****) mymalloc (sizeof (char ***) * (size_t) numpop);
        (*data)->numind = (long **) mymalloc (sizeof (long *) * (size_t) numpop);
        (*data)->numind[0] = (long *) mymalloc (sizeof (long) * (size_t) (numpop * loci));
        (*data)->numalleles = (long **) mymalloc (sizeof (long *) * (size_t) numpop);
        (*data)->numalleles[0] = (long *) mymalloc (sizeof (long) * (size_t) (numpop * loci));
	(*data)->locitypes = (char *) mycalloc(loci, sizeof(char));
	(*data)->locusname = (char **) mycalloc(loci, sizeof(char*));
        for (pop = 1; pop < numpop; pop++)
        {
	  (*data)->popnames[pop] = (*data)->popnames[0] + pop * STRSIZE;
	  (*data)->numind[pop] = (*data)->numind[0] + pop * loci;
	  (*data)->numalleles[pop] =  (*data)->numalleles[0] + pop * loci;
        }
        (*data)->seq[0]->sites = (long *) mycalloc (loci, sizeof (long));
        (*data)->seq[0]->links = (boolean *) mycalloc (loci, sizeof (boolean));
        (*data)->position = (long *) mycalloc (loci, sizeof (long));
	(*data)->datatype = (char *) mycalloc(loci, sizeof(char));
	(*data)->numdatatypealloc = loci;
	(*data)->numsublocialloc = loci;
	(*data)->subloci = (long *) mycalloc(loci,sizeof(long));
	//	for(pop=0;pop<loci;pop++)
	//  (*data)->subloci[pop] = 1;
	(*data)->totalsites = (long *) mycalloc((1+loci),sizeof(long));
	(*data)->numrepeatnumbers = (long *) mycalloc(loci, sizeof(long));
	(*data)->repeatnumbers = (long **) mycalloc(loci, sizeof(long *));
	(*data)->repeatlength = (long *) mycalloc(loci, sizeof(long));
	for(locus=0;locus<loci;locus++)
	  {
	    (*data)->numrepeatnumbers[locus] = 1;
	    (*data)->repeatnumbers[locus] = (long *) mycalloc(1, sizeof(long));
	  }
    }
    else
    {
        error ("Problem with initialization of data matrix yy\n");
    }
    (*data)->numregions = 1;
    (*data)->regions = (region_fmt *) mycalloc(loci, sizeof(region_fmt));
    if((*data)->haplotyping)
      {
	init_individuals(data,loci);
      }

}


void
init_data_structure2 (data_fmt ** data, option_fmt * options, world_fmt *world, long pop)
{
  //long sites;
  //long startsite;
  //long endsite;
  long ind, locus;
  long indalloc = -1;
  for(locus=0;locus<(*data)->loci;locus++)
    {
      if(indalloc < (*data)->numind[pop][locus])
	indalloc = (*data)->numind[pop][locus];
    }
  if (indalloc == 0)
    indalloc = 2;
#ifdef MPIDATAONDEMAND
  if(myID==MASTER)
    {
      (*data)->yy[pop] = (site_fmt ****) mymalloc (sizeof (site_fmt ***) * (size_t) indalloc);
    }
#else
  (*data)->yy[pop] = (site_fmt ****) mymalloc (sizeof (site_fmt ***) * (size_t) indalloc);
#endif
  (*data)->indnames[pop] = (char ***) mycalloc (1, sizeof (char **) * (size_t) indalloc);
  
  for (ind = 0; ind < indalloc; ind++)
    {
      (*data)->indnames[pop][ind] =
	(char **) mymalloc (sizeof (char *) * (size_t) (*data)->loci);
#ifdef MPIDATAONDEMAND    
  if(myID==MASTER)
    {
      (*data)->yy[pop][ind] =
	(site_fmt ***) mymalloc (sizeof (site_fmt **) * (size_t) (*data)->allsubloci);
    }
#else
      (*data)->yy[pop][ind] =
	(site_fmt ***) mymalloc (sizeof (site_fmt **) * (size_t) (*data)->allsubloci);
#endif
      for (locus = 0; locus < (*data)->loci; locus++)
        {
	  (*data)->indnames[pop][ind][locus] =
            (char *) mycalloc (1, sizeof (char) * (size_t) (1 + options->nmlength));
	}
#ifdef MPIDATAONDEMAND
      if (myID==MASTER)
	{
#endif
      long sublocus;
      long len;
      long site;
      long sites;

      for(sublocus=0; sublocus < (*data)->allsubloci;sublocus++)
	{
	  mutationmodel_fmt *s = & world->mutationmodels[sublocus];
	  sites = s->numsites;
	  (*data)->yy[pop][ind][sublocus] =
	    (site_fmt **) mycalloc (2, sizeof (site_fmt *));
	  
	  (*data)->yy[pop][ind][sublocus][0] =
	    (site_fmt *) mycalloc ((1+sites), sizeof (site_fmt));
	  if(!strchr(SEQUENCETYPES,options->datatype) && options->datatype!='@')
	    (*data)->yy[pop][ind][sublocus][1] =
	      (site_fmt *) mycalloc ((1+sites), sizeof (site_fmt));
	  
	  if(s->dataclass==SITEWORD)
	    {
	      len = (options->allelenmlength+1);
	    }
	  else
	    {
	      len = 2;
	    }
	  
	  for(site=0;site<sites;site++)
	    {
	      (*data)->yy[pop][ind][sublocus][0][site] = (site_fmt) mycalloc (len, sizeof (char));
	      if(!strchr(SEQUENCETYPES,options->datatype) && options->datatype!='@')
		(*data)->yy[pop][ind][sublocus][1][site] = (site_fmt) mycalloc (len, sizeof (char));
	    }
	}
#ifdef MPIDATAONDEMAND
	}
#endif /*MPIDATAONDEMAND*/
    }
}


short
findAllele (data_fmt * data, char s[], long locus)
{
    short found = 0;
    while ((strcmp (s, data->allele[locus][found])
            && data->allele[locus][found][0] != '\0'))
        found++;
    return found;
}

void
free_datapart (data_fmt * data, world_fmt * world , long locus)
{
  long ind, pop;
  //  long genomes = number_genomes (options->datatype);
  for (pop = 0; pop < data->numpop; pop++)
    {
      for (ind = 0; ind < data->numind[pop][locus]; ind++)
        {
	  const long sublocistart = world->sublocistarts[locus];
	  const long sublociend   = world->sublocistarts[locus+1];
	  long sublocus;
	  boolean xxxx = (boolean) (!strchr(SEQUENCETYPES,
					  world->options->datatype) && 
				world->options->datatype!='@');
	  for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
	    {
	      mutationmodel_fmt *s = &world->mutationmodels[sublocus];
	      long sites = s->numsites;
	      long site;
	      for(site=0;site<sites;site++)
		{
		  myfree(data->yy[pop][ind][sublocus][0][site]);
		  if(xxxx)
		    {
		      myfree(data->yy[pop][ind][sublocus][1][site]);
		    }
		}
	      myfree(data->yy[pop][ind][sublocus][0]);
	      if(xxxx)
		myfree(data->yy[pop][ind][sublocus][1]);
	    }
	}
    }
}

void
init_data_structure3 (data_fmt * data, option_fmt *options, world_fmt *world)
{
    long locus, pop, maxi;
    if (data->allele == NULL)
      data->allele = (char ***) mycalloc (1, sizeof (char **) * (size_t) data->allsubloci);
    if (data->subloci == NULL)
      data->subloci = (long *) mycalloc(data->loci,sizeof(long));
    for (locus = 0; locus < data->loci; locus++)
    {
      if(!data->oneliner)
	data->subloci[locus] = 1;
      maxi = 0;
      for (pop = 0; pop < data->numpop; pop++)
	maxi += data->numalleles[pop][locus];
      long sublocus;
      const long sublocistart = world->sublocistarts[locus];
      const long sublociend   = world->sublocistarts[locus+1];
      for(sublocus = sublocistart; sublocus < sublociend; sublocus++)
	{
	        data->allele[sublocus] =
        	    (char **) mycalloc (maxi, sizeof (char*));
        	for (pop = 0; pop < maxi; pop++)
	  	    data->allele[sublocus][pop] = (char *) mycalloc (options->allelenmlength+1, sizeof (char));
        }
    }
    data->maxalleles = (long *) mycalloc (data->allsubloci, sizeof (long));
    data->skiploci =
      (boolean *) mycalloc ((data->allsubloci + 1), sizeof (boolean));
    data->locusweight =
        (MYREAL *) mycalloc ((data->loci + 1), sizeof (MYREAL));
    for (locus=0;locus < data->loci; locus++)
      {
	data->locusweight[locus]=1.0;
      }
}

void print_utf_warning(void) 
{
  warning("Tried and failed to read a file that is not in ASCII TEXT format\n");
  warning("This is most commonly UTF-8 or UTF-16 formatted\n");
  warning("Try to convert the file into plain ASCII, enclodings like\n");
  warning("WESTERN (MAC ROMAN) or WESTERN (WINDOWS LATIN 1) or WESTERN ISO LATIN\n");
  usererror("Program aborts now\n");
}

void check_ascii(FILE *infile)
{
  wint_t wch, wch2;
  wch = fgetwc(infile);
  if (wch > 127)
    {
      print_utf_warning();
    }
  else
    {
      wch2 = fgetwc(infile);
      if(wch==0 && wch2>0)
	{
	  print_utf_warning();
	}
      else
	{
	  if (wch>0 && wch2==0)
	    {
	      print_utf_warning();
	    }
	}
      rewind(infile);
    }
}
///
/// read the first line of the data file
/// \param infile datafilename
/// \param data   data structure that holds all the data
/// \param options structure that contain all option information
/// retval none
void read_header (FILE * infile, data_fmt * data, option_fmt * options)
{
    char input[LINESIZE], *p;
    char title[LINESIZE];
    strcpy(title,"\0");
    //disabled for now: check_ascii(infile);
    input[0] = '#';
    while(input[0]=='#')
      {
	FGETS (input, sizeof (input), infile);
#ifdef DEBUG
	printf("@@>%s<@@\n",input);
#endif      
	if ((p = (char *) strpbrk (input, CRLF)) != NULL)
	  *p = '\0';
      }
    switch (lowercase (input[0]))
    {
    case 'a':
        sscanf (input, "%1s%ld%ld%[^\n]", &options->datatype, &(data->numpop),
                &(data->allsubloci), title);
        find_delimiter (title, &data->dlm);
	if(!(title[0] == '\0'))
	  strcpy(options->title,title);
        break;
    case 'b':
    case 'm':
        sscanf (input, "%1s%ld%ld%1s%[^\n]", &options->datatype,
                &(data->numpop), &(data->allsubloci), &data->dlm, title);
	if(!(title[0] == '\0'))
	  strcpy(options->title,title);
        break;
    case 's':
    case 'n':
    case 'h': //hapmap data
    case 'u':
    case 'f':
        sscanf (input, "%1s%ld%ld%[^\n]", &options->datatype, &(data->numpop),
                &(data->allsubloci), title);
	if(!(title[0] == '\0'))
	  strcpy(options->title,title);
        break;
    case 'g':   /* fall through if a menu change forces to analyze data
                               instead of using the already sampled genealogies */
        if (options->datatype == 'g')
            break;
        else
            memmove (input, input + 1, (strlen (input) - 1) * sizeof (char));
    default:
      if(input[0]== '<')
	{
	  usererror ("This data file may contain an XML or HTML tag,\nand cannot be read properly, check the data formatting section in the manual!");
	  //  exit(-1);
	}
        switch (options->datatype)
        {
        case 'a':
            sscanf (input, "%ld%ld%[^\n]", &(data->numpop), &(data->allsubloci),
                    title);
            find_delimiter (title, &data->dlm);
	    
	if(!(title[0] == '\0'))
	  strcpy(options->title,title);
            break;
        case 'b':
        case 'm':
            sscanf (input, "%ld%ld%1s%[^\n]", &(data->numpop), &(data->allsubloci),
                    &(data->dlm), title);
	    if(!(title[0] == '\0'))
	      strcpy(options->title,title);
            break;
        case 's':
        case 'n':
	case 'h': // hapmap data
        case 'u':
        case 'f':
            sscanf (input, "%ld%ld%[^\n]", &(data->numpop), &(data->allsubloci),
                    title);
	    if(!(title[0] == '\0'))
	      strcpy(options->title,title);
            break;
        default:
            usererror ("Datatype is wrong, please use a valid data type!");
        }
    }
    options->datatype = lowercase (options->datatype);
}

void
find_delimiter (char *title, char *dlm)
{
    char *p = title;
    size_t z = 0;
    while (*p == ' ')
    {
        p++;
        z++;
    }
    if (isalnum (*p))
      memmove (title, p, sizeof (char) * (strlen (title) - z));
    else
    {
        *dlm = *p;
        p++;
        while (*p == ' ')
        {
            p++;
            z++;
        }
        memmove (title, p, sizeof (char) * (strlen (title) - z));
    }
}

void set_datatype(char input, long locus, char ** locitypes, long *alloc)
{
  if(locus >= *alloc)
    {
      *alloc += SUBLOCICHUNKS;
      *locitypes = (char *) myrealloc(*locitypes, sizeof(char) * (size_t) (*alloc));
    }
  (*locitypes)[locus] = input;
}


///
/// setup all loci that are not sequences within the old framework
void setup_all_nonDNA_subloci(data_fmt * data, world_fmt *world)
{
  long locus;
  char val[2];
  val[0]='1';
  val[1]='\0';
  for(locus=0; locus < data->allsubloci; locus++)
    {
      init_mutationmodel_readsites(&world->mutationmodels[locus], data->datatype[0], val);
      data->subloci[locus]=1;
      world->sublocistarts[locus] = locus;
    }
  world->sublocistarts[locus] = locus;
}
///
/// setup all loci that are sequences within the old framework
void setup_all_oldDNA_subloci(data_fmt * data, world_fmt *world)
{
  long region;
  long locus;
  char *sites;
  char *input = (char *) mycalloc(LONGLINESIZE,sizeof(char));
  for(locus=0; locus < data->allsubloci; locus++)
    {
      read_word(data->infile, input, NULL);
      sites = input;//atol(input);
      init_mutationmodel_readsites(&world->mutationmodels[locus], data->datatype[0], sites);
      data->subloci[locus]=1;
      world->sublocistarts[locus] = locus;
      data->seq[0]->sites[locus] = world->mutationmodels[locus].numsites; 
      data->totalsites[locus] = world->mutationmodels[locus].numsites; 
    }
  world->sublocistarts[locus] = locus;
  data->numregions = data->allsubloci;
  for (region=0; region < data->numregions; region++)
    {
      data->regions[region].startlocus=region;
      data->regions[region].endlocus=region+1;
    }
  myfree(input);
}


void handle_bracket(char *c, char *input, char **newinput, long *allocbufsize, long *z, long *z2, long *pos)
{
  char *tmp;
  char *val;
  char *vallinked;
  char startloc='(';
  char stoploc=')';
  long count=0;
  long linking=1;
  long shortsites=0;
  char startsitetype;
  long repeats;
  long size;
  long startsite=0;
  long rest;
  long numsites=0;
  long newnumsites;
  long lastnumsites=0;
  //long laststartsite=0;
  long *random_startsites=NULL;
  long i;
  char datatype=' ';

  tmp = (char*) mycalloc(STRSIZE,  sizeof(char));
  while (*c!=']')
    {
      *c = input[*z];
      (*z)++;
      tmp[(*z2)++] = *c;
    }
  // tmp can be now:
  // 10 or 10o3 or 10r3 followed by total length of sequence
  // so that (1) sites/10 and last element is larger or smaller
  // or the subset is equally spaced
  // or the subset is randomly spaced
  // case 1: [10]
  // case 2: [10r3]
  // case 3: [10o3]
  // case 4 takes case 3 but rearranges the () to
  //        groups of n [example is 2 and 4 below]
  // case 4: [10o3_2]
  // case 4: [10o3_4]
  val = char_position("ro", tmp);// cases 1-3
  if (val==NULL)
    {
      shortsites = -1;
      startsitetype = '0';
      repeats = atol(tmp);
    }
  else
    {
      shortsites = atol(val+1); // if shortsites == 0 then use difference between start sites
      startsitetype = val[0];
      val = NULL;
      repeats = atol(tmp);
    }
  vallinked = char_position("_", tmp); // case 4
  if(vallinked!=NULL)
    {
      linking = atol(vallinked+1);
    }
  printf("repeats = %li\n",repeats);
  size = locate_char(input + *z,'(');
  *z += size;//location of '('
  size = locate_char(input+ *z,')');
  *z2=0;
  memset(tmp, '\0', STRSIZE *  sizeof(char));
  for(i=0;i<=size;i++)
    {
      tmp[*z2] = input[*z];
      (*z)++;
      (*z2)++;
    }
  tmp[*z2]='\0';
  // tmp is now (s1000)
  if (strchr("abmsn", tmp[1]))
    {
      datatype = tmp[1];
      numsites = atol(tmp+2); 
    }
  else
    error("Problem with assigning in handle_bracket()");
  if (startsitetype=='r')
    {
      assign_random_startsites(&random_startsites,numsites, shortsites,repeats);
    }
  newnumsites = numsites / repeats;
  printf("all loci except last = %li sites\n",newnumsites);
  rest = numsites % repeats;
  lastnumsites = newnumsites + rest;
  printf("last locus = %li sites\n",lastnumsites);
  if(startsitetype == '0')
    startsite = -newnumsites;
  else
    startsite = -newnumsites + rest/2; //sets start to rest/2
  //laststartsite = 0;
  for (i=0;i<repeats-1;i++)
    {
      switch(startsitetype)
	{
	case '0':
	  startsite += newnumsites;
	  //lastnumsites = newnumsites;
	  break;
	case 'o':
	  startsite += newnumsites;
	  lastnumsites = shortsites;
	  break;
	case 'r':
	  if (shortsites != 0)
	    {
	      startsite = random_startsites[i];
	      lastnumsites = shortsites;
	    }
	  else
	    {
	      startsite = random_startsites[i];
	      lastnumsites = random_startsites[i+1]-startsite;
	      //laststartsite = random_startsites[i];
	    }
	  break;
	default:
	  error("no startsitetype defined");
	}
      if(*pos+100 >= *allocbufsize) //enough for two large numbers and a char
	{
	  *allocbufsize += STRSIZE;
	  *newinput = (char *) realloc(*newinput, sizeof(char) * (size_t) (*allocbufsize));
	}
      if (startsitetype != '0')
	{
	  // if(startsitetype == 'o')
	  //  *pos += sprintf(*newinput+ *pos,"(%li%c%li) ", startsite, datatype, shortsites);
	  //else
	  //  *pos += sprintf(*newinput+ *pos,"(%li%c%li) ", startsite, datatype, newnumsites);
	  if(linking==1)
	    {
	      startloc='(';
	      stoploc=')';
	    }
	  else
	    {
	      if(count==0)
		{
		  startloc='(';
		  stoploc=' ';
		  count++;
		}
	      else
		{
		  if(count<linking-1)
		    {
		      startloc=' ';
		      stoploc=' ';
		      count++;		  
		    }
		  else
		    {
		      if(count==linking-1)
			{
			  startloc=' ';
			  stoploc=')';
			  count=0;
			}
		    }
		}
	    }
	  *pos += sprintf(*newinput+ *pos,"%c%li%c%li%c ", startloc, startsite, datatype, shortsites, stoploc);
	}
      else
	*pos += sprintf(*newinput+ *pos,"(%c%li) ", datatype, newnumsites);
    }
  switch(startsitetype)
    {
    case '0':
      startsite += newnumsites;
      //shortsites = newnumsites;
      break;
    case 'o':
      startsite += newnumsites;
      lastnumsites = shortsites;
      break;
    case 'r':
      if (shortsites != 0)
	{
	  startsite = random_startsites[i];
	  lastnumsites = shortsites;
	}
      else
	{
	  startsite = random_startsites[i];
	  lastnumsites = numsites - startsite;
	}
      break;
    default:
      error("no startsitetype defined");
    }
  if(*pos+100 >= *allocbufsize)
    {
      *allocbufsize += STRSIZE;
      *newinput = (char *) realloc(*newinput, sizeof(char) * (size_t) (*allocbufsize));
    }
  if (startsitetype != '0')
    {
      if(linking>1)
	{
	  if(count==linking-1)
	    {
	      startloc=' ';
	      stoploc=')';
	    }
	  if(count==0)
	    {
	      startloc='(';
	      stoploc=')';	  
	    }
	}
      else
	{
	  startloc='(';
	  stoploc=')';
	}
      *pos += sprintf(*newinput+*pos,"%c%li%c%li%c ", startloc,startsite, datatype, lastnumsites, stoploc);	  
    }
  else
    *pos += sprintf(*newinput+*pos,"(%c%li) ", datatype, lastnumsites);
  *c = input[*z++];
  myfree(tmp);
}


//==================================================================
/// read the number of sites for each locus in the dataset,
/// does not assume a fixed line length, but assumes that at the end of the line is either
//  a '\n' or '\r' or '\r\l' (similar to the sequence reader) to accommodate windows, mac and
/// unix line ends.
/// The generalized sites reader can also read link status using a parenthesis notation, 
/// and allows to give short names to the
/// loci, example: (ldh=s234 agdh=n3) gpi=s100  str1=b1 (str2=m1 a1)
/// the labels a, s, n etc are the same used for the datatype and will mark each locus
/// for the datatype=h it is assumed that all loci are the same and of type n
/// plain numbers are considered to be sequence loci with no name
//  examples:
//  (s10 s100) (s23)
//  (s10 s100), (s23)
//  (s10) (s100), (s23)
//  (5s10 345s100), (0s23) // start reading at position 4 for 10 sites then at position 345 for 100 sites, 
//                         // then go the a region and start at 0
//  or even
//  (ldh12=5s10 adh1=345s100), (idh=0s23)
//  10 100 23
//  implemented: [10] (s10000) leads to (s1000)(s1000)(s1000)(s1000)(s1000)(s1000)(s1000)(s1000)(s1000)(s1000)
//  it turns out that reading and working with contiguous >20 Mb is difficult, it needs tons of RAM
//  get n random or regular subsets of length x distributed over the genome could be done, syntax:
//  [10o3](s1000)--> (0s3) (99s3) (198s3) (297s3) ...
//  [10r3] (s1000) ---> (5s3) (10s3) (169s3) .....
//  [10o3_2](s1000)--> (0s3 99s3) (198s3 297s3) ...
void read_sites_new(data_fmt * data, world_fmt *world, option_fmt *options)
{
  long z=0;
  boolean oldsyntax;
  long l=0;
  long region=0;
  long size = 0;
  long loci = 0;
  long locus = 0;
  long sublocus;
  long subloci=0;
  char *regptr;
  // readline with whole sites string
  long allocbufsize = LINESIZE;
  long alloclongbufsize = LONGLINESIZE;
  char *input = (char*) mycalloc(allocbufsize,sizeof(char));
  //  char *tmp = (char*) mycalloc(allocbufsize,sizeof(char));
  char *tmp;
  char *word = (char*) mycalloc(allocbufsize,sizeof(char));
  char *inptr;
  char *regio = (char*) mycalloc(alloclongbufsize,sizeof(char));
  char *locistr = (char*) mycalloc(alloclongbufsize,sizeof(char));
  char *lociptr;
  long oldi=0;
  long len=0;
  long i=0;
  boolean adjust_msat_repeats=FALSE;
  if(options->datatype=='b' || options->datatype=='m')
    data->has_repeats = TRUE;
  FGETS2(&input,&allocbufsize,data->infile);
  //read_word(data->infile, input,";:\t\n");
  while (input[0]=='#')
    {
      if(input[1]=='@')
	{
	  // line is microsatellite instruction
	  if(input[2]=='M')
	    {
	      //FGETS2(&input,&allocbufsize,data->infile); // read the line
	      locus--;
	      oldi=0;
	      // long i;
	      len=3; // jump over #@M
	      data->has_repeats = FALSE;
	      adjust_msat_repeats=TRUE;
	      while(*(input+len) == ' ')
		len += 1;
	      for(i = 0; i < data->allsubloci; i++)
		{
		  len += read_word_delim(input+len,word,";,\t",FALSE);
		  if(word[0]!='\0')
		    {
		      data->repeatlength[i] = atol(word);
		      oldi=i;
		    }
		  else
		    {
		      data->repeatlength[i] = data->repeatlength[oldi];
		    }
		}
	      input[0]='\0';
	    }
	  else /*line starts with #@ buth no M read therefore this is a comment and discarded*/
	    {
	      FGETS2(&input,&allocbufsize,data->infile); // read the line 
	    }
	}
      else /*this is a plain comment*/
	{
	      FGETS2(&input,&allocbufsize,data->infile); // read the line 
	}
    }
  // count (): loci
  long repeatcount = count_char(input,'[');
  long locicount=0;
  if (repeatcount==0)
    {
      locicount = count_char(input,'(');
    }
  else
    {
      char c = input[0];
      z=1;
      long z2=0;
      long pos = 0;
      char *newinput;
      newinput = (char *) mycalloc(allocbufsize,sizeof(char));
      while (c != '\0')
	{
	  switch (c){
	    case '[':
	      z2=0;
	      handle_bracket(&c,input,&newinput, &allocbufsize, &z,&z2, &pos);
	      break;
	  default:
	    pos += sprintf(newinput+pos,"%c",input[z]);
	    c = input[z++];
	  }
	}
      printf("%s\n",newinput);
      myfree(input);
      input = newinput;
      locicount = count_char(input,'(');
    }
  // count ',': regions
  data->numregions = 1 + count_char(input,',');
  if (locicount == 0)
    {
      oldsyntax=TRUE;
      data->loci = count_words(input);
    }
  else
    {
      oldsyntax=FALSE;
      data->oneliner=TRUE;
      data->loci = locicount;
      options->datatype='@';
    }
  //printf("%li\n",locicount);
  data->oldsyntax = oldsyntax;
  // oldsyntax
  if (data->oldsyntax)
    {
      data->datatype[0] = options->datatype;
      for(sublocus=1;sublocus<data->allsubloci;sublocus++)
	{
	  data->datatype[sublocus] = options->datatype;
	}
      //printf("%i> =#=#=#=#=>[%c] %li %li \n",myID, data->datatype[0],data->loci,data->allsubloci);
      if (options->datatype == 'h' || strchr(ALLELETYPES, options->datatype)) //read data from hapmap allele frequency files
	{
	  data->loci = data->allsubloci;
	  if (!adjust_msat_repeats)
	    {
	      unread_word(data->infile, "\n");
	      unread_word(data->infile,input);
	    }
	  setup_all_nonDNA_subloci(data, world);// this includes hapmap data                                                         
	}
      else
	{
	      unread_word(data->infile,input);
	      setup_all_oldDNA_subloci(data, world);
	    }
	  myfree(input);
	  myfree(word);
	  myfree(regio);
	  myfree(locistr);
	  return;
	}
  //new syntax
  data->regions[0].startlocus = 0;
  data->regions[0].endlocus = 0;
  inptr = input;
  l = 0;
  z=0;
  for(region=0;region < data->numregions; region++)
    {
      size = locate_char(inptr,',');
      memcpy(regio, inptr,sizeof(char)* (size_t) size);
      inptr += size+1;
      loci = count_char(regio,'(');
      regptr = regio;
      for(locus=0; locus < loci; locus++)
	{
	  // assumes regio starts with '('
	  while(*regptr == '(' || *regptr == ' ')
	    regptr++;
	  size = locate_char(regptr,')');
	  memset(locistr,0,sizeof(char) * (size_t) alloclongbufsize);
	  memcpy(locistr, regptr,sizeof(char) * (size_t) size);
	  regptr += size + 1;
	  lociptr = locistr;
	  subloci = count_words(locistr);

	  for(sublocus=0; sublocus<subloci; sublocus++)
	    {	      
	      get_next_word(&lociptr," ",&tmp);
	      read_sublocus(tmp,&z,l,data,world);
	      world->sublocistarts[l+1] = z;	     
	    }
	  l++; //next locus
	}
      data->regions[region].endlocus += loci;
      if(region < data->numregions-1)
	{
	  data->regions[region+1].startlocus = data->regions[region].endlocus; 
	  data->regions[region+1].endlocus = data->regions[region].endlocus; 
	}
    }
  myfree(input);
  //myfree(tmp);
  myfree(word);
  myfree(regio);
  myfree(locistr);
}



long read_sublocus(char *input, long *z, long locus, data_fmt * data, world_fmt *world)
{
  char *sites;
  char *word;
  char *val;
  word = (char *) mycalloc(LONGLINESIZE ,sizeof(char));
  data->subloci[locus] += 1;
  read_word_delim(input,word,"=",FALSE);
  if(input[0]!='\0')
    {
      data->locusname[*z] = (char *) mycalloc(strlen(word)+2,sizeof(char));
      strcpy(data->locusname[*z],word+1);
      val=strchr((char*) "abmsnuh", input[0]);
      if(val!=NULL)
	{
	  sites = input; //atoi(input+1);
	  init_mutationmodel_readsites(&world->mutationmodels[*z],input[0],sites);
	  //	  set_datatype(input[0], *z, &data->locitypes, &data->numlocitypesalloc);
	}
      else
	{
	  // no datatype specified in the model therefore we assume there is a 
	  // unique data type specified in the parmfile or menu
	  // loci with explicit start sites appear also here and get a wrong datatype, but this gets
	  // fixed in the function
	  init_mutationmodel_readsites(&world->mutationmodels[*z], data->datatype[0], input);//atoi(input));
	  //	  set_datatype(data->datatype[0],locus, 
	  //	       &data->locitypes, &data->numlocitypesalloc);
	  //&world->mutationmodels[zz].datatype);
	}
    }
  else
    {
      // no locusname specified, we use the default
      data->locusname[*z] = (char *) mycalloc(1+strlen(input),sizeof(char));
      strcpy(data->locusname[*z],input);
      val=strchr((char *) "abmsnuh", input[0]);
      if(val!=NULL)
	{
	  init_mutationmodel_readsites(&world->mutationmodels[*z],input[0],input+1);//atoi(input+1));
	  //	  set_datatype(input[0], *z, &data->locitypes, &data->numlocitypesalloc);
	}
      else
	{
	  // no datatype specified in the model therefore we assume there is a 
	  // unique data type specified in the parmfile or menu
	  init_mutationmodel_readsites(&world->mutationmodels[*z], data->datatype[0],input);// atoi(input));
	  //      set_datatype(data->datatype[0],locus,
	  //		   &data->locitypes, &data->numlocitypesalloc);
	  //, &world->mutationmodels[zz].datatype);
	}
    }
  data->seq[0]->sites[*z] = world->mutationmodels[*z].numsites; 
  data->totalsites[data->loci] += world->mutationmodels[*z].numsites; 
  *z += 1;
  myfree(word);

  return *z;
}



//==================================================================
/// read the number of sites for each locus in the dataset,
/// does not assume a fixed line length, but assumes that at the end of the line is either
//   a '\n' or '\r' or '\r\l' (similar to the sequence reader) to accommodate windows, mac and
/// unix line ends.
/// The generalized sites reader can also read link status using a parenthesis notation, 
/// and allows to give short names to the
/// loci, example: (ldh=s234 agdh=n3) gpi=s100  str1=b1 (str2=m1 a1)
/// the labels a, s, n etc are the same used for the datatype and will mark each locus
/// for the datatype=h it is assumed that all loci are the same and of type n
/// plain numbers are considered to be sequence loci with no name
void
read_sites (data_fmt * data, world_fmt *world, option_fmt *options)
{
  char * val=NULL;
  boolean linked=TRUE;
  long allocbufsize=LINESIZE;
    long locus;
    char *input;
    char * word;
    long oldi;
    long i;
    long len;
    long z=0; //z is adding subloci and is reset with a closing ')'
    long zz=0; //zz is used for mutationmodels and is not reset at closing ')'
    input = (char *) mycalloc(allocbufsize ,sizeof(char));
    word = (char *) mycalloc(allocbufsize ,sizeof(char));
    // check whether the microsats are not repeats but fragmentlength  
    data->has_repeats = TRUE;
    for (locus = 0; locus < data->loci; locus++)
    {
      read_word(data->infile, input,NULL);
      switch(input[0])
	{
	case  '#':
	  if(input[1]=='@' || input[1]=='M')
	    {
	      if (input[1]=='M')
		{
		  warning("You should use #@M if you want to treat the input as fragment lengths in your infile\n");
		}
	      // line is microsatellite instruction
	      if(input[2]=='M' || input[1]=='M')
		{
		  FGETS2(&input,&allocbufsize,data->infile); // read the line 
		  locus--;
		  oldi=0;
		  // long i;
		  len=0;
		  data->has_repeats = FALSE;
		  for(i = 0; i < data->loci; i++)
		    {
		      len += read_word_delim(input+len,word,";,\t",FALSE);
		      if(word[0]!='\0')
			{
			  data->repeatlength[i] = atol(word);
			  oldi=i;
			}
		      else
			{
			  data->repeatlength[i] = data->repeatlength[oldi];
			}
		    }
		}
	      else
		{
		  FGETS2(&input,&allocbufsize,data->infile); // read the line 
		  locus--;
		}
	      myfree(input);
	      myfree(word);
	      return;
	    }
	  else
	    {
	      // line is a comment
	      FGETS2(&input,&allocbufsize,data->infile); // read the line 
	      locus--;
	      continue;
	    }
	  // break; will not be executed
	  
	  case '(':
	    // (locus1 locus2) (locus3 locus4 locus5) (locus6) on one line assumes that 
	    // there are 3 unlinked loci with 2 3 1 subloci. This would force the input to be strictly 
	    // on one long line alowing for genomic type input with braking up of the data within a run
	    // this is a change from the old format, dna is now treated the same way as other loci
	    // I am not sure whether this really works, but would make things simpler,
	    // still unclear how to deal with haplotype versus diplotype information
	    data->oneliner=TRUE;
	    options->oneliner = TRUE;
	    linked=TRUE;
	    if(z>=data->numsublocialloc)
	      {
		data->numsublocialloc += SUBLOCICHUNKS;
		data->subloci = (long *) myrealloc(data->subloci, data->numsublocialloc * sizeof(long));
	      }
	    data->subloci[z] += 1; //we opened a ( and therefore start a linkage group
	    if(input[1]!='\0' && isalpha(input[1]))
	      {
		read_word_delim(input,word,"=",FALSE);
		if(input[0]!='\0')
		  {
		    data->locusname[z] = (char *) mycalloc(strlen(word),sizeof(char));
		    strcpy(data->locusname[locus],word+1);
		    // what is this good for?
		    if(!strcmp(word,input))
		      val=strchr((char *) "bmsnuh", input[0]);
		    else
		      val=strchr((char *) "bmsnuh", input[1]);

		    if(val==NULL)
		      {
                        set_datatype((input[0]=='(' ? input[1] : input[0]),locus, &data->locitypes, &data->numlocitypesalloc);
			zz++;
		      }
		    else
		      {
			// no datatype specified in the model therewe assume there is a unique data type specified in the parmfile or menu
                        set_datatype(data->datatype[0],locus, &data->locitypes, &data->numlocitypesalloc);
		      }
		  }
		else
		  {
		    // no datatype specified in the model therefore we assume there is a unique data type specified in the parmfile or menu
                    set_datatype(data->datatype[0],locus, &data->locitypes, &data->numlocitypesalloc);
		  }
	      }
	    else
	      {
		if(strlen(input)==1 || input[1]==' ')
		  {
		    locus--;
		    continue;
		  }
	      }
	    break;
	  case ')':
	    z=0; // reset subloci for the next locus
	    break;
	  default:
	    if (!strchr (SEQUENCETYPES, options->datatype))
	      {
                unread_word(data->infile, input);
		myfree(input);
		myfree(word);
		return;
	      }
	    if(isalpha(input[0]))
	      {
		read_word_delim(input,word,"=",FALSE);
		data->locusname[locus] = (char *) mycalloc(strlen(word),sizeof(char));
		strcpy(data->locusname[locus],word);
	      }
	  }
	if(linked)
	  data->seq[0]->links[locus] = TRUE;
	else
	  data->seq[0]->links[locus] = FALSE;
	if(input[0]=='(')
	  {
	    if(strchr("bmsnuh",input[1]))
	      {
		val = input+2;
	      }
	    else
	      {
		val = input+1;
	      }
	  }
	else
	  {
	    if(strchr("bmsnuh",input[0]))
	      {
		val = input+1;
	      }
	    else
	      {
		val = input;
	      }
	  }
	init_mutationmodel_readsites(&world->mutationmodels[locus], data->datatype[0], val);//atoi(val));
	//        data->seq[0]->sites[locus] = atoi (val);
	val = NULL;
	if(strchr("n",options->datatype))
	  {
	    if (options->allelenmlength < data->seq[0]->sites[locus])
	      options->allelenmlength =  data->seq[0]->sites[locus];
	  }
    }
    //    FGETS2(&input,&allocbufsize,data->infile);
    //while(input[0]=='#')
    //  FGETS2(&input,&allocbufsize,data->infile);
    myfree(input);
    myfree(word);
}


///
/// read comments and check for haplotype specification and mutation model
///
void read_comments(option_fmt *options, data_fmt *data, world_fmt *world)
{
  char thesign='\0';
  long linesize = LINESIZE;
  char * input = (char *) mycalloc(linesize,sizeof(char));
  FGETS2(&input, &linesize, data->infile);
  while(input[0]=='#')
    {
      thesign = input[1];
      switch(thesign)
	{
	case '%': //haplotype specification
	  insert_individual_inlist(input+3, data, options);
	  break;
	case '$': // mutation model
	  read_mutationmodel_comments(input+2,data, world);
	  break;
	default:      // this is a plain comment
	  break;
	}
      FGETS2(&input, &linesize, data->infile);
    }
  unread_word(data->infile, "\n");
  unread_word(data->infile, input);
  myfree(input);
}

void
read_popheader (FILE * infile, data_fmt * data, world_fmt *world, long pop, long genomes)
{
    boolean havepopname = FALSE;
    long    minlength   = 0;
    long    lo;
    long    locus;
    char   *input;
    char   *word;
    long sublocus;
    long indregion;
    long region;
    long startlocus;
    long endlocus;
    long start;
    long end;
    //long lastchar=EOF;
    long linesize = LINESIZE;
    long readpos=0;
    long oldreadpos=0;
    input = (char *) mycalloc(linesize,sizeof(char));
    word = (char *) mycalloc(linesize,sizeof(char));
    

    // allows that sequence data can have different numbers of individuals for different loci
    // data syntax changes: #ind1 #ind2 #IND3 .... pop_name
    havepopname=FALSE;
    // with multiple regions we need to separate the numbers of individuals per regions 
    // from the population title, and also take care in case the number of 
    // individuals is not specified for all loci.
    if(data->numregions>1) 
    {
      FGETS2(&input,&linesize,data->infile);      //read whole line
      while(input[0]=='#')
	{
	  FGETS2(&input,&linesize,data->infile);      //read next line
	}
      readpos = read_word_delim(input,word," \t",TRUE);
      oldreadpos = readpos;
      indregion = atol(word);//set first numind value, this must be always present
      start = 0;
      end   = data->regions[0].endlocus;
      for(locus=start;locus<end;locus++)
	{
	  data->numind[pop][locus] = indregion;
	  for(sublocus=world->sublocistarts[locus]; sublocus < world->sublocistarts[locus+1]; sublocus++)
	    data->numalleles[pop][sublocus] = indregion * genomes;	
	}
      for(region=1; region < data->numregions; region++)
        {
	  oldreadpos = readpos;
	  readpos += read_word_delim(input+readpos, word," \t",TRUE);// if we encounter a # we treat the rest of the line as comment
	  if(isdigit(word[0]) && havepopname == FALSE )
	    {
	      indregion = atol(word);
	      startlocus = data->regions[region].startlocus;
	      endlocus   = data->regions[region].endlocus;
	      for(locus=startlocus; locus < endlocus; locus++)
		{
		  data->numind[pop][locus] = indregion;
		  for(sublocus=world->sublocistarts[locus]; sublocus < world->sublocistarts[locus+1]; sublocus++)
		    data->numalleles[pop][sublocus] = indregion * genomes;
		}
	    }
	  else
	    {
	      // encountered a letter and assume this is the population name
	      havepopname=TRUE;
	      if(input+oldreadpos != NULL)
		{
		  minlength = (long) strlen(input+oldreadpos);
		  minlength = MIN(minlength,80);
		  strncpy(data->popnames[pop],input+oldreadpos,(size_t) minlength);
		}
	      break;//this leaves the numind empty of not specified,  must be filled in later
	    }
        }
      if(!havepopname) // this is still false when all loci have numbers of ind specified
	  {
	    if(input+oldreadpos != NULL)
	      {
		minlength = (long) strlen(input+oldreadpos);
		minlength = MIN(minlength,80);
		sprintf(data->popnames[pop],"%-*s", (int) minlength, input+readpos);
	      }
	  }
	
      // fills numind for additional locus in case the numind was not specified
      start = data->subloci[locus-1];
      end = data->allsubloci;
      for(lo=locus; lo < data->loci; lo++)
	{
	  data->numind[pop][lo] = data->numind[pop][locus-1];
	  for(sublocus=start; sublocus < end; sublocus++)
	    data->numalleles[pop][sublocus] = data->numind[pop][lo] * genomes;
	}
    }
    else
      {
        // only one region so we can use old scheme and 
	// assume that all loci have the same number of individuals
        FGETS2(&input,&linesize,infile);
	while(input[0]=='#')
	  {
	    FGETS2(&input,&linesize,data->infile);
	  }
        sscanf (input, "%ld%[^\n]", &(data->numind[pop][0]), data->popnames[pop]);
	start = 1;
	end   = data->allsubloci;
	
	data->numalleles[pop][0] = data->numind[pop][0] * genomes;	
	for(sublocus=start; sublocus < end; sublocus++)
	  {
	    data->numalleles[pop][sublocus] = data->numind[pop][0] * genomes;	
	    data->numind[pop][sublocus] = data->numind[pop][0];
	  }
      }
    translate (data->popnames[pop], ' ', '_');
    translate (data->popnames[pop], '\t', '_');
    unpad(data->popnames[pop],"_");
    myfree(input);
    myfree(word);
}


void
read_indname (FILE * file, data_fmt * data, long pop, long ind, long locus, long nmlength)
{
    long i = 0;
    char ch;
    char input[LINESIZE];
    while (i < nmlength)
    {
      ch = (char) getc (file);
	while(ch =='#')
	  {
	    FGETS(input,LINESIZE,data->infile);
	    ch = (char) getc (file);
	  }
	// assign uncovered this proble that name have \t attached
	// weird that no other problem ever occured with this
        if(strchr("\t",ch))
            break;
        if(!strchr("\r\n",ch))
            data->indnames[pop][ind][locus][i++] = ch;
    }
    data->indnames[pop][ind][locus][nmlength] = '\0';
#ifdef DEBUG
    printf("indname: %s\n", data->indnames[pop][ind][locus]);
#endif
}


/// read a line from the datafile 
///
void read_data_line (FILE *infile, data_fmt *data, option_fmt * options, world_fmt * world, long pop, long ind, long region)
{
  long i;
  long locus=0;
  const long startlocus = data->regions[region].startlocus;
  const long endlocus   = data->regions[region].endlocus;

  for(i=0;i<data->loci;i++)
    {
      if(startlocus==world->sublocistarts[i])
	locus = i;
      if(endlocus==world->sublocistarts[i])
	break;
    }
  const long start = locus;
  const long end   = i;
  char c;
  read_indname (infile, data, pop, ind, locus, options->nmlength);//region was locus here
  if(data->oldsyntax)
    {
      if(options->datatype=='h')
	{
	  if ((c=(char) getc(infile))!='r')
	    {
	      ungetc(c,infile);
	      read_hapmap(infile, data, options, pop);
	    }
	  else
	    read_hapmap_genotypes(infile, data, options, pop);
	}
      else
	{
	  switch (options->datatype)
	    {
	    case 'a':
	    case 'b':
	    case 'm':
	      if (data->dlm == '\0')
		read_alleles (infile, data, pop, ind);
	      else
		read_microalleles (infile, data, options, pop, ind);
	      break;
	    case 's':
	    case 'n':
	    case 'u':
	    case 'f':
	      read_ind_seq (infile, data, options, locus, pop, ind, 0);
	      break;
	    default:
	      usererror
		("Wrong datatype, only the types a, m, s, h @\n       (electrophoretic alleles, \n       microsatellite data,\n       sequence data,\n       SNP polymorphism)\n       Hapmap data format\n        are allowed.\n");
	      //      break;
	    }
	}
    }
  else
    {
      if(data->haplotyping)
	{
	  store_indname(data->indnames[pop][ind][start],pop, ind, start, data);
	  for(i=start+1;i<end;i++)
	    {
	      strcpy(data->indnames[pop][ind][i],data->indnames[pop][ind][start]);
	      store_indname(data->indnames[pop][ind][i],pop, ind, i, data);
	    }
	}
      read_ind_seq_oneliner (infile, data, options, world, pop, ind, 0, startlocus, endlocus);
    }
}

/// reading all possible data the datafile need a new structure 
/// in which the second line is mandatory. If the dataype is not specified that it is read from the parmfile
/// the read_sites() reads this line and create unlinked regions (loci on different lines), unlinked loci on the
/// the same same line, and linked loci on the same line.
/// regions can have different #ind, loci on the same line cannot, the data->numind is set by the population header
/// and must specify the #ind per region if  all have the same then one number is OK, see there.
void read_popdata (FILE * infile, data_fmt * data, long pop, option_fmt * options, world_fmt *world)
{
  long locus;
  long region;
  long ind;
  char c;
  if(options->datatype == 'h')
    {
      if ((c=(char) getc(infile))!='r')
	{
	  ungetc(c,infile);
	  read_hapmap(infile, data, options, pop);
	}
      else
	read_hapmap_genotypes(infile, data, options,  pop);
    }
  else
    {
      for(region=0; region < data->numregions; region++)
	{
	  locus = data->regions[region].startlocus;
	  for (ind = 0; ind < data->numind[pop][locus]; ind++)
	    {
	      read_data_line (infile, data, options, world, pop, ind, region);
	    }
	}
    }
}

long check_list(long la1, long nla1, long locus, data_fmt *data)
{
  if(la1 > data->numrepeatnumbers[locus])
    {
      data->repeatnumbers[locus] = (long *) myrealloc(data->repeatnumbers[locus], sizeof(long) * (size_t) (la1+1));
      memset(data->repeatnumbers[locus]+data->numrepeatnumbers[locus], 0, sizeof(long) * (size_t) (la1+1-data->numrepeatnumbers[locus]));
      data->numrepeatnumbers[locus]= la1+1;
      data->repeatnumbers[locus][la1]=nla1;
    }
  else
    {
      if(data->repeatnumbers[locus][la1]!=0)
	return data->repeatnumbers[locus][la1];
    }
  return la1;
}

void len2repeat(char *a1, long rlen)
{
  long la1 = (long) atol(a1);
  long a=0;
  if((a=la1 % rlen) != 0)
    {
      //      long correction = (((((float) rlen) / 2.) < a) ? (-a) : (rlen/2 == a ? (UNIF_RANDUM()<0.5 ? -a : a) : (rlen - a))); 
      long correction = ((((rlen) / 2.) < a) ? (-a) : (rlen/2 == a ? (UNIF_RANDUM()<0.5 ? -a : a) : (rlen - a))); 
      //long nla = check_list(la1,la1+correction, locus, data);
      long nla = (la1+correction)/rlen;
      sprintf(a1,"%li",nla);
    }
  else
    {
      sprintf(a1,"%li",la1/rlen);
    }
}


void
read_microalleles (FILE * infile, data_fmt * data, option_fmt *options, long pop, long ind)
{
    char *input, *isave, dlm[2], ddlm[2], *p, *a, *a1, *a2;
    long locus, i;
    input = (char *) mycalloc (1, sizeof (char) * (SUPERLINESIZE + 1));
    isave = input;
    a = (char *) mycalloc (1, sizeof (char) * LINESIZE);
    a1 = (char *) mycalloc (1, sizeof (char) * LINESIZE);
    a2 = (char *) mycalloc (1, sizeof (char) * LINESIZE);
    dlm[0] = data->dlm;
    dlm[1] = '\0';
    ddlm[0] = ' ';
    ddlm[1] = '\0';
    FGETS (input, SUPERLINESIZE, infile);
    if ((p = (char *) strpbrk (input, CRLF)) != NULL)
        *p = '\0';
    for (locus = 0; locus < data->allsubloci; locus++)
    {
        while (isspace ((int) *input))
            input++;
        if (input[0] == '\0')
            FGETS (input, SUPERLINESIZE, infile);
        i = 0;
        while (strchr(" \t",input[i])==NULL && input[i] != dlm[0])
        {
            a1[i] = input[i];
            i++;
        }
        a1[i] = '\0';
        input += i;
        i = 0;
        if (input[i] == dlm[0])
        {
            input++;
            while (strchr(" \t",input[i])==NULL && input[i] != '\0')
            {
                a2[i] = input[i];
                i++;
            }
            a2[i] = '\0';
            if (a2[0] == '\0')
            {
                strcpy (a2, a1);
            }
            input += i;
        }
        else
        {
            strcpy (a2, a1);
        }
        sprintf (data->yy[pop][ind][locus][0][0], "%-.*s",(int) options->allelenmlength,a1);
        sprintf (data->yy[pop][ind][locus][1][0], "%-.*s",(int) options->allelenmlength,a2);
    }
    myfree(a);
    myfree(a1);
    myfree(a2);
    myfree(isave);
}

void
read_alleles (FILE * infile, data_fmt * data, long pop, long ind)
{
    char *input, *isave, *p, *a;
    long locus;
    a = (char *) mycalloc (1, sizeof (char) * LINESIZE);

    input = (char *) mycalloc (1, sizeof (char) * SUPERLINESIZE);
    isave = input;
    FGETS (input, SUPERLINESIZE, infile);
    if ((p = (char *) strpbrk (input, CRLF)) != NULL)
        *p = '\0';
    for (locus = 0; locus < data->allsubloci; locus++)
    {
        while (isspace ((int) *input))
        {
            input++;
        }
        if (sscanf (input, "%s", a) == 1)
        {
            input += (long) strlen (a);
        }

        data->yy[pop][ind][locus][0][0][0] = a[0];
        data->yy[pop][ind][locus][0][0][1] = '\0';
        if (a[1] == '\0')
        {
            data->yy[pop][ind][locus][1][0][0] = a[0];
            data->yy[pop][ind][locus][1][0][1] = '\0';
        }
        else
        {
            data->yy[pop][ind][locus][1][0][0] = a[1];
            data->yy[pop][ind][locus][1][0][1] = '\0';
        }
    }
    myfree(a);
    myfree(isave);
}

///
/// reads the standard hapmap genotype data file format using this header
/// rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode individual# ....
void
read_hapmap_genotypes (FILE * infile, data_fmt * data, option_fmt * options, long pop)
{
  int error;
  long ind;
    char label1;
    char label2;
    long label1count;
    long label2count;
    long label12count;
    char *input;
    long locus;
    long genomes = number_genomes(options->datatype);
    input = (char *) mycalloc(LONGLINESIZE,sizeof(char));
    // read header
    FGETS(input,LONGLINESIZE,infile);
    for(locus=0;locus<data->loci;locus++)
      {
	//rs# or comment
	read_word(infile,input, NULL);
	if(input[0]=='#')
	  {
	    fprintf(stdout,"%s\n",input);
	    FGETS(input,LONGLINESIZE,infile);
	    locus--;
	    continue;
	  }
	// discard rs#
	// read allele and record
	read_word(infile,input, NULL);
	label1 = input[0];// format is like this G/C
	label2 = input[2];
	// read chromosome and discard
	read_word(infile,input, NULL);
	// read position
	read_word(infile,input, NULL);
	data->position[locus] = atol(input);
	if(options->printdata)
	  {
	    fprintf(stdout,"%20li|",data->position[locus]);
	  }
        // read strand must be positive but not checked
	read_word(infile,input, NULL);
	// read assembly# and discard
	read_word(infile,input, NULL);
	// read center and discard
	read_word(infile,input, NULL);
	// read protLSID and discard
	read_word(infile,input, NULL);
	// read assyLSID and discard
	read_word(infile,input, NULL);
	// read panelLSID and discard
	read_word(infile,input, NULL);
	// read QCODE and use as dummy
	read_word(infile,input, NULL);
 	data->seq[0]->sites[locus]=1;
	label1count=0;
	label2count=0;
	read_word(infile,input, NULL);
	while(strstr(input,"rs")==NULL && !isdigit(input[0]))
	  {
	    if (label1 == input[0])
	      label1count++;
	    else
	      label2count++;
	    if (label2 == input[1])
	      label2count++;
	    else
	      label1count++;
	    error = read_word(infile,input, NULL);
	    if (error==EOF)
	      break;
	  }
	unread_word(infile,input);
	label12count = label1count + label2count;
	if(options->printdata)
	  {
	    fprintf(stdout, " freq[%c]=%f freq[%c]=%f\n",label1,
		    (double) label1count/label12count,label2,  (double) label2count/label12count);
	  }
	//was messing with the reading, use fraction data->numind[pop][locus] = label12count;
	if(options->randomsubset > 0 && options->randomsubset < data->numind[pop][locus])
	  {
	    data->numalleles[pop][locus] = options->randomsubset * genomes;
	    data->numind[pop][locus] = options->randomsubset / genomes;
	  }
	else
	  {
	    data->numalleles[pop][locus] = data->numind[pop][locus] = label12count;
	  //	  data->numalleles[pop][locus] = data->numind[pop][locus] * genomes;
	  }
	label1count = (long) (floor( (((double) label1count/label12count * data->numalleles[pop][locus]))+0.5 ));
	for(ind=0;ind < label1count; ind++)
	  {
	    data->yy[pop][ind][locus][0][0][0] = label1;
	  }
	for(ind=label1count;ind < data->numalleles[pop][locus]; ind++)
	  {
	    data->yy[pop][ind][locus][0][0][0] = label2;
	  }
      }
    myfree(input);
}

///
/// my own hapmap format derived from the frequency data
void
read_hapmap (FILE * infile, data_fmt * data, option_fmt * options,  long pop)
{
  long ind;
    char label1;
    char label2;
    long label1count;
    long label2count;
    long label12count;
    char *input;
    long locus;
    long genomes = number_genomes(options->datatype);
    input = (char *) mycalloc(LINESIZE,sizeof(char));
    for(locus=0;locus<data->allsubloci;locus++)
      {
	read_word(infile,input,NULL);
	if(input[0]=='#')
	  {
	    fprintf(stdout,"%s\n",input);
	    FGETS(input,LINESIZE,infile);
	    locus--;
	    continue;
	  }
	data->position[locus] = atol(input);
	if(options->printdata)
	  {
	    fprintf(stdout,"%20li|",data->position[locus]);
	  }
	data->seq[0]->sites[locus]=1;
	read_word(infile,input,NULL);
	label1 = input[0];;
	read_word(infile,input, NULL);
	label1count = atol(input);
	
	read_word(infile,input,NULL);
	label2 = input[0];
	read_word(infile,input,NULL);
	label2count = atol(input);
	read_word(infile,input,NULL);
	label12count = atol(input);
	fprintf(stdout,"%20li|",data->position[locus]);
	fprintf(stdout,"label12count=%li, not equal to label1count+label2count=%li+%li\n",
		label12count, label1count, label2count);

	if (label12count != label1count + label2count)
	  {
	    fprintf(stderr,"label12count=%li, not equal to label1count+label2count=%li+%li\n",
		    label12count, label1count, label2count);
	    warning("hapmap format has a counting error in locus %li and population %li\n",locus+1,pop+1);
	    error("The error need to be fixed");
	  }
	if(options->printdata)
	  {
	    fprintf(stdout, " freq[%c]=%f freq[%c]=%f\n",label1,
		    (double) label1count/label12count,label2,  (double) label2count/label12count);
	  }
	//was messing with the reading, use fraction data->numind[pop][locus] = label12count;
	if(options->randomsubset > 0 && options->randomsubset < data->numind[pop][locus])
	  {
	    data->numalleles[pop][locus] = options->randomsubset * genomes;
	    data->numind[pop][locus] = options->randomsubset / genomes;
	  }
	else
	  {
	    data->numalleles[pop][locus] = label12count;
	    data->numind[pop][locus] = label12count / genomes;
	  }
	// data->numind[pop][locus] * genomes;
	label1count = (long) (floor( (((double) label1count/label12count * data->numalleles[pop][locus]))+0.5 ));
	label12count = data->numalleles[pop][locus];
	for(ind=0;ind < label1count; ind++)
	  {
	    data->yy[pop][ind][locus][0][0][0] = label1;
	  }
	for(ind=label1count;ind < label12count; ind++)
	  {
	    data->yy[pop][ind][locus][0][0][0] = label2;
	  }
      }
    myfree(input);
}

long
read_ind_seq (FILE * infile, data_fmt * data, option_fmt * options,
              long locus, long pop, long ind, long baseread)
{
    long j;
#ifdef DEBUG
    printf("IN READ IND SEQ\n");
#endif
    char charstate;
    j = (options->interleaved) ? baseread : 0;
    charstate = (char) getc (infile);
    ungetc ((int) charstate, infile);
    if(options->printdata)
      {
    	fprintf(stdout,"%s|",data->indnames[pop][ind][locus]);
      }
    while (j < data->seq[0]->sites[locus]
            && !(options->interleaved && strchr(CRLF,charstate)))
    {
      charstate = (char) getc (infile);
        if (strchr(CRLF,charstate))
        {
#ifdef INTERLEAVED
            if (options->interleaved)
            {
                while(strchr(CRLF,charstate=getc(infile)))
                   ;
                ungetc ((int) charstate, infile);
                return j;
            }
            else
#endif
                charstate = ' ';
        }
        if (charstate == ' '
                || (charstate >= '0' && charstate <= '9') || charstate == '\\')
            continue;
        charstate = uppercase (charstate);
	if(options->printdata)
	{
	  fprintf(stdout, "%c",charstate);
	}
        if ((strchr ("ABCDGHKMNRSTUVWXY?O-", (int) charstate)) == NULL)
        {
            printf
            ("ERROR: BAD BASE: %c AT POSITION %5ld OF INDIVIDUUM %3li in POPULATION %ld\n",
             charstate, j, ind+1, pop+1);
            printf
            ("Last complete sequences was in population %li, individual %li and locus %li:\n%s",
             pop + 1, ind, locus+1, data->indnames[pop][ind - 1][locus]);
            for (j = 0; j < data->seq[0]->sites[locus]; j++)
                printf ("%c", data->yy[pop][ind - 1][locus][0][j][0]);
            exit (EXIT_FAILURE);
        }
        data->yy[pop][ind][locus][0][j++][0] = charstate;
    }
    charstate = (char) getc (infile); /* swallow the \n or \r */
    while(charstate == '\r' || charstate == '\n' || charstate == '\t' || charstate == ' ' || charstate == ';')
      {
	charstate = (char) getc(infile);
      }
    if(charstate!='\n')
      ungetc((int) charstate, infile);

    if(options->printdata)
     {
       fprintf(stdout,"\n");
     }
    return j;
}

void read_site(FILE *infile, char *site, long locus, long sublocus, long pop, long ind, world_fmt *world, data_fmt *data, option_fmt *options)
{
  long j;
  enum siteclass_enum dataclass = world->mutationmodels[sublocus].dataclass;
  if(dataclass==SITECHARACTER)
    {
      boolean done=FALSE;
      char charstate = '\0'; 
      while(!done)
	{
	  charstate = (char) getc (infile);
	  if (strchr(CRLF,charstate))
	    {
	      charstate = ' ';
	    }
	  if (charstate == ' '
	      || (charstate >= '0' && charstate <= '9') || charstate == '\\')
	    {
	      charstate = ' ';
	      continue;
	    }
	  else
	    {
	      done=TRUE;
	    }
	}
      charstate = uppercase (charstate);
      if(options->printdata)
	{
	  fprintf(stdout, "%c",charstate);
	}
      if ((strchr ("ABCDGHKMNRSTUVWXY?O-", (int) charstate)) == NULL)
	{
	  if(ind>0)
	    {
	      printf
		("Error: BAD BASE:%c. Last complete sequences was in population %li, individual %li and locus %li:\n%s",
		 charstate, pop + 1, ind, locus+1, data->indnames[pop][ind - 1][locus]);
	      for (j = 0; j < world->mutationmodels[/*world->sublocistarts[locus]+*/sublocus].numsites;j++)
		printf ("%c", data->yy[pop][ind - 1][locus][0][j][0]);
	    }
	  else
	    {
	      printf
		("Error: BAD BASE:%c. failed in first individual in population %li, individual %li and locus %li:\n%s",
		 charstate, pop + 1, ind, locus+1, data->indnames[pop][ind][locus]);
	      for (j = 0; j < world->mutationmodels[/*world->sublocistarts[locus]+*/sublocus].numsites;j++)
		printf ("%c", data->yy[pop][ind][locus][0][j][0]);
	    }
	  exit (EXIT_FAILURE);
	}
      site[0]=charstate;
    }
  else
    {
      read_word(data->infile,site,NULL);
    }
}


///
/// reads sequence data that comes where all loci are (1) on one line and (2) are not interleaved
/// this function allows a mixture of different datatypes using the mutationmodel to guide the reading and 
/// distribution into the different subloci and loci
long read_ind_seq_oneliner (FILE * infile, data_fmt * data, option_fmt * options, world_fmt *world,
			    long pop, long ind, long baseread, long startlocus, long endlocus)
{
#ifdef DEBUG
    fprintf(stdout,"IN READ IND SEQ ONELINER\n");
    fflush(stdout);
#endif
  (void) baseread;
  long j = 0;
  long locus=0;
  long sublocus;
  char charstate;
  char *site;
  long oldendsite = 0;
  long allocsite=LINESIZE;
  site = (char *) mycalloc(allocsite,sizeof(char));
  charstate = (char) getc (infile);
  ungetc ((int) charstate, infile);
  if(options->printdata)
    {
      fprintf(stdout,"%s",data->indnames[pop][ind][locus]);
    }
  //  for(locus = 0; locus < data->loci; locus++)
  for(locus = startlocus; locus < endlocus; locus++)
    {
      j=0;
      if(options->printdata)
	{
	  fprintf(stdout,"|");
	}
      //	while (j < data->seq[0]->sites[locus])
      
      const long sublocistart = world->sublocistarts[locus];
      const long sublociend   = world->sublocistarts[locus+1];
      for(sublocus = sublocistart; sublocus < sublociend; sublocus++)
	{
	  j   = oldendsite;
	  long numsites = world->mutationmodels[sublocus].numsites;
	  // random etc reading from data
	  // this translates the (s100) or (s10 4)
	  // into the yy matrix allowing for small fractions of data used from large genomic chunks
	  long startsite = world->mutationmodels[sublocus].startsite;
#ifdef DEBUG
	  printf("\nsublocus=%li numsites=%li startlocus=%li endlocus=%li\n",sublocus,numsites,startlocus, endlocus);
#endif
	  numsites += startsite;
	  if (startsite == 0)
	    {
	      j=0; // to supercede oldendstate
	    }
#ifdef DEBUG
	  printf("%li %li %li %li\n@",locus, sublocus, startsite, oldendsite);
#endif
	  while(j < numsites)
	    {	     
	      read_site(data->infile, site, locus, sublocus, pop, ind, world, data, options);
	      if (j>=startsite)
		{
		  //printf(".");
		  memcpy(data->yy[pop][ind][sublocus][0][j-startsite], site, sizeof(char) * strlen(site));
		}
	      j++;
	    }
	  oldendsite = j;
	  //printf("\n");
	  if(options->printdata)
	    {
	      fprintf(stdout,"|");
	    }
	}
    }
  //if (isspace(*site))
    FGETS2(&site,&allocsite,infile);
    //else
    //ungetc(*site,infile);
  if(options->printdata)
    {
      fprintf(stdout,"\n");
    }
  myfree(site);
  return j;
}

void
read_distance_fromfile (FILE * dfile, long tips, long nmlength, MYREAL **m)
{
    char input[SUPERLINESIZE];
    long i, j;
    //int retval;
    if (dfile != NULL)
    {
        // assumes that the dfile is in PHYLIP format
        // and that all n x n cells are filled.
        FGETS (input, LONGLINESIZE, dfile); //reads header line with
        for (i = 0; i < tips; i++) // of individuals
        {
            //reads the population label
            FGETS (input, nmlength + 1, dfile);
            for (j = 0; j < tips; j++)
            {
#ifdef USE_MYREAL_FLOAT
                fscanf (dfile, "%f", &m[i][j]);
#else
                fscanf (dfile, "%lf", &m[i][j]);
#endif
		if((i!=j) && (m[i][j] < EPSILON))
		  {
		    warning("Reading geofile: adjusting dist[%li][%li]=%f to %f\n",i,j,m[i][j],EPSILON);
		    m[i][j] = EPSILON;
		  }
            }
            // reads the last \n from the
            // data matrix
            FGETS (input, LONGLINESIZE, dfile);
        }
    }
}

#ifdef UEP
// uep function

void
read_uep_fromfile (FILE * uepfile, long tips, long nmlength, int **uep,
                   long *uepsites, long datatype)
{
    char input[LINESIZE];
    long i, j;
    long thistips;
    if (uepfile != NULL)
    {
        // assumes that the uepfile is in PHYLIP format
        // and that all n cells are filled.
        FGETS (input, LINESIZE, uepfile); //reads header line with
        // of individuals and uep sites
        sscanf (input, "%li%li", &thistips, uepsites);
        if (thistips != tips)
            error ("UEP datafile and infile are inconsistent!");
        if (strchr (SEQUENCETYPES, datatype))
        {
            for (i = 0; i < tips; i++)
            {
                uep[i] = (int *) mycalloc (*uepsites, sizeof (int));
                FGETS (input, nmlength + 1, uepfile); //reads each line
                for (j = 0; j < *uepsites; ++j)
                    fscanf (uepfile, "%i", &uep[i][j]);
                // reads the last \n from the data matrix
                FGETS (input, LINESIZE, uepfile);
            }
        }
        else
        {
            for (i = 0; i < tips; i++)
            {
                uep[i] = (int *) mycalloc (*uepsites, sizeof (int));
                uep[i + tips] = (int *) mycalloc (*uepsites, sizeof (int));
                FGETS (input, nmlength + 1, uepfile); //reads each line
                for (j = 0; j < *uepsites; ++j)
                    fscanf (uepfile, "%i", &uep[i][j]);
                // finished reading first allele, no onto the second
                for (j = 0; j < *uepsites; ++j)
                    fscanf (uepfile, "%i", &uep[i + tips][j]);
                // reads the last \n from the data matrix
                FGETS (input, LINESIZE, uepfile);
            }
        }
    }
}
#endif

void
finish_read_seq (FILE * infile, data_fmt * data, option_fmt * options,
                 long pop, long baseread)
{

    long ind, baseread2 = 0, locus = 0;
    if (options->interleaved)
    {
        while (baseread < data->seq[0]->sites[0])
        {
            for (ind = 0; ind < data->numind[pop][0]; ind++)
            {
                baseread2 =
                    read_ind_seq (infile, data, options, locus, pop, ind,
                                  baseread);
            }
            baseread = baseread2;
        }
    }
    for (locus = 1; locus < data->loci; locus++)
    {
        baseread = 0;
        for (ind = 0; ind < data->numind[pop][locus]; ind++)
        {
	  read_indname (infile, data, pop, ind, locus, options->nmlength);
            baseread = read_ind_seq (infile, data, options, locus, pop, ind, 0);
        }
        if (options->interleaved)
        {
            while (baseread < data->seq[0]->sites[locus])
            {
                for (ind = 0; ind < data->numind[pop][locus]; ind++)
                {
                    baseread2 =
                        read_ind_seq (infile, data, options, locus, pop, ind,
                                      baseread);
                }
                baseread = baseread2;
            }
        }
    }
}

long find_missing(data_fmt *data, long pop, long locus)
{
    long missing = 0;
    long ind;
    for(ind=0; ind < data->numind[pop][locus]; ind++)
    {
        if(data->yy[pop][ind][locus][0][0][0]=='?')
            missing++;
        if(data->yy[pop][ind][locus][1][0][0]=='?')
            missing++;
    }
    return missing;
}


///
/// Data set was subsampled and used a random sample of size
void print_random_subset(FILE * file, data_fmt * data, option_fmt *options)
{
  long locus;
  long pop;
  char *name;
  long maxnum;
  long count;
  long ind;
  long length;
  long index;
  if(options->randomsubset > 0)
    {
      fprintf (file, "\nData set was subsampled and used a random sample of size: %li\n\n", options->randomsubset);
      fprintf (file, "Locus Population Individuals\n");
      fprintf (file, "----- ---------- ---------------------------------------------------------------------\n");
      name = (char*) mycalloc(SMALLBUFSIZE,sizeof(char));
      for (locus=0; locus < data->loci; locus++)
	{
	  fprintf(file,"%5li ",locus+1);
	  for (pop=0; pop < data->numpop; pop++)
	    {
	      fprintf(file, "%-10.10s ", data->popnames[pop]);
	      maxnum = options->randomsubset < data->numind[pop][locus] ? options->randomsubset : data->numind[pop][locus];
	      count = 18;
	      for(ind=0;ind<maxnum;ind++)
		{
		  index = data->shuffled[pop][locus][ind];
		  memset(name,0,sizeof(char)*SMALLBUFSIZE);
		  if(data->indnames[pop][index][locus][0]=='\0')
		    memcpy(name,data->indnames[pop][index][0],sizeof(char) * (size_t) options->nmlength);
		  else
		    memcpy(name,data->indnames[pop][index][locus],sizeof(char) * (size_t) options->nmlength);
		  remove_trailing_blanks(&name);
		  if (options->has_datefile)
		    {
		      sprintf(name,"%s (%f) ",name, data->sampledates[pop][locus][ind].date);
		    }
		  length = (long) strlen(name);
		  if (count+length < LINELENGTH)
		    {
		      fprintf(file,"%s ",name);
		      count += length+1;
		    }
		  else
		    {
		      fprintf(file,"\n");
		      fprintf(file,"                 ");
		      fprintf(file,"%s ",name);
		      count=17+length+1;
		    }
		}
	      if(pop<data->numpop-1)
		fprintf(file,"\n      ");
	      else
		fprintf(file,"\n");
	    }
	}
      myfree(name);
    }
  fprintf(file,"\n\n");
}

void
print_ratetbl (FILE * outfile, world_fmt * world, option_fmt * options,
           long locus, char header)
{
  (void) options;
  boolean doprint=TRUE;
  long i;
  mutationmodel_fmt *s, *sp;
  long sublocus;
  const long sublocistart = world->sublocistarts[locus];
  const long sublociend   = world->sublocistarts[locus+1];
  if(header=='A')
    {
      FPRINTF (outfile, "\nLocus Sublocus Region type     Rate of change    Probability  Patch size\n");
      FPRINTF (outfile, "--------------------------------------------------------------------------\n");
    }
  if(header=='C')
    {
      FPRINTF (outfile, "\nLocus Sublocus Region type     Rate of change    Probability  Patch size\n");
      FPRINTF (outfile, "--------------------------------------------------------------------------\n");
      FPRINTF (outfile, "[compressed - only loci that are different than the one before are shown]\n");
    }
  for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
    {
      s = &world->mutationmodels[sublocus];
      if(header=='G')
	{
	  doprint=FALSE;
	  sp = &world->mutationmodels[sublocus-1];
	  for (i = 0; i < s->numsiterates; i++)
	    {
	      if (s->siterates[i] < sp->siterates[i])
		{
		  if (s->siterates[i] > sp->siterates[i])
		    {
		      continue;
		    }
		  else
		    {
		      doprint=TRUE;
		      break;
		    }
		}
	    }
	}
      if(doprint)
	{
	  for (i = 0; i < s->numsiterates; i++)
	    FPRINTF (outfile, "%4ld%8ld%9ld%16.3f%17.3f%17.3f\n", locus + 1, sublocus + 1- sublocistart, i + 1, s->siterates[i],
		     s->siteprobs[i], 1.0 / s->lambda);	  
	  FPRINTF (outfile,"\n");
	}
    }
}

void
print_categtbl (FILE * outfile, world_fmt * world, option_fmt * options,
           long locus, char header)
{
  (void) header;
      long i;
    mutationmodel_fmt *s;
    long sublocus;
    const long sublocistart = world->sublocistarts[locus];
    const long sublociend   = world->sublocistarts[locus+1];

    if (options->categs > 1)
    {
        FPRINTF (outfile, "Locus Sublocus Site category   Rate of change\n");
        FPRINTF (outfile, "---------------------------------------------\n");
	for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
	  {
	    s = &world->mutationmodels[sublocus];
	    for (i = 0; i < s->numcategs; i++)
	      FPRINTF (outfile, "%4ld%8ld%9ld%16.3f\n", locus + 1, sublocus + 1, i + 1, s->rate[i]);
	  }
    }
    FPRINTF (outfile, "\n");
}

void set_datatype_string(char datatype, char * dstring)
{
  switch (datatype)
    {
    case 'a':
      strcpy (dstring, "Allelic data");
      break;
    case 'f':
    case 's':
      strcpy (dstring, "Sequence data");
      break;
    case 'b':
      strcpy (dstring, "Microsatellite data [Brownian]");
      break;
    case 'm':
      strcpy (dstring, "Microsatellite data [Stepwise]");
      break;
    case 'n':
    case 'u':
      strcpy (dstring, "SNP data");
      break;
    case 'h':
      strcpy (dstring, "SNP data (Hapmap data)");
      break;
    case '@':
      strcpy (dstring, "Haplotype data");
      break;
    default:
      strcpy (dstring, "Unknown data [ERROR]");
      break;
    }
}

void
print_data_summary (FILE * file, world_fmt * world, option_fmt * options,
                    data_fmt * data)
{
  int maxlength = 24; // length of the the total string , for popnames
  int length;
  long locus;
  long pop;
  long numind;
  long nummiss;
  char dstring[LINESIZE];
  char modelname[LINESIZE];
  char modelparam[LINESIZE];
  long *total;
  long *totalmiss;
  boolean terse = options->tersepdf;
  total = (long *) mycalloc(data->loci,sizeof(long));
  totalmiss = (long *) mycalloc(data->loci,sizeof(long));
  set_datatype_string(options->datatype, dstring);
  fprintf (file, "\nSummary of data:\n");
  if(options->prioralone)
    {
      fprintf(file,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      fprintf(file,"Program is using NO DATA -- running WITHOUT DATA\n");
      fprintf(file,"option: NODATA=yes is set, to run real data remove this option\n");
      fprintf(file,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    }
  fprintf(file, "Title:%54.54s\n",options->title);
  fprintf (file, "Data file:  %48.48s\n", options->infilename);
  fprintf (file, "Datatype:   %48.48s\n", dstring);
  //  if (!strchr (MSATTYPES, options->datatype) && options->datatype!='@')
  if (dstring[0]=='M')
    {
      if(data->has_repeats==FALSE)
	{
	  fprintf (file, "  [Fragment length is translated to repeats]\n");
	}
      else
	{
	  fprintf (file, "  [Data was used as repeat-length information]\n");
	}
    }
  if(options->has_datefile)
    {
      fprintf (file, "Sample dates:%46.46s\n", options->datefilename);
      fprintf (file, "Generations per year:                        %13.8f\n", options->generation_year);
      fprintf (file, "Mutationrate per year: %.10g", options->mutationrate_year[0]);
      for(locus=1; locus < options->mutationrate_year_numalloc; locus++)
	{
	  if(locus % 4 == 0)
	    fprintf(file,"                      ");
	  fprintf(file,", %f", options->mutationrate_year[locus]); 
	} 
      fprintf(file,"\n"); 
    }
  fprintf (file, "Number of loci:                         %20li\nMutationmodel:\n",
	   data->loci);
  
  fprintf(file," Locus  Sublocus  Mutationmodel   Mutationmodel parameter\n");
  fprintf(file,"-----------------------------------------------------------------\n");
  long printloci = data->loci;
  if(terse)
    {
      printloci = TEN;
      fprintf(file,"Only the first 10 loci are shown\n");
    }
  for (locus=0; locus < printloci; locus++)
    {
      long   sublocus;
      long   sublocistart = world->sublocistarts[locus];
      long   sublociend = world->sublocistarts[locus+1];
      for(sublocus=sublocistart;sublocus<sublociend;sublocus++)
	{	  
	  get_mutationmodel_nameparam(modelname,modelparam, &world->mutationmodels[sublocus]);
	  fprintf(file,"%6li  %8li %-15.15s %s\n",locus+1,sublocus-sublocistart+1,modelname,modelparam);
	}
    }
  if(options->totalsites>0 && strchr(SNPTYPES,options->datatype))
    {
      fprintf(file,"[On average there are %li invariant sites per locus]\n",options->totalsites);
    }
    
  // 
  print_random_subset(file,data,options);
  if(file!=stdout)
    {
      pdf_print_random_subset(data,options);
    }
  for (pop = 0; pop < data->numpop; pop++)
    {
      length = (int) strlen(data->popnames[pop]);
      if(maxlength < length)
	 maxlength = length;
    }
    if(maxlength > 40)
      maxlength = 40;



    //print used sites per locus
    if (!(!strchr (SEQUENCETYPES, options->datatype) && options->datatype!='@'))
      {
	boolean compressed = FALSE;
	boolean print_siterates = FALSE;
	fprintf (file,"Sites per locus\n");
	fprintf (file,"---------------\n");
	if (data->loci > TEN)
	  compressed=TRUE;
	if (compressed)
	  {
	    long mini = MAXLONG;
	    long maxi = 0;
	    for(locus=0; locus< data->loci; locus++)
	      {
		long   m = 0;
		long   sublocus;
		long   sublocistart = world->sublocistarts[locus];
		long   sublociend = world->sublocistarts[locus+1];
		for(sublocus=sublocistart;sublocus<sublociend;sublocus++)
		  {
		    m = world->mutationmodels[sublocus].numsites;
		    if (world->mutationmodels[sublocus].numsiterates > 1)
		      print_siterates = TRUE;
		    if (m<mini)
		      mini = m;
		    if (m > maxi)
		      maxi = m;
		  }
	      }
	    fprintf(file,"%li loci with minimal of %li sites and maximal number of %li sites per sublocus\n",
		    data->loci, mini, maxi);

	    if (print_siterates)
	      {
		fprintf(file,"\n");
		fprintf (file,"Site Rate variation per locus\n");
		fprintf (file,"-----------------------------\n");
		print_ratetbl (file, world, options,0,'C');
		for(locus=1; locus< data->loci; locus++)
		  {
		    print_ratetbl (file, world, options,locus,'G');//followup but compressed
		  }
	      }
	    fprintf(file,"\n");

	  }
	else
	  {
	    fprintf (file,"Locus    Sites\n");
	    for(locus=0; locus< data->loci; locus++)
	      {
		fprintf(file,"%6li    ",locus+1);
		long   sublocus;
		long   sublocistart = world->sublocistarts[locus];
		long   sublociend = world->sublocistarts[locus+1];
		for(sublocus=sublocistart;sublocus<sublociend;sublocus++)
		  {
		    if (world->mutationmodels[sublocus].numsiterates > 1)
		      print_siterates = TRUE;
		    fprintf(file," %li",world->mutationmodels[sublocus].numsites);
		  }
		fprintf(file,"\n");
	      }
	    if (print_siterates)
	      {
		fprintf(file,"\n");
		fprintf (file,"Site Rate variation per locus\n");
		fprintf (file,"-----------------------------\n");
		print_ratetbl (file, world, options,0,'A');
		for(locus=1; locus< data->loci; locus++)
		  {
		    print_ratetbl (file, world, options,locus,'F');
		  }
	      }
	    fprintf(file,"\n");
	    if (!strchr (SEQUENCETYPES, options->datatype) && options->datatype!='@')
	      {
		fprintf (file,"%-*.*s     Locus   Gene copies    \n",maxlength,maxlength,"Population");
		fprintf (file,"%-*.*s             ---------------\n",maxlength, maxlength, " ");
		fprintf (file,"%-*.*s             data  (missing)\n",maxlength, maxlength, " ");
	      }
	    else
	      {
		fprintf (file,"%-*.*s     Locus   Gene copies    \n",maxlength,maxlength,"Population");
	      }
	    fprintf (file,"----%-*.*s------------------------\n",maxlength,maxlength,"------------------------------------------------------------------");
	    for (pop = 0; pop < data->numpop; pop++)
	      {
		for(locus=0; locus< data->loci; locus++)
		  {
		    long   sublocus;
		    long   sublocistart = world->sublocistarts[locus];
		    long   sublociend = world->sublocistarts[locus+1];
		    for(sublocus=sublocistart;sublocus<sublociend;sublocus++)
		      {
			if (!strchr (SEQUENCETYPES, options->datatype) && options->datatype!='@')
			  {
			    nummiss = find_missing(data,pop,sublocus);
			    numind = data->numalleles[pop][sublocus] - nummiss;
			    fprintf (file, "%3li %-*.*s %5li %6li (%li)\n", options->newpops[pop], maxlength,maxlength,(locus==0 ? data->popnames[pop] : " " ), locus+1 , numind, nummiss);
			  }
			else
			  {
			    nummiss = 0;
			    numind = data->numind[pop][locus];
			    fprintf (file, "%3li %-*.*s %5li    %6li\n", options->newpops[pop], maxlength,maxlength,(locus==0 ? data->popnames[pop] : " "), locus+1 , numind);
			  }
			total[locus] += numind;
			totalmiss[locus] += nummiss;
		      }
		  }
	      }
	    if (!strchr (SEQUENCETYPES, options->datatype) && options->datatype != '@')
	      {
		for(locus=0; locus< data->loci; locus++)
		  {
		    fprintf (file,"    %-*.*s %5li %6li (%li)\n",maxlength,maxlength, 
			     (locus == 0 ? "Total of all populations" : " "), locus+1, total[locus], totalmiss[locus]);
		  }
	      }
	    else
	      {
		for(locus=0; locus< data->loci; locus++)
		  {
		    fprintf (file,"    %-*.*s %5li    %6li\n",maxlength,maxlength, 
			     (locus == 0 ? "Total of all populations" : " "), locus+1, total[locus]);
		  }
	      }    
	  }
      }
    fprintf(file,"\n");
    myfree(total);
    myfree(totalmiss);
    fflush (file);
}

void
print_data (world_fmt * world, option_fmt * options, data_fmt * data)
{
    if (options->printdata)
    {
        switch (options->datatype)
        {
        case 'a':
	  //	  if(options->dlm=='\0')
          //  print_alleledata (world, data, options);
	  //else
	  print_microdata (world, data, options);
	  break;
        case 'b':
        case 'm':
            print_microdata (world, data, options);
            break;
        case 's':
        case 'n':
	case 'h':
        case 'u':
        case 'f':
            print_seqdata (world, options, data);
            break;
	case '@':
	  printf ("not implemented yet");
        }
    }
}

///
/// calculate allele frequency spectrum and print
/// allele population1 .. populationN All
/// taking into account population labeling and random_subsets
/// worker only calculate!!!! only master also prints
void print_spectra(world_fmt * world, option_fmt * options,data_fmt * data)
{
  const double one = 1.0;
  const double two = 2.0;
  long found;
  long locus;
  long a;
  long pop;
  long pop1;
  long ind;
  MYREAL **total;
  MYREAL *grandtotal;
  MYREAL allfreq;
  long nummaxalleles;
  MYREAL f;
  MYREAL homo = 0.0;
  MYREAL *avghet;
  MYREAL v;
  MYREAL avghet1;
  //MYREAL avghetall = 0.0;
  long *maxalleles;
  long *maxallelepop;
  MYREAL ***freq;
  char *thisallele;
  char *thatallele;
  MYREAL fx;
  MYREAL general_homo;
  FILE *outfile = world->outfile;
  long loctotal;
  long maxstates;
  boolean second = (strchr(ALLELETYPES,options->datatype)!=NULL);
  if(options->datatype=='@')
    return;
  if (strchr (DNASEQUENCETYPES, options->datatype))
    return; /*we do not calculate allele frequencies for DNA sequences*/
  // find total number of alleles
  maxalleles = (long *) mycalloc(data->loci, sizeof(long));
  maxallelepop = (long *) mycalloc(data->numpop, sizeof(long));
  avghet = (MYREAL *) mycalloc(data->numpop, sizeof(MYREAL));
  for (locus = 0; locus < data->loci; locus++)
    {
      const long sublocistart = world->sublocistarts[locus];
      const long sublociend   = world->sublocistarts[locus+1];
      long sublocus;
      for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
	{
	  mutationmodel_fmt *s = &world->mutationmodels[sublocus];
	  long smallest;
	  long biggest;
	  // this simple counts the alleles in data->alleles -- I hope
	  switch(s->datatype)
	    {
	    case 'a':
	    case 'b':
	    case 'm':
	      nummaxalleles = findAllele(data,"\0",sublocus);
	      if(nummaxalleles == 0)
		nummaxalleles = get_states(s, data, locus); 
	      s->numstates = nummaxalleles;
	      if (s->datatype=='a')
		s->maxalleles = nummaxalleles+1;
	      else
		s->maxalleles = nummaxalleles;
	      s->freq = 1.0 / s->maxalleles;
	      maxalleles[locus] = nummaxalleles;
	      break;
	    case 'x': // this needs to be fixed 
	      find_minmax_msat_allele (world, data, locus, &smallest, &biggest);
	      s->numstates  = biggest - smallest; //get_states(s,world); //biggest - smallest;
	      s->micro_threshold = options->micro_threshold;
	      if(s->micro_threshold > biggest-smallest)
		s->micro_threshold = biggest-smallest;//s->numstates;
	      if(s->micro_threshold % 2 != 0)
		s->micro_threshold += 1; // to make it even
	      s->microstart = smallest - 1 - s->micro_threshold;
	      if (s->microstart < 0)
		s->microstart = 0;
	      s->maxalleles = biggest - smallest + 2 * s->micro_threshold; //this is the scaled max
	      s->microrange = s->maxalleles;
	      s->basefreqs  = (double *) mycalloc((s->numstates+BASEFREQLENGTH-4),sizeof(double));
	      //printf("%i> micro-threshold=%li\n", myID, s->micro_threshold);
	      doublevec2d(&s->steps,s->micro_threshold,s->micro_threshold);
	      break;
	    case 'h':
	    case 'n':
	      maxalleles[locus] = s->maxalleles;
	    default:
	      s->numstates = get_states(s,data,locus);
	      break;
	    }
	  if(s->basefreqs==NULL)
	    s->basefreqs = (double *) mycalloc((s->numstates+BASEFREQLENGTH-4),sizeof(double));
	}
    } 
  // create bins for each population and all
  grandtotal = (MYREAL *) mycalloc(data->allsubloci, sizeof(MYREAL));
  freq = (MYREAL ***) mycalloc(data->numpop, sizeof(MYREAL **));
  doublevec2d(&total,data->numpop,data->allsubloci);
  for (pop = 0; pop < data->numpop; pop++)
    {
      freq[pop] = (MYREAL **) mycalloc(data->allsubloci + 1, sizeof(MYREAL *));
      maxstates = 0;
      for (locus = 0; locus < data->allsubloci; locus++)
	{
	  freq[pop][locus] = (MYREAL *) mycalloc(world->mutationmodels[locus].numstates+1, sizeof(MYREAL));
	  //if (maxstates < world->mutationmodels[locus].numstates)
	  //  maxstates =  world->mutationmodels[locus].numstates;	   
	}
      // over all loci
      freq[pop][locus] = (MYREAL *) mycalloc(maxstates, sizeof(MYREAL));
    }
  
  // calculate spectrum
  for (pop1 = 0; pop1 < data->numpop; pop1++)
    {
      //pop = options->newpops[pop1]-1;
      for (ind = 0; ind < data->numind[pop1][0]; ind++)
        {
	  for (locus = 0; locus < data->allsubloci; locus++)
	    {
	      thisallele = data->yy[pop1][ind][locus][0][0];
	      found = findAllele(data, thisallele, locus);
	      if(second)
		{
		  thatallele = data->yy[pop1][ind][locus][1][0];
		  if(!strcmp(thisallele,thatallele))
		    {
		      if(!strchr(thisallele,'?'))
			{
			  freq[pop1][locus][found] += two;
			  //total[pop][locus] += two;
			}
		    }
		  else	    
		    {
		      if(!strchr(thisallele,'?'))
			{
			  freq[pop1][locus][found] += one;
			  //total[pop][locus] += one;
			}
		      found = findAllele(data, thatallele, locus);
		      if(!strchr(thatallele,'?'))
			{
			  freq[pop1][locus][found] += one;
			  //total[pop][locus] += one;
			}
		    }
		}
	      else
		{
		  if(!strchr(thisallele,'?'))
		    {
		      freq[pop1][locus][found] += one;
		      //total[pop][locus] += one;
		    }
		}
	    }
	}
    }
  // now we calculate the poptotals and popfreqs if present
  // this will change freq and total therefore transform the
  // location score to population scores that are repeated
  for (pop1 = 0; pop1 < data->numpop; pop1++)
    {
      pop = options->newpops[pop1]-1;
      if (pop == pop1)
	continue;
      for (locus = 0; locus < data->allsubloci; locus++)
	{
	  for(a=0; a < world->mutationmodels[locus].numstates; a++)
	    {
	      freq[pop][locus][a] += freq[pop1][locus][a];
	    }
	  //replace all locations frequencies with the population freq
	  for(a=0; a < maxalleles[locus]; a++)
	    {
	      freq[pop1][locus][a] = freq[pop][locus][a];
	    }
	}
    }
  // now we have collected all numbers of alleles for each locus and each location
  // we calculate now all the totals for each locality, 
  for (pop = 0; pop < data->numpop; pop++)
    {
      for (locus = 0; locus < data->allsubloci; locus++)
	{
	  loctotal = 0;
	  for(a=0; a < maxalleles[locus]; a++)
	    {
	      loctotal += (long) freq[pop][locus][a];
	    }
	  total[pop][locus] = loctotal;
	  grandtotal[locus] += loctotal;
	}
    }

  // print in ascii
  //fprintf(stdout,"Allele frequency spectra\n");
  if (myID == MASTER)
    {
      fprintf(outfile,"Allele frequency spectra\n");
      fprintf(outfile,"========================\n\n");
      avghet1=0.0;
      for (locus = 0; locus < data->allsubloci; locus++)
	{
	  mutationmodel_fmt *s = &world->mutationmodels[locus];
	  fprintf(outfile,"Locus %li\n", locus + 1);
	  fprintf(outfile,"Allele  ");
	  general_homo = 0.0;
	  //@@ the migrate-old is different here
	  for (pop1 = 0; pop1 < data->numpop; pop1++)
	    {
	      maxallelepop[pop1] = 0;
	      pop = options->newpops[pop1]-1;
	      fprintf(outfile,"Pop%-2li  ",pop+1);
	    }
	  fprintf(outfile,"All\n-------");
	  for (pop = 0; pop < data->numpop+1; pop++)
	    {
	      fprintf(outfile,"-------");
	    }
	  fprintf(outfile, "\n");
	  //	  for(a=0; a < s->maxalleles; a++)
	  //for(a=0; a < s->numstates; a++)
	  for(a=0; a < s->maxalleles; a++)
	    {
	      //printf("%i> states=%li data->allele[%li][%li]=%s\n",myID, s->numstates, locus, a, data->allele[locus][a]);
	      allfreq = 0.0;
	      fprintf(outfile,"%6s ",data->allele[locus][a]);
	      for (pop1 = 0; pop1 < data->numpop; pop1++)
		{
		  pop = options->newpops[pop1]-1;
		  if(freq[pop][locus][a]>0.0)
		    {
		      maxallelepop[pop] += 1;
		      fprintf(outfile," %1.3f ",freq[pop][locus][a]/total[pop][locus]);
		      allfreq += freq[pop][locus][a];
		    }
		  else
		    {
		      fprintf(outfile,"   -   ");
		    }
		}
	      fx = allfreq/grandtotal[locus];
	      fprintf(outfile," %1.3f\n", (MYREAL) fx);	      
	      general_homo += fx * fx;	   
	    }
	  fprintf(outfile,"Alleles");
	  for (pop1 = 0; pop1 < data->numpop; pop1++)
	    {
	      //pop = options->newpops[pop1]-1;
	      fprintf(outfile,"%5li  ",maxallelepop[pop1]);
	    }
	  fprintf(outfile,"%5li\n", s->maxalleles);
	  //fprintf(outfile,"%5li\n", s->numstates);

	  fprintf(outfile,"Samples");
	  for (pop1 = 0; pop1 < data->numpop; pop1++)
	    {
	      fprintf(outfile,"%5li  ",  (long) total[pop1][locus]);
	    }
	  fprintf(outfile,"%5li\n", (long) grandtotal[locus]);
	  

	  if(maxalleles[locus]<=1 && strchr(SNPTYPES,options->datatype))
	    {
	      data->skiploci[locus] = TRUE;
	    }
	  fprintf(outfile,"H_exp  ");
	  for (pop1 = 0; pop1 < data->numpop; pop1++)
	    {
	      pop = options->newpops[pop1]-1;
	      homo = 0.0;	  
	      for(a=0;a<maxalleles[locus];a++)
		{
		  f = freq[pop][locus][a]/total[pop][locus];
		  homo += f*f;
		}
	      v = 1.0 - homo;
	      avghet[pop1] += v;
	      fprintf(outfile," %5.3f ",v);
	    }
	  fprintf(outfile," %5.3f\n\n",1.0-general_homo);
	  avghet1 += 1.0 - general_homo;
	  //      fprintf(outfile," %5.3f\n\n",avghet1/data->numpop);
	}
      fprintf(outfile,"Average expected heterozygosity\n");
      for (pop1 = 0; pop1 < data->numpop; pop1++)
	{
	  pop = options->newpops[pop1]-1;
	  fprintf(outfile,"Pop%-2li  ",pop+1);
	}
      fprintf(outfile,"All\n");
      for (pop = 0; pop < data->numpop+1; pop++)
	fprintf(outfile,"-------");
      fprintf(outfile, "\n");
      //      fprintf(outfile," %5li\n", maxalleles[locus]);
      //fprintf(outfile,"H_exp      ");
      for (pop1 = 0; pop1 < data->numpop; pop1++)
	{
	  //pop = options->newpops[pop1]-1;
	  fprintf(outfile,"%5.3f  ",avghet[pop1] / data->loci);
	}
      fprintf(outfile,"%5.3f\n\n",avghet1/data->loci);
      // printd in PDF
      
#ifdef PRETTY
      pdf_print_spectra(world, data, options, freq, total, grandtotal, avghet, avghet1, maxalleles);
#endif
      fflush(outfile);
    }
  // cleanup
  for (pop = 0; pop < data->numpop; pop++)
    {
      for (locus = 0; locus < data->allsubloci; locus++)
	{
	  myfree(freq[pop][locus]);
	}
      myfree(freq[pop]);
    }
#ifdef DEBUG
  printf("finished with printspectra()");
#endif
  myfree(freq);
  // myfree(maxalleles);
  myfree(maxallelepop);
  myfree(avghet);
  free_doublevec2d(total);
  myfree(grandtotal);
}

void
print_alleledata (world_fmt * world, data_fmt * data, option_fmt * options)
{
    long i, pop, ind, locus, mult80;
    for (pop = 0; pop < data->numpop; pop++)
    {
        print_header (world->outfile, pop, world, options, data);
        for (ind = 0; ind < data->numind[pop][0]; ind++)
        {
            fprintf (world->outfile, "%-*.*s ", (int) options->nmlength,
                     (int) options->nmlength, data->indnames[pop][ind][0]);
            mult80 = options->nmlength;
            for (locus = 0; locus < data->allsubloci; locus++)
            {
                mult80 +=
                    1 + (long) (strlen (data->yy[pop][ind][locus][0][0]));
                if (mult80 >= 80)
                {
                    mult80 = 0;
                    fprintf (world->outfile, "\n");
                    for (i = 0; i < options->nmlength; i++)
                        FPRINTF(world->outfile, " ");
                }
                fprintf (world->outfile, " %c.%c",
                         data->yy[pop][ind][locus][0][0][0],
                         data->yy[pop][ind][locus][0][1][0]);
            }
            fprintf (world->outfile, "\n");
        }
        fprintf (world->outfile, "\n");
    }
    fprintf (world->outfile, "\n\n");
    fflush (world->outfile);
}

void
print_microdata (world_fmt * world, data_fmt * data, option_fmt * options)
{
    long i, pop, ind, locus, mult80;
    for (pop = 0; pop < data->numpop; pop++)
    {
        print_header (world->outfile, pop, world, options, data);
        for (ind = 0; ind < data->numind[pop][0]; ind++)
        {
            fprintf (world->outfile, "%-*.*s ", (int) options->nmlength,
                     (int) options->nmlength, data->indnames[pop][ind][0]);
            mult80 = options->nmlength;
            for (locus = 0; locus < data->loci; locus++)
            {
                mult80 +=
                    1 + (long) (strlen (data->yy[pop][ind][locus][0][0]) +
                    strlen (data->yy[pop][ind][locus][1][0]));
                if (mult80 >= 80)
                {
                    mult80 = 0;
                    fprintf (world->outfile, "\n");
                    for (i = 0; i < options->nmlength; i++)
                        FPRINTF(world->outfile, " ");
                }
                fprintf (world->outfile, " %s.%-s",
                         data->yy[pop][ind][locus][0][0],
                         data->yy[pop][ind][locus][1][0]);
            }
            fprintf (world->outfile, "\n");
        }
        fprintf (world->outfile, "\n");
    }
    fprintf (world->outfile, "\n\n");
    fflush (world->outfile);
}

void
print_seqdata (world_fmt * world, option_fmt * options, data_fmt * data)
{
    long pop, locus;
    for (pop = 0; pop < data->numpop; pop++)
    {
        print_header (world->outfile, pop, world, options, data);
        for (locus = 0; locus < data->loci; locus++)
        {
            print_locus_head (locus, world, options, data);
            print_seq_pop (locus, pop, world, options, data);
        }
    }
    fflush (world->outfile);
}

void
print_header (FILE * outfile, long pop, world_fmt * world,
              option_fmt * options, data_fmt * data)
{
    long i;
    long locus, mult80 = 80;
    fprintf (outfile, "\n%-s", data->popnames[pop]);
    for (i = 0; i < (long) (80 - (long) strlen (data->popnames[pop])); i++)
         fprintf(world->outfile, "-");
    fprintf (outfile, "\n\n");
    if (!strchr (SEQUENCETYPES, options->datatype))
    {
        fprintf (outfile, "%-s  ", (data->loci == 1 ? "locus" : "loci "));
        for (i = 0; i < (options->nmlength - 6); i++)
            fprintf(world->outfile, " ");
        for (locus = 0; locus < data->loci; locus++)
        {
            if (locus * 4 + options->nmlength > mult80)
            {
                mult80 += 80;
                fprintf (outfile, "\n");
                for (i = 0; i < options->nmlength; i++)
                    fprintf (outfile, " ");
            }
            fprintf (outfile, "  %2li", locus + 1);
        }
        fprintf (outfile, "\n%-*.*s\n", (int) options->nmlength, (int) options->nmlength, "indiv.");
    }
}


MYREAL findleastsquare(MYREAL *rawdata, long total, long repeatlength, long shift, MYREAL *startvalue)
{
  long i;
  long j;
  long high;
  long low;
  MYREAL lowvalue=MYREAL_MAX;
  MYREAL highvalue=0.;
  MYREAL sum = 0.0;
  MYREAL value;
  MYREAL oldvalue;
  for(i=0;i<total;i++)
    {
      if(rawdata[i] < lowvalue)
	lowvalue = rawdata[i];
      if(rawdata[i]> highvalue)
	highvalue = rawdata[i];
    }
  high = (long) (highvalue+repeatlength+shift);
  low  = (long) (lowvalue -repeatlength+shift);
  if(low<0){
    low = 0;
  }
  *startvalue = low;
  for(i=0;i<total;i++)
    {
      oldvalue = MYREAL_MAX;
      for(j=low; j < high; j += repeatlength)
	{
	  value = rawdata[i] - j;
	  value *= value;
	  if(value < oldvalue)
	    {
	      oldvalue = value;
	    }
	}
      sum += oldvalue;
    }
  return sum;
}

///
/// calculates and sets the microsatellite or brownian repeatlengths from fragmentlength data
/// this is executed when data->has_repeats == FALSE
void find_allele_repeatlength(data_fmt *data, option_fmt *options, long locus)
{
  (void) options;
  long z=0;
  long zz=0;
  long pop;
  long ind;
  long r;
  MYREAL *rawdata;
  //char *a1;
  //char *a2;
  long total;
  MYREAL minLS=MYREAL_MAX;
  MYREAL value;
  long intval;
  MYREAL startvalue= 0.0;
  MYREAL keepstartvalue= 0.0;
  MYREAL leastsquare;
  MYREAL nonfloored;
  MYREAL floored;
  MYREAL diff;
  for (pop = 0; pop < data->numpop; pop++)
    {
      zz += data->numind[pop][locus];
    }
  rawdata = (MYREAL *) mycalloc((1+ 2 * zz), sizeof(MYREAL));

  for (pop = 0; pop < data->numpop; pop++)
    {
      for (ind = 0; ind < data->numind[pop][locus]; ind++)
	{
	  char *a1 = data->yy[pop][ind][locus][0][0];
	  char *a2 = data->yy[pop][ind][locus][1][0];
	  if (a1[0]!='?')
	    {
	      rawdata[z++] = atof(a1);
	    }
	  if (a2[0]!='?')
	    {
	      rawdata[z++] = atof(a2);
	    }
	}
    }

  total=z;
  minLS=MYREAL_MAX;
  startvalue= 0.0;
  keepstartvalue= 0.0;

  for(r=0;r<data->repeatlength[locus];r++)
    {
      //MYREAL endvalue  = 0.0;
      leastsquare=findleastsquare(rawdata,total,data->repeatlength[locus],r,&startvalue);
      if(leastsquare < minLS)
	{
	  minLS = leastsquare;
	 // keepr = r;
	  keepstartvalue = startvalue;
	}
    }
  for (pop = 0; pop < data->numpop; pop++)
    {
      for (ind = 0; ind < data->numind[pop][locus]; ind++)
	{
	  char *a1 = data->yy[pop][ind][locus][0][0];
	  char *a2 = data->yy[pop][ind][locus][1][0];
	  if (a1[0]!='?')
	    {
	      value = atof(a1);
	      //Floor[(279.47 - 259)/4 + If[(Random[]) < 1 - FractionalPart[(279.47 - 259)/4 ], 0, 1]], {1000}] // Tally
	      nonfloored = (value - keepstartvalue)/data->repeatlength[locus]; 
	      floored = floor(nonfloored);
	      diff = nonfloored - floored; 
	      intval = MSAT_OFFSET + (long) floored + (RANDUM() < (1.0-diff) ? 0 : 1); 
	      sprintf(data->yy[pop][ind][locus][0][0],"%li",intval);
	    }
	  if (a2[0]!='?')
	    {
	      value = atof(a2);
	      nonfloored = (value - keepstartvalue)/data->repeatlength[locus]; 
	      floored = floor(nonfloored);
	      diff = nonfloored - floored; 
	      intval =  MSAT_OFFSET + (long) floored + (RANDUM() < (1.0-diff) ? 0 : 1); 
	      sprintf(data->yy[pop][ind][locus][1][0],"%li",intval);
	    }
	}
    }
}


MYREAL
create_alleles (world_fmt *world, data_fmt * data, option_fmt *options)
{
  long n=0;
  MYREAL  mean=0.;
  MYREAL  delta=0.;
  long locus, pop, ind;
  long z;
  char *a1;//DEFAULT_ALLELENMLENGTH];
  char *a2;//DEFAULT_ALLELENMLENGTH];
  boolean second = (strchr(ALLELETYPES,options->datatype)!=NULL);
  long lena1 = LINESIZE;
  long lena2 = LINESIZE;
  a1 = (char *) malloc(sizeof(char)*LINESIZE);
  a2 = (char *) malloc(sizeof(char)*LINESIZE);
  for (locus = 0; locus < data->loci; locus++)
    {
      long sublocus;
      const long sublocistart = world->sublocistarts[locus];
      const long sublociend   = world->sublocistarts[locus+1];
      for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
	{
	  mutationmodel_fmt *s= &world->mutationmodels[sublocus];
	  z = 0;
	  if(data->repeatlength[sublocus]!=0)
	    find_allele_repeatlength(data, options,sublocus);
	  for (pop = 0; pop < data->numpop; pop++)
	    {
	      for (ind = 0; ind < data->numind[pop][locus]; ind++)
		{

		  char *data1 = data->yy[pop][ind][locus][0][0];
		  long lendata1 = (long) strlen(data1);
		  if (lena1 < lendata1)
		    {
		      a1 = (char *) myrealloc(a1,sizeof(char)* (size_t) (lendata1+1));
		      lena1 = lendata1 + 1;
		    }
		  strcpy (a1, data1);
		  if(second)
		    {
		      char *data2 = data->yy[pop][ind][locus][1][0];
		      long lendata2 = (long) strlen(data2);
		      if (lena2 < lendata2)
			{
			  a2 = (char *) myrealloc(a2,sizeof(char)* (size_t) (lendata2+1));
			  lena2 = lendata2 + 1;
			}
		      strcpy (a2, data2);
		      if (strcmp (a1, a2))
			{
			  addAllele (data, a1, locus, &z);
			  addAllele (data, a2, locus, &z);
			}
		      else
			{
			  addAllele (data, a1, locus, &z);
			}
		    }
		  else
		    {
		      addAllele (data, a1, locus, &z);
		    }
		}
            }
	  if(z==0)
	    {
	      data->skiploci[locus] = TRUE;
	      continue;
	    }	 
	  if (data->datatype[sublocus]=='a')
 	    {
	      data->maxalleles[sublocus] = z + 1;
	      s->maxalleles = z + 1;
	    }
	  else
	    {
	      data->maxalleles[sublocus] = z;
	      s->maxalleles = z;
	    }
	  /* + 1: for all the unencountered alleles */
	  if(options->murates_fromdata)
	    {
	      if(options->mu_rates==NULL)
		{
		  options->mu_rates = (MYREAL * ) mycalloc(data->loci,sizeof(MYREAL));
		}
	      options->mu_rates[locus] = z+1;
	      n = n + 1;
	      delta = options->mu_rates[locus] - mean;
	      mean += delta/n;
	    }
	}
    }
  myfree(a1);
  myfree(a2);
  if(options->murates_fromdata)
      {
	options->muloci = data->loci;
	for (locus=0; locus < data->allsubloci; locus++)
	  {
	    if(!data->skiploci[locus])
	      options->mu_rates[locus] /= mean;
	  }
      }
    return mean;
}


MYREAL
create_mixed_data (data_fmt * data, option_fmt *options)
{
  long n=0;
  MYREAL  mean=0.;
  MYREAL  delta=0.;
  long locus, pop, ind;
  long z;
  boolean second = (strchr(ALLELETYPES,options->datatype)!=NULL);
  char a1[DEFAULT_ALLELENMLENGTH];
  char a2[DEFAULT_ALLELENMLENGTH];
  for (locus = 0; locus < data->allsubloci; locus++)
    {
      if(data->repeatlength[locus]!=0)
	find_allele_repeatlength(data, options,locus);
      z = 0;
      for (pop = 0; pop < data->numpop; pop++)
        {
            for (ind = 0; ind < data->numind[pop][locus]; ind++)
            {
                strcpy (a1, data->yy[pop][ind][locus][0][0]);
		if(second)
		  {
		    strcpy (a2, data->yy[pop][ind][locus][1][0]);
		    if (!strcmp (a1, a2))
		      {
			addAllele (data, a1, locus, &z);
		      }
		    else
		      {
			addAllele (data, a1, locus, &z);
			addAllele (data, a2, locus, &z);
		      }
		  }
		else
		  {
		    addAllele (data, a1, locus, &z);
		  }
	    }
        }
	if(z==0)
	  {
	    data->skiploci[locus] = TRUE;// this skips loci with _no_ data
	    continue;
	  }

        data->maxalleles[locus] = z + 1;
        /* + 1: for all the unencountered alleles */
	if(options->murates_fromdata)
	  {
	    if(options->mu_rates==NULL)
	      {
		options->mu_rates = (MYREAL * ) mycalloc(data->loci,sizeof(MYREAL));
		//		printf("%i> data.c: 1557 murate size %li\n",myID,data->loci * sizeof (MYREAL));
	      }
	    options->mu_rates[locus] = z+1;
	    n = n + 1;
	    delta = options->mu_rates[locus] - mean;
	    mean += delta/n;
	  }
    }
  // warning("unfinished in data.c:3840");
  //  return mean;
    if(options->murates_fromdata)
      {
	options->muloci = data->loci;
	for (locus=0; locus < data->loci; locus++)
	  {
	    if(!data->skiploci[locus])
	      options->mu_rates[locus] /= mean;
	  }
      }
    return mean;
}

void
addAllele (data_fmt * data, char s[], long locus, long *z)
{
    long found = 0;
    if(!strcmp("?",s))
      return;
    while ((data->allele[locus][found++][0] != '\0')
            && (strcmp (s, data->allele[locus][found - 1])))
       ;
    if (found > (*z))
    {
        strcpy (data->allele[locus][*z], s);
        (*z)++;
    }
}

void
set_numind (data_fmt * data)
{
    long locus, pop;
    for (locus = 1; locus < data->loci; locus++)
    {
        for (pop = 0; pop < data->numpop; pop++)
        {
            data->numind[pop][locus] = data->numind[pop][0];
            data->numalleles[pop][locus] = data->numalleles[pop][0];
        }
    }
}


void
print_seq_pop (long locus, long pop, world_fmt * world, option_fmt * options,
               data_fmt * data)
{
    long ind;
    for (ind = 0; ind < data->numalleles[pop][locus]; ind++)
    {
        print_seq_ind (locus, pop, ind, world, options, data);
    }
}

void
print_seq_ind (long locus, long pop, long ind, world_fmt * world,
               option_fmt * options, data_fmt * data)
{
    long site;
    char blank[2] = " ";
    fprintf (world->outfile, "%-*.*s", (int) options->nmlength,
             (int) options->nmlength, data->indnames[pop][ind][0]);
    fprintf (world->outfile, " %c", data->yy[pop][ind][locus][0][0][0]);
    for (site = 1; site < data->seq[0]->sites[locus]; site++)
    {
        if ((site) % 60 == 0)
        {
            fprintf (world->outfile, "\n%-*.*s %c", (int) options->nmlength,
                     (int) options->nmlength, blank,
                     data->yy[pop][ind][locus][0][site][0]);
        }
        else
        {
            if ((site) % 10 == 0)
            {
                fprintf (world->outfile, " ");
            }
            fprintf (world->outfile, "%c", data->yy[pop][ind][locus][0][site][0]);
        }
    }
    fprintf (world->outfile, "\n");
}


void
print_locus_head (long locus, world_fmt * world, option_fmt * options,
                  data_fmt * data)
{
  (void) data;
    char *head;
    head = (char *) mycalloc (1, sizeof (char) * (size_t) (MAX (10, options->nmlength)));
    sprintf (head, "Locus %li", locus);
    fprintf (world->outfile, "%-*.*s --------10 --------20 --------30",
             (int) options->nmlength, (int) options->nmlength, head);
    fprintf (world->outfile, " --------40 --------50 --------60\n");

    myfree(head);
}

void
read_geofile (data_fmt * data, option_fmt * options, long numpop)
{
    long i, j, pop;
    long numpop2 = numpop * numpop;
    data->geo = (MYREAL *) mycalloc (1, sizeof (MYREAL) * (size_t) numpop2);
    data->lgeo = (MYREAL *) mycalloc (1, sizeof (MYREAL) * (size_t) numpop2);
    if (!options->geo)
    {
        for (i = 0; i < numpop2; i++)
            data->geo[i] = 1.0;
    }
    else
    {
      data->ogeo = (MYREAL **) mycalloc (1, sizeof (MYREAL *) * (size_t) numpop);
      data->ogeo[0] = (MYREAL *) mycalloc (1, sizeof (MYREAL) * (size_t) numpop2);
        for (pop = 1; pop < numpop; pop++)
            data->ogeo[pop] = data->ogeo[0] + numpop * pop;
        read_distance_fromfile (data->geofile, numpop, options->nmlength,
                                data->ogeo);
        for (i = 0; i < numpop; i++)
        {
            for (j = 0; j < numpop; j++)
            {
                if(i!=j)
                {
                    data->geo[mm2m (i, j, numpop)] =   1. / data->ogeo[i][j];
                    data->lgeo[mm2m (i, j, numpop)] =  data->ogeo[i][j] > 0.0 ?
                                                       LOG (1. / data->ogeo[i][j]) : -MYREAL_MAX;
                }
            }
        }
    }
}

///
/// read the file with the tip dates and returns the oldest date
MYREAL read_date_fromfile (FILE * datefile, data_fmt *data, option_fmt *options, long nmlength)
{
    char input[LINESIZE];
    long pop;
    long locus;
    //long l;
    long ind;
    MYREAL oldest = 0.;
    MYREAL youngest = DBL_MAX;
    MYREAL temp;
    char *name;
    boolean backward = TRUE;
    name = (char *) mycalloc(LINESIZE,sizeof(char));
    if (datefile != NULL)
    {
        FGETS (input, LINESIZE, datefile); //title line
	while(input[0]=='#')
	  FGETS(input,LINESIZE,datefile);
	if(input[0]=='F') //this checks the direction of time
	  backward=FALSE; 
	fprintf(stdout,"\n%i> Tip dates from file %s\n----------------------------------------\n", myID, options->datefilename);	
	fprintf(stdout,"Generations per year:            %f\n", options->generation_year);
	for(pop=0;pop<data->numpop;pop++)
	  {
	    FGETS (input, LINESIZE, datefile); //first populations 
	    while(input[0]=='#')
	      FGETS(input,LINESIZE,datefile);
	    for(ind=0; ind < data->numind[pop][0]; ind++)
	      {
		FGETS (input, LINESIZE, datefile);
		while(input[0]=='#')
		  FGETS(input,LINESIZE,datefile);
#ifdef USE_MYREAL_FLOAT
		sscanf (input, "%s",name);
		sscanf (input+nmlength, "%f",&temp);
#else
		sscanf (input, "%s",name);
		sscanf (input+nmlength, "%lf",&temp);
#endif
		for (locus=0; locus < data->loci; locus++)
		  {
		    //l = (locus >= options->mutationrate_year_numalloc ? options->mutationrate_year_numalloc -1 : locus);
		    //fprintf(stdout,"Mutation rate of Locus %li: %g\n", locus, options->mutationrate_year[l]);
		    
		    //fprintf(stdout,"Locus %li: Tipdate %*.*s %f %f\n", locus, (int) nmlength, (int) nmlength, input, temp, temp * options->mutationrate_year[l] / options->generation_year);
		    data->sampledates[pop][locus][ind].date = temp;
		    data->sampledates[pop][locus][ind].name = (char *) mycalloc(strlen(name)+1,sizeof(char));
		    strcpy(data->sampledates[pop][locus][ind].name, name);
		    unpad(data->sampledates[pop][locus][ind].name," ");
		    translate(data->sampledates[pop][locus][ind].name,' ', '_');
		  }
		if(oldest < temp)
		  oldest = temp;
		if(youngest > temp)
		  youngest = temp;
	      }
	  }
	// are the times forward or backward?
	for(pop=0;pop<data->numpop;pop++)
	  {
	    for (locus=0; locus < data->loci; locus++)
	      {
		for(ind=0; ind < data->numind[pop][locus]; ind++)
		  {
		    if(backward)
		      {
			data->sampledates[pop][locus][ind].date =
			  data->sampledates[pop][locus][ind].date - youngest;
		      }
		    else
		      {
			data->sampledates[pop][locus][ind].date = 
			  (oldest - data->sampledates[pop][locus][ind].date) 
			  + (youngest - oldest);
		      }
		  }
	      }
	  }
    }
    free(name);
    return oldest;
}

///
/// read the file with the tip dates
void
read_datefile (data_fmt * data, option_fmt * options, long numpop)
{
  (void) numpop;
  long locus, pop;
  data->sampledates = (tipdate_fmt ***) mycalloc (data->numpop, sizeof (tipdate_fmt **));
  data->sampledates[0] = (tipdate_fmt **) mycalloc ((data->numpop * data->loci), sizeof (tipdate_fmt *));
  for (pop = 1; pop < data->numpop; pop++)
    {
      data->sampledates[pop] = data->sampledates[0] + data->loci * pop;
    }
  for(locus=0;locus<data->loci;locus++)
    {
      for(pop=0;pop < data->numpop; pop++)
	{
	  data->sampledates[pop][locus] = (tipdate_fmt*) mycalloc(data->numind[pop][locus],sizeof(tipdate_fmt));
	}
    }

  if (options->has_datefile)
    {
      data->maxsampledate = read_date_fromfile (data->datefile, data, options, options->nmlength);
      //      printf("%i> in data section maxsampledate=%f\n",myID, data->maxsampledate);
    }
  else
    {
      data->maxsampledate=0.0;
    }
}

#ifdef NEWSTUFF
void  read_header_fasta(FILE **infile, fasta_fmt * f)
{
}
void  read_fasta(FILE **infile, fasta_fmt * f)
{
}
void  read_fastq(FILE **infile, fasta_fmt * f)
{
}

long get_new_fastq(FILE *infile, data_fmt * data, options_fmt *options, world_fmt * world)
{
  long locus=0;
  size_t i=0;
  while (*infile != NULL)
    {
      read_header_fasta(&infile,data->fasta);
      read_fasta(&infile,data->fasta);
      read_fastq(&infile,data->fasta);
      i++;
    }
  return locus;
}

#endif


#ifdef UEP
void
read_uepfile (data_fmt * data, option_fmt * options, long numpop)
{
    long i;
    long sumtips = 0;

    if (!options->uep)
        return;

    for (i = 0; i < numpop; ++i)
        sumtips += data->numind[i][0];   //Assumes that UEP has the same number of individuals as
    // locus 1 (Is this OK? most dataset with UEP will have 1 locus?????)
    data->uep = (int **) mycalloc (number_genomes (options->datatype) * sumtips,
                                 sizeof (int *));
    read_uep_fromfile (data->uepfile, sumtips, options->nmlength, data->uep,
                       &data->uepsites, options->datatype);
}

#endif
