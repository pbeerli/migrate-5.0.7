/* \file tree.c */
/*------------------------------------------------------
Maximum likelihood estimation
of migration rate  and effectice population size
using a Metropolis-Hastings Monte Carlo algorithm
-------------------------------------------------------
	T R E E B U I L D I N G   R O U T I N E S

	Peter Beerli 1996, Seattle
	beerli@fsu.edu

	Copyright 1997-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
	Copyright 2003-2016 Peter Beerli, Tallahassee FL

	some code in this file are successors of code in dnaml in the PHYLIP
        package of Joseph Felsenstein. Several changes were made to the conditionl
        likelihood methods to improve speed, but the original design idea is Joe's.

	A new model that includes treatment of gaps in the alignment is in the works
        and may superceede all other models because it seems that this model
        is capable to handle all other sequence models. <THIS IS NOY WORKING YET>
 
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


$Id: tree.c 2169 2013-08-24 19:02:04Z beerli $

-------------------------------------------------------*/
#ifndef HAS_STRDUP
#ifndef _XOPEN_SOURCE
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#define _XOPEN_SOURCE 600
#pragma clang diagnostic pop
#endif
#endif
#include <stdlib.h>
#include <math.h>
#include "migration.h"
#include "sighandler.h"
#include "assignment.h"
#include "random.h"
#include "options.h"
#include "data.h"
#include "sequence.h"
#include "world.h"
#include "tools.h"
#include "migrate_mpi.h"
#include "mcmc.h"
#include "mutationmodel.h"
#include "haplotype.h"
#include "speciate.h"
#include "tree.h"
#ifdef UEP
#include "uep.h"
#endif
#include "speciate.h"
#include "tree.h"
#include "skyparam.h"
#ifdef BEAGLE
#include "calculator.h"
#endif
#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif
#include <assert.h>

#define NOTIPS 0
#define WITHTIPS 1

extern long unique_id_global;
extern int simulator;

MYREAL treelike_anc (mutationmodel_fmt *s, long xs, world_fmt * world, long locus);
void treeout (FILE * treefile, node * joint, node * p, long s);
void allocate_tree (world_fmt * world,  data_fmt * data,
                    long locus);
void allocateinterior (world_fmt * world, option_fmt *options, data_fmt * data, long locus);
void allocatepoproot (world_fmt * world);
void allocate_tip (world_fmt * world, option_fmt * options, node ** p,
                   long pop, long locus, long a, long ind, char **tipnames);
void alloc_seqx (world_fmt * world, node * theNode, long locus);
void allocate_xseq(xarray_fmt *x, long sites, long categs);

/* first tree material (upgma, distance) */
void set_tree (world_fmt * world, option_fmt * options, data_fmt * data,
               long locus);
// sets migration events into a tree read from the user (dna only)
void set_migrations (world_fmt * world, long locus);
void distance_EP (char **data, long tips, MYREAL **m);
void distance_micro (char **data, long tips, MYREAL **m);
void distance_sequence (data_fmt * data, option_fmt *options, long locus, long tips, long sites,
                        long nmlength, MYREAL **m);
void distance_allele (world_fmt * world, option_fmt * options, long tips, MYREAL **distm);
void randomize_distm(MYREAL **distm, long tips, char *custm);
void constrain_distance_zeromig (MYREAL **m, option_fmt *options, data_fmt * data, long locus,
                                 long tips, char *custm);

void makevalues (world_fmt * world, option_fmt * options, data_fmt * data,
                 long locus);

void upgma (world_fmt * world, MYREAL **x, long tips, node ** nodep);
void set_top (world_fmt * world, node * p, long locus);
void set_v (node * p);
void calc_sancost (MYREAL **cost, world_fmt *world);


void free_treetimes (world_fmt * world, long size);
void traverseNodes (node * theNode, timelist_fmt ** timevector, long *slice, world_fmt *world, long *tips);
void increase_timelist (timelist_fmt ** timevector);
void increase_timelist2 (timelist_fmt ** timevector,long extender);
void allocate_lineages (timelist_fmt **timevector, const long offset, const long numpop);
void smooth (const node * root, node * p, world_fmt * world,
             const long locus);
void which_nuview (world_fmt *world, boolean fastlike, boolean use_gaps, int watkins);
void which_pseudonuview (char datatype, boolean fastlike, boolean use_gaps, int watkins);
void nuview_allele (mutationmodel_fmt *s, long sublocus, long xs, node * mother, world_fmt * world, const long locus);
void nuview_micro (mutationmodel_fmt *s, long sublocus, long xs, node * mother, world_fmt * world, const long locus);
void nuview_brownian (mutationmodel_fmt *s, long sublocus, long xs, node * mother, world_fmt * world, const long locus);
void nuview_sequence (mutationmodel_fmt *s, long sublocus, long xs, node * mother, world_fmt * world, const long locus);
void nuview_sequence_slow (mutationmodel_fmt *s, long sublocus, long xs, node * mother, world_fmt * world,
                           const long locus);
void nuview_ancestral (mutationmodel_fmt *s, long sublocus, long xs, node * mother, world_fmt * world, const long locus);
void adjustroot (node * r);
MYREAL pseudo_tl_seq (mutationmodel_fmt *s, long xs, phenotype xx1, phenotype xx2, MYREAL v1, MYREAL v2,
                      proposal_fmt * proposal, world_fmt * world);
 MYREAL pseudo_tl_snp (mutationmodel_fmt *s, long xs, phenotype xx1, phenotype xx2, MYREAL v1, MYREAL v2,
                       proposal_fmt * proposal, world_fmt * world);
 MYREAL pseudo_tl_snp_unlinked (mutationmodel_fmt *s, long xs, phenotype xx1, phenotype xx2, MYREAL v1,
                                MYREAL v2, proposal_fmt * proposal,
                                world_fmt * world);
MYREAL pseudo_tl_anc (mutationmodel_fmt *s, phenotype xx1, phenotype xx2, MYREAL v1, MYREAL v2,
                      proposal_fmt * proposal, world_fmt * world);
void pseudonu_allele    (mutationmodel_fmt *s, proposal_fmt * proposal, xarray_fmt * xxx1, MYREAL * lx1,
			 MYREAL v1, xarray_fmt * xxx2, MYREAL * lx2, MYREAL v2, long xs);
void pseudonu_micro     (mutationmodel_fmt *s, proposal_fmt * proposal, xarray_fmt * xxx1, MYREAL * lx1,
			 MYREAL v1, xarray_fmt * xxx2, MYREAL * lx2, MYREAL v2, long xs);
void pseudonu_brownian  (mutationmodel_fmt *s, proposal_fmt * proposal, xarray_fmt * xxx1, MYREAL * lx1,
			 MYREAL v1, xarray_fmt * xxx2, MYREAL * lx2, MYREAL v2, long xs);

void pseudonu_seq       (mutationmodel_fmt *s, proposal_fmt * proposal, xarray_fmt * xxx1, MYREAL * lx1,
			 MYREAL v1, xarray_fmt * xxx2, MYREAL * lx2, MYREAL v2, long xs);
void pseudonu_seq_slow  (mutationmodel_fmt *s, proposal_fmt * proposal, xarray_fmt * xxx1, MYREAL * lx1,
			 MYREAL v1, xarray_fmt * xxx2, MYREAL * lx2, MYREAL v2, long xs);
void pseudonu_anc       (mutationmodel_fmt *s, proposal_fmt * proposal, xarray_fmt * xxx1, MYREAL * lx1,
			 MYREAL v1, xarray_fmt * xxx2, MYREAL * lx2, MYREAL v2, long xs);
void calculate_steps (world_fmt * world);
MYREAL logfac (long n);


boolean treereader (world_fmt * world,  option_fmt *options,  data_fmt * data, long locus);
void length_to_times (node * p);
boolean treeread (FILE * file, world_fmt * world, option_fmt *options, node ** pp, node * q);
char processlength (FILE * file, node ** p);
node *allocate_nodelet (world_fmt *world, long num, char type);
void find_tips (node * p, node ** nodelist, long *z);
node *add_migration (world_fmt *world, node * p, char event, long from, long to, MYREAL utime);
node *create_interior_node (world_fmt * world, node ** q);
node *create_root_node (world_fmt *world, node ** q);
node *create_tip_node (FILE * file, world_fmt * world, option_fmt *options, node ** q, char *ch);
boolean processbracket (FILE * file, world_fmt *world, node ** p, char *ch);
void set_tree_pop (node * p, long *pop);
void allocate_x (node * p, world_fmt * world, long locus,
                 boolean withtips);
long find_firstpop (node * p);

void sankoff (world_fmt * world);
MYREAL minimum (MYREAL *vec1, MYREAL *vec2, long n);
void santraverse (world_fmt *world, node * theNode, MYREAL **cost, long numpop);
long ranbest (MYREAL *array, long tie, MYREAL best, long n);
void jumble (long *s, long n);
void jumble_ownseed (long *s, long n);
long number_genomes (int datatype);

/* copy whole tree */
void copy_tree (world_fmt * original, world_fmt * kopie);
node *copy_node (world_fmt * original, node * o, world_fmt * kopie,
                 node * last);
void copy_node_content (world_fmt * original, world_fmt * kopie, node * o,
                        node * t);
void swap_tree (world_fmt * tthis, world_fmt * tthat);
void swap (contribarr *a, contribarr *b);

void free_tree (node * p, world_fmt * world);
void free_tipnodelet (node * p, world_fmt * world);
void free_mignodelet (node * p, world_fmt * world);
void free_nodelet (node * p, long num, world_fmt * world);
void free_nodedata (node * p, world_fmt * world);

MYREAL inverse_logprob_noevent (world_fmt * world, long interval);
MYREAL sum_migprob (world_fmt * world, long pop, long interval);

MYREAL prob_micro_watkins (MYREAL t, long diff, world_fmt * world, mutationmodel_fmt *s, pair *helper);
MYREAL prob_micro_singlestep (MYREAL t, long diff, world_fmt * world, mutationmodel_fmt *s, pair *helper);

void make_sequences(long xs, mutationmodel_fmt *s, node *theNode, site_fmt **datapart);
void make_invarsites(long sublocus, mutationmodel_fmt *s, node *theNode);
long make_alleles(long xs, mutationmodel_fmt *s, node *theNode, option_fmt *options,
		  data_fmt *data, site_fmt *datapart2);
long make_brownian(long xs, mutationmodel_fmt *s, node *theNode, option_fmt *options, site_fmt *datapart);
long make_microsatellites (long sublocus, mutationmodel_fmt *s, node *theNode,
			   option_fmt *options, site_fmt *datapart2);
void debugtreeout (FILE * file, node * joint, node * p, long s);

void my_random_tree (world_fmt * world, long tips);
void print_sequences_stats(FILE *file, world_fmt * world, option_fmt * options, long locus);
void timeslices (timelist_fmt ** timevector);
void print_dead_tree(long numpop, timelist_fmt * timevector, world_fmt *world);
void pseudonu_allele (mutationmodel_fmt *s, proposal_fmt *proposal,
		      xarray_fmt *xxx1, MYREAL *lx1, MYREAL v1, xarray_fmt *xxx2, MYREAL *lx2, MYREAL v2, long xs);
node * get_random_sibling(node *origin, node **nodelist, long *nodenum, node *interior, long *lineages, double utime);
void print_lineages(long *l, long n);
void my_start_eventtime(double age, world_fmt *world, node **nodelist, long *lineages,long simtips, 
			double *shortt, long *shorti, char * shorte, long *to, long *from);
void zero_xseq(xarray_fmt *x, world_fmt *world);
void treeout_string (char ** file, long *filesize, long *pos, node * joint, node * p, long s);
MYREAL calc_pseudotreelength (proposal_fmt * proposal, MYREAL treelen);
void free_mignodelet (node * p, world_fmt * world);
void debugline(node *up);

void myfree_datapart(site_fmt ***datapart,  long numsites);

//##

/** global variable NUVIEW points to function nuview_datatype() */
//typedef int (*pt2Function)(float, char, char);
nuview_function * nuview;
//static void (*ppseudonuview)(proposal_fmt *, xarray_fmt , MYREAL *, MYREAL , xarray_fmt , MYREAL *, MYREAL);

prob_micro_function * prob_micro;
pseudonuview_function * pseudonuv;



/* ======================================================= */

void print_sequences_stats(FILE *file, world_fmt * world, option_fmt * options, long locus)
{
    long sublocus;
    if (strchr (SEQUENCETYPES, options->datatype))
      {
	const long sublocistart = world->sublocistarts[locus];
	const long sublociend   = world->sublocistarts[locus+1];
	for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
	  {
	    print_seqfreqs (file, world, options);
	    print_ratetbl (file, world, options, sublocus, 'A');
	    print_weights (file, world, options, sublocus);
	  }
      }
}

/*
 * Creates a start-genealogy using a coalescence approach
 *
 * - set NUVIEW according to datatype
 * initializes tree structure - fills tree with data - set_tree():
 * upgma-tree, adjust for times, sankoff() for migration events
 */
/// \brief builds the starting tree
///
/// Builds the starting tree. 
// \callgraph
void
buildtree (world_fmt * world, option_fmt * options, data_fmt * data,
           long locus)
{
    long pop;
    long genomes = number_genomes (options->datatype);

    if (world->data->skiploci[locus])
      return;

    free_tree(world->root, world);
    world->sumtips = 0;
    world->migration_counts = 0;
    for (pop = 0; pop < data->numpop; pop++)
      {
	if(options->randomsubset > 0 && options->randomsubset < data->numind[pop][locus])
	  world->sumtips += options->randomsubset * genomes; 
	else
	  world->sumtips += data->numalleles[pop][locus];
      }
    which_nuview (world, options->fastlike, FALSE, options->msat_option);

    switch (options->datatype)
      {
      case 's':
      case 'n':
      case 'h':
      case 'u':
      case 'f':
	break;
      case 'b':
	world->data->seq[0]->endsite = 1;
	data->freq = -10000000000000.;
	break;
      case 'm':
	world->data->seq[0]->endsite = 1;
	break;
      case 'a':
	world->data->seq[0]->endsite = 1;
	world->data->freq = 1. / (data->maxalleles[locus]);
	world->data->freqlast = 1. - world->data->freq;
      }

    if (options->usertree) // this seems not used by anyone at all, delete (20161030)
      {
        options->usertreewithmig = treereader (world, options, data, locus);
        makevalues (world, options, data, locus);//insert values into tips
      }
    else
      {
        allocate_tree (world, data, locus); //allocate nodep
        makevalues (world, options, data, locus); //allocate/insert data into tips
        allocateinterior(world,options, data,locus); //allocate nodep guts
        allocatepoproot (world); //allocate bottom parts
      }
    set_subloci_frequencies(world, options, data,locus);  
    init_tbl (world, locus);

    if (world->options->progress)
      {
	if (world->cold &&  world->replicate == 0)
	  {
#ifndef MPI
	    print_sequences_stats(stdout, world, options, locus);
	       //		print_seqfreqs (stdout, world, options);
	       //print_ratetbl (stdout, world, options, locus, 'A');
	       //print_weights (stdout, world, options, locus);
#endif
		if (world->options->writelog)
		  {
		    print_sequences_stats(world->options->logfile, world, options, locus);
		    //		    print_seqfreqs (world->options->logfile, world, options);
		    //print_ratetbl (world->options->logfile, world, options, locus,'A');
		    //print_weights (world->options->logfile, world, options,
		    //		   locus);
		  }
	      }
	  }
    if(!options->usertree)
      set_tree (world, options, data, locus);
    else
      {
        if(!options->usertreewithmig)
	  set_migrations (world, locus);
      }
    if (world->options->datatype == 'b')
      world->data->maxalleles[locus] = XBROWN_SIZE;
    
    // calculate the migration counts for the first tree
    world->migration_counts = 0;
    count_migrations (world->root->next->back, &world->migration_counts);

    if(simulator==TRUE)
      {
	treeout (stdout, crawlback (world->root->next),
		 crawlback (world->root->next), 0);      
	exit(0);
      }
}

/*
 * creates the timelist which represents all time intervals on a tree. The
 * timelist is an array of pointers and not a linked list.
 *
 * - allocates memory using an arbitrary value this will be later adjusted if a
 * longer list is needed - construct
 */
void
create_treetimelist (world_fmt * world,  timelist_fmt ** ltl)
{
    if ((*ltl)->tl == NULL)
    {
      (*ltl)->allocT =  TIMELIST_GUESS;
      (*ltl)->tl = (vtlist *) mycalloc ((*ltl)->allocT, sizeof (vtlist));
        allocate_lineages (ltl, 0, world->numpop);
    }
    (*ltl)->copies = 0;
    construct_tymelist (world, (*ltl));
    traverse_check(crawlback (world->root->next));
}


void
allocate_lineages (timelist_fmt **timevector, const long offset, const long numpop)
{
    long i;
    const long allocT = (*timevector)->allocT;
    vtlist * tl;
    long *lineages;

    // speedup and simplification for allocating and freeing
    // accessing is still done the old way through tl[i].lineages
    // but instead having memory scattered the allocation is one long string
    // freeing time for lineages should be reduced this way [I hope]
    if((*timevector)->lineages != NULL)
      {
	(*timevector)->lineages = (long *) myrealloc((*timevector)->lineages,sizeof(long)*(size_t)(numpop*allocT));
	if(offset != allocT)
	  {
	    memset((*timevector)->lineages+offset,0,sizeof(long)*(size_t) (numpop*(allocT-offset)));
	  }
      }
    else
      (*timevector)->lineages = (long *) mycalloc((numpop*allocT),sizeof(long));
   
    tl = (*timevector)->tl;
    lineages = (*timevector)->lineages;

    for (i = offset; i < allocT; i++)
      {
	tl[i].lineages = lineages + i * numpop;	
      }
}

/*
 * start first pass through the tree to calculate the tree-likleihood
 */
void
first_smooth (world_fmt * world, long locus)
{
    smooth (world->root->next, crawlback (world->root->next), world, locus);
}

/*
 * Marks a node, so that TREELIKELIHOOD() will recalulated values in node
 */
void
set_dirty (node * p)
{
    p->dirty = TRUE;
}
///
/// inserts the from and to into the timelist from the tree
void timeslices (timelist_fmt ** timevector)
{
    long z;
    vtlist * tl = (*timevector)->tl;
    vtlist * tlz;
    node * eventnode;

    long timeslice=0;

    const long T = (*timevector)->T;
 
    for (z = 0; z < T; z++)
    {
      tlz        = &(tl[z]);
      eventnode  = tlz->eventnode;
      if(eventnode!=NULL)
	{
	  tlz->from  = eventnode->pop;
	  tlz->to    = eventnode->actualpop;
	  tlz->timeslice = timeslice;
	}
      else
	{
	  if(z!=0)
	    {
	      tlz->from = tl[z-1].from;
	      tlz->to   = tlz[z-1].from;
	      timeslice++;
	      tlz->timeslice = timeslice;
	    }
	  else
	    {
	      error("confused with time slices");
	    }
	}
      tlz->slice = z;
    }
}


void print_dead_tree(long numpop, timelist_fmt * timevector, world_fmt *world)
{
  fprintf(stderr,"> %i ///////////////////\n",myID);
  long T = timevector->T-2;
  long ii, pop;
  for (ii = T; ii >= 0; ii--)
    {
      for (pop = 0; pop < numpop; pop++)
	fprintf(stderr,"%li ",timevector->tl[ii].lineages[pop]);
      fprintf(stderr,"%20.10f %20.10f | %f %c %li -> %li (%li)(%li)\n",timevector->tl[ii].age,
	      timevector->tl[ii].eventnode->tyme, timevector->tl[ii].eventnode->tyme,
	      timevector->tl[ii].eventnode->type, timevector->tl[ii].from,
	      timevector->tl[ii].to, timevector->tl[ii].eventnode->id,
	      showtop(timevector->tl[ii].eventnode->back)->id);
    }
  fprintf(stderr,"> %i /////////////////////////////////////////////\n", myID);
  treeout (stderr, crawlback (world->root->next),
	   crawlback (world->root->next), 0);      
  return;
}

boolean
add_partlineages426bottomup (long numpop, timelist_fmt ** timevector, world_fmt* world);
boolean
  add_partlineages426 (long numpop, timelist_fmt ** timevector, world_fmt* world);

void
add_partlineages (long numpop, timelist_fmt ** timevector, world_fmt* world)
{
  long k,pop;
  boolean dead = add_partlineages426bottomup (numpop,timevector,world);
  if (dead)
    {
      for(k=0;k<(*timevector)->T; k++)
	{
	  for(pop=0;pop<numpop;pop++)
	    {
	      if((*timevector)->tl[k].lineages[pop] <0)
		{
		  printf("%i> FAILED: pop=%li %li\n",myID, pop, (*timevector)->tl[k].lineages[pop]);
		}
	    }
	}
      error("extracting lineages from timelist failed");
    }
}

boolean
add_partlineages426bottomup (long numpop, timelist_fmt ** timevector, world_fmt* world)
{
  // this points to the MRCA
  long T;
  long i, pop;
  vtlist *tl = (*timevector)->tl;
  vtlist *tli;
  vtlist *tli1;
  long from;
  long to;
  long *lineages;
  char type;
  long tips=0;
  boolean dead=FALSE;
  T = (*timevector)->T;
  T = T - 1;
  memset(tl[0].lineages, 0, (size_t) numpop * sizeof(long));  
  memset(tl[1].lineages, 0, (size_t) numpop * sizeof(long));  
  //wrongtl[0].lineages[tl[0].eventnode->pop] = 1;
  for (i = 1; i < T; i++)
    {
      tli1 = &tl[i-1];
      memset(tl[i].lineages, 0, (size_t) numpop * sizeof(long));  
      tli = &tl[i];
      lineages = tli->lineages;
      for(pop=0;pop<numpop;pop++)
	lineages[pop] = tli1->lineages[pop];
      //memcpy(tli->lineages,tli1->lineages, (size_t) numpop * sizeof(long));
      from = tli1->from;
      to = tli1->to;
      if(tli1->eventnode != NULL)
      type = tli1->eventnode->type;
      else
	type = tli1->type; //'b'; // boundary

      switch(type)
	{
	case 't':
	  tips += 1;
	  lineages[to] += 1;
	  break;
	case 'i':
	  lineages[to] -= 1;
	  if (lineages[to] < 0)
	    {
	      print_dead_tree(numpop, *timevector, world);
	      dead=TRUE;
	      //error("lineages < 0 in addpartlineages");
	    }
	  break;
	case 'm':
	case 'd':
	  lineages[to] -= 1;
	  lineages[from] += 1;
	  if (lineages[to] < 0)
	    {
	      print_dead_tree(numpop, *timevector, world);
	      dead=TRUE;
	      //error("lineages < 0 in addpartlineages");
	    }
	  break;
	case 'r':
	  //case 'b':
	  warning("Root node in tymelist, this should not happen");
	default:
	  {
	    printf("%i> failed in constructing timelist: with type=%i\n",myID,type);
	    print_dead_tree(numpop, *timevector, world);
	    error("funny node received and died in add_partlineages");
	  }
	  //  break;
	}
    }
  //#ifdef DEBUGXX
  //fprintf(stderr,"%i> sumtips=%li,tips=%li\n",myID, world->sumtips, tips);
  //print_dead_tree(numpop, *timevector, world);
  //#endif
  return dead;
}


boolean
add_partlineages426 (long numpop, timelist_fmt ** timevector, world_fmt* world)
{
  // this points to the MRCA
  long T;
  long i, pop;
  vtlist *tl = (*timevector)->tl;
  vtlist *tli;
  vtlist *tli1;
  long from;
  long to;
  long *lineages;
  char type;
  boolean dead=FALSE;
  T = (*timevector)->T;
  T = T - 2;
  memset(tl[T+1].lineages, 0, (size_t) numpop * sizeof(long));  
  memset(tl[T].lineages, 0, (size_t) numpop * sizeof(long));  
  tl[T+1].lineages[tl[T].eventnode->pop] = 1;
  //tl[T].lineages[tl[T].eventnode->pop] = 2;
  for (i = T; i >= 0; i--)
    {
      tli1 = &tl[i+1];
      tli = &tl[i];
      lineages = tli->lineages;
      for(pop=0;pop<numpop;pop++)
	lineages[pop] = tli1->lineages[pop];
      //memcpy(tli->lineages,tli1->lineages, (size_t) numpop * sizeof(long));
      from = tli->from;
      to = tli->to;
      //if(tli->eventnode != NULL)
      type = tli->eventnode->type;
	//else
	//type = 'b'; // boundary

      switch(type)
	{
	case 't':
	  lineages[to] -= 1;
	  if (lineages[to] < 0)
	    {
	      print_dead_tree(numpop, *timevector, world);
	      dead=TRUE;
	      //error("lineages < 0 in addpartlineages");
	    }
	  break;
	case 'i':
	  lineages[to] += 1;
	  break;
	case 'm':
	case 'd':
	  lineages[to] += 1;
	  lineages[from] -= 1;
	  if (lineages[from] < 0)
	    {
	      print_dead_tree(numpop, *timevector, world);
	      dead=TRUE;
	      //error("lineages < 0 in addpartlineages");
	    }
	  break;
	case 'r':
	  //case 'b':
	  warning("Root node in tymelist, this should not happen");
	default:
	  {
	    printf("%i> failed in constructing timelist: with type=%i\n",myID,type);
	    print_dead_tree(numpop, *timevector, world);
	    error("funny node received and died in add_partlineages");
	  }
	  //  break;
	}
    }
  return dead;
}


 /*313
void add_partlineages (long numpop, timelist_fmt ** timevector, world_fmt *world)
{
  // this points to the MRCA
  const long T = (*timevector)->T - 2;

  long i, pop;
  vtlist *tl = (*timevector)->tl;
  vtlist *tli;
  vtlist *tli1;
  long from;
  long to;
  long *lineages;
  char type;
  // this should add a lineages for the MRCA
  memset(tl[T+1].lineages, 0, (size_t) numpop * sizeof(long));  
  tl[T+1].lineages[tl[T].eventnode->pop] = 1;
  for (i = T; i >= 0; i--)
    {
      tli1 = &tl[i+1];
      tli = &tl[i];
      lineages = tli->lineages;
      from = tli1->from;
      to = tli1->to;
      type = tli1->eventnode->type;
      memset(lineages, 0, (size_t) numpop * sizeof(long));  
      // if the node is an internode (a coalescent node) then add an additional line
      // if it is a tip reduce one
      if(type == 't')
	{
	  lineages[to] -= 1;
	  if (lineages[to] < 0)
	    error("313: lineages < 0 in addpartlineages");

	}
      else
	{
	  if(type != 'r')
	    {
	      if (from == to)
		{
		  lineages[to] += 1;
		}
	      else
		{
		  lineages[to] += 1;
		  lineages[from] -= 1;
		  if (lineages[to] < 0)
		    error("313 lineages < 0 in addpartlineages");

		}
	    }
	}
      // this copies the content from the last timeinterval (tli1) to the next (tli)
      //      printf("%i> lineages %5li:", myID, i);
      for (pop = 0; pop < numpop; pop++)
	{
	  lineages[pop] += tli1->lineages[pop];
	  //  printf(" %li",lineages[pop]);
	}
      //      printf(" %f %c %li-->%li %s\n",tli->eventnode->tyme, tli->eventnode->type, 
      //                 tli->from, tli->to, tli->eventnode->nayme);
    }
}
*/

/*
static void add_partlineages_312 (long numpop, timelist_fmt ** timevector)
{
  // this points to the MRCA
  const long T = (*timevector)->T - 2;

  long i, pop;
  vtlist *tl = (*timevector)->tl;
  vtlist *tli;
  vtlist *tli1;
  long from;
  long to;
  long *lineages;
  char type;
  // this should add a lineages for the MRCA
  memset(tl[T+1].lineages, 0, (size_t) numpop * sizeof(long));  
    tl[T].lineages[tl[T].eventnode->pop] = 1;// the tMRCA has then 1 lineages looking to the root
  tl[T+1].lineages[tl[T].eventnode->pop] = 0;// the tMRCA has then 1 lineages looking to the root
  for (i = T; i > 0; i--)
    {
      tli1 = &tl[i-1];
      tli = &tl[i];
      lineages = tli1->lineages;
      from = tli->from;
      to = tli->to;
      type = tli->eventnode->type;
      memset(lineages, 0, (size_t) numpop * sizeof(long));  
      // if the node is an internode (a coalescent node) then add an additional line
      // if it is a tip reduce one
      if(type == 't')
	{
	    lineages[to] -= 1;
	  if (lineages[to] < 0)
	    error("312 lineages < 0 in addpartlineages");

	}
      else
	{
	  if(type != 'r')
	    {
	      if (from == to)
		{
		  lineages[to] += 1;
		}
	      else
		{
		  lineages[to] += 1;
		  lineages[from] -= 1;
		  if (lineages[to] < 0)
		    error("312 lineages < 0 in addpartlineages");

		}
	    }
	}
      // this copies the content from the last timeinterval (tli1) to the next (tli)
      //printf("%i> lineages %5li:", myID, i);
      for (pop = 0; pop < numpop; pop++)
	{
	  lineages[pop] += tli->lineages[pop];
	  //printf(" %li",lineages[pop]);
	}
      //      printf(" %f %c %li-->%li %s lin=%li\n",tli->eventnode->tyme, tli->eventnode->type, tli->from, tli->to, tli->eventnode->nayme, lineages[0]);
    }
}
*/

/*
 calculates the tree-likelihood according to datatype a, m, b, s,u
 */
MYREAL
treelikelihood (world_fmt * world)
{
    long a;
    MYREAL aterm = 0.0;
    MYREAL term = 0.0;
    node *nn = crawlback (world->root->next);
#ifdef BEAGLE
    term = calcLnL(world, world->beagle->scalingIndex1);
    //    printf("Start LnL from Beagle: %f\n",term);
    return term;
#endif
    set_dirty (nn);
    smooth (world->root->next, crawlback(world->root->next), world, world->locus);
    adjustroot (world->root);

    mutationmodel_fmt *s;
    long sublocus;
    long locus = world->locus;
    const long sublocistart = world->sublocistarts[locus];
    const long sublociend   = world->sublocistarts[locus+1];
    term = 0.0;
    MYREAL localterm;
    for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
      {
	s = &world->mutationmodels[sublocus];
	const long xs = sublocus - sublocistart;
	switch (s->datatype)
	  {
	  case 's':
	    localterm = treelike_seq (s, xs, world, world->locus);
	    //printf("localterm=%f ",localterm);fflush(stdout);
	    term += localterm;
	    //printf("term=%f\n",term);fflush(stdout);
	    break;
	  case 'n':
	  case 'h':
	    term += treelike_snp (s, xs, world, world->locus);
	    break;
	  case 'u':
	    term += treelike_snp_unlinked (s, xs, world, world->locus);
	    break;
	  case 'a':
	    aterm = 0.0;
	    for (a = 0; a < s->maxalleles - 1; a++)
	      {
		aterm += (s->freq * nn->x[xs].a[a]);
	      }
	    aterm += (s->freqlast * nn->x[xs].a[a]);
	    term += (aterm != 0.0) ? (LOG (aterm) + nn->scale[xs][0]) : -MYREAL_MAX;
	    break;
	  case 'm':
	    aterm = 0.0;
	    for (a = 0; a < s->maxalleles; a++)
	      aterm += nn->x[xs].a[a];
	    term += (aterm != 0.0) ? (LOG (aterm) + nn->scale[xs][0]) : -MYREAL_MAX;
	    break;
	  case 'b':
	    term += nn->x[xs].a[2];
	    break;
	  case 'f':
	    term += treelike_anc (s, xs, world, world->locus);
	    break;
	  }
      }
#ifdef UEP
    if (world->options->uep)
      ueplikelihood (world);
#endif
#ifdef DEBUG
    //printf("@L(D=%li|G)=%f\n",world->locus,term);
#endif
#ifdef TREEDEBUG
    fprintf(stdout,"%i> [temperature:%f] ",myID, world->heat);
    int i;
    for(i=0;i<world->numpop2;i++)
      fprintf(stdout,"%f ",world->param0[i]);
    debugtreeout (stdout, crawlback (world->root->next), crawlback (world->root->next),0);
#endif

    return term;
}

/*
 * calculates tree-likelihood using only arrays DOES NOT CHANGE ARRAYS IN THE
 * TREE
 */
MYREAL
pseudotreelikelihood (world_fmt * world, proposal_fmt * proposal)
{
  long a, locus = world->locus;
  /* freq is not different between pop */
  MYREAL aterm = 0.0;
  MYREAL term = 0.0;
  //  printf("@");//DEBUG MPI 2010
#ifdef BEAGLE
      //      printf("PseudoLike: ");
      return calcLnL(world, world->beagle->scalingIndex2);
#endif
      mutationmodel_fmt *s;
      long sublocus;
      const long sublocistart = world->sublocistarts[locus];
      const long sublociend   = world->sublocistarts[locus+1];
      for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
	{
	  s = &world->mutationmodels[sublocus];
	  const long xs = sublocus - sublocistart;
	  xarray_fmt xxf = proposal->xf[xs];
	  xarray_fmt xxt = proposal->xt[xs];
	  switch (s->datatype)
	    {
	    case 's':
	      term += pseudo_tl_seq (s, xs, xxf.s, xxt.s, proposal->v,
				    proposal->vs, proposal, world);
	      
	      break;
	    case 'n':
	    case 'h':
	      term += pseudo_tl_snp (s, xs, xxf.s, xxt.s, proposal->v,
				    proposal->vs, proposal, world);
	      break;
	    case 'u':
	      term += pseudo_tl_snp_unlinked (s, xs, xxf.s, xxt.s,
					     proposal->v, proposal->vs, proposal,
					     world);
	      break;
	    case 'a':
	      aterm = 0.0;
	      for (a = 0; a < s->maxalleles - 1; a++)
		{
		  aterm += (s->freq * xxf.a[a]);
		  //printf("%f ",term);
		}
	      aterm += (s->freqlast * xxf.a[a]);
	      //printf("%f ",term);
	      if (aterm == 0.0)
		term = -MYREAL_MAX;
	      else
		term = (LOG (aterm) + proposal->mf[xs][0]);
	      //			printf("\n@@pseudolike=%f scale=%f\n",term, proposal->mf[0]);
	      break;
	    case 'b':
	      term += xxf.a[2];
	      break;
	    case 'm':
	      aterm = 0.0;
	      for (a = 0; a < s->maxalleles; a++)
		{
		  aterm += xxf.a[a];
		}
	      if (aterm == 0.0)
		term = -MYREAL_MAX;
	      else
		term = (LOG (aterm) + proposal->mf[xs][0]);
	      break;
	    case 'f':
	      term = pseudo_tl_anc (s, xxf.s, xxt.s, proposal->v,
				    proposal->vs, proposal, world);
	      break;
	    default:
	      warning("no datatype found for treelikelihood() calculation!");
	      term = -MYREAL_MAX;
	      break;
	    }
	}
#ifdef UEP
  if (proposal->world->options->uep)
    term += pseudo_tl_uep (&(proposal->uf), &proposal->ut, proposal->v,
			   proposal->vs, proposal, world);
#endif
  if(MYFINITE(((double) term)) == 0)
    term = -MYREAL_MAX;
#ifdef DEBUG
  //    printf("   proposed L(D|G)=%f\n",term);
#endif
  return term;
}


/*
 * Calculates the sub-likelihoods but does not change the arrays in the tree,
 * it uses the passed arrays and overwrites the xx1 array DOES NOT CHANGE THE
 * TREE
 */

//void pseudonuview (proposal_fmt * proposal, xarray_fmt *axx1, MYREAL **alx1, MYREAL v1,
//              xarray_fmt **axx2, MYREAL **alx2, MYREAL v2)
void pseudonuview (proposal_fmt * proposal, xarray_fmt *axx1, MYREAL **alx1, MYREAL v1, 
		   xarray_fmt *axx2, MYREAL **alx2, MYREAL v2)

{
  world_fmt *world  = proposal->world;
  const long locus        = world->locus;
  mutationmodel_fmt *s;
  long sublocus;
  const long sublocistart = world->sublocistarts[locus];
  const long sublociend   = world->sublocistarts[locus+1];
  MYREAL *lx1;
  MYREAL *lx2;
  for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
    {
      s = &world->mutationmodels[sublocus];
      const long xs = sublocus - sublocistart;
      //xarray_fmt  xx1   = axx1[xs];
      //xarray_fmt  xx2   = axx2[xs];
      lx1   = alx1[xs];
      lx2   = alx2[xs];
      (*pseudonuv[xs]) (s, proposal, axx1, lx1, v1, axx2, lx2, v2, xs);
      /*switch (s->datatype)
	{
	case 'a':
	  pseudonu_allele (s, proposal, &xx1.a, &(lx1[0]), v1, xx2.a, lx2[0], v2);
	  break;
	case 'b':
	  pseudonu_brownian (s, proposal, &xx1.a, lx1, v1, xx2.a, lx2[0], v2);
	  break;
	case 'm':
	  pseudonu_micro (s, proposal, &xx1.a, lx1, v1, xx2.a, lx2[0], v2, xn);
	  break;
	case 'u':
	case 'n':
	case 'h':
	  pseudonu_seq (s, proposal, xx1.s, v1, xx2.s, v2);
	  break;
	case 's':
	  if (proposal->world->options->fastlike)
	    pseudonu_seq (s, proposal, xx1.s, v1, xx2.s, v2);
	  else
	    pseudonu_seq_slow (s, proposal, xx1.s, lx1, v1, xx2.s, lx2, v2);
	  break;
	case 'f':
	  pseudonu_anc (s, proposal, xx1.s, v1, xx2.s, v2);
	  break;
	  }*/
    }
#ifdef  UEP
  if (proposal->world->options->uep)
    {
      pseudonu_twostate (proposal, &proposal->uf, proposal->umf, v1,
			 &proposal->ut, proposal->umt, v2);
    }
#endif
}

void pseudonu_allele (mutationmodel_fmt *s, proposal_fmt *proposal, xarray_fmt *xxx1, MYREAL *lx1, MYREAL v1, xarray_fmt *xxx2, MYREAL *lx2, MYREAL v2, long xs)
//(mutationmodel_fmt *s, proposal_fmt * proposal, MYREAL **xx1, MYREAL *lx1,
//                 MYREAL vv1, MYREAL *xx2, MYREAL lx2, MYREAL vv2)
{
  (void) proposal;
    long   a;
    long   aa;
    long   mal   = s->maxalleles; /* maxalleles */
    MYREAL freq  = s->freq;
    //MYREAL freqlast = s->freqlast; // same as 1 - k/(k+1) = 1/(k+1) = freq
    MYREAL w1 = 0.0; /* time variables */
    MYREAL w2 = 0.0;
    MYREAL vv1;
    MYREAL vv2;
    MYREAL v1freq;
    MYREAL v2freq;
    MYREAL pija1;/* summary of probabilities */
    MYREAL pija2;  
    MYREAL x3m = -MYREAL_MAX;
    MYREAL inv_x3m;
    MYREAL *xx3;
    MYREAL *xx1 = xxx1[xs].a;
    MYREAL *xx2 = xxx2[xs].a;
    xx3 = (MYREAL *) mymalloc(sizeof(MYREAL) * (size_t) mal);
    vv1 = 1.0 - EXP (-v1);
    vv2 = 1.0 - EXP (-v2);
    if (vv1 >= 1.)
    {
        w1 = 0.0;
        vv1 = 1.0;
    }
    else
    {
        w1 = 1.0 - vv1;
    }
    if (vv2 >= 1.)
    {
        w2 = 0.0;
        vv2 = 1.0;
    }
    else
    {
        w2 = 1.0 - vv2;
    }
    //    printf("@@pseudonu 1:{%f,%f,%f,%f}   2:(%f,%f,%f,%f}\n",(*xx1)[0],(*xx1)[1],*lx1,v1,
    //	   xx2[0], xx2[1],lx2,v2);
  
    v1freq = vv1 * freq;
    v2freq = vv2 * freq;

    for (aa = 0; aa < mal; aa++)
    {
        pija1 = pija2 = 0.0;
        for (a = 0; a < mal; a++)
        {
	  if(aa==a)
	    {
	      pija1 += (w1 + v1freq) * xx1[a];
	      pija2 += (w2 + v2freq) * xx2[a];
	    }
	  else
	    {
	      pija1 += v1freq * xx1[a];
	      pija2 += v2freq * xx2[a];
	    }
	}
        xx3[aa] = pija1 * pija2;
        if (xx3[aa] > x3m)
            x3m = xx3[aa];
    }
    inv_x3m = 1./ x3m;
    for (aa = 0; aa < mal; aa++)
    {
        xx1[aa] = xx3[aa] * inv_x3m;
    }
    //    printf("@@finish: %f %f\n",(*xx1)[0],(*xx1)[1]);
    lx1[0] = LOG (x3m) + lx2[0] + lx1[0];
    myfree(xx3);
}

void
pseudonu_micro (mutationmodel_fmt *s, proposal_fmt *proposal, xarray_fmt *xxx1,
		MYREAL *lx1, MYREAL v1, xarray_fmt *xxx2, MYREAL *lx2, MYREAL v2, long xs)
{
  long a, ss, diff;
  long aa1, aa2;
  long smax = s->maxalleles;
  long margin = s->micro_threshold;
  MYREAL pija1s, pija2s, vv1, vv2;
  MYREAL x3m = -MYREAL_MAX;
  world_fmt *world = proposal->world;
  MYREAL inv_x3m;
  MYREAL *xx3;
  MYREAL *pm1;
  MYREAL *pm2;
  pair *helper = &world->options->msat_tuning;

  MYREAL * xx1 = xxx1[xs].a;
  MYREAL * xx2 = xxx2[xs].a;
  xx3 = (MYREAL *) mymalloc(sizeof(MYREAL) * (size_t) smax);
  vv1 = v1;
  vv2 = v2;
  pm1 = (MYREAL *) mymalloc(sizeof(MYREAL) * (size_t) (2 * margin));
    pm2 = pm1 +  margin;

    for (diff = 0; diff < margin; diff++)
      {
	pm1[diff] = (*prob_micro[xs]) (vv1, diff, world, s, helper);
	pm2[diff] = (*prob_micro[xs]) (vv2, diff, world, s, helper);
	//printf("%li: pm1=%10.10f pm2=%10.10f\n", diff, pm1[diff], pm2[diff]);
      }
    
    for (ss = 0; ss < smax; ss++)
    {
        pija1s = pija2s = 0.0;
	aa1 = MAX (0, ss - margin);
	aa2 = MIN(ss + margin,smax);
	for (a = aa1; a < aa2; a++)
	  //for (a = 0; a < smax; a++)
        {
            diff = labs (ss - a);
	    if(diff >= margin)
	      continue;

            if (xx1[a] > 0)
            {
                pija1s += pm1[diff] * xx1[a];
            }
            if (xx2[a] > 0)
            {
                pija2s += pm2[diff] * xx2[a];
            }
        }
	//	printf("%li: pija1s=%10.10f pija2s=%10.10f\n", ss, pija1s,pija2s);
        xx3[ss] = pija1s * pija2s;
        if (xx3[ss] > x3m)
            x3m = xx3[ss];
    }
    inv_x3m = 1./ x3m;
    for (ss = 0; ss < smax; ss++)
    {
        xx1[ss]= xx3[ss] * inv_x3m;
    }
    lx1[0] += LOG (x3m) + lx2[0];
    myfree(xx3);
    myfree(pm1);
}

//================================================
// brownian motion calculation for data likelihood: ghost-version
// changed divisions to multiply inverse
// simplified check and include of boundary for vtot
void pseudonu_brownian (mutationmodel_fmt *s, proposal_fmt *proposal, xarray_fmt *xxx1, MYREAL *lx1, MYREAL v1, xarray_fmt *xxx2, MYREAL *lx2, MYREAL v2, long xs)
//(mutationmodel_fmt * s, proposal_fmt * proposal, MYREAL **xx1, MYREAL *lx1,
//                   MYREAL v1, MYREAL *xx2, MYREAL lx2, MYREAL v2)
{
  (void) s;
  (void) proposal;
  (void) lx1;
  (void) lx2; 
    MYREAL vtot, rvtot, c12;
    MYREAL mean1, mean2, mean, vv1, vv2, f1, f2, diff;
    MYREAL *x1 =  xxx1[xs].a;
    MYREAL *x2 =  xxx2[xs].a;
    mean1 = x1[0];
    mean2 = x2[0];
	
    vv1 = v1 + x1[1];
    vv2 = v2 + x2[1];
    vtot = vv1 + vv2;
    if (vtot > 0.0)
    {
        rvtot = 1./vtot;
        f1 = vv2 * rvtot;
        f2 = 1.0 - f1;
        mean = f1 * mean1 + f2 * mean2;
        diff = mean1 - mean2;
        c12 = diff * diff * rvtot;
        x1[2] = x1[2] + x2[2] + MIN (0, -0.5 * (LOG (vtot) + c12) + LOG2PIHALF);
        x1[1] = vv1 * f1;
        x1[0] = mean;
    }
    else
    {
        //xcode rvtot = HUGE;
        //xcode vtot = SMALL_VALUE;
        f1 = 0.5;
        x1[2] = (double) HUGE;
        x1[1] = vv1 * f1;
        x1[0] = f1 * (mean1 + mean2);
    }
}

/*
 adjust the variables POP and ACTUALPOP in interior nodes
 */
void
set_pop (node * theNode, long pop, long actualpop)
{
    switch (theNode->type)
    {
    case 'm':
    case 'd':
      theNode->pop = theNode->next->pop = pop;
      theNode->actualpop = theNode->next->actualpop = actualpop;
      break;
    case 'i':
    case 'r':
      theNode->pop = theNode->next->pop = theNode->next->next->pop = pop;
      theNode->actualpop = theNode->next->actualpop = actualpop;
      theNode->next->next->actualpop = actualpop;
      break;
    case 't':
      if (theNode->pop != pop)
	error ("Population designation scrambled");
      break;
    default:
      error ("Undefined node type?!");
      //break;
    }
}



/*
 * ======================================================= local functions
 */

void
allocate_tree(world_fmt * world, data_fmt * data,
              long locus)
{
    long nodenum = 0, pop, numpop = data->numpop;
    for (pop = 0; pop < numpop; pop++)
    {
        nodenum += data->numalleles[pop][locus] * 2;
    }
    //    world->nodep = (node **) mycalloc (nodenum+1, sizeof (node *));
    world->nodep = (node **) mycalloc (world->sumtips * 2, sizeof (node *));
}


void allocateinterior (world_fmt * world, option_fmt *options, data_fmt * data, long locus)
{
  (void) options;
    node *p;
    long i;
    long temp=0;
    long mini = 0, maxi = 0;
    long numpop = data->numpop;
    //long genomes = number_genomes (options->datatype);
    for (i = 0; i < numpop; i++)
    {
      //temp += genomes * max_shuffled_individuals(options, data, i, locus);
      temp += data->numalleles[i][locus];
    }
    mini = temp;
    maxi = temp * 2 - 1; 
    for (i = mini; i < maxi; i++)
    {
      p = allocate_nodelet (world, 3, 'i');
      p->top = TRUE;
      p->s = (MYREAL *) mycalloc (1, sizeof (MYREAL) * (size_t) world->numpop);
      p->id = i;
      world->nodep[i] = p;
    }
}

node *
allocate_nodelet (world_fmt *world, long num, char type)
{
  boolean isfirst = TRUE;
  long j;
  node *p, *q = NULL, *pfirst = NULL;
  for (j = 0; j < num; j++)
    {
#ifdef DISPENSER
      p = dispense_nodelet(world);
#else
      p = (node *) mycalloc(1,sizeof(node)); 
#endif
      if(p==NULL)
	{
	  warning( "%i> heat=%f Nodelet dispenser failed to supply nodelet\n",myID,1./world->heat);
	  p = (node *) mymalloc (sizeof (node));
	}
      p->tip = FALSE;
      p->visited = FALSE;
      p->number = unique_id_global;
      p->pop = -1;
      p->actualpop = -1;
      p->type = type;
      p->id = -1;
      p->top = FALSE;
      p->dirty = TRUE;
      p->next = q;
      p->scale = NULL;
      p->s = NULL;
      p->x = NULL;
      //p->x.a = NULL;
#ifdef UEP
      
      p->uep = NULL;
      p->ux.s = NULL;
      p->ux.a = NULL;
#endif
      p->back = NULL;
      p->nayme = NULL;
      p->truename = NULL;
      p->v = 0.0;
      p->tyme = 0.0;
      p->length = 0.0;
      if (isfirst)
        {
	  isfirst = FALSE;
	  pfirst = p;
        }
      q = p;
    }
    if(pfirst !=NULL)
	  pfirst->next = q;
  return q;
}

void
allocatepoproot (world_fmt * world)
{
  node *p, *q;
  q = NULL;
  p = allocate_nodelet (world, 3, 'r');
  p->top = TRUE;
  world->root = p;
}


void
allocate_tip (world_fmt * world, option_fmt * options, node ** p, long pop,
              long locus, long a, long ind, char **tipnames)
{
  (void) ind;
  long i;
  if(options->usertree)
    return; //do nothing because the tip is alread allocated
  char *tipname;
  unsigned long tiplen;
  unpad(tipnames[locus]," ");
  switch (tipnames[locus][0])
    {
    case '\0':
    case ' ':
      tipname = tipnames[0];
      break;
    default:
      tipname = tipnames[locus];
    }
#ifdef DEBUG  
  if ((*p) != NULL)
    warning("%i> node was not empty: %s\n",myID, (*p)->nayme);
  else
#else
    if((*p) == NULL)
#endif
    {
      (*p) = allocate_nodelet (world, 1, 't');
    }
  if((*p)->s == NULL)
    (*p)->s = (MYREAL *) mycalloc (1, sizeof (MYREAL) * (size_t) world->numpop);
  else
    (*p)->s = (MYREAL *) myrealloc ((*p)->s, sizeof (MYREAL) * (size_t) world->numpop);
  tiplen = 15 + strlen(tipname);
  if((*p)->nayme == NULL)
    (*p)->nayme = (char *) mycalloc((2*tiplen), sizeof(char));
  else
    (*p)->nayme = (char *) myrealloc((*p)->nayme, (size_t) (2*tiplen) * sizeof(char));

  (*p)->tip = TRUE;
  (*p)->top = TRUE;
  (*p)->id  = a;
  if(options->has_datefile)
    {
      (*p)->tyme = find_tipdate(tipname, pop, world);
    }
  else
    (*p)->tyme = 0.0;

  translate(tipname,' ', '_');
  unpad(tipname,"_"); 
  (*p)->pop = (*p)->actualpop = options->newpops[pop]-1;
  (*p)->truepop = (*p)->pop;
  for (i = 0; i < world->numpop; i++)
    (*p)->s[i] = MYREAL_MAX;
  (*p)->s[options->newpops[pop]-1] = 0;
  strncpy((*p)->nayme,tipname,tiplen);
  (*p)->truename = (*p)->nayme + tiplen;
  strncpy((*p)->truename,tipname,tiplen);

#ifdef DEBUG
      //      printf("%i> allocate_tip: type=%c truename=%s\n",myID, (*p)->type,(*p)->truename);
#endif
  set_unassigned(*p, world);
  alloc_seqx (world, (*p), locus);
#ifdef UEP
  if (world->options->uep)
    {
      if((*p)->uep == NULL)
	{
	  (*p)->uep = (int *) mycalloc (world->data->uepsites, sizeof (int));
	  (*p)->ux.s = (pair *) mycalloc (world->data->uepsites, sizeof (pair));
	}
    }
#endif
}

void myfree_datapart(site_fmt ***datapart, long numsites)
{
  long site;
  if (*datapart != NULL)
    {
      for(site=0;site<numsites;site++)
	{
	  myfree((*datapart)[0][site]);
	}
      myfree((*datapart)[0]);
      if ((*datapart)[1]!=NULL)
	{
	  for(site=0;site<numsites;site++)
	    {
	      myfree((*datapart)[1][site]);
	    }
	  myfree((*datapart)[1]);
	}
      myfree(*datapart);
      *datapart = NULL;
    }
}

void makevalues(world_fmt *world, option_fmt *options, data_fmt *data, long locus)
{
  long   zpop=0;
  long   pop;
  long   top;
  long   ind;
  long   z=0;
  long   ii;
  long invi;	
  node **treenode = world->nodep;
  long   sublocus;
  long   sublocistart = world->sublocistarts[locus];
  long   sublociend = world->sublocistarts[locus+1];
  mutationmodel_fmt *s;
  boolean old = (!data->oneliner);// && !strchr(SEQUENCETYPES, world->mutationmodels[sublocistart].datatype));
  long zminus=0;
  long oldsite;
  site_fmt **datapart = NULL;  // MPI data on demand will use this, the other compile version treat this as a
  // a reference to the data->yy section
  //empty_world_unassigned(world);
  for (pop = 0; pop < data->numpop; pop++)
    {
      top = max_shuffled_individuals(options, data, pop, locus);
      for (ii = 0; ii < top; ii++)
	{
	  ind = data->shuffled[pop][locus][ii];
	  
	  if (!options->usertree)
	    {
#ifndef MPIDATAONDEMAND
	      if(treenode[z]==NULL)
		{
		  allocate_tip (world, options, &treenode[z], pop, locus, z,
				ind, data->indnames[pop][ind]);
		  world->data->sampledates[pop][locus][ind].id = treenode[z]->id;
		}
	      if(options->haplotyping)
		link_individual_node(data->indnames[pop][ind][locus],treenode[z], locus, world);
#endif
	      for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
		{
		  s = &world->mutationmodels[sublocus];
		  if (s->datatype==0)
		    s->datatype = world->options->datatype;
		  const long xs = sublocus - sublocistart;
#ifdef DEBUG
		  fprintf(stderr,"%i>locus=%li (%li),substart=%li,substop=%li\n",myID,locus,sublocus,sublocistart,sublociend);
#endif		  
#ifdef MPI
#ifdef MPIDATAONDEMAND
		  if (myID != MASTER)
		    {
		      // the heated chains and replicates all need to 
		      // request the data and storing it in the data->datapart
		      // does not make sense because only one record is stored
		      // leading to problems with the heated chains,
		      // it is unclear how this works with threading because
		      // the heated chains are thread using Grand Central 
		      // or pthreads and thus may have access collisions on 
		      // systems that need semaphores or that may have issues
		      // with interaction among master and worker node.
		      //if(world->cold)
		      //	{			 
		      if (!strchr (SEQUENCETYPES, s->datatype))
			{
			  request_data(pop,ind,locus,sublocus,1, world, data, options, &datapart);
			}		
		      else
			{
			  request_data(pop,ind,locus,sublocus,0, world, data, options, &datapart);
			}
		      if(treenode[z]==NULL)
			{
			  allocate_tip (world, options, &treenode[z], pop, locus, z,
					ind, data->indnames[pop][ind]);
			  world->data->sampledates[pop][locus][ind].id = treenode[z]->id;
			}
		      if(options->haplotyping)
			link_individual_node(data->indnames[pop][ind][locus],treenode[z], locus, world);		      
		    }
		  else
		    {
		      error("Master is not supposed to usemakevalues()");
		    }
#else
		  datapart = data->yy[pop][ind][sublocus];
#endif /*MPIONDEMAND*/
#else
		  datapart = data->yy[pop][ind][sublocus];
#endif /*MPI*/
		  switch(s->datatype)
		    {
		      // sequences
		    case 's':
		    case 'f':
		      make_sequences(xs, s,treenode[z], datapart);
		      //z--;
		      break;
		      // SNPs
		    case 'h':
		    case 'n':
		      make_sequences(xs, s,treenode[z], datapart);
		      make_invarsites(xs, s,treenode[z]);
		      //z--;
		      break;
		      // allelic
		    case 'a':
		      zminus = make_alleles(xs, s,treenode[z], options, data, datapart[0]);
		      if(old)
			{
			  if(zminus == 1)
			    {
			      if(world->has_unassigned)
				remove_node_assigndb(world,treenode[z]);
			    }
			  else
			    z++;
			  allocate_tip (world, options, &treenode[z], pop, locus, z,
					ind, data->indnames[pop][ind]);
			  zminus = make_alleles(xs, s,treenode[z], options, data, datapart[1]);
			}
		      break;
		      // brownian
		    case 'b':
		      zminus = make_brownian(xs, s,treenode[z], options, datapart[0]);
		      if(old)
			{
			  if(zminus == 1)
			    {
			      if(world->has_unassigned)
				remove_node_assigndb(world,treenode[z]);
			    }
			  else
			    z++;
			  allocate_tip (world, options, &treenode[z], pop, locus, z,
					ind, data->indnames[pop][ind]);
			  zminus = make_brownian(xs, s,treenode[z], options, datapart[1]);
			}
		      break;
		      // microsatellites
		    case 'm':
		      zminus = make_microsatellites(xs, s,treenode[z], options, datapart[0]);
		      if(old)
			{
			  if(zminus == 1)
			    {
			      if(world->has_unassigned)
				remove_node_assigndb(world,treenode[z]);
			    }
			  else
			    z++;
			  allocate_tip (world, options, &treenode[z], pop, locus, z,
					ind, data->indnames[pop][ind]);
			  zminus = make_microsatellites(xs, s,treenode[z], options, datapart[1]);
			}
		      break;
		    default:
		      error ("Wrong datatype");
		    }
#ifdef MPIDATAONDEMAND
		  // is there a better solution, there were problems with reallocations
		  // this make sure that datapart=NULL
		  myfree_datapart(&datapart,s->numsites);
#ifdef DEBUG
		  printf("%i> datapart freed\n",myID);
#endif
#endif
		}
	      if(zminus == 1)
		{
		  if(world->has_unassigned)
		    remove_node_assigndb(world,treenode[z]);
		}
	      else
		z++;
	    }
	}
      data->numalleles[pop][locus]=z-zpop;
      zpop = z;
    }
  // snp and hapmap stuff:
  for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
    {
      s = &world->mutationmodels[sublocus];
      if (strchr(SNPTYPES,s->datatype))
	{
	  oldsite = s->numpatterns;
	  invi = s->numpatterns; //seq->endsite;
	  //s->numpatterns += 4; //seq->endsite += 4;
	  if (options->totalsites>0)
	    {
	      long invsites = options->totalsites - oldsite;
	      s->aliasweight[invi++] = invsites*s->basefreqs[NUC_A];
	      s->aliasweight[invi++] = invsites*s->basefreqs[NUC_C];
	      s->aliasweight[invi++] = invsites*s->basefreqs[NUC_G];
	      s->aliasweight[invi]   = invsites*s->basefreqs[NUC_T];
	      world->options->datatype = 's';
	    }
	}
    }
  world->sumtips = z;
  if(world->sumtips==0)
    {
      data->skiploci[locus] = TRUE;
      world->data->skiploci[locus] = TRUE;
      world->skipped += 1;
    }
#ifdef UEP
  if (world->options->uep)
    {
      make_uep_values (world, data, locus);
    }
#endif
  //  set_subloci_frequencies(world, options, data);
  //#ifdef MPIDATAONDEMANDXXx
  //if (!strchr (SEQUENCETYPES, s->datatype))
  //{
  //  myfree(datapart2[1]);
  //  myfree(datapart2[0]);
  //  myfree(datapart2);
  //}
  //else
  //{
  //  myfree(datapart[0]);
  //  myfree(datapart);
  //}
  //#endif
}



///
///
/// creates the branchlength and adds migration nodes to the
/// start tree
/// - creates rough genetic distance for upgma
/// - upgma
/// - find branches where we need to insert migrations
/// - insert migrations
/// - adjust time of all nodes using the coalescent with migration
void
set_tree (world_fmt * world, option_fmt * options, data_fmt * data,
          long locus)
{
  (void) options;
  (void) data;
  long     tips = world->sumtips;
  //MYREAL   **distm;
  //node     **topnodes;
  //if(options->randomtree)
  //  {
  //warning("Current Migrate 4.0 works only with random starting trees, may fail to construct them\n");
  //warning("if program crashes immediately try to restart a few time\n");
      my_random_tree(world,tips);
      world->root->tyme = world->root->next->tyme =
	world->root->next->next->tyme = world->root->next->back->tyme + 10000.;
      set_top (world, world->root->next->back, locus);
      set_v (world->root->next->back);
      allocate_x (world->root, world, locus, NOTIPS);

#ifdef TREEDEBUG
      printf("TREE DEBUG START RANDOM TREE\n");
      treeout (stdout, crawlback (world->root->next),
	       crawlback (world->root->next), 0);      
#endif

      return;
      //  }
#if 0 /* this stuff needs to go because it does not work with divergence*/
      topnodes = (node **) mycalloc (1, sizeof (node *) * tips);
  doublevec2d(&distm,tips, tips);
  if (!options->randomtree)
    {
      // create a crude distance matrix according to the datatype 
      switch (world->options->datatype)
	{
	case 'a':
	case 'b':
	case 'm':
	  distance_allele (world, options, tips, distm);
	  break;
	case 's':
	case 'n':
	case 'h':
	case 'u':
	case 'f':
	  distance_sequence (data, options, locus, tips,
			     world->data->seq[0]->sites[locus],
			     options->nmlength, distm);
	  break;
	}
    }
  else
    {
      randomize_distm (distm,  tips, world->options->custm);
      //printf("\n\nRANDOM TREE START\n\n");
    }
  constrain_distance_zeromig (distm, options, data, locus, tips,
			      world->options->custm);
  //    printf("distm[0][1]=%f distm[1][2]=%f\n",distm[0][1],distm[1][2]);
  //fflush(stdout);
  
#ifdef UEP
  if (options->uep)
    constrain_distance_uep (data->uep, world->data->uepsites, distm,
			    tips);
#endif
  if(!options->usertree)
    upgma (world, distm, tips, world->nodep);
  myfree(distm[0]);
  myfree(distm);
  //}
  world->root->tyme = world->root->next->tyme =
    world->root->next->next->tyme = world->root->next->back->tyme + 10000.;
  /* orient the tree up-down, set the length and v */
  set_top (world, world->root->next->back, locus);
  set_v (world->root->next->back);
#ifdef BEAGLE
  long bid=0;
  bid = set_branch_index (world->root->next->back, &bid);
#endif
  /*
   * insert migration nodes into the tree using the Slatkin and
   * Maddison approach (Fitch parsimony)
   */
  memcpy (topnodes, world->nodep, sizeof (node *) * tips);
  //zzz = 0;
  if(!options->usertree) //TRIAL
    allocate_x (world->root, world, locus, NOTIPS);
#ifdef UEP
  
  if (world->options->uep)
    {
      //      allocate_uep (world->root, world, world->options->datatype, NOTIPS);
      update_uep (world->root->next->back, world);
      check_uep_root (world->root->next->back, world);
    }
#endif
  //    debugtreeout (stdout, crawlback (world->root->next), crawlback (world->root->next),0);
  sankoff (world);
  //    debugtreeout (stdout, crawlback (world->root->next), crawlback (world->root->next),0);
  myfree(topnodes);
#endif /*this stuff needs to go */
}    /* set_tree */


/*
 creates the branchlength and adds migration nodes to the
 start tree
 - creates rough genetic distance for upgma
 - upgma
 - find branches where we need to insert migrations
 - insert migrations
 - adjust time of all nodes using the coalescent with migration
 */
void
set_migrations (world_fmt * world, long locus)
{
    long tips = world->sumtips;
	
    node **topnodes;
    topnodes = (node **) mycalloc (tips, sizeof (node *));
	// can we do this here ???#ifdef UEP
	//            if (options->uep)
	//            constrain_distance_uep (data->uep, world->data->uepsites, distm,  tips);
	//#endif
	world->root->tyme = world->root->next->tyme =
		world->root->next->next->tyme = world->root->next->back->tyme + 10000.;
	/* orient the tree up-down, set the length and v */
	set_top (world, world->root->next->back, locus);
	set_v (world->root->next->back);
	/*
	 * insert migration nodes into the tree using the Slatkin and
	 * Maddison approach (Fitch parsimony)
	 */
	memcpy (topnodes, world->nodep, sizeof (node *) * (size_t) tips);
#ifdef UEP
	if (world->options->uep)
	{
		//      allocate_uep (world->root, world, world->options->datatype, NOTIPS);
		update_uep (world->root->next->back, world);
		check_uep_root (world->root->next->back, world);
	}
#endif
	sankoff (world);
	myfree(topnodes);
}    /* set_migrations */


void randomize_distm(MYREAL **distm, long tips, char *custm)
{
  (void) custm;
  long i;
  long j;
  for(i=0;i<tips; i++)
    {
      distm[i][i] = 0.0;
      for(j=0;j<i;j++)
	{
	  distm[i][j] = distm[j][i] = RANDUM();
	  //printf("%f ",distm[i][j]);
	}
      //printf("\n");
    }
}


void
constrain_distance_zeromig (MYREAL **m, option_fmt *options, data_fmt * data, long locus,
                            long tips, char *custm)
{
  long top, ii;
    long pop;
  //  ind
   long   j, i = 0;
  long *pops;
  pops = (long *) mycalloc (tips, sizeof (long));
  for (pop = 0; pop < data->numpop; pop++)
    {
      if(options->randomsubset>0 && options->randomsubset <  data->numind[pop][locus])
	{	
	  top = options->randomsubset;
	}
      else
	{
	  top = data->numind[pop][locus];
	}
      for (ii = 0; ii < top; ii++)
	{
	  //ind = data->shuffled[pop][locus][ii];
	  pops[i++] = options->newpops[pop];
	}
    }
    for (i = 0; i < tips; i++)
    {
        for (j = 0; j < i; j++)
        {
            if (custm[pops[i] + pops[i] * pops[j]] == '0')
                m[i][j] = m[j][i] = 1000;
	    //printf("%f ",m[i][j]);
        }
	//printf("\n");
    }
    myfree(pops);
}

void
distance_EP (char **data, long tips, MYREAL **m)
{
    long i, j;
    for (i = 0; i < tips; i++)
    {
        for (j = 0; j < i; j++)
        {
            if (!strcmp (data[i], data[j]))
                m[i][j] = m[j][i] = fabs (rannor (1., 0.1));
            else
                m[i][j] = m[j][i] = fabs (rannor (0., 0.1));
        }
    }
}

void
distance_micro (char **data, long tips, MYREAL **m)
{
    long i, j;
    for (i = 0; i < tips; i++)
    {
        for (j = 0; j < i; j++)
        {
            m[i][j] = m[j][i] = pow (atof (data[i]) - atof (data[j]), 2.);
            m[i][j] = m[j][i] = fabs (rannor (m[i][j], 0.1));
        }
    }
}

/* calculate pairwise distances using allele disimilarity */
void
distance_allele (world_fmt * world, option_fmt * options, 
                 long tips, MYREAL **distm)
{
    char **mdata;
    long pop;
	
    mdata = (char **) mycalloc ((tips+1), sizeof (char *));
    for (pop = 0; pop < tips; pop++)
    {
        mdata[pop] =
		(char *) mycalloc (options->allelenmlength, sizeof (char));
        strcpy (mdata[pop], strrchr(world->nodep[pop]->nayme,'!')+1);
    }
    if (world->options->datatype == 'a')
        distance_EP (mdata, tips, distm);
    else
        distance_micro (mdata, tips, distm);
    for (pop = 0; pop < tips; pop++)
    {
        myfree(mdata[pop]);
    }
    myfree(mdata);
}

/* calculate  pairwise distances using sequence similarity */
void
distance_sequence (data_fmt * data, option_fmt *options, long locus, long tips, long sites,
                   long nmlength, MYREAL **m)
{
  long top;
  long ii;
  long i = 0, j, z, pop, ind;
  char **dat;
  if (data->distfile != NULL)
    {
        read_distance_fromfile (data->distfile, tips, nmlength, m);
    }
  else
    {
      charvec2d (&dat, tips, data->totalsites[locus]);
      for (pop = 0; pop < data->numpop; pop++)
	  {
	    if(options->randomsubset>0 && options->randomsubset <  data->numind[pop][locus])
	      {	
		top = options->randomsubset;
	      }
	    else
	      {
		top = data->numind[pop][locus];
	      }
	    for (ii = 0; ii < top; ii++)
	      {
		ind = data->shuffled[pop][locus][ii];
		for(j=0; j < data->totalsites[locus]; j++) 
		  dat[i][j] = data->yy[pop][ind][locus][0][j][0];
		i++;
            }
        }
        if (i != tips)
        {
            error ("Mistake in distance_sequence() tips is not equal sum(i)\n");
        }
        for (i = 0; i < tips; i++)
        {
			
            for (j = i + 1; j < tips; j++)
            {
                //to come m[i][j] = m[j][i] = make_ml_distance(dat[i], dat[j], i, j);
				
                for (z = 0; z < sites; z++)
                {
					
                    if (dat[i][z] != dat[j][z])
                    {
                        //m[i][j] = m[j][i] += fabs(rannor(1.0, 0.1));
                        m[i][j] = m[j][i] += 1.0;
                    }
                }
            }
        }
        free_charvec2d(dat);
    }
}

void
fix_times (world_fmt * world, option_fmt * options)
{
    long k;
    //MYREAL interval;
    MYREAL tipdate=0.0;
    vtlist *tl;  // = &world->treetimes[0].tl[world->treetimes[0].T - 1];
    node *theNode;
    MYREAL age = world->data->maxsampledate 
      * world->options->meanmu[world->locus]
      * world->options->generation_year  * world->options->mu_rates[world->locus];
    if (!options->usertree)
    {
#ifdef TREEDEBUG1
      printf("%i> %s\nage= %f\nmaxsamp=%f\nmut=%f\ngen=%f\nmu_rate=%g\n-----------------------------\n",
	     myID, "start timelist ---------------------", age,world->data->maxsampledate,
	     world->options->meanmu[world->locus],1./world->options->generation_year,
	     world->options->mu_rates[world->locus]);
#endif
      for (k = 0; k < world->treetimes[0].T - 1; k++)
	{
	  tl = &world->treetimes[0].tl[k];
	  theNode = tl->eventnode;
	  if(theNode!=NULL)
	    {
	      if(theNode->type == 't')
		{
		  //tipdate = find_tipdate(theNode->nayme, theNode->pop, world);
		  tipdate = theNode->tyme;
#ifdef TREEDEBUG1
		  printf("DEBUG TIPDATE %i> %s %g\n",myID, theNode->nayme, tipdate);
#endif
		  //theNode->tyme = tipdate;
		  tl->age = tipdate;	
		}
	      else
		{
		  //k here is never 0, so this k-1 seems correct
		  MYREAL interval = inverse_logprob_noevent (world, k-1);
		  if (interval >= ROOTLENGTH)
		      interval = age; //hack
		  age += interval;
		  tl->age = age;	
		  adjust_time (theNode, age);
#ifdef TREEDEBUG1
		  printf("%i> %10.10s %g\n",myID, " ", theNode->tyme);
#endif
		}

	    }
	}
      tl = &world->treetimes[0].tl[k];
      tl->age = age + 10000.; /* this is the root */
      adjust_time (tl->eventnode, age+10000);
#ifdef TREEDEBUG1
      printf("%i> %s %g\n",myID, "end timelist ---------------------", age+10000 );
#endif
    }
    set_v (world->root->next->back);
}

/* single addition tree
using the custm2 migration matrix and also the split times to generate the first tree
 */
node * get_random_sibling(node *origin, node **nodelist, long *nodenum, node *interior,
			  long *lineages, double utime)
{
  (void) lineages;
  node * sibling;
  long  z=0;
  long i;
  long origin_loc= -1;
  //long chosen[*nodenum];
  // lineages[origin->pop]
  long * chosen = (long *) mycalloc((*nodenum),sizeof(long));
  for (i=0; i < *nodenum; i++)
    {
      if (origin != nodelist[i])
	{
	  if(origin->pop == nodelist[i]->pop)
	    {
	      chosen[z] = i;
	      z++;
	    }
	}
      else
	{
	  origin_loc=i;
	}
    }
  z--;
  //printf("[%li] z=%li ",*nodenum, z);
  z  = RANDINT(0,z);
  //printf("(%li) ",z);
  sibling = nodelist[chosen[z]];
  interior->tyme = utime;
  interior->next->back = origin;
  origin->back = interior->next;
  interior->next->next->back = sibling;
  sibling->back = interior->next->next;
  interior->actualpop = interior->pop = origin->pop;
  nodelist[origin_loc] = interior;
  nodelist[chosen[z]] = nodelist[*nodenum-1];
  nodelist[*nodenum-1] = NULL;
  *nodenum -= 1;
  free(chosen);
  return interior;
}

void print_lineages(long *l, long n)
{
  int i;
  for(i=0;i<n;i++)
    {
      printf("%li ",l[i]);
    }
  printf("\n");
}


void my_start_eventtime(double age, world_fmt *world, node **nodelist, long *lineages,long simtips, 
			double *shortt, long *shorti, char * shorte, long *to, long *from)
{
  long i;
  double utime;
  long tox = -1;
  long fromx = -1;
  char event = ' ';
  *shortt = (double) HUGE;
  *shorte = '_';
  *shorti = -1;
  for (i=0;i<simtips;i++)
    {
      if(lineages[nodelist[i]->pop] == 0)
	{
	  printf("%li %li %li %li -- %li, %li of %li\n",lineages[0],lineages[1],lineages[2],
		 lineages[3],nodelist[i]->pop,i,simtips);
	  //error("node in pop without lineage");
	}
      // set lineages to conditional tree, aka we work with nodelist[i]
      // and want to connect or migrate/diverge with other lineages
      lineages[nodelist[i]->pop] -= 1;
      utime = eventtime_single(NULL, world, nodelist[i]->pop, 0, lineages, age, &event, &tox, &fromx);
      lineages[nodelist[i]->pop] += 1;
      //printf("%i: %f %f %li %li\n",i,utime,age,tox,fromx);
      // set lineages back to full tree
      if (utime <= *shortt)
	{
	  *shorti = i;
	  *shortt = utime;
	  *shorte = event;
	  *from = fromx;
	  *to = tox;
	}
    }
}

void
my_random_tree (world_fmt * world, long tips)
{
  double shortt;
  char shorte;
  long simtips;
  long maxsimtips;
  long i;
  long shorti;
  long z = tips;
  long zz;
    node *interior;//= world->nodep[z];
    node ** nodelist = (node **) mycalloc(tips,sizeof(node *));
    vtlist * datelist = (vtlist *) mycalloc(tips+1,sizeof(vtlist));
    long * lineages = (long *) mycalloc(world->numpop,sizeof(long)); 
  node *origin;
  long from;
  long to;
  double age=0.0;
  //double otime = 0.0;
  double starttime = (double) HUGE;
  precalc_world(world);
  for (i=0;i<tips;i++)
    {
      datelist[i].eventnode = world->nodep[i];
      datelist[i].age = world->nodep[i]->tyme;
    }
  datelist[i].eventnode = NULL;
  datelist[i].age = (double)HUGE;
  qsort(datelist,(size_t) (tips+1), sizeof(vtlist),agecmp);
  starttime = datelist[0].age;
  simtips = 0;
  for (i=0; i<tips;i++)
    {
      origin = world->nodep[i]; // all tips are loaded into nodelist	 
      if(origin->tyme <= starttime) // and lineages is set up [without dates, n1 n2 ... nn
	{
	  nodelist[simtips] = origin;
	  lineages[origin->pop] += 1;
	  simtips++;
	}
    }
  zz = simtips;
  maxsimtips = simtips;
  while(maxsimtips < tips || simtips > 1)
    {
      while (zz < tips && (fabs(age - datelist[zz].age) < DBL_EPSILON))
	{
	  nodelist[simtips] = datelist[zz].eventnode;
	  lineages[nodelist[simtips]->pop] += 1;
	  simtips++;
	  maxsimtips++;
	  zz++;
	}
      //printf("myrandomtree: lineages: %li %li %li ",lineages[0],lineages[1],lineages[2]);
      my_start_eventtime(age, world,nodelist,lineages,simtips,&shortt,&shorti,&shorte, &to, &from);
      assert(shorte != ' ');
      //printf(" |%c|\n",shorte);
      if (age < datelist[zz].age && age + shortt > datelist[zz].age && datelist[zz].age < (double) HUGE)
	{
	  age = datelist[zz].age;
	  shorte= 't';
	  nodelist[simtips] = datelist[zz].eventnode;
	  lineages[nodelist[simtips]->pop] += 1;
	  simtips++;
	  maxsimtips++;
	  if(zz<tips)
	    {
	      zz++;
	    }
	  continue;
	}
      else
	{
	  if(shortt < (double) HUGE)	  
	    age += shortt;
	  else
	    {
	      if (maxsimtips < tips)
		{
		  age = datelist[zz].age;
		  if(zz<tips)
		    zz++;
		}
	      continue;
	    }
	}
      origin = nodelist[shorti];
      switch(shorte)
	{
	case 'c':
	  interior = world->nodep[z]; //gets first nodep element after all the tips
	  origin = get_random_sibling(origin, nodelist, &simtips, interior, lineages, age);
	  nodelist[shorti] = origin;
	  lineages[origin->pop] -= 1;
	  //printf("c:%f %li [pop=%li:lin=%li]", age,z,origin->pop,lineages[origin->pop]);
	  //print_lineages(lineages,world->numpop);
	  z++;
	  break;
	case 'm':
	  origin = add_migration(world,origin,'m',from, to,shortt);
	  lineages[origin->actualpop] -= 1;
	  lineages[origin->pop] += 1;
	  nodelist[shorti] = origin;
	  //printf("m:%f ",age);
	  //print_lineages(lineages,world->numpop);
	  break;
	case 'd':
	  origin = add_migration(world,origin,'d',from, to,shortt);
	  origin->tyme = age;
	  lineages[origin->actualpop] -= 1;
	  lineages[origin->pop] += 1;
	  nodelist[shorti] = origin;
	  //printf("d:%f [pop=%li/%li, lin=%li/%li] ",age, origin->pop,origin->actualpop,lineages[origin->pop],lineages[origin->actualpop]);
	  //print_lineages(lineages,world->numpop);
	  break;
	case 't':
	  //added a new tip
	  break;
	default :
	  fprintf(stderr,"Unknown type: |%c|\n",shorte);
	  error("Wrong Type");
	}
      // printf("+");
    }
  if (nodelist[0]!=NULL)
    {
      world->root->next->back = nodelist[0];
      nodelist[0]->back = world->root->next;
    }
  else
    {
      error("problem in assigning nodelist[0] value");
    }
  myfree(nodelist);
  myfree(datelist);
  myfree(lineages);
  //#ifdef TREEDEBUG
  //debugtreeout (stdout, crawlback (world->root->next), crawlback (world->root->next),0);
  //#endif
}

/*
 * creates a UPGMA tree: x     = distance matrix which will be destroyed
 * through the process, tips  = # of sequences/alleles, nodep = treenodes
 * have to be allocated for ALL nodes
 *
 * This code is stripped neighbor-joining code out of phylip v3.6. Only the
 * upgma option is present.
 */
void
upgma (world_fmt * world, MYREAL **x, long tips, node ** nodep)
{
    long nc, nextnode, mini = -900, minj = -900, i, j, jj, ia, iaa, ja, jaa; 
    MYREAL zz = (world->data->maxsampledate 
		 * world->options->meanmu[world->locus]
		 * world->options->generation_year  * world->options->mu_rates[world->locus]+ 1);
    MYREAL total, tmin, bi, bj, /* ti, tj, */ da;
    MYREAL *av;
    long *oc;
    node **cluster;
    long *enterorder;
	MYREAL denom=1.;
    /* First initialization */
    enterorder = (long *) mycalloc (tips, sizeof (long));
    for (ia = 0; ia < tips; ia++)
        enterorder[ia] = ia;
    jumble (enterorder, tips);
    nextnode = tips;
    av = (MYREAL *) mycalloc (tips, sizeof (MYREAL));
    oc = (long *) mymalloc (tips * sizeof (long));
    cluster = (node **) mycalloc (tips, sizeof (node *));
    for (i = 0; i < tips; i++)
        oc[i] = 1;
    for (i = 0; i < tips; i++)
        cluster[i] = nodep[i];
    /* Enter the main cycle */
    for (nc = 0; nc < tips - 1; nc++)
    {
        tmin = 99999.0;
        /* Compute sij and minimize */
        for (jaa = 1; jaa < tips; jaa++)
        {
            ja = enterorder[jaa];
            if (cluster[ja] != NULL)
            {
                for (iaa = 0; iaa < jaa; iaa++)
                {
                    ia = enterorder[iaa];
                    if (cluster[ia] != NULL)
                    {
                        total = x[ia][ja];
			//printf("ia=%li ja=%li x[ia][ja]=%f\n",ia,ja,x[ia][ja]);
                        if (total < tmin)
                        {
                            tmin = total;
                            mini = ia;
                            minj = ja;
                        }
                    }
                }
            }
        }   /* compute lengths and print */
        bi = x[mini][minj] / 2.0 - av[mini];
        bj = x[mini][minj] / 2.0 - av[minj];
        av[mini] += bi;
        nodep[nextnode]->next->back = cluster[mini];
      	//xcode
        if(cluster[mini]!=NULL)
        {
            if(nodep[nextnode]->next!=NULL)
            {
                cluster[mini]->back = nodep[nextnode]->next;
            }
            nodep[nextnode]->next->next->back = cluster[minj];
            cluster[minj]->back = nodep[nextnode]->next->next;
            cluster[mini]->back->v = cluster[mini]->v = bi;
            cluster[minj]->back->v = cluster[minj]->v = bj;
            cluster[mini] = nodep[nextnode];
        }
        adjust_time (nodep[nextnode], (MYREAL) zz++);
        cluster[minj] = NULL;
        nextnode++;
        /* re-initialization */
        denom = (oc[mini] + oc[minj]);
        for (jj = 0; jj < tips; jj++)
        {
            if (cluster[jj] != NULL)
            {
                da = (x[mini][jj] * oc[mini] + x[minj][jj] * oc[minj])/denom;
                x[mini][jj] = da;
                x[jj][mini] = da;
            }
        }
        for (j = 0; j < tips; j++)
        {
            x[minj][j] = x[j][minj] = 0.0;
        }
        oc[mini] += oc[minj];
    }
    /* the last cycle */
    for (i = 0; i < tips; i++)
    {
        if (cluster[i] != NULL)
            break;
    }
    world->root->next->back = cluster[i];
    cluster[i]->back = world->root->next;
    myfree(av);
    myfree(oc);
    myfree(cluster);
    myfree(enterorder);
}

void
set_top (world_fmt * world, node * p, long locus)
{
	
	
    if (p->type == 't')
    {
        p->top = TRUE;
        //p->tyme = 0.0;
        return;
    }
    p->top = TRUE;
    p->next->top = FALSE;
    if (!(p->type == 'm' || p->type == 'd'))
    {
        p->next->next->top = FALSE;
    }
    set_top (world, p->next->back, locus);
    if (!(p->type == 'm' || p->type == 'd'))
    {
        set_top (world, p->next->next->back, locus);
    }
    if (p == crawlback (world->root->next))
    {
        p->back->top = FALSE;
        p->back->tyme = ROOTLENGTH;
    }
}    /* set_top */


void
set_v (node * p)
{
    if (p->type == 't')
    {
        p->v = p->length = lengthof (p);
        return;
    }
    ltov (p);
    set_v (crawlback (p->next));
    set_v (crawlback (p->next->next));
}    /* set_v */

void
ltov (node * p)
{
    p->v = lengthof (p);
}    /* ltov */

/*
 * cost matrix COST is for the ISLAND and MATRIX model all 1 with a 0
 * diagonal, this will be changable through options, but perhaps this is not
 * so important
 */
void
calc_sancost (MYREAL **cost, world_fmt *world)
{
    long i, j;
    long z=0;
    long numpop = world->numpop;
    char *custm2 = world->options->custm2;
    for (i = 0; i < numpop; i++)
    {
        for (j = 0; j < numpop; j++)
        {
	  z = mm2m(i,j,numpop);
	  if (z<numpop)
	    cost[i][i] = 0.0;
	  switch(custm2[z])
	    {
	    case '0':
	      cost[i][j] = 1000.0;
	      break;
	    case 'd':
	      cost[i][j] = 100.0;
	      break;
	    case 'm':
	    default:
	      cost[i][j] = 1.0;
	      break;
	    }
        }
    }
}

void
sankoff (world_fmt * world)
{
    MYREAL **cost;
    long i;
    cost = (MYREAL **) mymalloc (sizeof (MYREAL *) * (size_t) world->numpop);
    cost[0] =
      (MYREAL *) mymalloc (sizeof (MYREAL) * (size_t) (world->numpop * world->numpop));
    for (i = 1; i < world->numpop; i++)
    {
        cost[i] = cost[0] + world->numpop * i;
    }
    calc_sancost (cost, world);
    santraverse (world, crawlback (world->root->next), cost, world->numpop);
    myfree(cost[0]);
    myfree(cost);
}

void
jumble (long *s, long n)
{
    long *temp, i, rr, tn = n;

    if (n>0)
      {
	temp = (long *) mycalloc (1, sizeof (long) * (size_t) n);
	memcpy (temp, s, sizeof (long) * (size_t) n);
	for (i = 0; i < n && tn > 0; i++)
	  {
	    s[i] = temp[rr = RANDINT (0, tn - 1)];
	    temp[rr] = temp[tn - 1];
	    tn--;
	  }
	myfree(temp);
      }
}


#ifdef WIN32
#define SUBSETRANDOM (((double) rand())/((double) RAND_MAX))
#else
#define SUBSETRANDOM (drand48())
#endif

// this uses a different stream than the usual random number 
// and is used to guarantee that the selection of individual subsets can be 
// maintained between different runs
void
jumble_ownseed (long *s, long n)
{
    long *temp, i, rr, tn = n;
	
    temp = (long *) mycalloc (n, sizeof (long));
    memcpy (temp, s, sizeof (long) * (size_t) n);
    for (i = 0; i < n && tn > 0; i++)
    {
      s[i] = temp[rr = (long) ( SUBSETRANDOM * (tn - 1))];     
      temp[rr] = temp[tn - 1];
      tn--;
    }
    myfree(temp);
}

void
santraverse (world_fmt *world, node * theNode, MYREAL **cost, long numpop)
{
    long i, ii, tie, which = 0;
    node *p = NULL, *q = NULL, *tmp = NULL, *left = NULL, *right = NULL;
    MYREAL best;
    long *poplist;
    //char type;
    poplist = (long *) mycalloc (numpop, sizeof (long));
    if (theNode->type != 't')
    {
        if (RANDUM () > 0.5)
        {
            left = theNode->next;
            right = theNode->next->next;
        }
        else
        {
            left = theNode->next->next;
            right = theNode->next;
        }
        if (left->back != NULL)
        {
          p = crawlback (left);
	  santraverse (world, p, cost, numpop);
        }
        if (right->back != NULL)
        {
          q = crawlback (right);
	  santraverse (world, q, cost, numpop);
        }
        best = MYREAL_MAX;
        tie = 0;
        for (i = 0; i < numpop; i++)
            poplist[i] = i;
        jumble (poplist, numpop);
        for (ii = 0; ii < numpop; ii++)
        {
            i = poplist[ii];
          	if(p!=NULL && q!=NULL)
            {              
            	theNode->s[i] =
                	minimum (cost[i], p->s, numpop) + minimum (cost[i], q->s, numpop);
            }
            if (theNode->s[i] < best)
            {
                best = theNode->s[i];
                which = i;
                tie = 0;
            }
            else
            {
	      if (fabs(theNode->s[i] - best) <= DBL_EPSILON)
                {
                    tie++;
                    which = i;
                }
            }
        }
        if (tie != 0)
        {
            theNode->pop = theNode->actualpop = which =
			ranbest (theNode->s, tie, best, numpop);
            theNode->s[which] -= SANKOFF_DELTA;
        }
        else
        {
            theNode->pop = theNode->actualpop = which;
        }
        if (p!=NULL && p->actualpop != which)
        {
	  //	  printf("%i> p: theNode->tyme=%f, p->tyme=%f mignode=%f\n",myID, theNode->tyme,p->tyme,p->tyme + (theNode->tyme - p->tyme) / 2.);
	  tmp = set_type2(world, p,theNode, world->options->custm2);
	  //tmp = set_type(world, p->actualpop,theNode->actualpop, world->options->custm2, world->numpop);
	  //tmp = add_migration (world, p, type, theNode->actualpop, p->actualpop,
	  //		     (MYREAL) RANDDOUBLE(0.0, theNode->tyme - p->tyme));
            left->back = tmp;
            tmp->back = left;
            //debug FPRINTF(startfile, "%li %li\n", theNode->actualpop, p->actualpop);
        }
        if (q!=NULL && q->actualpop != which)
        {
	  tmp = set_type2(world, q,theNode, world->options->custm2);
	  //type = set_type(world, q->actualpop,theNode->actualpop, world->options->custm2, world->numpop);

	  //tmp = add_migration (world, q, type, theNode->actualpop, q->actualpop,
	  //	     (MYREAL) RANDDOUBLE(0.0, theNode->tyme - q->tyme));

	    //printf("%i> q: theNode->tyme=%f, q->tyme=%f mignode=%f\n",myID, theNode->tyme,q->tyme,q->tyme + (theNode->tyme - q->tyme) / 2.);
            right->back = tmp;
            tmp->back = right;
            //debug FPRINTF(startfile, "%li %li\n", theNode->actualpop, q->actualpop);
        }
    }
    myfree(poplist);
}

/* returns minimum for sankoff routine */
MYREAL
minimum (MYREAL *vec1, MYREAL *vec2, long n)
{
    long j;
    MYREAL summ, min = MYREAL_MAX;
    for (j = 0; j < n; j++)
    {
        if (vec2[j] < MYREAL_MAX)
        {
            if ((summ = vec1[j] + vec2[j]) < min)
                min = summ;
        }
    }
    return min;
}

long
ranbest (MYREAL *array, long tie, MYREAL best, long n)
{
    long i;
    long which = RANDINT (0, tie-1);
    for (i = 0; i < n; i++)
    {
        if (fabs (array[i] - best) < EPSILON)
        {
            if (which == 0)
                return i;
            else
                which--;
        }
    }
    return -2;
}

///
/// allocate memory for data in nodes for all linked loci
void
alloc_seqx (world_fmt * world, node * theNode, long locus)
{
  mutationmodel_fmt *s;
  long endsite;
  long sublocus;
  long sublocistart = world->sublocistarts[locus];
  long sublociend = world->sublocistarts[locus+1];
  //  printf("%i> locus:%li subloci:%li\n",myID, locus, world->numsubloci[locus]);
  if(theNode->x==NULL)
    {
      theNode->x = (xarray_fmt *) mycalloc(world->numsubloci[locus],sizeof(xarray_fmt));
      //printf("C");fflush(stdout);
    }
  else
    {
      theNode->x = (xarray_fmt *) myrealloc(theNode->x, (size_t) world->numsubloci[locus] * sizeof(xarray_fmt));
      //printf("R");fflush(stdout);
    }
  theNode->scale = (MYREAL **) mycalloc(world->numsubloci[locus],sizeof(MYREAL*));
  for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
    {
      s = &world->mutationmodels[sublocus];
      endsite = s->numpatterns;
      const long xs = sublocus - sublocistart;
      if (strchr ("u", s->datatype))
	{
	  endsite = endsite * (s->addon + 1) + s->addon;
	}
      if (strchr ("nh", s->datatype))
	{
	  endsite += s->addon;
	}
      if (strchr (SEQUENCETYPES, s->datatype))
	{	 
	  allocate_xseq(&theNode->x[xs], endsite, s->numsiterates);
    	  theNode->scale[xs] = (MYREAL *) mycalloc (endsite, sizeof (MYREAL));
	}
      else
	{
	  if (s->datatype=='b' && s->maxalleles != XBROWN_SIZE)
	    {
	      // on initializing we calculate maxalleles but we do not need this for brownian
	      s->maxalleles = XBROWN_SIZE;
	    }
	  theNode->x[xs].a =
	    (MYREAL *) mycalloc (s->maxalleles + 1 , sizeof (MYREAL));
	  theNode->scale[xs] =
	    (MYREAL *) mycalloc (s->maxalleles + 1 , sizeof (MYREAL));
	}
    }
}

///
/// allocate sequence data type
void     allocate_xseq(xarray_fmt *x, long sites, long categs)
{
    long j;
    (*x).s = (phenotype) mycalloc (sites, sizeof (ratelike *));
#ifdef VARMUT
    (*x).s[0] = (ratelike) mycalloc (linkedloci * categs, sizeof (MYREAL) * sites[datamodeltype] * sitelikesize[datamodeltype]);
#else
    (*x).s[0] = (ratelike) mycalloc ((categs * sites), sizeof (sitelike));
#endif /*VARMUT*/
    for (j = 1; j < sites; j++)
      (*x).s[j] = (*x).s[0] + categs * j;
}

///
/// zeroes sequence data vector
void     zero_xseq(xarray_fmt *x, world_fmt *world)
{
  long rcategs, sites;
    long j;
    mutationmodel_fmt *s;
    long locus = world->locus;
    const long sublocistart = world->sublocistarts[locus];
    const long sublociend   = world->sublocistarts[locus+1];
    long sublocus;
    for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
      {
	s = &world->mutationmodels[sublocus];
	rcategs = s->numsiterates;
	sites   = s->numpatterns;
	const long xs = sublocus - sublocistart;
	memset(x[xs].s[0], 0, (size_t) (rcategs * sites) * sizeof (sitelike));
	for (j = 1; j < sites; j++)
	  x[xs].s[j] = x[xs].s[0] + rcategs * j;
      }
}


void make_sequences(long sublocus, mutationmodel_fmt *s, node *theNode, site_fmt **datapart)
{
  long k;
  long j;
  long l;
  long numpatterns = s->numpatterns;
  //MYREAL * scoring_error = s->scoring_error;
  MYREAL scoring_error[4];
  scoring_error[0]=0.0;
  scoring_error[1]=0.0;
  scoring_error[2]=0.0;  
  scoring_error[3]=0.0;
  if (s->estimateseqerror)
    {
      if (theNode->sequence==NULL)
	theNode->sequence = (char *) mycalloc(numpatterns,sizeof(char));
      else
	theNode->sequence = (char *) myrealloc(theNode->sequence, (size_t) numpatterns * sizeof(char));
    }
  for (k = 0; k < numpatterns; k++)
    {
      // once I figured out how to send only the condensed data this needs to be activated TODO
      //#ifdef MPI
      //j = k;
      //#else
      j = s->alias[k] - 1;
      //#endif
      //
      //if (world->options->has_seqerror)
      // ordering sequence the same way as the x array using k and not j
      if (s->estimateseqerror)
	{
	  theNode->sequence[k] = datapart[0][j][0];
	}
      for (l = 0; l < s->numsiterates; l++)
        {
          set_nucleotide(theNode->x[sublocus].s[k][l],datapart[0][j][0], scoring_error);
        }
      //printf("%c",datapart[0][j][0]);
    }
  //  printf("@@@@#@\n");
}

void make_invarsites(long sublocus, mutationmodel_fmt *s, node *theNode)
{
  long k;
  //long j;
  long l;
  long numpatterns = s->numpatterns;
  MYREAL * scoring_error = s->scoring_error;
  // add four monomorphic sites to the list -- SNP
  //long added = numpatterns + s->addon;
  //  for (k = numpatterns; k < added; k++)
  //  {
  //    long thesite = k - numpatterns;
  k = numpatterns;
  for (l = 0; l < s->numsiterates; l++)
    {
      set_nucleotide(theNode->x[sublocus].s[k][l],'A', scoring_error);
      set_nucleotide(theNode->x[sublocus].s[k+1][l],'C', scoring_error);
      set_nucleotide(theNode->x[sublocus].s[k+2][l],'G', scoring_error);
      set_nucleotide(theNode->x[sublocus].s[k+3][l],'T', scoring_error);
      //          theNode->x[sublocus].s[k][l][thesite] = 1.0 - scoring_error;
    }
}


long make_alleles(long sublocus, mutationmodel_fmt *s, node *theNode, option_fmt *options, data_fmt *data, site_fmt *datapart2)
{
  //long ii,top;
  //long pop, ind;
  //long zz=0;
  //long zpop;
  long iu;
  char a1[DEFAULT_ALLELENMLENGTH];
  
  strcpy (a1, datapart2[0]);
  if (strcmp (a1, "?"))
    {
      sprintf(theNode->nayme,"A!%s",a1);
      theNode->x[sublocus].a[findAllele (data, a1, sublocus)] = 1.0;
      return 0;
    }
  else
    {
      if(options->include_unknown)
	{
	  sprintf(theNode->nayme,"A!%s",a1);
	  for(iu=0;iu < s->maxalleles; iu++)
	    theNode->x[sublocus].a[iu] = 1.0;
	  return 0;
	}
    }
  return 1;
}

void
find_minmax_msat_allele (world_fmt * world, data_fmt * data, long locus, 
                         long *smallest, long *biggest)
{
    long pop, ind, tmp;
    *biggest = 0;
    *smallest = LONG_MAX;
    for (pop = 0; pop < data->numpop; pop++)
    {
        for (ind = 0; ind < data->numind[pop][locus]; ind++)
        {
	  const long sublocistart = world->sublocistarts[locus];
	  const long sublociend   = world->sublocistarts[locus+1];
	  long sublocus;
	  for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
	    {
	      if (data->yy[pop][ind][sublocus][0][0][0] != '?')
		{
		  if ((tmp = atoi (data->yy[pop][ind][sublocus][0][0])) > *biggest)
                    *biggest = tmp;
		  if (tmp < *smallest)
                    *smallest = tmp;
		}
	      if(world->options->datatype!='@')
		{
		  if (data->yy[pop][ind][sublocus][1][0][0] != '?')
		    {
		      if ((tmp = atoi (data->yy[pop][ind][sublocus][1][0])) > *biggest)
			*biggest = tmp;
		      if (tmp < *smallest)
			*smallest = tmp;
		    }
		}
	    }
	}
    }
}




///
/// make microsatellites using the mutation model syntax
long make_microsatellites (long sublocus, mutationmodel_fmt *s, node *theNode, option_fmt *options, site_fmt *datapart2)
{
  long iu;
  const long smax = s->microrange;
  const long smallest = s->microstart;
  char a1[DEFAULT_ALLELENMLENGTH];
  strcpy (a1, datapart2[0]);
  if (strcmp (a1, "?"))
    {
      theNode->x[sublocus].a = (MYREAL *) myrealloc (theNode->x[sublocus].a, sizeof (MYREAL) * (size_t) smax);
      theNode->x[sublocus].a[atoi (a1) - smallest] = 1.0;
      sprintf(theNode->nayme,"A!%s",a1);
      return 0;
    }
  else
    {
      if(options->include_unknown)
	{
	  for(iu=0;iu<s->maxalleles; iu++)
	    theNode->x[sublocus].a[iu] = 1.0;
	  sprintf(theNode->nayme,"A!%s",a1);
	  return 0;
	}
    }	  
  return 1;
}

/// NEW
///
long make_brownian(long sublocus, mutationmodel_fmt *s, node *theNode, option_fmt *options, site_fmt *datapart)
{
  //long len;
  char a1[DEFAULT_ALLELENMLENGTH];
  MYREAL bsum = s->browniandefault;
  strcpy (a1, datapart[0]);
  if (strcmp (a1, "?"))
    {
      sprintf(theNode->nayme,"A!%s",a1);
      theNode->x[sublocus].a[0] = atof (a1);
      return 0;
    }
  else
    {
      if(options->include_unknown)
	{
	  sprintf(theNode->nayme,"A!%s",a1);
	  theNode->x[sublocus].a[0] = bsum; //average over other values
	  // this might not be very sensible at all.
	  return 0;
	}
    }
  return 1;
}


/*---------------------------------------------
free_treetimes frees the ptr_array of timeslist
*/
void
free_treetimes (world_fmt * world, long size)
{
  //  long i;
  while (size >= 0)
    {
      //      for(i=world->treetimes[size].allocT-1; i>= 0; i--)
      //	myfree(world->treetimes[size].tl[i].lineages);
      myfree(world->treetimes[size].lineages);
      myfree(world->treetimes[size].tl);
      size--;
    }
}

///
/// construct timelist
void
construct_tymelist (world_fmt * world, timelist_fmt * timevector)
{
    long z = 0;
    long tips = 0;
    //long ii;
    //long T;
    //    long pop;
    MYREAL tmp;
    timevector->numpop = world->numpop;
    traverseNodes (world->root, &timevector, &z, world, &tips);
    timevector->T = z;
#ifdef TREEDEBUG1
    printf("timevector->T=%li tips=%li\n",timevector->T,tips);
#endif
    //    insert_time_boundaries(timevector,world);
    qsort ((void *) timevector->tl, (size_t) timevector->T, sizeof (vtlist), agecmp);
    if ((*timevector).tl[(*timevector).T - 1].eventnode->type != 'r')
    {
        z = 0;
        while ((*timevector).tl[z].eventnode->type != 'r')
            z++;
		
	tmp = 	(*timevector).tl[(*timevector).T - 1].eventnode->tyme + 10000.;
	(*timevector).tl[z].eventnode->tyme = tmp;
	(*timevector).tl[z].eventnode->next->tyme = tmp;
	(*timevector).tl[z].eventnode->next->next->tyme = tmp;
	
        (*timevector).tl[z].age = (*timevector).tl[z].eventnode->tyme;
        qsort ((void *) timevector->tl, (size_t) timevector->T, sizeof (vtlist), agecmp);
        warning("construct_tymelist root moved: new time = %f\n",(*timevector).tl[z].eventnode->tyme);
    }
    timeslices (&timevector);
#ifdef DEBUGXX
    printf("%i> in construct_tymelist(): %li\n", myID, timevector->T);
#endif
    add_partlineages (world->numpop, &timevector, world);
#ifdef TREEDEBUG1
    fprintf(stderr,"> %i @///////////////////\n",myID);
    long T = timevector->T-2;
    long ii, pop;
    for (ii = T; ii >= 0; ii--)
      {
	for (pop = 0; pop < world->numpop; pop++)
	  fprintf(stderr,"%li ",timevector->tl[ii].lineages[pop]);
	fprintf(stderr,"%20.10f %20.10f | %f %c %li -> %li (%li)(%li)\n",timevector->tl[ii].age,timevector->tl[ii].eventnode->tyme, timevector->tl[ii].eventnode->tyme, timevector->tl[ii].eventnode->type, timevector->tl[ii].from, timevector->tl[ii].to, timevector->tl[ii].eventnode->id, showtop(timevector->tl[ii].eventnode->back)->id);
      }
    fprintf(stderr,"> %i @/////////////////////////////////////////////\n", myID);
#endif
}



/// find a tipdate
MYREAL 
find_tipdate(char * id, long pop, world_fmt *world)
{
  MYREAL date;
  long ind;
  long slen;
  long locus = world->locus;
  tipdate_fmt *sampledates = world->data->sampledates[pop][locus];

  for(ind=0;ind < world->data->numind[pop][locus]; ind++)
    {
      if(sampledates[ind].name != NULL)
	{
	  slen = (long) strlen(sampledates[ind].name);
	  if(!strncmp(sampledates[ind].name,id, (size_t) slen))
	    {
	      date = (sampledates[ind].date     
		      * world->options->generation_year 
		      * world->options->meanmu[locus]) * world->options->mu_rates[locus];
	      //	      printf("%i> id='%s' name='%s' realdate=%f date=%f locus=%li pop=%li\n",myID, id,sampledates[ind].name, sampledates[ind].date, date, locus, pop);
	      //fflush(stdout);
	      return date;
	    }
	}
    }
  return 0.0;
}

///
/// traverse the tree and writes node-information into the real timelist also
/// takes care that the size of timelist is increased accordingly
void
traverseNodes (node * theNode, timelist_fmt ** timevector, long *slice, world_fmt *world, long *tips)
{
    //#ifdef DEBUG
    // MYREAL tipdate = 0.0;
    //#endif
    if (theNode != NULL)
      {  
          if (theNode->type != 't')
            {
              if (theNode->next->back != NULL)
                {
                  traverseNodes (theNode->next->back, timevector, slice, world, tips);
                }
              if (theNode->type != 'm' && theNode->type != 'd' && theNode->next->next->back != NULL)
                {
                  traverseNodes (theNode->next->next->back, timevector, slice, world, tips);
                }
              if (theNode->top)
                {
                  /*
                   * Here we are on the save side if we increase the
                   * timelist so never a fence-write can happen
                   */
                  if (*slice >= (*timevector)->allocT-1)
                    {
                      increase_timelist (timevector);
                    }
                  /*
                   * (*timevector)->tl[*slice].pop =
                   * theNode->actualpop;
                   */
                  (*timevector)->tl[*slice].age = theNode->tyme;
                  //mark the visited node this will not work with recycled node pointers
                  //does it work at all?
                  //if(theNode != (*timevector)->tl[*slice].eventnode)
                  //  (*timevector)->tl[*slice].visited = TRUE;
                  //else
                  //  (*timevector)->tl[*slice].visited = FALSE;
                  //
                  (*timevector)->tl[*slice].eventnode = theNode;
                  (*timevector)->tl[*slice].slice = *slice;
                  (*timevector)->tl[*slice].from = theNode->pop;
                  (*timevector)->tl[*slice].to = theNode->actualpop;
                  (*slice) += 1;
                }
              else
                {
                  error("traverseNodes expects to look only at TOP nodes, but received another one and died\n");
                }
            }
          else
            {
              //tipdate = find_tipdate(theNode->nayme, theNode->pop, world);
              //theNode->tyme = tipdate;
              if (*slice >= (*timevector)->allocT-1)
                {
                  increase_timelist (timevector);
                }
              (*timevector)->tl[*slice].age = theNode->tyme;
              (*timevector)->tl[*slice].eventnode = theNode;
              (*timevector)->tl[*slice].slice = *slice;
              (*timevector)->tl[*slice].from = theNode->pop;
              (*timevector)->tl[*slice].to = theNode->actualpop;
              (*tips) += 1;
              (*slice) += 1;
            }
      }
    else
        error("no node?????????");
}

void
increase_timelist2 (timelist_fmt ** timevector,long extender)
{
  (*timevector)->oldT = (*timevector)->allocT;
  (*timevector)->allocT += extender;
  (*timevector)->tl = (vtlist *) myrealloc ((*timevector)->tl,
					    sizeof (vtlist) * (size_t) ((*timevector)->allocT + 1));
  memset ((*timevector)->tl + (*timevector)->oldT, 0,
	  (size_t) ((*timevector)->allocT - (*timevector)->oldT) * sizeof (vtlist));
  allocate_lineages (timevector, 0, (*timevector)->numpop);
}

void
increase_timelist (timelist_fmt ** timevector)
{
  (*timevector)->oldT = (*timevector)->allocT;
  increase_timelist2(timevector, (*timevector)->allocT / 4); /* increase timelist by 25%*/
}

void set_all_dirty (const node * root, node * p, world_fmt * world, const long locus)
{
//  node *left;
//  node *right;
  if (p->type == 'm' || p->type == 'd')
    error ("MIGRATION/SPECIATION NODE IN SMOOTH FOUND, PLEASE REPORT!!!!\n");
  if (p->type == 'i')
    {
        set_all_dirty (root, /*right=*/crawlback (p->next), world, locus);
        set_all_dirty (root, /*left=*/crawlback (p->next->next), world, locus);
        p->dirty = TRUE;
#ifdef BEAGLE_DEBUG
	printf("(INT %li, %li, %f)-->%li\n",p->id, p->bid, p->v, showtop(crawlback(p))->id);
#endif
    }
  else
    {
#ifdef BEAGLE_DEBUG
	printf("(TIP %li, %li, %f)-->%li\n",p->id, p->bid, p->v, showtop(crawlback(p))->id);
#endif
    }
}    /* set_all_dirty */

// multibranch safe [next->next]
void smooth (const node * root, node * p, world_fmt * world, const long locus)
{
  node *left;
  node *right;
  if (!p->dirty)
    return;
  
  if (p->type == 'm' || p->type == 'd')
    error ("MIGRATION/SPECIATION NODE IN SMOOTH FOUND, PLEASE REPORT!!!!\n");
  if (p->type == 'i')
    {
      if (!p->top)
	error("smooth: does not start at top nodelet");
      right = p->next;
      while(!right->top)
	{
	  left = crawlback (right);
	  smooth (root, left, world, locus);
	  right = right->next;
	}
#ifdef BEAGLE
      //      prepare_beagle_instances(p,left, right, world->beagle);
      error("smooth does not work");
      //	printf("(PINT %li, %li, %f)-->%li\n",p->id, p->bid, p->v, showtop(crawlback(p))->id);
#else
	const long sublocistart = world->sublocistarts[locus];
	const long sublociend   = world->sublocistarts[locus+1];
	long sublocus;
	for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
	  {
	    mutationmodel_fmt *s = &world->mutationmodels[sublocus];
	    const long xs = sublocus - sublocistart;
	    (*nuview[xs]) (s, sublocus, xs, p, world, locus);
	  }
#endif
#ifdef UEP		
        if(world->options->uep)
            twostate_nuview (p, world, locus);
#endif		
        p->dirty = FALSE;
    }
    if (p->type == 't')
      p->dirty = FALSE;
}    /* smooth */


/// sets the nuview machinery so that depending on datatype the correct
/// conditional likelihood calculator is used
void
which_nuview (world_fmt *world, boolean fastlike, boolean use_gaps, int watkins)
{  
  long sublocus;
  long locus = world->locus;
  mutationmodel_fmt *s;
  const long sublocistart = world->sublocistarts[locus];
  const long sublociend   = world->sublocistarts[locus+1];
#ifdef DEBUG
  //printf("Check nuview allocation: numsubloci[locus=%li]=%li\n",world->locus, world->numsubloci[world->locus]);
#endif
  if(nuview==NULL)
    {
      nuview = (nuview_function *) mycalloc(world->numsubloci[world->locus], sizeof(nuview_function));
      pseudonuv = (pseudonuview_function *) mycalloc(world->numsubloci[world->locus], sizeof(pseudonuview_function));
    }
  else
    {
      nuview = (nuview_function *)myrealloc(nuview, (size_t) world->numsubloci[world->locus] *  sizeof(nuview_function));
      pseudonuv = (pseudonuview_function *)realloc(pseudonuv, (size_t) world->numsubloci[world->locus] * sizeof(pseudonuview_function));
    }
  if(prob_micro==NULL)
    {
      prob_micro = (prob_micro_function *) mycalloc(world->numsubloci[world->locus], sizeof(prob_micro_function)); 
    }
  else
    {
      prob_micro = (prob_micro_function *) myrealloc(prob_micro, (size_t) world->numsubloci[world->locus] * sizeof(prob_micro_function)); 
    }
  
  for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
    {
      s = &world->mutationmodels[sublocus];
      const long xs = sublocus - sublocistart;
      switch (s->datatype)
	{
	case 'a':
	  nuview[xs] = (void (*)(mutationmodel_fmt *, long, long, node *, world_fmt *, long)) nuview_allele;
	  pseudonuv[xs] = (void (*)(mutationmodel_fmt *, proposal_fmt *, xarray_fmt *, MYREAL *, MYREAL , xarray_fmt *, MYREAL *, MYREAL, long)) pseudonu_allele; 
	  break;
	case 'b':
	  nuview[xs] = (void (*)(mutationmodel_fmt *, long,long, node *, world_fmt *, long)) nuview_brownian;
	  pseudonuv[xs] = (void (*)(mutationmodel_fmt *, proposal_fmt *, xarray_fmt *, MYREAL *, MYREAL , xarray_fmt *, MYREAL *, MYREAL, long)) pseudonu_brownian; 
	  break;
	case 'm':
	  switch(watkins)
	    {
	    case MULTISTEP:
	      prob_micro[xs] = (double (*)(MYREAL, long, world_fmt *, mutationmodel_fmt *, pair *)) prob_micro_watkins;
	      break;
	    case SINGLESTEP:  
	    default:
	      prob_micro[xs] = (double (*)(MYREAL, long, world_fmt *, mutationmodel_fmt *, pair *)) prob_micro_singlestep;
	    }
	  nuview[xs] = (void (*)(mutationmodel_fmt *, long,long, node *, world_fmt *, long)) nuview_micro;
	  pseudonuv[xs] = (void (*)(mutationmodel_fmt *, proposal_fmt *, xarray_fmt *, MYREAL *, MYREAL , xarray_fmt *, MYREAL *, MYREAL, long)) pseudonu_micro; 
	  break;
	case 's':
	case 'n':
	case 'h':
	case 'u':
	  switch(s->model)
	    {
	    case JC69:
	      force_basefreqs(&s->basefreqs,0.25,0.25,0.25);
	      s->parameters[0]=1.0;
	      s->parameters[1]=1.0;
	      s->parameters[2]=1.0;
#ifdef DEBUG
	      printf("JC: 0:%f 1:%f 2:%f \n",s->parameters[0],s->parameters[1],s->parameters[2]);
#endif	      

	      nuview[xs] = (void (*)(mutationmodel_fmt *, long, long, node *, world_fmt *, long)) nuview_tn93;
	      pseudonuv[xs] = (void (*)(mutationmodel_fmt *, proposal_fmt *, xarray_fmt *, MYREAL *, MYREAL , xarray_fmt *, MYREAL *, MYREAL, long)) pseudonu_tn93; 
	      break;
	    case K2P:
	      force_basefreqs(&s->basefreqs,0.25,0.25,0.25);
	      //s->parameters[0]=s->ttratio;
	      //s->parameters[1]=1.0;
	      //s->parameters[2]=1.0;
#ifdef DEBUG
	      printf("K2P: 0:%f 1:%f 2:%f \n",s->parameters[0],s->parameters[1],s->parameters[2]);
#endif	      
	      nuview[xs] = (void (*)(mutationmodel_fmt *, long, long, node *, world_fmt *, long)) nuview_tn93;
	      pseudonuv[xs] = (void (*)(mutationmodel_fmt *, proposal_fmt *, xarray_fmt *, MYREAL *, MYREAL , xarray_fmt *, MYREAL *, MYREAL, long)) pseudonu_tn93; 
	      break;
	    case F81:
	      s->parameters[0]=1.0;
	      s->parameters[1]=1.0;
	      s->parameters[2]=1.0;
#ifdef DEBUG
	      printf("F81: 0:%f 1:%f 2:%f \n",s->parameters[0],s->parameters[1],s->parameters[2]);
#endif	      
	      nuview[xs] = (void (*)(mutationmodel_fmt *, long, long, node *, world_fmt *, long)) nuview_tn93;
	      pseudonuv[xs] = (void (*)(mutationmodel_fmt *, proposal_fmt *, xarray_fmt *, MYREAL *, MYREAL , xarray_fmt *, MYREAL *, MYREAL, long)) pseudonu_tn93; 
	      break;
	    case F84:
	      //s->parameters[0]=1.0;
	      //s->parameters[1]=1.0;
	      // done in model s->parameters[2]=1.0;
#ifdef DEBUG
	      printf("F84: 0:%f 1:%f 2:%f \n",s->parameters[0],s->parameters[1],s->parameters[2]);
#endif	      
	      nuview[xs] = (void (*)(mutationmodel_fmt *, long, long, node *, world_fmt *, long)) nuview_tn93;
	      pseudonuv[xs] = (void (*)(mutationmodel_fmt *, proposal_fmt *, xarray_fmt *, MYREAL *, MYREAL , xarray_fmt *, MYREAL *, MYREAL, long)) pseudonu_tn93; 
	      break;
	    case HKY:
	      //s->parameters[0]=1.0;
	      //s->parameters[1]=1.0;
	      //s->parameters[2]=1.0;
#ifdef DEBUG
	      printf("HKY: 0:%f 1:%f 2:%f \n",s->parameters[0],s->parameters[1],s->parameters[2]);
#endif	      
	      nuview[xs] = (nuview_function) nuview_tn93;
	      pseudonuv[xs] = (pseudonuview_function) pseudonu_tn93; 
	      break;
	    case TN:
#ifdef DEBUG
	      printf("TN: 0:%f 1:%f 2:%f \n",s->parameters[0],s->parameters[1],s->parameters[2]);
#endif	      
	      nuview[xs] = (nuview_function) nuview_tn93;
	      pseudonuv[xs] = (pseudonuview_function) pseudonu_tn93; 
	      break;
	    default:
	      if (fastlike && !use_gaps)
		{
		  nuview[xs] = (nuview_function) nuview_sequence;
		  pseudonuv[xs] = (pseudonuview_function) pseudonu_seq; 
		}
	      else
		{
		  if(fastlike)
		    {
		      nuview[xs] = (nuview_function) nuview_sequence_slow;
		      pseudonuv[xs] = (pseudonuview_function) pseudonu_seq_slow;
		    }
		  else /*use_gaps*/
		    {
		      nuview[xs] = (nuview_function) nuview_sequence_slow; //nuview_gaps;	
		      pseudonuv[xs] = (pseudonuview_function) pseudonu_seq_slow; 
		    }
		}
	      break;
	    case 'f':   /* fitch, reconstruction of ancestral state
			 * method */
	      nuview[xs] = (nuview_function) nuview_ancestral;
	      pseudonuv[xs] = (pseudonuview_function) pseudonu_anc; 
	    }
	}
    }
}

void
nuview_allele (mutationmodel_fmt *s, long sublocus, long xs, node * mother, world_fmt * world, const long locus)
{
  (void) s;
  (void) sublocus;
  long a;
  long aa;
  long mal = world->data->maxalleles[locus];

  MYREAL freq = world->data->freq;
  //  MYREAL freqlast = world->data->freqlast; see pseudonuview
  MYREAL w1;
  MYREAL w2;
  MYREAL v1;
  MYREAL v2;
  MYREAL v1freq;
  MYREAL v2freq;
  MYREAL pija1;
  MYREAL pija2;
  MYREAL lx1;
  MYREAL lx2;

  MYREAL x3m = -MYREAL_MAX;

  MYREAL *xx1, *xx2;
  MYREAL *xx3;

  node *d1 = NULL; 
  node *d2 = NULL;

  //  MYREAL test = 0.0;
  
  children (mother, &d1, &d2);
  xx1 = d1->x[xs].a;
  xx2 = d2->x[xs].a;
  xx3 = mother->x[xs].a;
  lx1 = d1->scale[xs][0];
  lx2 = d2->scale[xs][0];
  v1 = 1 - EXP (-d1->v);
  v2 = 1 - EXP (-d2->v);
    if (v1 >= 1.)
    {
        w1 = 0.0;
        v1 = 1.0;
    }
    else
    {
        w1 = 1.0 - v1;
    }
    if (v2 >= 1.)
    {
        w2 = 0.0;
        v2 = 1.0;
    }
    else
    {
        w2 = 1.0 - v2;
    }
    //    printf("@@nuview 1:{%f,%f, %f,%f}   2:(%f,%f, %f,%f}\n",xx1[0],xx1[1],lx1,v1,xx2[0],xx2[1],lx2,v2);

    v1freq = v1 * freq;
    v2freq = v2 * freq;
    
    for (aa = 0; aa < mal; aa++)
    {
        pija1 = pija2 = 0.0;
        for (a = 0; a < mal; a++)
        {
	  if(aa==a)
	    {
	      pija1 += (w1 + v1freq) * xx1[a];
	      pija2 += (w2 + v2freq) * xx2[a];
	    }
	  else
	    {
	      pija1 += v1 * freq * xx1[a];
	      pija2 += v2 * freq * xx2[a];
	    }
	}
        xx3[aa] = pija1 * pija2;
        //test += xx3[aa];
        if (xx3[aa] > x3m)
            x3m = xx3[aa];
		
    }
    //if (test <= 0.0)
    //    error ("xx3 is 0 or garbage!");
    for (aa = 0; aa < mal; aa++)
    {
        xx3[aa] /= x3m;
    }
    //    printf("@@finish: (%f=%f) (%f=%f)\n",xx3[0],mother->x.a[0],xx3[1],mother->x.a[1]);
    mother->scale[xs][0] = LOG (x3m) + lx2 + lx1;
}

void
nuview_brownian (mutationmodel_fmt *s, long sublocus, long xs, node * mother, world_fmt * world, const long locus)
{
  (void) s;
  (void) sublocus;
  (void) world;
  (void) locus;
    node *d1 = NULL, *d2 = NULL;
    MYREAL rvtot = 0.;
    MYREAL xx1, xx2, c12, diff;
    MYREAL mean1, mean2, mean, v1, v2, vtot, f1, f2;
    xarray_fmt x1,x2;
    children (mother, &d1, &d2);
    x1 = d1->x[xs];
    x2 = d2->x[xs];
    mean1 = x1.a[0];
    mean2 = x2.a[0];
    xx1 = x1.a[2];
    xx2 = x2.a[2];
	
    v1 = d1->v + x1.a[1]; /* di->v == length of branch time(n1) - * time(n2) */
    v2 = d2->v + x2.a[1]; /* x.a[1] contains the deltav */
    vtot = v1 + v2;
    diff = mean1 - mean2;
    if (vtot > 0.0)
      {
	rvtot = 1./ vtot;
        f1 = v2 * rvtot;
	c12 = diff * diff * rvtot;
      }
    else
      {
        f1 = 0.5;
	c12 = (double) HUGE;
	vtot = SMALL_VALUE;
      }
    f2 = 1.0 - f1;
    mean = f1 * mean1 + f2 * mean2;


    mother->x[xs].a[2] =
        xx1 + xx2 + MIN (0.0, -0.5 * (LOG (vtot) + c12) + LOG2PIHALF);
    /*
     * printf("L=%f , L1=%f, L2=%f, log(vtot=%f)=%f,
     * c12=%f\n",mother->x.a[2], xx1, xx2,vtot,log(vtot),c12);
     */
    mother->x[xs].a[1] = v1 * f1;
    mother->x[xs].a[0] = mean;
	
}


void
nuview_micro (mutationmodel_fmt *s, long sublocus, long xs, node * mother, world_fmt * world, const long locus)
{
  (void) sublocus;
  (void) locus;

    node *d1 = NULL;
    node *d2 = NULL;
    long ss;
    long a;
    long aa1, aa2;
    long diff;
    long margin = s->micro_threshold;
    MYREAL vv1, vv2, lx1, lx2;
    MYREAL x3m = -MYREAL_MAX;
    MYREAL pija1s, pija2s;
    MYREAL *xx1 = NULL, *xx2 = NULL;

    MYREAL *pm1;
    MYREAL *pm2;

    long smax = s->maxalleles;
    MYREAL *xx3 = NULL;
    // needed for watkins' microsatellite model
    pair *helper = &world->options->msat_tuning;

    children (mother, &d1, &d2);
    vv1 = d1->v;
    vv2 = d2->v;
    xx1 = d1->x[xs].a;
    xx2 = d2->x[xs].a;
    xx3 = mother->x[xs].a;
    lx1 = d1->scale[xs][0];
    lx2 = d2->scale[xs][0];

    pm1 = (MYREAL *) mymalloc(sizeof(MYREAL) * (size_t) (4 * margin));
    pm2 = pm1 +  2 * margin;

    for (diff = 0; diff < margin; diff++)
      {
	pm1[diff] = (*prob_micro[xs]) (vv1, diff, world, s, helper);
	pm2[diff] = (*prob_micro[xs]) (vv2, diff, world, s, helper);
	//fprintf(stdout,"pm[diff=%li]=(%g %g)\n",diff,pm1[diff],pm2[diff]);
      }

    for (ss = 0; ss < smax; ss++)
    {
        pija1s = pija2s = 0.0;
	aa1 = MAX (0, ss - margin);
	aa2 = MIN(ss + margin,smax);
	for (a = aa1; a < aa2; a++)
	  {
            diff = labs (ss - a);
	    if(diff>=margin)
	      continue;
	    //	    fprintf(stdout,"***pm[diff=%li]=(%g %g) xx[a=%li]=(%f %f)\n",diff,pm1[diff],pm2[diff],a,xx1[a],xx2[a]);
            if (xx1[a] > 0)
            {
                pija1s += pm1[diff] * xx1[a];
            }
            if (xx2[a] > 0)
            {
                pija2s += pm2[diff] * xx2[a];
            }
        }
        xx3[ss] = pija1s * pija2s;
        if (xx3[ss] > x3m)
            x3m = xx3[ss];
    }
    if (x3m == 0.0)
    {
        mother->scale[xs][0] = -MYREAL_MAX;
    }
    else
    {
        for (ss = 0; ss < smax; ss++)
        {
            xx3[ss] /= x3m;
        }
        mother->scale[xs][0] = LOG (x3m) + lx1 + lx2;
    }
    myfree(pm1);
}


///
/// calculates probability of a mutation of size diff in time t. used for msat calculations.
/// i=diff
/// prob[i_, t_] := Exp[-t] Sum[(t/2)^(i + 2 k)/((i + k)! k!), {k, 0, Infinity}]
/// calculated as Sum[Exp[(-t + ((log(t)-log(2))*i) + (2*(log(t)-log(2)) * k) - steps(i,k))
MYREAL
prob_micro_singlestep(MYREAL t, long diff, world_fmt * world, mutationmodel_fmt *s, pair *helper)
{
  (void) world;
  (void) helper;
  // helper is a dummy so that the call is the same as for prob_micro_watkins
  const long stepnum = s->micro_threshold;
  //const long locus = world->locus;
  const MYREAL *steps = s->steps[diff];
  long k;
  long k2;

  MYREAL temp1;
  MYREAL temp2;
  MYREAL temp3;
  MYREAL temp4;
  MYREAL const_part;
  MYREAL summ = 0.0;
  // MYREAL oldsum = 0.0;
  const MYREAL logt = LOG (t) - LOG2;
  const MYREAL logt2 = 2 * logt;
  if (diff >= stepnum)
    return summ;
  // linking to precalculated "denominator" Log[(i+k)! k!]
  // was precalculated in calculate_steps()
  //  steps = world->options->steps[locus][diff];
  const_part = -t + logt * diff;
  // loop unrolling to remove possible stalls
  // assumes stepnum is even
  for (k = 0; k < stepnum; k += 2)
    {
      k2 = k + 1;
      temp1 = const_part + logt2 * k - steps[k];
      temp2 = const_part + logt2 * k2 - steps[k2];
      temp3 = EXP(temp1);
      temp4 = EXP(temp2);
      summ += temp3 + temp4;
      //      if (fabs (oldsum - summ) < eps)
      //	break;
      //oldsum = summ;
    }
    return summ;
}

///
/// Calculates probability of a change of repeat number in time. Used for msat calculations.
/// based on Joe Watkins' 2007 in theoretical population biology, chapter 4.2
/// prob[diff,t,tune, upchance]
/// diff: is negative or positive depending on direction
/// t: is the time
/// tune:    0.0 with the single step mutation model
///          1.0 with the infinite allele model
///          all values in between allow multiple steps   
/// upchance: is the probability that the repeat number increases (upchance < 2/3!)
/// helper is the container for the additional variables
MYREAL
prob_micro_watkins (MYREAL t, long diff, world_fmt * world, mutationmodel_fmt *s, pair *helper)
{
  (void) world;
  (void) s;
  const MYREAL tune = (*helper)[0];
  const MYREAL upchance = (*helper)[1]; 
  MYREAL x;
  MYREAL summ = 0.0;
  const MYREAL delta = 2. * PI / 100.;
  const MYREAL expt = EXP(-t);
  MYREAL oneplustunesq;
  MYREAL tunep;
  MYREAL oneminustunet; 
  MYREAL sinx;
  MYREAL cosx;
  MYREAL invdenom;
  MYREAL first;
  MYREAL second;

  //if (diff >= stepnum)
  //  return summ;

  oneplustunesq = 1 + tune * tune;
  tunep = (1.0 - tune) * (2.0 * upchance - 1.0) * t; 
  oneminustunet = (1.0 - tune) * t;
  for (x = -PI; x < PI; x += delta)
    {
      cosx = cos(x);
      sinx = sin(x);
      invdenom = (oneplustunesq - 2. * tune * cosx);
      if(invdenom < SMALL_VALUE)
	invdenom = 1. / SMALL_VALUE;
      first = (diff * x + tunep * sinx) * invdenom;
      second = oneminustunet * cosx * invdenom;
      summ += cos(first) * EXP(second);
    }
#ifdef DEBUG
  if(summ == 0.0)
    error("underflow in msat prob calc\n");
#endif
  return expt * INV2PI * delta * summ;
}

///
/// test implementation of Watkins' calculations
//MYREAL
//prob_micro (MYREAL t, long diff, world_fmt * world)
//{
//  return prob_micro_watkins (t, diff, 0.2 , 0.5, world);
//}
//MYINLINE MYREAL
//prob_micro (MYREAL t, long diff, world_fmt * world)
//{
//  return prob_micro_singlestep (t, diff, world);
//}


///
/// conditional likelihoods
void
nuview_sequence (mutationmodel_fmt *s, long sublocus, long xs, node * mother, world_fmt * world, const long locus)
{
  (void) sublocus;
  (void) world;
  (void) locus;
  long endsite     = s->numpatterns + s->addon;
  long i, j;
#ifndef GRANDCENTRAL_notworking /*replace SNOWLEOPARD*/
    long k;
    register sitelike *xx1, *xx2, *xx3;
    register MYREAL xx1t0, xx1t1,xx1t2, xx1t3;
    register MYREAL xx2t0, xx2t1, xx2t2, xx2t3;
    register MYREAL vzsumr1, vzsumr2, vzsumy1, vzsumy2, sum1, sum2, sumr1, sumr2,
      sumy1, sumy2; 
    register valrec *tbljk;
#endif
    register MYREAL ww1, ww2, zz1, zz2, lw1, lw2, yy1, yy2, ww1zz1, vv1zz1, ww2zz2, vv2zz2;
    register MYREAL freqa, freqc, freqg, freqt;// freqr, freqy;
    register MYREAL freqar, freqcy, freqgr, freqty;
    register node *q, *r;
    register long rcategs = s->numsiterates;
    register long categs = s->numcategs;
    register tbl_fmt tbl = s->tbl;

    register valrec *tbl00;
    register valrec *tblij;
    register MYREAL fracchange = s->fracchange;
    q = crawlback (mother->next);
    r = crawlback (mother->next->next);

    freqa = s->basefreqs[NUC_A];
    freqc = s->basefreqs[NUC_C];
    freqg = s->basefreqs[NUC_G];
    freqt = s->basefreqs[NUC_T];
//    freqr = s->basefreqs[NUC_R];
//    freqy = s->basefreqs[NUC_Y];
    freqar = s->basefreqs[NUC_AR];
    freqcy = s->basefreqs[NUC_CY];
    freqgr = s->basefreqs[NUC_GR];
    freqty = s->basefreqs[NUC_TY];
    
    lw1 = -q->v * s->fracchange;
    
    if ((rcategs | categs) == 1)
    {
        tbl00 = tbl[0][0];
        ww1 = EXP (tbl00->ratxi * lw1);
        zz1 = EXP (tbl00->ratxv * lw1);
        ww1zz1 = ww1 * zz1;
        vv1zz1 = (1.0 - ww1) * zz1;
        lw2 = -r->v * fracchange;
        ww2 = EXP (tbl00->ratxi * lw2);
        zz2 = EXP (tbl00->ratxv * lw2);
        ww2zz2 = ww2 * zz2;
        vv2zz2 = (1.0 - ww2) * zz2;
        yy1 = 1.0 - zz1;
        yy2 = 1.0 - zz2;
#ifdef GRANDCENTRAL_notworking /*replace SNOWLEOPARD*/
	dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
	dispatch_apply(endsite, queue,
		       ^(unsigned long i) {
			 register sitelike *xx1, *xx2, *xx3;
			 register MYREAL xx1t0, xx1t1,xx1t2, xx1t3;
			 register MYREAL xx2t0, xx2t1, xx2t2, xx2t3;
			 register MYREAL vzsumr1, vzsumr2, vzsumy1, vzsumy2, sum1, sum2, sumr1, sumr2,
			   sumy1, sumy2;
#else
	for (i = 0; i < endsite; i++)
	  {
#endif
            xx1 = &(q->x[xs].s[i][0]);
            xx2 = &(r->x[xs].s[i][0]);
            xx3 = &(mother->x[xs].s[i][0]);
			
            xx1t0 = (*xx1)[0];
            xx1t1 = (*xx1)[1];
            xx1t2 = (*xx1)[2];
            xx1t3 = (*xx1)[3];
			
            xx2t0 = (*xx2)[0];
            xx2t1 = (*xx2)[1];
            xx2t2 = (*xx2)[2];
            xx2t3 = (*xx2)[3];
			
            sum1 = yy1 * (freqa * xx1t0 + freqc * xx1t1 + freqg * xx1t2 + freqt * xx1t3);
            sum2 = yy2 * (freqa * xx2t0 + freqc * xx2t1 + freqg * xx2t2 + freqt * xx2t3);
            
            sumr1 = freqar * xx1t0 + freqgr * xx1t2;
            sumr2 = freqar * xx2t0 + freqgr * xx2t2;
            sumy1 = freqcy * xx1t1 + freqty * xx1t3;
            sumy2 = freqcy * xx2t1 + freqty * xx2t3;
            
            vzsumr1 = vv1zz1 * sumr1;
            vzsumr2 = vv2zz2 * sumr2;
            vzsumy1 = vv1zz1 * sumy1;
            vzsumy2 = vv2zz2 * sumy2;
            (*xx3)[0] =
                (sum1 + ww1zz1 * xx1t0 + vzsumr1) * (sum2 +
													 ww2zz2 * xx2t0 +
													 vzsumr2);
            (*xx3)[1] =
                (sum1 + ww1zz1 * xx1t1 + vzsumy1) * (sum2 +
													 ww2zz2 * xx2t1 +
													 vzsumy2);
            (*xx3)[2] =
                (sum1 + ww1zz1 * xx1t2 + vzsumr1) * (sum2 +
													 ww2zz2 * xx2t2 +
													 vzsumr2);
            (*xx3)[3] =
                (sum1 + ww1zz1 * xx1t3 + vzsumy1) * (sum2 +
													 ww2zz2 * xx2t3 +
													 vzsumy2);
	  }
#ifdef GRANDCENTRAL_notworking /*replace SNOWLEOPARD*/
	); // end of block
#endif
    }
    else
    {
        for (i = 0; i < rcategs; i++)
            for (j = 0; j < categs; j++)
            {
                tblij = tbl[i][j];
                tblij->ww1 = EXP (tblij->ratxi * lw1);
                tblij->zz1 = EXP (tblij->ratxv * lw1);
                tblij->ww1zz1 = tblij->ww1 * tblij->zz1;
                tblij->vv1zz1 = (1.0 - tblij->ww1) * tblij->zz1;
            }
				lw2 = -r->v * s->fracchange;
        for (i = 0; i < rcategs; i++)
            for (j = 0; j < categs; j++)
            {
                tblij = tbl[i][j];
                tblij->ww2 = EXP (tblij->ratxi * lw2);
                tblij->zz2 = EXP (tblij->ratxv * lw2);
                tblij->ww2zz2 = tblij->ww2 * tblij->zz2;
                tblij->vv2zz2 = (1.0 - tblij->ww2) * tblij->zz2;
            }
#ifdef GRANDCENTRAL_notworking /*replace SNOWLEOPARD*/
	dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
	dispatch_apply(endsite, queue,
		       ^(unsigned long i) {
			 long j, k;
			 register sitelike *xx1, *xx2, *xx3;
			 register MYREAL xx1t0, xx1t1,xx1t2, xx1t3;
			 register MYREAL xx2t0, xx2t1, xx2t2, xx2t3;
			 register MYREAL vzsumr1, vzsumr2, vzsumy1, vzsumy2, sum1, sum2, sumr1, sumr2,
			   sumy1, sumy2;
			 valrec * tbljk;
			 MYREAL ww1zz1, vv1zz1,yy1,ww2zz2,vv2zz2,yy2;
#else
	for (i = 0; i < endsite; i++)
	  {
#endif
	    k = s->category[s->alias[i] - 1] - 1;
	    for (j = 0; j < rcategs; j++)
	      {
		tbljk = tbl[j][k];
		ww1zz1 = tbljk->ww1zz1;
		vv1zz1 = tbljk->vv1zz1;
		yy1 = 1.0 - tbljk->zz1;
		ww2zz2 = tbljk->ww2zz2;
		vv2zz2 = tbljk->vv2zz2;
		yy2 = 1.0 - tbljk->zz2;
		xx1 = &(q->x[xs].s[i][j]);
		xx2 = &(r->x[xs].s[i][j]);
		xx3 = &(mother->x[xs].s[i][j]);
		
		xx1t0 = (*xx1)[0];
		xx1t1 = (*xx1)[1];
		xx1t2 = (*xx1)[2];
		xx1t3 = (*xx1)[3];
		
		xx2t0 = (*xx2)[0];
		xx2t1 = (*xx2)[1];
		xx2t2 = (*xx2)[2];
		xx2t3 = (*xx2)[3];
		
						
		sum1 =
		  yy1 * (freqa * xx1t0 + freqc * xx1t1 +
			 freqg * xx1t2 + freqt * xx1t3);
		sum2 =
		  yy2 * (freqa * xx2t0 + freqc * xx2t1 +
			 freqg * xx2t2 + freqt * xx2t3);
		sumr1 = freqar * xx1t0 + freqgr * xx1t2;
		sumr2 = freqar * xx2t0 + freqgr * xx2t2;
		sumy1 = freqcy * xx1t1 + freqty * xx1t3;
		sumy2 = freqcy * xx2t1 + freqty * xx2t3;
		vzsumr1 = vv1zz1 * sumr1;
		vzsumr2 = vv2zz2 * sumr2;
		vzsumy1 = vv1zz1 * sumy1;
		vzsumy2 = vv2zz2 * sumy2;
		(*xx3)[0] =
		  (sum1 + ww1zz1 * xx1t0 + vzsumr1) * (sum2 +
						       ww2zz2 * xx2t0 +
						       vzsumr2);
		(*xx3)[1] =
		  (sum1 + ww1zz1 * xx1t1 + vzsumy1) * (sum2 +
						       ww2zz2 * xx2t1 +
						       vzsumy2);
		(*xx3)[2] =
		  (sum1 + ww1zz1 * xx1t2 + vzsumr1) * (sum2 +
						       ww2zz2 * xx2t2 +
						       vzsumr2);
		(*xx3)[3] =
		  (sum1 + ww1zz1 * xx1t3 + vzsumy1) * (sum2 +
						       ww2zz2 * xx2t3 +
						       vzsumy2);
	      }
	  }
#ifdef GRANDCENTRAL_notworking /*replace SNOWLEOPARD*/
	);
#endif
    }
}    /* nuview */

///
/// accurate version of conditional likleihood calculator using a scaler
#ifdef GAP
void nuview_sequence_slow (node * mother, world_fmt * world, const long locus)
{
  nuview_f84gap_slow (mother, world, locus);
}
#else
 void nuview_sequence_slow (mutationmodel_fmt *s, long sublocus, long xs, node * mother, world_fmt * world, const long locus)
{
  (void) locus;
  (void) sublocus;
    //static long count = 0;
    const long endsite = s->numpatterns + s->addon;
    long i, j;
    //#ifndef GRANDCENTRAL_notworking /*replace SNOWLEOPARD*/
    long k;
    register sitelike *xx1, *xx2, *xx3;
    register MYREAL xx1t0, xx1t1,xx1t2, xx1t3;
    register MYREAL xx2t0, xx2t1, xx2t2, xx2t3;
    register MYREAL vzsumr1, vzsumr2, vzsumy1, vzsumy2, sum1, sum2, sumr1, sumr2,
      sumy1, sumy2; 
    MYREAL sxx3m, tempsxx3m, invsxx3m;
    valrec *tbljk;
    //#endif
    register MYREAL lw1, lw2, yy1, yy2, ww1zz1, vv1zz1, ww2zz2, vv2zz2, ww1, ww2, zz1, zz2;
    MYREAL *sxx1 = NULL, *sxx2 = NULL;

    node *q, *r;

    long rcategs = world->options->rcategs;
    long categs = world->options->categs;
    tbl_fmt tbl = world->tbl;
    register MYREAL freqa, freqc, freqg, freqt;//, freqr, freqy;
    register MYREAL freqar, freqcy, freqgr, freqty;

    valrec *tbl00;
    valrec *tblij;
    freqa = s->basefreqs[NUC_A];
    freqc = s->basefreqs[NUC_C];
    freqg = s->basefreqs[NUC_G];
    freqt = s->basefreqs[NUC_T];
//    freqr = s->basefreqs[NUC_R];
//    freqy = s->basefreqs[NUC_Y];
    freqar = s->basefreqs[NUC_AR];
    freqcy = s->basefreqs[NUC_CY];
    freqgr = s->basefreqs[NUC_GR];
    freqty = s->basefreqs[NUC_TY];
    MYREAL scalesum = 0.0;
    boolean newscale=FALSE;;

    q = crawlback (mother->next);
    r = crawlback (mother->next->next);
    lw1 = -q->v * s->fracchange;
    sxx1 = q->scale[xs];
    sxx2 = r->scale[xs];
    tbl00 = tbl[0][0];
    if ((rcategs | categs) == 1)
    {
        ww1 = EXP (tbl00->ratxi * lw1);
        zz1 = EXP (tbl00->ratxv * lw1);
        ww1zz1 = ww1 * zz1;
        vv1zz1 = (1.0 - ww1) * zz1;
        lw2 = -r->v * s->fracchange;
        ww2 = EXP (tbl00->ratxi * lw2);
        zz2 = EXP (tbl00->ratxv * lw2);
        ww2zz2 = ww2 * zz2;
        vv2zz2 = (1.0 - ww2) * zz2;
        yy1 = 1.0 - zz1;
        yy2 = 1.0 - zz2;
	//#ifdef GRANDCENTRAL_notworking /*replace SNOWLEOPARD*/
	//dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
	//dispatch_apply(endsite, queue,
	//	       ^(unsigned long i) {
	//		 register sitelike *xx1, *xx2, *xx3;
	//		 register MYREAL xx1t0, xx1t1,xx1t2, xx1t3;
	//		 register MYREAL xx2t0, xx2t1, xx2t2, xx2t3;
	//		 register MYREAL vzsumr1, vzsumr2, vzsumy1, vzsumy2, sum1, sum2, sumr1, sumr2,
	//		   sumy1, sumy2;
	//#else
        for (i = 0; i < endsite; i++)
        {
	  //#endif
	  xx1 = &(q->x[xs].s[i][0]);
	  xx2 = &(r->x[xs].s[i][0]);
	  xx3 = &(mother->x[xs].s[i][0]);
	  
	  xx1t0 = (*xx1)[0];
	  xx1t1 = (*xx1)[1];
	  xx1t2 = (*xx1)[2];
	  xx1t3 = (*xx1)[3];
          
	  xx2t0 = (*xx2)[0];
	  xx2t1 = (*xx2)[1];
	  xx2t2 = (*xx2)[2];
	  xx2t3 = (*xx2)[3];
          
            sum1 = yy1 * (freqa * xx1t0 + freqc * xx1t1 + freqg * xx1t2 + freqt * xx1t3);
            sum2 = yy2 * (freqa * xx2t0 + freqc * xx2t1 + freqg * xx2t2 + freqt * xx2t3);
            sumr1 = freqar * xx1t0 + freqgr * xx1t2;
            sumr2 = freqar * xx2t0 + freqgr * xx2t2;
            sumy1 = freqcy * xx1t1 + freqty * xx1t3;
            sumy2 = freqcy * xx2t1 + freqty * xx2t3;
            vzsumr1 = vv1zz1 * sumr1;
            vzsumr2 = vv2zz2 * sumr2;
            vzsumy1 = vv1zz1 * sumy1;
            vzsumy2 = vv2zz2 * sumy2;
            (*xx3)[0] = (sum1 + ww1zz1 * xx1t0 + vzsumr1) * (sum2 + ww2zz2 * xx2t0 + vzsumr2);
	    scalesum = (*xx3)[0];
            (*xx3)[1] =
                (sum1 + ww1zz1 * xx1t1 + vzsumy1) * (sum2 +
                                                     ww2zz2 * xx2t1 +
                                                     vzsumy2);
	    scalesum += (*xx3)[1];
            (*xx3)[2] =
                (sum1 + ww1zz1 * xx1t2 + vzsumr1) * (sum2 +
                                                     ww2zz2 * xx2t2 +
                                                     vzsumr2);
	    scalesum += (*xx3)[2];
            (*xx3)[3] =
                (sum1 + ww1zz1 * xx1t3 + vzsumy1) * (sum2 +
                                                     ww2zz2 * xx2t3 +
                                                     vzsumy2);
	    scalesum += (*xx3)[3];
            mother->scale[xs][i] = sxx1[i] + sxx2[i];
	    if(!(scalesum > EPSILON))
              {
                sxx3m = MAX ((*xx3)[0], (*xx3)[1]);
                sxx3m = MAX (sxx3m, (*xx3)[2]);
                sxx3m = MAX (sxx3m, (*xx3)[3]);
		invsxx3m = 1.0/sxx3m;
                (*xx3)[0] *= invsxx3m;
		(*xx3)[1] *= invsxx3m;
		(*xx3)[2] *= invsxx3m;
		(*xx3)[3] *= invsxx3m;
                mother->scale[xs][i] += LOG (sxx3m);
	      }
#ifdef DEBUG
	    //	printf("[(%f,%f,%f,%f) %f]",(*xx3)[0],(*xx3)[1],(*xx3)[2],(*xx3)[3],mother->scale[i]);
#endif
	    //#ifdef GRANDCENTRAL_notworking /*replace SNOWLEOPARD*/
	    //);
	    //#endif
        }
#ifdef DEBUG
	//   printf("id=%li\n",mother->id);
#endif
    }
    else
    {
      lw2 = -r->v * s->fracchange;
      for (i = 0; i < rcategs; i++)
	for (j = 0; j < categs; j++)
	  {
	    tblij = tbl[i][j];
	    tblij->ww1 = EXP (tblij->ratxi * lw1);
	    tblij->zz1 = EXP (tblij->ratxv * lw1);
	    tblij->ww1zz1 = tblij->ww1 * tblij->zz1;
	    tblij->vv1zz1 = (1.0 - tblij->ww1) * tblij->zz1;
	  }
      for (i = 0; i < rcategs; i++)
	{
	  for (j = 0; j < categs; j++)
	    {
	      tblij = tbl[i][j];
	      tblij->ww2 = EXP (tblij->ratxi * lw2);
	      tblij->zz2 = EXP (tblij->ratxv * lw2);
	      tblij->ww2zz2 = tblij->ww2 * tblij->zz2;
	      tblij->vv2zz2 = (1.0 - tblij->ww2) * tblij->zz2;
	    }
	}
      for (i = 0; i < endsite; i++)
	{
	  newscale = FALSE;
	  //sxx3m = -MYREAL_MAX;
	  k = s->category[s->alias[i] - 1] - 1;
	  for (j = 0; j < rcategs; j++)
	    {
	      tbljk = tbl[j][k];
		ww1zz1 = tbljk->ww1zz1;
		vv1zz1 = tbljk->vv1zz1;
		yy1 = 1.0 - tbljk->zz1;
		ww2zz2 = tbljk->ww2zz2;
		vv2zz2 = tbljk->vv2zz2;
		yy2 = 1.0 - tbljk->zz2;
		xx1 = &(q->x[xs].s[i][j]);
		xx2 = &(r->x[xs].s[i][j]);
		xx3 = &(mother->x[xs].s[i][j]);
                
                
		xx1t0 = (*xx1)[0];
		xx1t1 = (*xx1)[1];
		xx1t2 = (*xx1)[2];
		xx1t3 = (*xx1)[3];
                
		xx2t0 = (*xx2)[0];
		xx2t1 = (*xx2)[1];
		xx2t2 = (*xx2)[2];
		xx2t3 = (*xx2)[3];
                
                
                
		sum1 =
		  yy1 * (freqa * xx1t0 + freqc * xx1t1 +
			 freqg * xx1t2 + freqt * xx1t3);
		sum2 =
		  yy2 * (freqa * xx2t0 + freqc * xx2t1 +
			 freqg * xx2t2 + freqt * xx2t3);
		sumr1 = freqar * xx1t0 + freqgr * xx1t2;
		sumr2 = freqar * xx2t0 + freqgr * xx2t2;
		sumy1 = freqcy * xx1t1 + freqty * xx1t3;
		sumy2 = freqcy * xx2t1 + freqty * xx2t3;
		vzsumr1 = vv1zz1 * sumr1;
		vzsumr2 = vv2zz2 * sumr2;
		vzsumy1 = vv1zz1 * sumy1;
		vzsumy2 = vv2zz2 * sumy2;
		(*xx3)[0] =
		  (sum1 + ww1zz1 * xx1t0 + vzsumr1) * (sum2 +
						       ww2zz2 * xx2t0 +
						       vzsumr2);
		(*xx3)[1] =
		  (sum1 + ww1zz1 * xx1t1 + vzsumy1) * (sum2 +
						       ww2zz2 * xx2t1 +
						       vzsumy2);
		(*xx3)[2] =
		  (sum1 + ww1zz1 * xx1t2 + vzsumr1) * (sum2 +
						       ww2zz2 * xx2t2 +
						       vzsumr2);
		(*xx3)[3] =
		  (sum1 + ww1zz1 * xx1t3 + vzsumy1) * (sum2 +
						       ww2zz2 * xx2t3 +
						       vzsumy2);
		if(!((*xx3)[0]+ (*xx3)[1]+ (*xx3)[2]+ (*xx3)[3]>EPSILON))
		  newscale=TRUE;
	      }
	    if(newscale)
	      {
		sxx3m = -MYREAL_MAX;
                for (j = 0; j < rcategs; j++)
		  {
                    xx3 = &(mother->x[xs].s[i][j]);
                    tempsxx3m = MAX ((*xx3)[0], (*xx3)[1]);
                    tempsxx3m = MAX (tempsxx3m, (*xx3)[2]);
                    tempsxx3m = MAX (tempsxx3m, (*xx3)[3]);
                    if (tempsxx3m > sxx3m)
                        sxx3m = tempsxx3m;
		  }
		invsxx3m = 1.0 / sxx3m;
                for (j = 0; j < rcategs; j++)
                {
		  xx3 = &(mother->x[xs].s[i][j]);
		  (*xx3)[0] *= invsxx3m;
		  (*xx3)[1] *= invsxx3m;
		  (*xx3)[2] *= invsxx3m;
		  (*xx3)[3] *= invsxx3m;
                }
                mother->scale[xs][i] = LOG (sxx3m)  + sxx1[i] + sxx2[i];
	      }
	    else
	      {
		mother->scale[xs][i] = sxx1[i] + sxx2[i];
	      }
	  }
    }
}    /* nuview */
#endif /*GAP*/

///
/// conditional likleihood calcculator using a simplified scheme called ancestral
/// with two lineages the one with the higher value will win and is the ancestor
/// non-Altivec version
void
  nuview_ancestral (mutationmodel_fmt *s, long sublocus, long xs, node * mother, world_fmt * world, const long locus)
{
  (void) sublocus;
  (void) world;
  (void) locus;
  const long endsite     = s->numpatterns;
  long i;
  MYREAL xx1t0, xx1t1,xx1t2, xx1t3;
  MYREAL xx2t0, xx2t1, xx2t2, xx2t3;
  
  MYREAL lw1, lw2, ratio1, yy1, yy2, sum1, sum2;
  node *q, *r;
  sitelike *xx1, *xx2;
  register MYREAL freqa, freqc, freqg, freqt;
  //register MYREAL freqar, freqcy, freqgr, freqty;
 
  freqa = s->basefreqs[NUC_A];
  freqc = s->basefreqs[NUC_C];
  freqg = s->basefreqs[NUC_G];
  freqt = s->basefreqs[NUC_T];
  //    freqr = s->basefreqs[NUC_R];
  //freqy = s->basefreqs[NUC_Y];
  //freqar = s->basefreqs[NUC_AR];
  //freqcy = s->basefreqs[NUC_CY];
  //freqgr = s->basefreqs[NUC_GR];
  //freqty = s->basefreqs[NUC_TY];
  
  q = crawlback (mother->next);
  r = crawlback (mother->next->next);
  lw1 = -q->v * s->fracchange;
  lw2 = -r->v * s->fracchange;
  ratio1 = lw1 / (lw1 + lw2);
  yy1 = (1. - ratio1);
  yy2 = ratio1;
  //printf("%f ", q->tyme);
  //    for (i = 0; i < endsite; i++)
  //printf("(%f %f %f %f)", q->x.s[i][0][0], q->x.s[i][0][1], q->x.s[i][0][2], q->x.s[i][0][3]);
  //printf("\n%f ", r->tyme);
  //for (i = 0; i < endsite; i++)
  //printf("(%f %f %f %f)", r->x.s[i][0][0], r->x.s[i][0][1], r->x.s[i][0][2], r->x.s[i][0][3]);
  //printf("\n");
  for (i = 0; i < endsite; i++)
    {
      xx1 = &(q->x[xs].s[i][0]);
      xx2 = &(r->x[xs].s[i][0]);
      
      
      xx1t0 = (*xx1)[0];
      xx1t1 = (*xx1)[1];
      xx1t2 = (*xx1)[2];
      xx1t3 = (*xx1)[3];
      
      xx2t0 = (*xx2)[0];
      xx2t1 = (*xx2)[1];
      xx2t2 = (*xx2)[2];
      xx2t3 = (*xx2)[3];
      
      
      
      sum1 =
	yy1 * (freqa * xx1t0 + freqc * xx1t1 +
	       freqg * xx1t2 + freqt * xx1t3);
      sum2 =
	yy2 * (freqa * xx2t0 + freqc * xx2t1 +
	       freqg * xx2t2 + freqt * xx2t3);
      if (fabs(sum1-sum2) <= DBL_EPSILON)
	sum1 += RANDUM () > 0.5 ? -1. : 1.;
      if (sum1 > sum2)
	memcpy (mother->x[xs].s[i][0], xx1, sizeof (sitelike));
      else
	memcpy (mother->x[xs].s[i][0], xx2, sizeof (sitelike));
    }
}    /* nuview_ancestral */

	
void
adjustroot (node * r)
{
  r->next->tyme = r->tyme;
  r->next->length = r->length;
  r->next->v = r->v;
  r->next->next->tyme = r->tyme;
  r->next->next->length = r->length;
  r->next->next->v = r->v;
}


/// \brief Calculation of conditional likelihood for new genealogy
///
/// Calculation of conditional likelihood for new genealogy
/// most of the code is unrolled for efficiency and may local variables are
/// declared also for speed
/// non-Altivec version 
void
  pseudonu_seq (mutationmodel_fmt *s, proposal_fmt * proposal, xarray_fmt * xxx1,  MYREAL *lx1, MYREAL v1, xarray_fmt * xxx2, MYREAL * lx2, MYREAL v2, long xs)
{
  (void) proposal;
  (void) lx1;
  (void) lx2;
  const long endsite     = s->numpatterns + s->addon;
  const long rcategs = s->numsiterates;
  const long categs = s->numcategs;
  const tbl_fmt tbl = s->tbl;
  phenotype xxxx1 = xxx1[xs].s;
  phenotype xxxx2 = xxx2[xs].s;
  long i, j;
#ifndef GRANDCENTRAL_notworking /*replace SNOWLEOPARD*/
  long k;
  register sitelike *xx1, *xx2;//, *xx3;
  register MYREAL xx1t0, xx1t1,xx1t2, xx1t3;
  register MYREAL xx2t0, xx2t1, xx2t2, xx2t3;
  register MYREAL vzsumr1, vzsumr2, vzsumy1, vzsumy2, sum1, sum2, sumr1, sumr2,
    sumy1, sumy2; 
  register sitelike *xx1copy, *xx2copy;
  
  register MYREAL ta, tc,  tg, tt;
  register MYREAL tta, ttc,  ttg, ttt;
  register valrec *tbljk;
  register MYREAL xxsumr1, xxsumr2, xxsumy1, xxsumy2;
#endif
  register MYREAL ww1, ww2, zz1, zz2, lw1, lw2, yy1, yy2, ww1zz1, vv1zz1, ww2zz2, vv2zz2;
  register MYREAL freqa, freqc, freqg, freqt, freqar, freqgr, freqcy, freqty;
  
  register valrec *tbl00;
  register valrec *tblij;
  
  freqa = s->basefreqs[NUC_A];
  freqc = s->basefreqs[NUC_C];
  freqg = s->basefreqs[NUC_G];
  freqt = s->basefreqs[NUC_T];
  freqar = s->basefreqs[NUC_AR];
  freqcy = s->basefreqs[NUC_CY];
  freqgr = s->basefreqs[NUC_GR];
  freqty = s->basefreqs[NUC_TY];
  
  lw1 = -v1 * s->fracchange;
  // use shortcut for dataset that do not use mutiple categories
  if ((rcategs | categs ) == 1)
    {
      tbl00 = tbl[0][0];
      lw2 = -v2 * s->fracchange;
      
      ww1 = EXP (tbl00->ratxi * lw1);
      zz1 = EXP (tbl00->ratxv * lw1);
      ww2 = EXP (tbl00->ratxi * lw2);
      zz2 = EXP (tbl00->ratxv * lw2);
      
      ww1zz1 = ww1 * zz1;
      vv1zz1 = (1.0 - ww1) * zz1;
      
      ww2zz2 = ww2 * zz2;
      vv2zz2 = (1.0 - ww2) * zz2;
      
      yy1 = 1.0 - zz1;
      yy2 = 1.0 - zz2;
#ifdef GRANDCENTRAL_notworking /*replace SNOWLEOPARD*/
      dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
      dispatch_apply(endsite, queue,
		     ^(unsigned long i) {
		       register sitelike *xx1, *xx2;
		       register MYREAL xx1t0, xx1t1,xx1t2, xx1t3;
		       register MYREAL xx2t0, xx2t1, xx2t2, xx2t3,  ta, tc,tg, tt, tta, ttc, ttg, ttt,xxsumr1, xxsumr2,xxsumy1,xxsumy2;
		       register MYREAL vzsumr1, vzsumr2, vzsumy1, vzsumy2, sum1, sum2, sumr1, sumr2,
			 sumy1, sumy2;
		       sitelike *xx1copy, *xx2copy;
#else
		       for (i = 0; i < endsite; i++)
			 {
#endif
			   xx1 = xx1copy = &(xxxx1[i][0]);
			   xx2 = xx2copy = &(xxxx2[i][0]);
			   
			   xx1t0 = (*xx1)[0];
			   ta = freqa * xx1t0; 
			   xx1t1 = (*xx1copy)[1];
			   xx2t1 = (*xx2)[1];
			   tc = freqc * xx1t1; 
			   xx1t2 = (*xx1)[2];
			   xx2t2 = (*xx2copy)[2];
			   tg = freqg * xx1t2; 
			   xx1t3 = (*xx1copy)[3];
			   xx2t3 = (*xx2)[3];
			   xx2t0 = (*xx2copy)[0];
			   tt = freqt * xx1t3;
			   
			   tta = freqa * xx2t0;
			   ttc = freqc * xx2t1;
			   ttg = freqg * xx2t2;
			   ttt = freqt * xx2t3;
			   
			   //sum1 = freqa * xx1t0 + freqc * xx1t1 + freqg * xx1t2 + freqt * xx1t3;
			   sum1 = ta + tc + tg + tt;
			   sum2 = tta + ttc + ttg + ttt;
			   sum1 *= yy1;
			   sum2 *= yy2;
			   
			   sumr1 = freqar * xx1t0 + freqgr * xx1t2;
			   sumr2 = freqar * xx2t0 + freqgr * xx2t2;
			   sumy1 = freqcy * xx1t1 + freqty * xx1t3;
			   sumy2 = freqcy * xx2t1 + freqty * xx2t3;
			   
			   vzsumr1 = vv1zz1 * sumr1;
			   vzsumr2 = vv2zz2 * sumr2;
			   vzsumy1 = vv1zz1 * sumy1;
			   vzsumy2 = vv2zz2 * sumy2;
			   
			   xxsumr1 = sum1 + vzsumr1;
			   xxsumr2 = sum2 + vzsumr2;
			   
			   (*xx1)[0] = (xxsumr1 + ww1zz1 * xx1t0) * (xxsumr2 + ww2zz2 * xx2t0);
			   (*xx1)[2] = (xxsumr1 + ww1zz1 * xx1t2) * (xxsumr2 + ww2zz2 * xx2t2);
			   xxsumy1 = sum1 + vzsumy1;
			   xxsumy2 = sum2 + vzsumy2;
			   (*xx1)[1] = (xxsumy1 + ww1zz1 * xx1t1) * (xxsumy2 + ww2zz2 * xx2t1);
			   (*xx1)[3] = (xxsumy1 + ww1zz1 * xx1t3) * (xxsumy2 + ww2zz2 * xx2t3);
			 }
#ifdef GRANDCENTRAL_notworking /*replace SNOWLEOPARD*/
		       );
#endif
		     }
      else
	{
	  for (i = 0; i < rcategs; i++)
	    for (j = 0; j < categs; j++)
	      {
		tblij = tbl[i][j];
		tblij->ww1 = EXP (tblij->ratxi * lw1);
		tblij->zz1 = EXP (tblij->ratxv * lw1);
		tblij->ww1zz1 = tblij->ww1 * tblij->zz1;
		tblij->vv1zz1 = (1.0 - tblij->ww1) * tblij->zz1;
	      }
	  lw2 = -v2 * s->fracchange;
	  for (i = 0; i < rcategs; i++)
	    for (j = 0; j < categs; j++)
	      {
		tblij = tbl[i][j];
		tblij->ww2 = EXP (tblij->ratxi * lw2);
		tblij->zz2 = EXP (tblij->ratxv * lw2);
		tblij->ww2zz2 = tblij->ww2 * tblij->zz2;
		tblij->vv2zz2 = (1.0 - tblij->ww2) * tblij->zz2;
	      }
#ifdef GRANDCENTRAL_notworking /*replace SNOWLEOPARD*/
	  dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
	  dispatch_apply(endsite, queue,
			 ^(unsigned long i) {
			   long j, k;
			   register sitelike *xx1, *xx2;
			   register MYREAL xx1t0, xx1t1,xx1t2, xx1t3;
			   register MYREAL xx2t0, xx2t1, xx2t2, xx2t3;
			   register MYREAL vzsumr1, vzsumr2, vzsumy1, vzsumy2, sum1, sum2, sumr1, sumr2,
			     sumy1, sumy2;
			   valrec * tbljk;
			   MYREAL ww1zz1, vv1zz1,yy1,ww2zz2,vv2zz2,yy2;
#else
			   for (i = 0; i < endsite; i++)
			     {
#endif
			       k = s->category[s->alias[i] - 1] - 1;
			       for (j = 0; j < rcategs; j++)
				 {
				   tbljk = tbl[j][k];
				   ww1zz1 = tbljk->ww1zz1;
				   vv1zz1 = tbljk->vv1zz1;
				   yy1 = 1.0 - tbljk->zz1;
				   ww2zz2 = tbljk->ww2zz2;
				   vv2zz2 = tbljk->vv2zz2;
				   yy2 = 1.0 - tbljk->zz2;
				   xx1 = &(xxxx1[i][j]);
				   xx2 = &(xxxx2[i][j]);
				   
				   xx1t0 = (*xx1)[0];
				   xx1t1 = (*xx1)[1];
				   xx1t2 = (*xx1)[2];
				   xx1t3 = (*xx1)[3];
				   
				   xx2t0 = (*xx2)[0];
				   xx2t1 = (*xx2)[1];
				   xx2t2 = (*xx2)[2];
				   xx2t3 = (*xx2)[3];
				   
				   sum1 = yy1 * (freqa * xx1t0 + freqc * xx1t1 +
						 freqg * xx1t2 + freqt * xx1t3);
				   sum2 = yy2 * (freqa * xx2t0 + freqc * xx2t1 +
						 freqg * xx2t2 + freqt * xx2t3);
				   sumr1 = freqar * xx1t0 + freqgr * xx1t2;
				   sumr2 = freqar * xx2t0 + freqgr * xx2t2;
				   sumy1 = freqcy * xx1t1 + freqty * xx1t3;
				   sumy2 = freqcy * xx2t1 + freqty * xx2t3;
				   vzsumr1 = vv1zz1 * sumr1;
				   vzsumr2 = vv2zz2 * sumr2;
				   /* xx3[j][0] */
				   (*xx1)[0] =
				     (sum1 + ww1zz1 * xx1t0 + vzsumr1) * (sum2 +
									  ww2zz2 * xx2t0 +
									  vzsumr2);
				   /* xx3[j][2] */
				   (*xx1)[2] =
				     (sum1 + ww1zz1 * xx1t2 + vzsumr1) * (sum2 +
									  ww2zz2 * xx2t2 +
									  vzsumr2);
				   vzsumy1 = vv1zz1 * sumy1;
				   vzsumy2 = vv2zz2 * sumy2;
				   /* xx3[j][1] */
				   (*xx1)[1] =
				     (sum1 + ww1zz1 * xx1t1 + vzsumy1) * (sum2 +
									  ww2zz2 * xx2t1 +
									  vzsumy2);
				   /* xx3[j][3] */
				   (*xx1)[3] =
				     (sum1 + ww1zz1 * xx1t3 + vzsumy1) * (sum2 +
									  ww2zz2 * xx2t3 +
									  vzsumy2);
				 }
			     }
#ifdef GRANDCENTRAL_notworking /*replace SNOWLEOPARD*/
			   );
#endif
			 }    /* pseudonu_seq */
	
}
      /// formatting problems?
      //}}
///
/// calculates conditiional likleihood on a fake tree before acceptance/rejection scheme
/// non-altivec version
void  pseudonu_seq_slow (mutationmodel_fmt *s, proposal_fmt * proposal, xarray_fmt * xxx1,  MYREAL *lx1, MYREAL v1, xarray_fmt * xxx2, MYREAL * lx2, MYREAL v2, long xs) 
{
  (void) proposal;
  const long endsite = s->numpatterns + s->addon;
  boolean newscale=FALSE;
  long i, j, k;
  //static long count = 0;
  phenotype xxxx1 = xxx1[xs].s;
  phenotype xxxx2 = xxx2[xs].s;
  MYREAL *sxx1 = lx1;
  MYREAL *sxx2 = lx2;
  register MYREAL xx1t0, xx1t1,xx1t2, xx1t3;
  register MYREAL xx2t0, xx2t1, xx2t2, xx2t3;
  register MYREAL ta, tc,  tg, tt;
  register MYREAL tta, ttc,  ttg, ttt;

  register MYREAL lw1, lw2, yy1, yy2, ww1zz1, vv1zz1, ww2zz2, vv2zz2, vzsumr1,
    vzsumr2, vzsumy1, vzsumy2, sum1, sum2, sumr1, sumr2,
    sumy1, sumy2, ww1, ww2, zz1, zz2;
  register MYREAL freqa, freqc, freqg, freqt, freqar, freqgr, freqcy, freqty;
  
  register MYREAL xxsumr1, xxsumr2, xxsumy1, xxsumy2;
  register  MYREAL sxx3m, tempsxx3m, invsxx3m;
  register  sitelike *xx1, *xx2, *xx1copy, *xx2copy;
  register  valrec *tblij;
  register  valrec *tbljk;
  register  valrec *tbl00;
  const  long rcategs = s->numsiterates;
  const  long categs = s->numcategs;
  register  tbl_fmt tbl = s->tbl;

  freqa = s->basefreqs[NUC_A];
  freqc = s->basefreqs[NUC_C];
  freqg = s->basefreqs[NUC_G];
  freqt = s->basefreqs[NUC_T];
  //freqr = s->basefreqs[NUC_R];
  //freqy = s->basefreqs[NUC_Y];
  freqar = s->basefreqs[NUC_AR];
  freqcy = s->basefreqs[NUC_CY];
  freqgr = s->basefreqs[NUC_GR];
  freqty = s->basefreqs[NUC_TY];
  lw1 = -v1 * s->fracchange;
  if ((rcategs | categs) == 1)
    {
        tbl00 = tbl[0][0];
        lw2 = -v2 * s->fracchange;

        ww1 = EXP (tbl00->ratxi * lw1);
        zz1 = EXP (tbl00->ratxv * lw1);
        ww2 = EXP (tbl00->ratxi * lw2);
        zz2 = EXP (tbl00->ratxv * lw2);

        ww1zz1 = ww1 * zz1;
        vv1zz1 = (1.0 - ww1) * zz1;

        ww2zz2 = ww2 * zz2;
        vv2zz2 = (1.0 - ww2) * zz2;

        yy1 = 1.0 - zz1;
        yy2 = 1.0 - zz2;

        for (i = 0; i < endsite; i++)
        {
            xx1 = xx1copy = &(xxxx1[i][0]);
            xx2 = xx2copy = &(xxxx2[i][0]);
 
	    xx1t0 = (*xx1)[0];
	    ta = freqa * xx1t0; 
	    xx1t1 = (*xx1copy)[1];
	    xx2t1 = (*xx2)[1];
	    tc = freqc * xx1t1; 
	    xx1t2 = (*xx1)[2];
	    xx2t2 = (*xx2copy)[2];
	    tg = freqg * xx1t2; 
	    xx1t3 = (*xx1copy)[3];
	    xx2t3 = (*xx2)[3];
	    xx2t0 = (*xx2copy)[0];
	    tt = freqt * xx1t3;
	    
	    tta = freqa * xx2t0;
	    ttc = freqc * xx2t1;
	    ttg = freqg * xx2t2;
	    ttt = freqt * xx2t3;
	  
	    sum1 = ta + tc + tg + tt;
	    sum2 = tta + ttc + ttg + ttt;
	    sum1 *= yy1;
	    sum2 *= yy2;
	    
            sumr1 = freqar * xx1t0 + freqgr * xx1t2;
            sumr2 = freqar * xx2t0 + freqgr * xx2t2;
            sumy1 = freqcy * xx1t1 + freqty * xx1t3;
            sumy2 = freqcy * xx2t1 + freqty * xx2t3;

            vzsumr1 = vv1zz1 * sumr1;
            vzsumr2 = vv2zz2 * sumr2;
            vzsumy1 = vv1zz1 * sumy1;
            vzsumy2 = vv2zz2 * sumy2;

	    xxsumr1 = sum1 + vzsumr1;
	    xxsumr2 = sum2 + vzsumr2;
	    
	    (*xx1)[0] = (xxsumr1 + ww1zz1 * xx1t0) * (xxsumr2 + ww2zz2 * xx2t0);
	    (*xx1)[2] = (xxsumr1 + ww1zz1 * xx1t2) * (xxsumr2 + ww2zz2 * xx2t2);
	    xxsumy1 = sum1 + vzsumy1;
	    xxsumy2 = sum2 + vzsumy2;
	    (*xx1)[1] = (xxsumy1 + ww1zz1 * xx1t1) * (xxsumy2 + ww2zz2 * xx2t1);
	    (*xx1)[3] = (xxsumy1 + ww1zz1 * xx1t3) * (xxsumy2 + ww2zz2 * xx2t3);
            sxx1[i] += sxx2[i];
	    if(!((*xx1)[0]+ (*xx1)[1]+ (*xx1)[2]+ (*xx1)[3]>EPSILON))
	      {
		//		if(!((*xx1)[0]+ (*xx1)[1]+ (*xx1)[2]+ (*xx1)[3]>0.0))
		//  {
		//    warning("pseudolike problem\n");
		//  }
		sxx3m = -MYREAL_MAX;
		tempsxx3m = MAX ((*xx1)[0], (*xx1)[1]);
		tempsxx3m = MAX (tempsxx3m, (*xx1)[2]);
		tempsxx3m = MAX (tempsxx3m, (*xx1)[3]);
		if (tempsxx3m > sxx3m)
		  sxx3m = tempsxx3m;
		invsxx3m = 1. / sxx3m;
		(*xx1)[0] *= invsxx3m;
		(*xx1)[1] *= invsxx3m;
		(*xx1)[2] *= invsxx3m;
		(*xx1)[3] *= invsxx3m;
		sxx1[i] += LOG (sxx3m);// + sxx2[i];
	      }
	}
    }
  else
    {
        lw2 = -v2 * s->fracchange;
        for (i = 0; i < rcategs; i++)
	  {
            for (j = 0; j < categs; j++)
	      {
                tblij = tbl[i][j];
                tblij->ww1 = EXP (tblij->ratxi * lw1);
                tblij->zz1 = EXP (tblij->ratxv * lw1);
                tblij->ww1zz1 = tblij->ww1 * tblij->zz1;
                tblij->vv1zz1 = (1.0 - tblij->ww1) * tblij->zz1;
                tblij->ww2 = EXP (tblij->ratxi * lw2);
                tblij->zz2 = EXP (tblij->ratxv * lw2);
                tblij->ww2zz2 = tblij->ww2 * tblij->zz2;
                tblij->vv2zz2 = (1.0 - tblij->ww2) * tblij->zz2;
	      }
	  }
	for (i = 0; i < endsite; i++)
	  {
	    k = s->category[s->alias[i] - 1] - 1;
	    for (j = 0; j < rcategs; j++)
	      {
		tbljk = tbl[j][k];
		ww1zz1 = tbljk->ww1zz1;
		vv1zz1 = tbljk->vv1zz1;
		yy1 = 1.0 - tbljk->zz1;
		ww2zz2 = tbljk->ww2zz2;
		vv2zz2 = tbljk->vv2zz2;
		yy2 = 1.0 - tbljk->zz2;
		xx1 = xx1copy = &(xxxx1[i][j]);
		xx2 = xx2copy = &(xxxx2[i][j]);                
                
		xx1t0 = (*xx1)[0];
		ta = freqa * xx1t0; 
		xx1t1 = (*xx1copy)[1];
		xx2t1 = (*xx2)[1];
		tc = freqc * xx1t1; 
		xx1t2 = (*xx1)[2];
		xx2t2 = (*xx2copy)[2];
		tg = freqg * xx1t2; 
		xx1t3 = (*xx1copy)[3];
		xx2t3 = (*xx2)[3];
		xx2t0 = (*xx2copy)[0];
		tt = freqt * xx1t3;
		
		tta = freqa * xx2t0;
		ttc = freqc * xx2t1;
		ttg = freqg * xx2t2;
		ttt = freqt * xx2t3;
		
		sum1 = ta + tc + tg + tt;
		sum2 = tta + ttc + ttg + ttt;
		sum1 *= yy1;
		sum2 *= yy2;

		sumr1 = freqar * xx1t0 + freqgr * xx1t2;
		sumr2 = freqar * xx2t0 + freqgr * xx2t2;
		sumy1 = freqcy * xx1t1 + freqty * xx1t3;
		sumy2 = freqcy * xx2t1 + freqty * xx2t3;
                
		vzsumr1 = vv1zz1 * sumr1;
		vzsumr2 = vv2zz2 * sumr2;
		vzsumy1 = vv1zz1 * sumy1;
		vzsumy2 = vv2zz2 * sumy2;
                
		xxsumr1 = sum1 + vzsumr1;
		xxsumr2 = sum2 + vzsumr2;
		xxsumy1 = sum1 + vzsumy1;
		xxsumy2 = sum2 + vzsumy2;
		
		(*xx1)[0] = (xxsumr1 + ww1zz1 * xx1t0) * (xxsumr2 + ww2zz2 * xx2t0);
		(*xx1)[1] = (xxsumy1 + ww1zz1 * xx1t1) * (xxsumy2 + ww2zz2 * xx2t1);
	        (*xx1)[2] = (xxsumr1 + ww1zz1 * xx1t2) * (xxsumr2 + ww2zz2 * xx2t2);
		(*xx1)[3] = (xxsumy1 + ww1zz1 * xx1t3) * (xxsumy2 + ww2zz2 * xx2t3);
		if(!((*xx1)[0]+ (*xx1)[1]+ (*xx1)[2]+ (*xx1)[3]>EPSILON))
		  newscale=TRUE;
	      }
	    if(newscale)
	      {
		sxx3m = -MYREAL_MAX;
		for (j = 0; j < rcategs; j++)
		  {
		    xx1 = &(xxxx1[i][j]);
		    tempsxx3m = MAX ((*xx1)[0], (*xx1)[1]);
		    tempsxx3m = MAX (tempsxx3m, (*xx1)[2]);
		    tempsxx3m = MAX (tempsxx3m, (*xx1)[3]);
		    if (tempsxx3m > sxx3m)
		      sxx3m = tempsxx3m;
		  }
		invsxx3m = 1. / sxx3m;
		for (j = 0; j < rcategs; j++)
		  {
		    xx1 = &(xxxx1[i][j]);
		    (*xx1)[0] *= invsxx3m;
		    (*xx1)[1] *= invsxx3m;
		    (*xx1)[2] *= invsxx3m;
		    (*xx1)[3] *= invsxx3m;
		  }
		sxx1[i] += LOG (sxx3m) + sxx2[i];
	      }
	    else
	      {
		sxx1[i] += sxx2[i];
	      }
	  }
	/*count++;
        if (count == SCALEINTERVAL)
        {
            count = 0;
            for (i = 0; i < endsite; i++)
            {
	      sxx3m = -MYREAL_MAX;
	      for (j = 0; j < rcategs; j++)
                {
		  xx1 = &(xxxx1[i][j]);
		  tempsxx3m = MAX ((*xx1)[0], (*xx1)[1]);
		  tempsxx3m = MAX (tempsxx3m, (*xx1)[2]);
		  tempsxx3m = MAX (tempsxx3m, (*xx1)[3]);
		  if (tempsxx3m > sxx3m)
		    sxx3m = tempsxx3m;
                }
	      invsxx3m = 1. / sxx3m;
	      for (j = 0; j < rcategs; j++)
                {
		  xx1 = &(xxxx1[i][j]);
		  (*xx1)[0] *= invsxx3m, (*xx1)[1] *= invsxx3m;
		  (*xx1)[2] *= invsxx3m, (*xx1)[3] *= invsxx3m;
                }
	      sxx1[i] += LOG (sxx3m);
            }
	    }*/
    }
}    /* pseudonu_seq */
 
///
/// calculates ancestral likelihood before acc/reject scheme
/// non-atlivec version
void  pseudonu_anc (mutationmodel_fmt *s, proposal_fmt * proposal, xarray_fmt * xxx1,  MYREAL *lx1, MYREAL v1, xarray_fmt * xxx2, MYREAL * lx2, MYREAL v2, long xs) 
{
  (void) proposal;
  (void) lx1;
  (void) lx2;
    long i;
    phenotype xxxx1 = xxx1[xs].s;
    phenotype xxxx2 = xxx2[xs].s;
    const long endsite = s->numpatterns;
    MYREAL xx1t0, xx1t1,xx1t2, xx1t3;
    MYREAL xx2t0, xx2t1, xx2t2, xx2t3;
	
    MYREAL lw1, lw2, ratio1, yy1, yy2, sum1, sum2;
    MYREAL freqa, freqc, freqg, freqt;
    //seqmodel_fmt *seq;
    sitelike *xx1, *xx2;
    //seq = proposal->world->data->seq[0];
    freqa = s->basefreqs[NUC_A];
    freqc = s->basefreqs[NUC_C];
    freqg = s->basefreqs[NUC_G];
    freqt = s->basefreqs[NUC_T];
    //    freqr = s->basefreqs[NUC_R];
    //freqy = s->basefreqs[NUC_Y];
    //freqar = s->basefreqs[NUC_AR];
    //freqcy = s->basefreqs[NUC_CY];
    //freqgr = s->basefreqs[NUC_GR];
    //freqty = s->basefreqs[NUC_TY];

    lw1 = -v1 * s->fracchange;
    lw2 = -v2 * s->fracchange;
    ratio1 = lw1 / (lw1 + lw2);
    yy1 = (1. - ratio1);
    yy2 = ratio1;
    //for (i = 0; i < endsite; i++)
    //printf("(%f %f %f %f)", xxx1[i][0][0], xxx1[i][0][1], xxx1[i][0][2], xxx1[i][0][3]);
    //printf("\n");
    //for (i = 0; i < endsite; i++)
    //printf("(%f %f %f %f)", xxx2[i][0][0], xxx2[i][0][1], xxx2[i][0][2], xxx2[i][0][3]);
    //printf("\n");
    for (i = 0; i < endsite; i++)
    {
        xx1 = &(xxxx1[i][0]);
        xx2 = &(xxxx2[i][0]);
		
        
		xx1t0 = (*xx1)[0];
		xx1t1 = (*xx1)[1];
		xx1t2 = (*xx1)[2];
		xx1t3 = (*xx1)[3];
		
		xx2t0 = (*xx2)[0];
		xx2t1 = (*xx2)[1];
		xx2t2 = (*xx2)[2];
		xx2t3 = (*xx2)[3];
		
		
        sum1 =
            yy1 * (freqa * xx1t0 + freqc * xx1t1 + freqg * xx1t2 +
                   freqt * xx1t3);
        sum2 =
            yy2 * (freqa * xx2t0 + freqc * xx2t1 + freqg * xx2t2 +
                   freqt * xx2t3);
        if (fabs(sum1-sum2) <= DBL_EPSILON)
            sum1 += RANDUM () > 0.5 ? -1. : 1.;
        if (sum1 > sum2)
            memcpy (xxxx1[i][0], *xx1, sizeof (sitelike));
        else
            memcpy (xxxx1[i][0], *xx2, sizeof (sitelike));
    }
}    /* pseudo_nu_anc */

///
/// ancestral conditional method, calculates ancestral tree likelihood
MYINLINE MYREAL
  pseudo_tl_anc (mutationmodel_fmt *s, phenotype xx1, phenotype xx2, MYREAL v1, MYREAL v2,
               proposal_fmt * proposal, world_fmt * world)
{
  (void) xx2;
  (void) v1;
  (void) v2;
  (void) proposal;
  (void) world;
  const long endsite = s->numpatterns;
    contribarr tterm;
    MYREAL summ;
    long i;
    sitelike *x1;
	
    register MYREAL freqa, freqc, freqg, freqt;
    //seqmodel_fmt *seq;
    //seq = proposal->world->data->seq[0];

    freqa = s->basefreqs[NUC_A];
    freqc = s->basefreqs[NUC_C];
    freqg = s->basefreqs[NUC_G];
    freqt = s->basefreqs[NUC_T];

    summ = 0.0;
    for (i = 0; i < endsite; i++)
    {
        x1 = &(xx1[i][0]);
        tterm[0] =
            freqa * (*x1)[0] + freqc * (*x1)[1] + freqg * (*x1)[2] +
            freqt * (*x1)[3];
        summ += s->aliasweight[i] * LOG (tterm[0]);
        //printf("pseudo %3li> %f %f \n", i, tterm[0], summ);
    }
    return summ;
} /*anc*/


///
/// calculates the conditional likelihood on a tree
MYINLINE MYREAL
pseudo_tl_seq (mutationmodel_fmt *s, long xs, phenotype xx1, phenotype xx2, MYREAL v1, MYREAL v2,
               proposal_fmt * proposal, world_fmt * world)
{
  (void) xx2;
  (void) v1;
  (void) v2;
  (void) world;
  const MYREAL freqa = s->basefreqs[NUC_A];
  const MYREAL freqc = s->basefreqs[NUC_C];
  const MYREAL freqg = s->basefreqs[NUC_G];
  const MYREAL freqt = s->basefreqs[NUC_T];
  const long numpatterns= s->numpatterns;
  const long rcategs = s->numsiterates;
  const long categs  = s->numcategs;
  //const long numsites = s->numsites;
  MYREAL *probcat    = s->siteprobs;
  contribarr tterm;
  contribarr clai;
  contribarr like;
  contribarr nulike;
  //long          size = sizeof(MYREAL) * world->options->rcategs;
  MYREAL summ    = 0.0;
  MYREAL sum2    = 0.0;
  MYREAL sumc    = 0.0;
  MYREAL sumterm = 0.0;
  MYREAL lterm;
  MYREAL scale;
  long i; 
  long j; 
  long k;
  long lai;
  //worldoption_fmt *opt;
  sitelike *x1;
  
  //opt = world->options;
  summ = 0.0;
  if (rcategs == 1 && categs == 1)
    {
        for (i = 0; i < numpatterns; i++)
        {
            x1 = &(xx1[i][0]);
            tterm[0] = 
                freqa * (*x1)[0] + freqc * (*x1)[1] + freqg * (*x1)[2] +
                freqt * (*x1)[3];
            summ += s->aliasweight[i] * (LOG (tterm[0]) + proposal->mf[xs][i]);
        }
    }
    else
    {
      for (i = 0; i < numpatterns; i++)
	{
	  sumterm = 0.0;
	  scale = proposal->mf[xs][i];
	  for (j = 0; j < rcategs; j++)
            {
                x1 = &(xx1[i][j]);
                tterm[j] = freqa * (*x1)[0] + freqc * (*x1)[1] + freqg * (*x1)[2] + freqt * (*x1)[3];
                sumterm += probcat[j] * tterm[j];
	    }
	  lterm = LOG (sumterm) + scale;
	  for (j = 0; j < rcategs; j++)
	    clai[j] = tterm[j] / sumterm;
	  swap (&clai, &s->contribution[i]);
            //memcpy(s->contribution[i], clai, size);
            summ += s->aliasweight[i] * lterm;
        }
        for (j = 0; j < rcategs; j++)
            like[j] = 1.0;
        for (i = 0; i < s->numsites; i++)
        {
            sumc = 0.0;
            for (k = 0; k < rcategs; k++)
                sumc += probcat[k] * like[k];
            sumc *= s->lambda;
            if ((s->ally[i] > 0) && (s->location[s->ally[i] - 1] > 0))
            {
                lai = s->location[s->ally[i] - 1];
                swap (&s->contribution[lai - 1], &clai);
                //memcpy(clai, s->contribution[lai - 1], size);
                for (j = 0; j < rcategs; j++)
                    nulike[j] = ((1.0 - s->lambda) * like[j] + sumc) * clai[j];
            }
            else
            {
                for (j = 0; j < rcategs; j++)
                    nulike[j] = ((1.0 - s->lambda) * like[j] + sumc);
            }
            swap (&nulike, &like);
            //memcpy(like, nulike, size);
        }
        sum2 = 0.0;
        for (i = 0; i < rcategs; i++)
            sum2 += probcat[i] * like[i];
        summ += LOG (sum2);
    }
    //    printf("summ=%f\n",summ);
    return summ;
}

MYINLINE 
MYREAL
pseudo_tl_snp (mutationmodel_fmt *s, long xs, phenotype xx1, phenotype xx2, MYREAL v1, MYREAL v2,
               proposal_fmt * proposal, world_fmt * world)
{
  (void) xx2;
  (void) v1;
  (void) v2;
  (void) world;

  const long endsite = s->numpatterns;
  const long rcategs = s->numsiterates;
  const long categs  = s->numcategs;
  contribarr tterm, invariants;
  contribarr like, nulike, clai;
  //long          size = sizeof(MYREAL) * world->options->rcategs;
    MYREAL summ, sum2, sumc, sumterm, lterm;
    long i, j, k, lai;
    //worldoption_fmt *opt = world->options;
    register MYREAL freqa, freqc, freqg, freqt;//, freqr, freqy;
 
    sitelike *x1;
    assert(rcategs>0);
    freqa = s->basefreqs[NUC_A];
    freqc = s->basefreqs[NUC_C];
    freqg = s->basefreqs[NUC_G];
    freqt = s->basefreqs[NUC_T];
//    freqr = s->basefreqs[NUC_R];
//    freqy = s->basefreqs[NUC_Y];

    summ = 0.0;
    memset (invariants, 0, sizeof (contribarr));
    snp_invariants (s, xs, invariants, world, world->locus, xx1, proposal->mf[xs]);
    if ((rcategs | categs) == 1)
    {
        for (i = 0; i < endsite; i++)
        {
            x1 = &(xx1[i][0]);
            tterm[0] =
	      ((freqa * (*x1)[0] + freqc * (*x1)[1] +
		freqg * (*x1)[2] + freqt * (*x1)[3]));
	    //
	    //
	    
            if (tterm[0] == 0.0)
	      {
#ifdef TESTDEBUG
		fprintf(stderr,"Freq/condL=(%f,%f | %f,%f | %f,%f | %f,%f) Invars=%f\n", 
			freqa , (*x1)[0] , freqc , (*x1)[1] ,
			freqg , (*x1)[2] , freqt , (*x1)[3] , invariants[0]);
       
                warning ("Tree incompatible with data\n");
#endif
		return (double) -HUGE;
	      }
            lterm = LOG (tterm[0]) + proposal->mf[xs][i] - log(invariants[0]);
	    //
	    if (lterm>0.0)
	      {
		//fprintf(stderr,"L%5li: %f Freq/condL=(%f,%f | %f,%f | %f,%f | %f,%f) Invars=%f\n", 
		//	i, lterm, freqa , (*x1)[0] , freqc , (*x1)[1] ,
		//	freqg , (*x1)[2] , freqt , (*x1)[3] , invariants[0]);
		lterm = 0.0;
	      }
	    //
	    summ += s->aliasweight[i] * lterm;			
        }
        like[0] = 1.0;
        for (i = 0; i < endsite; i++)
        {
            sumc = s->lambda * like[0];
            nulike[0] = ((1.0 - s->lambda) * like[0] + sumc);
            //memcpy(like, nulike, size);
            swap (&nulike, &like);
        }
        summ += LOG (like[0]);
        return summ;
    }
    else
    {
        for (i = 0; i < endsite - 4; i++)
        {
            //check against JF code: k = s->category[s->alias[i] - 1] - 1;
            for (j = 0; j < rcategs; j++)
            {
                x1 = &(xx1[i][j]);
                tterm[j] =
                    (freqa * (*x1)[0] + freqc * (*x1)[1] +
                     freqg * (*x1)[2] +
                     freqt * (*x1)[3]) / invariants[j];
		//fixing positive logL
		if (tterm[j] > 1.0)
		  tterm[j] = 1.0;
            }
            sumterm = 0.0;
            for (j = 0; j < rcategs; j++)
                sumterm += s->siteprobs[j] * tterm[j];
            lterm = LOG (sumterm) + proposal->mf[xs][i];
            for (j = 0; j < rcategs; j++)
                clai[j] = tterm[j] / sumterm;
            swap (&clai, &s->contribution[i]);
            summ += s->aliasweight[i] * lterm;
        }
        for (j = 0; j < rcategs; j++)
            like[j] = 1.0;
        for (i = 0; i < endsite; i++)
        {
            sumc = 0.0;
            for (k = 0; k < rcategs; k++)
                sumc += s->siteprobs[k] * like[k];
            sumc *= s->lambda;
            if ((s->ally[i] > 0) && (s->location[s->ally[i] - 1] > 0))
            {
                lai = s->location[s->ally[i] - 1];
                swap (&s->contribution[lai - 1], &clai);
                for (j = 0; j < rcategs; j++)
                    nulike[j] = ((1.0 - s->lambda) * like[j] + sumc) * clai[j];
            }
            else
            {
                for (j = 0; j < rcategs; j++)
                    nulike[j] = ((1.0 - s->lambda) * like[j] + sumc);
            }
            swap (&nulike, &like);
            //memcpy(like, nulike, size);
        }
        sum2 = 0.0;
        for (i = 0; i < rcategs; i++)
            sum2 += s->siteprobs[i] * like[i];
        summ += LOG (sum2);
        return summ;
    }
}

MYINLINE
MYREAL
pseudo_tl_snp_unlinked (mutationmodel_fmt *s, long xs, phenotype xx1, phenotype xx2, MYREAL v1, MYREAL v2,
                        proposal_fmt * proposal, world_fmt * world)
{
  (void) xx2;
  (void) v1;
  (void) v2;
  (void) world;

  const long endsite = s->numpatterns;
  //const long rcategs = s->numsiterates;
  
    contribarr tterm, invariants;
    MYREAL summ, datasum = 0, lterm, result = 0;
    long i;
    //worldoption_fmt *opt;
    MYREAL freqa, freqc, freqg, freqt;//, freqr, freqy;

    sitelike *x1;
    freqa = s->basefreqs[NUC_A];
    freqc = s->basefreqs[NUC_C];
    freqg = s->basefreqs[NUC_G];
    freqt = s->basefreqs[NUC_T];
//    freqr = s->basefreqs[NUC_R];
//    freqy = s->basefreqs[NUC_Y];
	
    //opt = world->options;

    summ = 0.0;
    memset (invariants, 0, sizeof (contribarr));
    snp_invariants (s, xs, invariants, world, world->locus,  xx1, proposal->mf[xs]);
    for (i = 0; i < endsite; i++)
    {
        x1 = &(xx1[i][0]);
        tterm[0] =
            (freqa * (*x1)[0] + freqc * (*x1)[1] +
             freqg * (*x1)[2] + freqt * (*x1)[3]);
        if (tterm[0] == 0.0)
            error ("Tree incompatible with data\n");
		
        if (i % 5 == 0)
        {
            lterm = LOG (tterm[0]) + proposal->mf[xs][i];
            summ = 0;
            datasum = s->aliasweight[i / 5] * lterm;
        }
        else
            summ += pow (tterm[0], (MYREAL) s->aliasweight[i / 5]);
        if (((i + 1) % 5) == 0 && i != 0)
            result +=
                datasum + LOG ((1 - EXP (LOG (summ) - datasum)) / invariants[0]);
    }
    //EXP (sum) is the prob(xa | g)
    //              EXP (datasum) is prob(? a | g)
    //              panelsum = invariants is prob(x ? |g)
    // (datasum - sum) / invariants
    // ++some small number business
    return result;
}




MYREAL
treelike_anc (mutationmodel_fmt *s, long xs, world_fmt * world, long locus)
{
  (void) locus;
  const long endsite = s->numpatterns;
  contribarr tterm;
  
  MYREAL summ;
  long i;
  node *p;
  register MYREAL freqa, freqc,freqg,freqt;
  sitelike *x1;
  freqa = s->basefreqs[NUC_A];
  freqc = s->basefreqs[NUC_C];
  freqg = s->basefreqs[NUC_G];
  freqt = s->basefreqs[NUC_T];

    p = crawlback (world->root->next);
    summ = 0.0;
    for (i = 0; i < endsite; i++)
    {
        x1 = &(p->x[xs].s[i][0]);
        tterm[0] =
            freqa * (*x1)[0] + freqc * (*x1)[1] +
            freqg * (*x1)[2] + freqt * (*x1)[3];
        summ += s->aliasweight[i] * LOG (tterm[0]);
        //printf("real  %3li> %f %f \n", i, tterm[0], summ);
    }
    return summ;
}    /* treelike_anc */


///
/// write a tree to a diskfile, compile setting will allow to do this as a NEXUS file
void
treeout (FILE * file, node * joint, node * p, long s)
{
    /* write out file with representation of final tree */
    MYREAL x;
    char migstring[100];
    if (p->type == 't')
    {
#ifdef UEP	
        if(p->uep!=NULL)
            FPRINTF (file, "%s [%c]", p->nayme, p->uep[0]);
        else
            FPRINTF (file, "%s ", p->nayme);
#else
        FPRINTF (file, "%s ", p->nayme);
#endif
#ifdef TREECOMMENTS
	FPRINTF(file,"[& t %li:%.10f]",p->actualpop, p->tyme);
#endif
    }
    else
    {
        FPRINTF (file, "(");
        treeout (file, joint, crawlback (p->next), s);
        FPRINTF (file, ",");
        treeout (file, joint, crawlback (p->next->next), s);
        FPRINTF (file, ")");
#ifdef TREECOMMENTS
	FPRINTF(file,"[& c %li:%.10f]",p->actualpop, p->tyme);
#endif
    }
    if (p != joint)
    {
      x = showtop(crawlback (p))->tyme - p->tyme;
#ifdef UEP	
	if(p->uep!=NULL)
	  FPRINTF (file, ":%.10f [%c]", x, p->uep[0]);
	else
	  FPRINTF (file, ":%.10f ",  x);
#else
	FPRINTF (file, ":%.10f ", x);
#endif
        p = showtop (p->back);
        while (p->type == 'm' || p->type == 'd')
	  {
	    if(p->type == 'd')
	      sprintf (migstring, " [&D %li %li:%g]", p->pop, p->actualpop,
		       p->tyme - showtop (p->next->back)->tyme);
	    else
	      sprintf (migstring, " [&M %li %li:%g]", p->pop, p->actualpop,
		       p->tyme - showtop (p->next->back)->tyme);
#ifdef TREECOMMENTS
	    FPRINTF(file,"%s[& %c %li %li:%.10f]",migstring, p->type, p->pop, p->actualpop, p->tyme);
#else
            FPRINTF (file, "%s", migstring);
#endif
            p = showtop (p->back);
	  }
    }
    else
    {
#ifdef NEXUSTREE
	FPRINTF (file, ";\n");
#else
	FPRINTF (file, ":0;\n");
#endif
    }
}    /* treeout */

///
/// write a tree to a diskfile, compile setting will allow to do this as a NEXUS file
void
debugtreeout (FILE * file, node * joint, node * p, long s)
{
    /* write out file with representation of final tree */
  //  MYREAL x;
    char migstring[100];
    if (p->type == 't')
    {
      FPRINTF (file, "[%i]%li ", myID,p->id);
    }
    else
    {
        FPRINTF (file, "(");
        debugtreeout (file, joint, crawlback (p->next), s);
        FPRINTF (file, ",");
        debugtreeout (file, joint, crawlback (p->next->next), s);
        FPRINTF (file, ")<%li>",p->id);
    }
    if (p != joint)
    {
        //x = crawlback (p)->tyme - p->tyme;
        //	FPRINTF (file, ":%.10f ", x);
        p = showtop (p->back);
        while (p->type == 'm' || p->type == 'd')
	  {
            sprintf (migstring, " [%c %li-%li]", p->type, p->actualpop, p->pop);
            FPRINTF (file, "%s", migstring);
            p = showtop (p->back);
	  }
    }
    else
    {
	FPRINTF (file, ":0\n");

    }
}    /* treeout */

///
/// writes tree in newick format to a string
void treeout_string (char ** file, long *filesize, long *pos, node * joint, node * p, long s)
{
  char *tmp;
  MYREAL x;
  char migstring[100];
  tmp = (char *) mycalloc(1024,sizeof(char));
  /* write out file with representation of final tree */
  //    long w;
  if (p->type == 't')
    {
#ifdef UEP
        if(p->uep!=NULL)
	  print_to_buffer(file, filesize, tmp, pos, "%s [%c]", p->nayme, p->uep[0]);
        else
	  print_to_buffer(file, filesize, tmp, pos, "%s ", p->nayme);
#else
        print_to_buffer(file, filesize, tmp, pos, "%s ", p->nayme);
#endif		
#ifdef TREECOMMENTS
	print_to_buffer(file, filesize, tmp, pos,"[& t %li:%.10f]",p->actualpop,p->tyme);
#endif
    }
    else
    {
      print_to_buffer(file, filesize, tmp, pos, "(");
      treeout_string (file, filesize, pos, joint, crawlback (p->next), s);
      print_to_buffer(file, filesize, tmp, pos, ",");
      treeout_string (file, filesize, pos, joint, crawlback (p->next->next), s);
      print_to_buffer(file, filesize, tmp, pos, ")");
#ifdef TREECOMMENTS
      print_to_buffer(file, filesize, tmp, pos,"[& c %li:%.10f]",p->actualpop,p->tyme);
#endif

    }
    if (p == joint)
    {
        x = 0.0;
    }
    else
    {
        x = crawlback (p)->tyme - p->tyme;
    }
#ifdef UEP
    if(p->uep!=NULL)
      print_to_buffer(file, filesize, tmp, pos, ":%.10f [%c]", x, p->uep[0]);
    else
      print_to_buffer(file, filesize, tmp, pos, ":%.10f ", x);
#else
    print_to_buffer(file, filesize, tmp, pos, ":%.10f ", x);
#endif
    if (p != joint)
    {
        p = showtop (p->back);
        while (p->type == 'm'  || p->type == 'd')
        {
	  if(p->type == 'd')
	    sprintf (migstring, " [&D %li %li:%g]", p->pop, p->actualpop,
                     p->tyme - showtop (p->next->back)->tyme);
	  else
            sprintf (migstring, " [&M %li %li:%g]", p->pop, p->actualpop,
                     p->tyme - showtop (p->next->back)->tyme);
#ifdef TREECOMMENTS
	  sprintf(migstring,"%s[& %c %li %li:%.10f]",migstring, p->type, p->pop, p->actualpop, p->tyme);
#endif
	  print_to_buffer(file, filesize, tmp, pos, "%s", migstring);
	  p = showtop (p->back);
        }
    }
    else
    {
      print_to_buffer(file, filesize, tmp, pos, ";\n\0");
    }
    myfree(tmp);
}    /* treeout_string */


void
print_tree (world_fmt * world, long g, long *filepos)
{
#ifdef NEXUSTREE
  static long counter = 0;
#endif
  static long count   = 1;
  long pos            = 0;
  long allocval       = 0;
  char *tmp;
  tmp = (char *) calloc(1024,sizeof(char));

  switch (world->options->treeprint)
    {
    case BEST:
      if(world->options->treeinmemory)
	{
	  *filepos = 0;
	  pos = 0;
	  allocval = world->treespacealloc[world->locus];
	  //	  if (world->likelihood[g] >= world->besttreelike[world->locus])
	  if (world->param_like >= world->besttreelike[world->locus])
	    {
	      //printf("%i> print_tree(BEST) at %p  locus=%li bestlike=%f like=%f\n",myID, &world->treespace[world->locus], world->locus,world->besttreelike[world->locus], world->likelihood[g]);
	      world->besttreelike[world->locus] = world->param_like;
	      // world->besttreelike[world->locus] = world->likelihood[g];
	      world->treespace[world->locus] = (char *) myrealloc(world->treespace[world->locus],
								  sizeof(char) * LONGLINESIZE);
	      allocval = LONGLINESIZE;
	      //speed problem	      memset(world->treespace[world->locus],0, sizeof(char) * LONGLINESIZE);
	      world->treespace[world->locus][0]='\0';
	      pos = sprintf (world->treespace[world->locus], "\n[& Locus %li, best ln(L) = %f (c=coalescent node, t=tipnode) ]\n",
			     world->locus + 1, world->likelihood[g]);
#ifdef NEXUSTREE
	      pos = sprintf (world->treespace[world->locus], 
			     "\n[& Locus %li, best ln(L) = %f (c=coalescent node, t=tipnode) ]\ntree repl.%li = [&R] ",
			     world->locus + 1, world->likelihood[g],counter++);
#else
	      pos = sprintf (world->treespace[world->locus], "\n[& Locus %li, best ln(L) = %f (c=coalescent node, t=tipnode) ]\n",
			     world->locus + 1, world->likelihood[g]);
#endif
	      treeout_string (&(world->treespace[world->locus]), 
			      &allocval,&pos, 
			      crawlback (world->root->next),
			      crawlback (world->root->next), 0);
	      world->treespacealloc[world->locus] = allocval;
	      world->treespacenum[world->locus] = pos;
	    }
	  //	  else
	  //  {
	  //    printf("%i> *** print_tree(BEST) locus=%li bestlike=%f like=%f\n", myID, world->locus, world->besttreelike[world->locus], world->likelihood[g]);
	  //  }
	}
      else
	{
	  if (world->likelihood[g] > world->allikemax)
	    {
	      if (world->allikemax <= -MYREAL_MAX)
		{
		  *filepos = ftell (world->treefile);
		}
	      else
		{
		  fseek (world->treefile, *filepos, SEEK_SET);
		}
	      world->allikemax = world->likelihood[g];
	      FPRINTF (world->treefile,
			     "\n[& Locus %li, best ln(L) = %f (c=coalescent node, t=tipnode) ]\n",
		       world->locus + 1, world->likelihood[g]);
#ifdef NEXUSTREE
		  FPRINTF (world->treefile,"tree repl.%li = [&R] ",counter++);
#endif
	      treeout (world->treefile, crawlback (world->root->next),
		       crawlback (world->root->next), 0);
	    }
	}
      break;
    case ALL:
    case LASTCHAIN:
      if((count++ % world->options->treeinc) == 0)
	{
	  if (world->in_last_chain)
	    {
	      if(world->options->treeinmemory)
		{
		  *filepos = 0;
		  pos = world->treespacenum[world->locus];
		  allocval = world->treespacealloc[world->locus];
#ifdef NEXUSTREE
		  pos = print_to_buffer(&(world->treespace[world->locus]), &allocval, tmp, &pos,
					"\n[& Locus %li, best ln(L) = %f (c=coalescent node, t=tipnode) ]\ntree repl.%li = [&R] ",
					world->locus + 1, world->likelihood[g], counter++);
#else
		  pos = print_to_buffer(&(world->treespace[world->locus]), &allocval, tmp, &pos,
					"\n[& Locus %li, best ln(L) = %f (c=coalescent node, t=tipnode) ]\n",
					world->locus + 1, world->likelihood[g]);
#endif
		  treeout_string (&(world->treespace[world->locus]), 
				  &allocval,&pos, 
				  crawlback (world->root->next),
				  crawlback (world->root->next), 0);
		  world->treespacealloc[world->locus] = allocval;
		  world->treespacenum[world->locus] = pos;
		}
	      else
		{
		  // writing direct to file
		  FPRINTF (world->treefile,
			   "\n[& Locus %li, ln(L) = %f (c=coalescent node, t=tipnode) ]\n",
			   world->locus + 1, world->likelihood[g]);
#ifdef NEXUSTREE
		  FPRINTF (world->treefile,"tree repl.%li = [&R] ",counter++);
#endif
		  treeout (world->treefile, crawlback (world->root->next),
			   crawlback (world->root->next), 0);
		}
	    }
	}
      break;
    case myNONE:
      break;
    default:
      break;
    }
  fflush(world->treefile);
  myfree(tmp);
}


boolean
treereader (world_fmt * world,  option_fmt *options, data_fmt * data, long locus)
{
    /*
     * read a  tree with or without migration events from the usertree and set up nodes
     * and pointers
     */
    boolean has_migration=FALSE;
    node **nodelist;
    char *nayme;
    char *temp, *temp2;
    long pop, w, zz, z = 0, zzz = 0;
    world->nodep = (node **) mycalloc ((world->sumtips+(world->sumtips+1)),sizeof (node *));
    temp = (char *) mymalloc (LINESIZE * sizeof (char));
    temp2 = (char *) mymalloc (LINESIZE * sizeof (char));
    has_migration = treeread (data->utreefile, world, options, &(world->root), NULL);
    length_to_times (world->root->next->back);
    nodelist = (node **) mycalloc ((world->sumtips + 1), sizeof (node *));
    pop = -1;
    set_tree_pop (world->root, &pop);
    allocate_x (world->root, world, locus, WITHTIPS);
    find_tips (world->root, nodelist, &z);
    for (pop = 0; pop < data->numpop; pop++)
    {
        for (w = 0; w < data->numind[pop][world->locus]; w++)
        {
            strcpy (temp2, data->indnames[pop][w][world->locus]);
	    unpad(temp2," ");
	    //            temp2[strcspn (temp2, " ")] = '\0';
            sprintf (temp, "%li%s", pop, temp2);
            for (zz = 0; zz < z; zz++)
            {
                nayme = nodelist[zz]->nayme;
		unpad(nayme," _");
                if (!strcmp (temp, nayme) || !strcmp (temp2, nayme))
                {
                    world->nodep[zzz++] = nodelist[zz];
                    break;
                }
            }
        }
    }
    myfree(nodelist);
    myfree(temp);
    myfree(temp2);
    return has_migration;
}


char
processlength (FILE * file, node ** p)
{
    char ch;
    long digit, ordzero;
    MYREAL valyew, divisor;
    boolean pointread, minusread;
    assert(*p);
    ordzero = '0';
    pointread = FALSE;
    minusread = FALSE;
    valyew = 0.0;
    divisor = 1.0;
    ch = (char) getc (file);
    digit = ch - ordzero;
    while (((unsigned long) digit <= 9) | (ch == '.') || (ch == '-'))
    {
        if (ch == '.')
            pointread = TRUE;
        else if (ch == '-')
            minusread = TRUE;
        else
        {
            valyew = valyew * 10.0 + digit;
            if (pointread)
                divisor *= 10.0;
        }
        ch = (char) getc (file);
        digit = ch - ordzero;
    }
    if (!minusread)
        (*p)->length = valyew / divisor;
    else
        (*p)->length = 0.0;
    return ch;
}

boolean
treeread (FILE * file, world_fmt * world, option_fmt *options, node ** pp, node * q)
{
    node *p=NULL;
    boolean has_migration=FALSE;
    char ch = (char) getc (file);
    //int retval;
    while (ch != ';')
    {
      switch (ch)
        {
	case '(':
	  p = create_interior_node (world, &q);
	  q = p->next;
	  ch = (char) getc (file);
	  break;
	case ',':
	  assert(q->next);			   
	  q = q->next;
	  if (q->top)
	    {
	      usererror ("Multifurcation handling not yet installed");
	    }
	  ch = (char) getc (file);
	  break;
	case ')':
	  p = showtop (q);
	  q = p->back;
	  ch = (char) getc (file);
	  break;
	case ' ':
	case '\n':
	case '\t':
	  ch = (char) getc (file);
	  break;
	case ':':
	  ch = processlength (file, &p);
	  break;
	case '[':
	  has_migration = processbracket (file, world, &p, &ch);
	  if(has_migration>0)
	    {
	      q->back = p;
	      p->back = q;
	    }
	  break;
	default:
	  p = create_tip_node (file, world, options, &q, &ch);
	  break;
        }
    }
    if(p!=NULL)
      {
	    p->length = 10000.;
    	(*pp) = showtop (p->back);
      }
    fscanf (file, "%*[^\n]");
    getc (file);
    return has_migration;
}

void
length_to_times (node * p)
{
    node *q;
    if (p->type != 't')
    {
        length_to_times ((p)->next->back);
        if ((p)->type == 'i')
            length_to_times ((p)->next->next->back);
    }
    q = showtop ((p)->back);
    q->tyme = q->next->tyme = q->next->next->tyme = (p)->tyme + (p)->length;
}

void
find_tips (node * p, node ** nodelist, long *z)
{
    if (p->type == 't')
    {
        nodelist[(*z)++] = p;
    }
    else
    {
        if (p->next->back != NULL)
            find_tips (crawlback (p->next), nodelist, z);
        if (p->next->next->back != NULL)
            find_tips (crawlback (p->next->next), nodelist, z);
    }
}

long
find_firstpop (node * p)
{
    static boolean found = FALSE;
    static long pop = -1;
    if (p->type == 'm'  || p->type == 'd')
    {
        found = TRUE;
        pop = p->pop;
    }
    else
    {
        if (p->next->back != NULL)
        {
            find_firstpop (p->next->back);
            if (found)
                return pop;
        }
        if (p->next->next->back != NULL)
            find_firstpop (p->next->next->back);
    }
    return pop;
}

/* touches only coalescent nodes! migration nodes are already set */
/*
static void set_tree_pop_old (node * p, long *pop)
{
    if (p->type != 'r')
    {
		
        (*pop) =
		(showtop (p->back)->actualpop !=
		 *pop) ? showtop (p->back)->actualpop : *pop;
    }
    p->actualpop = p->pop = *pop;
    if (p->type != 't')
    {
        if (p->next->back != NULL)
        {
            set_tree_pop (crawlback (p->next), pop);
        }
        if (p->type != 'm' && p->type != 'd' && p->next->next->back != NULL)
        {
            set_tree_pop (crawlback (p->next->next), pop);
        }
    }
}
*/

/// sets the actualpop and pop values deducting from the mgiration nodes that
/// are already set
void
set_tree_pop (node * p, long *pop)
{
  static boolean done=FALSE;
  switch(p->type)
    {
    case 'r':
	  set_tree_pop(p->next->back, pop);
	  *pop = p->actualpop = p->pop = p->next->back->pop;
      if(!done)
	{
	  done=TRUE;
	  if(*pop == -1)
	    *pop = 0;
	  set_tree_pop(p,pop);
	}
      break;
    case 't':
      if (*pop != -1)
	p->actualpop = p->pop = *pop;
      break;
    case 'i':
      set_tree_pop(p->next->back,pop);
      if(*pop != -1)
	p->actualpop = p->pop = p->next->back->pop;
      set_tree_pop(p->next->next->back,pop);
      if(*pop != -1)
	{
	  p->actualpop = p->pop = p->next->next->back->pop;
	  *pop = p->pop;
	}
      break;
    case 'm':
    case 'd':
      *pop = p->actualpop;
      set_tree_pop(p->next->back, pop);
      *pop = p->pop;
      break;
    }
}


node *
create_interior_node (world_fmt * world, node ** q)
{
    node *p;
    p = allocate_nodelet (world, 3, 'i');
    p->top = TRUE;
    p->s = (MYREAL *) mycalloc (world->numpop, sizeof (MYREAL));
    p->back = *q;
    if ((*q) == NULL)
      create_root_node (world, &p);
    else
        (*q)->back = p;
    return p;
}

node *
create_root_node (world_fmt *world, node ** q)
{
    node *p;
    p = allocate_nodelet (world, 3, 'r');
    p->top = TRUE;
    p->next->back = *q;
    (*q)->back = p->next;
    return p;
}


node *
create_tip_node (FILE * file, world_fmt * world, option_fmt *options, node ** q, char *ch)
{
    long pop;
    node *p;
    char c;
    char *nayme;
    long i = 1;
    nayme = (char *) mycalloc (options->nmlength+1, sizeof (char));
    nayme[0] = (*ch);
    while (strchr ("[):;,\t\n\r", (int) (c = (char) getc (file))) == NULL)
        nayme[i++] = c;
    nayme[i] = '\0';
    p = allocate_nodelet (world, 1, 't');
    p->nayme = (char *) mycalloc ((2*(options->nmlength+15)), sizeof (char));
    p->truename = p->nayme + 15 + options->nmlength;
    p->top = TRUE;
    p->tip = TRUE;
    p->s = (MYREAL *) mycalloc (world->numpop, sizeof (MYREAL)  );
    for (i = 0; i < world->numpop; i++)
      {
        p->s[i] = MYREAL_MAX;    
      }

    sprintf(p->nayme, "%-*s",(int) options->nmlength,  nayme);
    unpad(p->nayme," _");
    translate(p->nayme,' ', '_');
    strcpy(p->truename,p->nayme);
    sscanf(nayme,"%li ",&pop);		       
    p->s[pop] = 0;
    p->back = *q;
    (*q)->back = p;
    myfree(nayme);
    (*ch) = c;
    return p;
}

boolean
processbracket (FILE * file, world_fmt *world, node ** p, char * ch)
{
    boolean migfound=FALSE;
    long pop1, pop2;
    MYREAL utime;
    char c;
    //int retval;
    c = (char) getc (file);
    if (c == '&')
    {
      c = (char) getc (file);
        switch (c)
        {
			case 'M':
#ifdef USE_MYREAL_FLOAT
				fscanf (file, "%li %li:%f", &pop1, &pop2, &utime);
#else
				fscanf (file, "%li %li:%lf", &pop1, &pop2, &utime);
#endif
				/*c=*/getc (file);
				(*p) = add_migration (world, *p, 'm', pop1, pop2, utime);
				migfound = TRUE;
				break;
			default:
				while (c != ']')
				  c = (char) getc (file);
				break;
        }
    }
    else
    {
        while (c != ']')
	  c = (char) getc (file);
    }
    *ch = (char) getc (file);
    return migfound;
}


node *
add_migration (world_fmt *world, node * p, char event, long from, long to, MYREAL utime)
{
    node *tmp;
    tmp = allocate_nodelet (world, 2, event);
    tmp->top = TRUE;
    tmp->next->back = p;
    p->back = tmp->next;
    tmp->length = p->length - utime;
    p->length = utime;
    tmp->tyme = p->tyme + utime;
    tmp->pop = tmp->next->pop = from;
    tmp->actualpop = tmp->next->actualpop = to;
    //if (event == 'd')
    //  {
    //	fprintf(stderr,"@@%i> divergence: %li->%li %f\n",myID, tmp->pop, tmp->actualpop,tmp->tyme);
    //  }
    return tmp;
}

void
allocate_x (node * p, world_fmt * world, long locus, boolean withtips)
{
    if (p->type != 't')
    {
        if (p->next->back != NULL)
            allocate_x (crawlback (p->next), world, locus, withtips);
        if (p->next->next->back != NULL)
            allocate_x (crawlback (p->next->next), world, locus, withtips);
	alloc_seqx (world, p, locus);
#ifdef UEP
        if(world->options->uep)
        {
            p->uep = (int *) mycalloc (world->data->uepsites, sizeof (int));
            p->ux.s = (pair *) mycalloc (world->data->uepsites, sizeof (pair));
        }
#endif
		
    }
    else
    {
      if (withtips)
        {
	  alloc_seqx (world, p, locus);
        }
      //#ifdef UEP
      //      if(world->options->uep)
      // {
      //   p->uep = (long *) mycalloc (world->data->uepsites, sizeof (long));
      //   p->ux.s = (pair *) mycalloc (world->data->uepsites, sizeof (pair));
      // }
      //#endif
      
    }
}

long
number_genomes (int datatype)
{
    switch (datatype)
    {
    case 'a':
    case 'b':
    case 'm':
      return 2;
    case 's':
    case 'n':
    case 'h':
    case 'u':
    case 'f':
    case '@':
      return 1;
    default:
      error ("Wrong data type");
    }
    //return 0;
	
}


void
copy_tree (world_fmt * original, world_fmt * kopie)
{
    kopie->root = copy_node (original, original->root, kopie, NULL);
}

node *
copy_node (world_fmt * original, node * o, world_fmt * kopie, node * last)
{
    static long z = 0;
	
    node *t = NULL, *t2, *t3;
    if (o == NULL)
        return NULL;
    if (!o->top)
        error ("copy_tree messed up");
	
    switch (o->type)
    {
    case 'r':
      z = 0;
    case 'i':
      t = (node *) mycalloc (1, sizeof (node));
      t2 = (node *) mycalloc (1, sizeof (node));
      t3 = (node *) mycalloc (1, sizeof (node));
      t->next = t2;
      t2->next = t3;
      t3->next = t;
      copy_node_content (original, kopie, o, t);
      copy_node_content (original, kopie, o->next, t2);
      copy_node_content (original, kopie, o->next->next, t3);
      if (o->next->back != NULL)
	t2->back = copy_node (original, o->next->back, kopie, t2);
      if (o->next->next->back != NULL)
	t3->back = copy_node (original, o->next->next->back, kopie, t3);
      t->back = last;
      break;
    case 'd':
    case 'm':
      t = (node *) mycalloc (1, sizeof (node));
      t2 = (node *) mycalloc (1, sizeof (node));
      t->next = t2;
      t2->next = t;
      copy_node_content (original, kopie, o, t);
      copy_node_content (original, kopie, o->next, t2);
      t2->back = copy_node (original, o->next->back, kopie, t2);
      t->back = last;
      break;
    case 't':
      t = (node *) mycalloc (1, sizeof (node));
      //kopie->nodep[z++] = t;
      t->next = t;
      copy_node_content (original, kopie, o, t);
      t->back = last;
      break;
    }
    return t;
}

///
/// copy the guts of a node in the tree (but to where?)
void
copy_node_content (world_fmt * original, world_fmt * kopie, node * o,
                   node * t)
{
    long i;
//    long j;
    long endsite;
    long rcategs;
    t->type = o->type;
    t->number = o->number;
    t->pop = o->pop;
    t->actualpop = o->actualpop;
    t->id = o->id;
    t->top = o->top;
    t->dirty = o->dirty;
    t->v = o->v;
    t->tyme = o->tyme;
    t->length = o->length;
	
    if (t->top && t->type != 'm' && t->type != 'd')
    {
#ifdef UEP
		
        if (original->options->uep)
        {
            t->uep = (int *) mycalloc (original->data->uepsites, sizeof (int));
            memcpy (t->uep, o->uep, sizeof (long) * original->data->uepsites);
            t->ux.s = (pair *) mycalloc (original->data->uepsites, sizeof (pair));
            for (i = 0; i < original->data->uepsites; i++)
            {
                memcpy (t->ux.s[i], o->ux.s[i], sizeof (pair));
            }
        }
#endif
	mutationmodel_fmt *s;
	long sublocus;
	long locus = original->locus;
	const long sublocistart = original->sublocistarts[locus];
	const long sublociend   = original->sublocistarts[locus+1];
	for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
	  {
	    s = &original->mutationmodels[sublocus];
	    endsite = s->numpatterns + s->addon;
	    rcategs = s->numsiterates;
	    const long xs = sublocus - sublocistart;
	    if (strchr (SEQUENCETYPES, s->datatype))
	      {
		alloc_seqx (kopie, t, locus);
		memcpy (t->scale[xs], o->scale[xs], sizeof (MYREAL) * (size_t) endsite);
		for (i = 0; i < endsite; i++)
		  {
		    memcpy (t->x[xs].s[i], o->x[xs].s[i], (size_t) rcategs * sizeof (sitelike));
		  }
	      }
	    else
	      {
		t->x[xs].a =
		  (MYREAL *) mycalloc (1,
				       s->maxalleles *
				       sizeof (MYREAL));
		memcpy (t->x[xs].a, o->x[xs].a,
			sizeof (MYREAL) *
			(size_t) s->maxalleles);
	      }
	  }
    }
    //if (o->s != NULL)
    ///{
    //t->s = (MYREAL *) mycalloc(1, sizeof(MYREAL) * original->numpop);
    //memcpy] (t->s, o->s, sizeof(MYREAL) * original->numpop);
    //
    //      }
    //
    //else
    //    t->s = NULL;
    if (o->type == 't')
      {
	if (o->nayme != NULL)
	  {
	    if (t->nayme == NULL)
	      t->nayme = (char*) mycalloc(2*(80 + 15),sizeof(char));
	    strcpy(t->nayme,o->nayme);
	  }
	else
	  t->nayme = NULL;
	if (o->truename != NULL && o->nayme != NULL)
	  {
	    t->truename = t->nayme + 15 + 80;
	    strcpy(t->truename,o->truename);
	    //printf("%i> copy_node_content: type=%c truename=%s\n",myID, t->type,t->truename);
	  }
      }
    else
      t->truename = NULL;
}

void
swap_tree (world_fmt * tthis, world_fmt * tthat)
{
    node *tmp;
    node **nodetemps;
    tmp = tthis->root;
    tthis->root = tthat->root;
    tthat->root = tmp;
    nodetemps = tthis->nodep;
    tthis->nodep = tthat->nodep;
    tthat->nodep = nodetemps;

}

void
calc_treelength (node * p, MYREAL *treelen)
{
    node *pn, *pnn;
    switch (p->type)
    {
    case 't':
      break;
    case 'd':
    case 'm':
      error ("yelp\n");
      //break;
    case 'i':
      pn = crawlback (p->next);
      calc_treelength (pn, treelen);
      pnn = crawlback (p->next->next);
      calc_treelength (pnn, treelen);
      break;
    default:
      error ("default reached");
    }
    pn = showtop (crawlback (p));
    if (pn->type != 'r')
        *treelen += pn->tyme - p->tyme;
}

MYREAL calc_pseudotreelength (proposal_fmt * proposal, MYREAL treelen)
{
    MYREAL len = 0.0;
    MYREAL ot = proposal->origin->tyme;
    MYREAL obt = proposal->oback->tyme;
    MYREAL tt = proposal->target->tyme;
    MYREAL rt = proposal->world->root->next->back->tyme;
    //target is not root
    if (proposal->target != proposal->world->root)
    {
        //oback is not root
        if (proposal->oback != proposal->world->root)
        {
            len = treelen - (obt - ot) + (proposal->time - ot);
            //printf("pseudo_treelen: ob!=r t!=r %f\n", len);
        }
        else
        {
            //oback is root
            len = treelen - (obt - ot) - (rt - tt) +
			(proposal->time - ot) + (proposal->time - tt);
            //printf("pseudo_treelen: ob=r t!=r %f\n", len);
        }
    }
    else
        //target is root
    {
        //oback is not root
        if (proposal->oback != proposal->world->root)
        {
            len = treelen - (obt - ot) + (proposal->time - ot) - (rt - tt) +
			(proposal->time - tt);
            //printf("pseudo_treelen: ob!=r t=r %f\n", len);
        }
        else
        {
            //oback is root
            len = treelen - (obt - ot) - (obt - tt) + (proposal->time - ot)
			+ (proposal->time - tt);
            //printf("pseudo_treelen: ob=r t=r %f\n", len);
        }
    }
    return len;
}


void swap (contribarr *a, contribarr *b)
{
    contribarr *t;
    t = a;
    a = b;
    b = t;//clang static analyzer claims this assignment is not doing anything, but it does!
}


void
free_tree (node * p, world_fmt * world)
{
    if (p != NULL)
    {
        if (p->type != 't')
        {
            if (p->next->back != NULL)
            {
                free_tree (p->next->back, world);
            }
            if (p->type != 'm' && p->type != 'd' && p->next->next->back != NULL)
            {
                free_tree (p->next->next->back, world);
            }
        }
        switch (p->type)
        {
	case 'd':
	case 'm':
	  free_mignodelet (p, world);
	  break;
	case 't':
	  free_tipnodelet (p, world);
	  break;
	case 'i':
	  free_nodelet (p, 3, world);
	  break;
	case 'r':
	  free_nodelet (p, 3, world);
	  world->root = NULL;
	  break;
	default:
	  error("error in freeing nodes, a node with no type found");
	  //  break;
        }
    }
}

///
/// delete tipnode-nodelet and its content data
void
free_tipnodelet (node * p, world_fmt * world)
{
  free_nodedata (p, world);
#ifdef DISPENSER
  collect_nodelet(world, p);
#else
  myfree(p); 
#endif
  //p->id = -p->id;
}

///
/// delete migratenode-nodelets
void free_mignodelet (node * p, world_fmt * world)
{
  (void) world;
  node *q = p->next;
#ifdef DISPENSER
  collect_nodelet(world, q);
  collect_nodelet(world, p);
#else
  myfree(p);
  myfree(q);
#endif
  //  myfree(q);
  //myfree(p);
}

///
/// delete interior node and its data
void
free_nodelet (node * p, long num, world_fmt * world)
{
    long i;
    node *q , *r;
    switch(num)
      {
      case 1:
	if(p->top == TRUE)
	  free_nodedata (p, world);
#ifdef DISPENSER
	collect_nodelet(world, p);
#else
	myfree(p); 
#endif
	break;
      case 2:	
	q = p->next;
	if(p->top == TRUE)
	  free_nodedata (p, world);
#ifdef DISPENSER
	collect_nodelet(world, p);
#else
	myfree(p); 
#endif
	if(q!=NULL)
	  {
	    if(q->top == TRUE)
	      free_nodedata (q, world);
#ifdef DISPENSER
	    collect_nodelet(world, q);
#else
	    myfree(q);
#endif
	  }
	break;
      case 3:	
	q = p->next;
	r = p->next->next;
	if(p->top == TRUE)
	  free_nodedata (p, world);
#ifdef DISPENSER
	collect_nodelet(world, p);
#else
	myfree(p);
#endif
	//p->id = -p->id;
	if(q!=NULL)
	  {
	    if(q->top == TRUE)
	      free_nodedata (q, world);
#ifdef DISPENSER
	    collect_nodelet(world, q);
#else
	    myfree(q);
#endif
	    //q->id = -q->id;
	  }
	if(r!=NULL)
	  {
	    if(r->top == TRUE)
	      free_nodedata (r, world);
#ifdef DISPENSER
	    collect_nodelet(world, r);
#else
	    myfree(r);
#endif
	    //r->id = -r->id;
	  }
	break;
      default:
	for (i = 0; i < num; i++)
	  {
	    q = p->next;
	    if(p->top == TRUE)
	      free_nodedata (p, world);
#ifdef DISPENSER
	    collect_nodelet(world, p);
#else
	    myfree(p);
#endif
	    p = q;
	  }
      }
}

void
free_nodedata (node * p, world_fmt * world)
{
  mutationmodel_fmt *s;
  long xs;
  long locus = world->locus;
  const long start = world->sublocistarts[locus];
  const long sublociend   = world->sublocistarts[locus+1] - start;
  //boolean type;
  if (p->nayme != NULL)
    myfree(p->nayme);

  if (p->x !=NULL)
    {
      for(xs=0;xs<sublociend;xs++)
	{
	  s = &world->mutationmodels[start+xs];
	  if (strchr (SEQUENCETYPES, s->datatype))
	    {
	      if(p->x[xs].s != NULL)
		{
		  myfree(p->x[xs].s[0]);
		  myfree(p->x[xs].s);
		}
	    }
	  else
	    {
	      if(p->x[xs].a != NULL)
		{
		  myfree(p->x[xs].a);
		}
	    }
	}
      myfree(p->x);
    }
  if (p->scale != NULL)
    {
      for(xs=0;xs<sublociend;xs++)
	{
	  myfree(p->scale[xs]);
	}
      myfree(p->scale);
    }
  if(p->s!=NULL)
    myfree(p->s);
#ifdef UEP
  if (world->options->uep)
    myfree(p->uep);
#endif
}

///
/// calculates the 1/probability that no event happens
/// is safeguarded against a problem on two tip trees 
/// that will return a result of zero and so a inverse of INF
/// INF is replaced by MYREAL_MAX, but a large value 100000.0 because this function
/// ignores the speciation time, since this is only used to construct a first tree
/// for the moment we ignore this.
MYREAL
inverse_logprob_noevent (world_fmt * world, long interval)
{
    long pop, k;
    MYREAL result = 0.0;
    for (pop = 0; pop < world->numpop; pop++)
    {
        k = world->treetimes[0].tl[interval].lineages[pop];
        result +=
            (k * (k - 1) / world->param0[pop]) + sum_migprob (world, pop, interval);
    }
    if(result > 0.0)
      return 1./result;
    else
      {
	if(world->treetimes[0].tl[interval].eventnode!=NULL && world->treetimes[0].tl[interval].eventnode->type=='t')
	  {
	    return world->treetimes[0].tl[interval].age;
	  }
	else
	  return MYREAL_MAX;
      }
}

/// calculates the sum of all immigrations into a specific population
MYREAL
sum_migprob (world_fmt * world, long pop, long interval)
{
    long i;
    MYREAL result = 0.0;
    long *lineages = world->treetimes[0].tl[interval].lineages;
    long msta = world->mstart[pop];
    long msto = world->mend[pop];
    boolean usem = world->options->usem;
    MYREAL pk;
    for (i = msta; i < msto; i++)
    {
      pk = usem ? world->param0[i] : world->param0[i]/world->param0[pop];
      result += pk;
    }
    return result * lineages[pop];
}

void debugline(node *up)
{
  node *thenode = up;
  node * down = crawlback(showtop(up));
  while (thenode != down)
    {
      printf("<%li:%li %c>\n",thenode->id, showtop(thenode)->id, thenode->type);
      thenode = showtop(thenode)->back;
    }
}

#ifdef NEXUSTREE
void nexus_treesheader(world_fmt *world, option_fmt *options, data_fmt *data)
{
  long ii;
  long ind;
  long pop;
  long start;
  long stop;
  FPRINTF(world->treefile,"#nexus\n\n");
  FPRINTF(world->treefile,"begin taxa;\n");
  FPRINTF(world->treefile,"dimensions ntax=%li;\n",world->sumtips);
  FPRINTF(world->treefile,"taxlabels\n");
  for (pop=0; pop< world->numpop; pop++)
    {
      long top = max_shuffled_individuals(options, data, pop, 0);
      for (ii = 0; ii < top; ii++)
	{
	  ind = data->shuffled[pop][0][ii];
	  FPRINTF(world->treefile,"        %10.10s\n", data->indnames[pop][ind][0]); //assumed locus 0 	
	}
    }
  FPRINTF(world->treefile,";\nend;\n\nbegin sets;");
  start = 1;
  for (pop=0;pop<world->numpop;pop++)
    {
      stop = start-1 + max_shuffled_individuals(options, data, pop, 0); 
      FPRINTF(world->treefile,"taxset deme%li = %li-%li\n",pop,start, stop);
      start = stop+1;
    }
  FPRINTF(world->treefile,"end;\n\nbegin trees;\n");
}
#endif
