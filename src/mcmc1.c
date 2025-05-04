/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effective population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 M C M C   R O U T I N E S 
 
 Markov Monte Carlo stuff: treechange, acceptance
 
 
 Peter Beerli 1996, Seattle
 beerli@fsu.edu
 
 Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
 Copyright 2003-2013 Peter Beerli, Tallahassee FL
 
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
 
 $Id: mcmc1.c 2169 2013-08-24 19:02:04Z beerli $
 -------------------------------------------------------*/
/* \file mcmc1.c
 
 Tree changer and acceptance rejection scheme
 
 */
#include "migration.h"
#include "sighandler.h"
#include "tools.h"
#include "random.h"
#include "tree.h"
#include "mcmc2.h"
#include "speciate.h"
#include "world.h"
#include "mittag_leffler.h"
#ifdef UEP
#include "uep.h"
#endif

#include "bayes.h"
#include "assignment.h"

#ifdef BEAGLE
#include "calculator.h"
#endif
#include <assert.h>
#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif

//#define MIGRATION_AIR (boolean) 1
//#define MIGRATION_IN_TREE (boolean) 0
#define NO_MIGR_NODES 0
#define WITH_MIGR_NODES 1

//boolean debugtip;

extern int myID;
extern MYREAL probg_treetimesX(world_fmt *world, vtlist *tl, long T);

extern void     zero_xseq(xarray_fmt *x, long sites, long categs);

/* prototypes ------------------------------------------- */
// metropolize over trees
long tree_update (world_fmt * world, long g, boolean assign);
/* private functions */
void new_localtimelist (timelist_fmt ** ntl, timelist_fmt * otl, long numpop);
void new_proposal (proposal_fmt ** proposal, timelist_fmt * tl,
                   world_fmt * world);
//CLANG_ANALYZER_NORETURN
void set_new_proposal (proposal_fmt ** proposal, timelist_fmt * tl, world_fmt * world);
void chooseOrigin (proposal_fmt * proposal);
void construct_localtimelist (timelist_fmt * timevector,
                              proposal_fmt * proposal);
void traverseAllNodes (node * theNode, node *** nodelist, long *node_elem,
                       long *oldnode_elem, int include_migration, div_fmt **divlist);
void chooseTarget (proposal_fmt * proposal, timelist_fmt * timevector,
                   node ** bordernodes, long *bordernum);
void findbordernodes (node * theNode, proposal_fmt * proposal, long pop,
                      node *** bordernodes, long *allocsize, long *bordernum, vtlist ** tyme,
                      long gte);
void free_masterproposal (proposal_fmt * proposal);
void free_timevector (timelist_fmt * timevector);
void prune_timelist (timelist_fmt * oldtv, timelist_fmt * newtv, register node ** __restrict ptr, proposal_fmt * proposal, long numpop);

long xor (node ** ptrl1, node ** ptrl2);
long rmigrcount (proposal_fmt * proposal);

int migrate (proposal_fmt * proposal, node * up, char event, long from, long to);
int migrateb (proposal_fmt * proposal, node * up, char event, long from, long to);

//int pre_population (proposal_fmt * proposal, vtlist * ltime, long gte,
//                    long *slider); this was replaced by beyond_last_node() in speciate.c

boolean acceptlike (world_fmt * world, proposal_fmt * proposal, long g,
                    timelist_fmt * tyme, boolean assign);
MYREAL eventtime (proposal_fmt * proposal, long pop, vtlist * tentry,
                  char *event);
node *showsister (node * theNode);
void count_migrations (node * p, long *count);
long migration_from (long to, proposal_fmt * proposal);

MYREAL prob_tree (world_fmt * world, timelist_fmt * tyme);
void traverse_check(node *theNode);
void reset_proposal (proposal_fmt ** proposal, world_fmt *world);
void set_tree_dirty (node * p);
void set_tree_down_dirty (node * p);
void new_localtimelist_new (timelist_fmt ** ntl, timelist_fmt * otl, long numpop);
void jumblenodes (node ** s, long n);
void traverse_tagnew(node *theNode, node *origin);
void record_above_origin_divergence(div_fmt **divlist, long from, long to);
void allocate_nodelist(node ***nodelist, long *oldnode_elem, long elements);
void test_traverse_error(node * theNode);
void findbordernodes_recursive (node * theNode, proposal_fmt * proposal, long pop,
                 node *** bordernodes, long *allocsize, long *bordernum, vtlist ** tyme,
                 long gte);
void add_node_border(node *tmp, node ***bordernodes,long *bordernum, long *allocsize);
boolean is_in_bracket(node *tmp, node*back, MYREAL thetime, long pop);
void increase_stack(node ***stack, long *allocsize,  long newallocsize);
void free_timevector_new (timelist_fmt * timevector);

void add_proposal_to_tl(proposal_fmt *proposal, timelist_fmt ** timevector);

//##
/* Functions implementation ++++++++++++++++++++++++++++++++++++++++++++++++*/
void set_tree_dirty (node * p)
{
  node *circle=NULL;
    switch (p->type)
    {
    case 'd':
    case 'm':
      set_dirty (p);
      set_tree_dirty (p->next->back);
      break;
    case 't':
      break;
    case 'i':
      set_dirty (p);
      circle = p->next;
      while(!circle->top)
	{
	  set_tree_dirty(circle->back);
	  circle = circle->next;
	}
      // set_tree_dirty (p->next->back);
      //set_tree_dirty (p->next->next->back);
      break;
    case 'r':
      set_dirty (p);
      set_tree_dirty (p->next->back);
      break;
    }
}

void set_tree_down_dirty (node * p)
{
  node *q=NULL;
  switch (p->type)
    {
    case 'm':
    case 'd':
      q = showtop(p);
      assert(q!=NULL);
      set_tree_down_dirty (q->back);
      break;
    case 't':
      set_tree_down_dirty (p->back);
      break;
    case 'i':
      q=showtop(p);
      set_dirty (q);
      set_tree_down_dirty (q->back);
      break;
    case 'r':
      q = showtop(p);
      set_dirty (q);
      break;
    default:
      error("p->type is an undefined  node type");
    }
}

/*=======================================================*/
long
tree_update (world_fmt * world, long g, boolean assign)
{    
/*
return 1 if tree was accepted, 0 otherwise 
assign is set for assigning individuals, caller is responsible to reset the 
origin back to the original state when fail
*/
    if (assign && world->unassignednum<2)
      return 0;
    return newtree_update(world,g,assign);
}


#ifdef SLATKIN_IMPORTANCE
long slatkin_importance(world_mt *world, long g)
{
    long slice=0;
    // create timelist with times drawn using the parameters
    new_localtimelist (&timevector, &world->treetimes[0], world->numpop);
    while (slice < sumlineages)
    {
        actualpop =
        (proposal->migr_table_counter >
         0) ? proposal->migr_table[proposal->migr_table_counter -
                                   1].from : proposal->origin->pop;
        age = (*timevector).tl[slice].age;
        (*timevector).tl[slice].age = age + eventtime (proposal, actualpop, tentry, &event);
    }
    // calculate which two datapoints to join give the distance
}
#endif

/*=======================================================*/
///
/// allocate the timelist, contains  times and lineages at that time
/// lineages have a large storage device that is accessed from the 
/// tl[i].lineages.
void
new_localtimelist (timelist_fmt ** ntl, timelist_fmt * otl, 
long numpop)
{
    (*ntl) = (timelist_fmt *) mycalloc (1, sizeof (timelist_fmt));
    (*ntl)->tl = (vtlist *) mymalloc ((*otl).allocT * sizeof (vtlist));
    (*ntl)->allocT = otl->allocT;
    //(*ntl)->allocT = otl->T+2;
    //printf("%i> alloc new_localtimelist: %li\n",myID,(*ntl)->allocT);
    (*ntl)->T = otl->T;
    (*ntl)->oldT = otl->oldT;
    //memcpy ((*ntl)->tl, otl->tl, otl->allocT * sizeof (vtlist));
    allocate_lineages (ntl, 0, numpop);
    //memcpy ((*ntl)->lineages, otl->lineages, (*ntl)->allocT * numpop * sizeof (long));
}

///
/// reuse a global timelist to store values, this should save 
/// malloc/free calls
void new_localtimelist_new (timelist_fmt ** ntl, timelist_fmt * otl, long numpop)
{
    // UNFINISHED
    //    long i;
    //    (*ntl) = (timelist_fmt *) mycalloc (1, sizeof (timelist_fmt));
    if((*ntl)->allocT < otl->allocT)
    {
        (*ntl)->tl = (vtlist *) myrealloc ((*ntl)->tl,(*otl).allocT * sizeof (vtlist));
        //      memset((*ntl)->tl+((*ntl)->allocT),0,sizeof(vtlist)*((*otl)->allocT-(*ntl)->allocT));
        (*ntl)->allocT = otl->allocT;
    }
    (*ntl)->T = otl->T;
    memcpy ((*ntl)->tl, otl->tl, (size_t) otl->allocT * sizeof (vtlist));
    allocate_lineages (ntl, 0, numpop);
    memcpy ((*ntl)->lineages, otl->lineages, (size_t) ((*ntl)->allocT * numpop) * sizeof (long));
}



void
new_proposal (proposal_fmt ** proposal, timelist_fmt * tl, world_fmt * world)
{
  (void) tl;
  const long np = world->numpop2 + world->bayes->mu + world->species_model_size * 2 + world->grownum;
#ifdef UEP    
    long j;
#endif
    long listsize = 2*(world->sumtips + 2);
    long sumtips = world->sumtips;
    mutationmodel_fmt *s;
    // allocate the scratchpad (contains a shadow of the tree)
    (*proposal) = (proposal_fmt *) mycalloc (1, sizeof (proposal_fmt));
    (*proposal)->likelihood = (double) -HUGE;
    // pointers and values to outside structures
    (*proposal)->world = world;
    (*proposal)->sumtips = world->sumtips;
    (*proposal)->numpop = world->numpop;
    (*proposal)->param0 = world->param0;
    (*proposal)->root = world->root;
    (*proposal)->migration_model = world->options->migration_model;
    // precalculated values
    (*proposal)->mig0list = world->mig0list;
    (*proposal)->design0list = world->design0list;
    
    (*proposal)->aboveorigin_allocsize = listsize;    
    // nodes above the picked node + migration nodes
    (*proposal)->aboveorigin = (node **) mycalloc (listsize, sizeof (node *));
    // list to record divergences above picked node
    (*proposal)->divlist = (div_fmt *) mycalloc (1, sizeof (div_fmt));
    (*proposal)->divlist->divlist = (longpair *) mycalloc (listsize, sizeof (longpair));
    (*proposal)->divlist->div_allocsize = listsize;
    (*proposal)->param0save = (MYREAL *) mycalloc (np, sizeof (MYREAL));
    // node data holding vector for bordernodes, line_f, line_t
    //(*proposal)->nodedata =
    //(node **) mycalloc (3 * listsize, sizeof (node *));
    // adjacent nodes of the picked nodes
    (*proposal)->bordernodes_allocsize = listsize;    
    (*proposal)->bordernodes = (node **) mycalloc (listsize, sizeof (node *)); 
    (*proposal)->line_f_allocsize = listsize;    
    (*proposal)->line_f =      (node **) mycalloc (listsize, sizeof (node *));
    (*proposal)->line_t_allocsize = listsize;    
    (*proposal)->line_t =      (node **) mycalloc (listsize, sizeof (node *));
    (*proposal)->mf = (MYREAL **) mycalloc(world->numsubloci[world->locus],sizeof(MYREAL*));
    (*proposal)->mt = (MYREAL **) mycalloc(world->numsubloci[world->locus],sizeof(MYREAL*));
    long sublocus;
    long endsite;
    //long rcategs;
    long sublocistart = world->sublocistarts[world->locus];
    long sublociend = world->sublocistarts[world->locus+1];
    (*proposal)->xf = (xarray_fmt *) mycalloc(world->numsubloci[world->locus],sizeof(xarray_fmt));
    (*proposal)->xt = (xarray_fmt *) mycalloc(world->numsubloci[world->locus],sizeof(xarray_fmt));
    for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
      {
	s = &world->mutationmodels[sublocus];
	const long xs = sublocus - sublocistart;
	if (strchr (SEQUENCETYPES, s->datatype))
	  {
	    endsite = s->numpatterns + s->addon;
	    (*proposal)->mf[xs] = (MYREAL *) mycalloc(endsite,sizeof(MYREAL));
	    (*proposal)->mt[xs] = (MYREAL *) mycalloc(endsite,sizeof(MYREAL));
	    allocate_xseq(&(*proposal)->xf[xs], endsite, s->numsiterates);
	    allocate_xseq(&(*proposal)->xt[xs], endsite, s->numsiterates);
	  }
	else
	  {
	    long mal = s->maxalleles + 1;
	    (*proposal)->mf[xs] = (MYREAL *) mycalloc(mal,sizeof(MYREAL));
	    (*proposal)->mt[xs] = (MYREAL *) mycalloc(mal,sizeof(MYREAL));
	    (*proposal)->xf[xs].a = (MYREAL *) mycalloc (mal, sizeof (MYREAL));
	    (*proposal)->xt[xs].a = (MYREAL *) mycalloc (mal, sizeof (MYREAL));
	  }
      }
    (*proposal)->old_migr_table_counter = 4 * sumtips /* 100 */;
    (*proposal)->old_migr_table_counter2 = 4 * sumtips /* 100 */;
    (*proposal)->migr_table =
    (migr_table_fmt *) mycalloc ((*proposal)->old_migr_table_counter,
                                 sizeof (migr_table_fmt));
                                 
    (*proposal)->migr_table2 = (migr_table_fmt *) mycalloc ((*proposal)->old_migr_table_counter2,
							    sizeof (migr_table_fmt));
    (*proposal)->migr_table_counter = 0;
    (*proposal)->migr_table_counter2 = 0;
#ifdef UEP
    if (world->options->uep)
    {
        (*proposal)->ueplike =
        (MYREAL **) mycalloc (world->data->uepsites, sizeof (MYREAL *));
        (*proposal)->ueplike[0] =
        (MYREAL *) mycalloc (world->numpop * world->data->uepsites,
                             sizeof (MYREAL));
        for (j = 1; j < world->data->uepsites; ++j)
            (*proposal)->ueplike[j] = (*proposal)->ueplike[0] + j * world->numpop;
        
        (*proposal)->ut.s = (pair *) mycalloc (world->data->uepsites, sizeof (pair));
        (*proposal)->uf.s = (pair *) mycalloc (world->data->uepsites, sizeof (pair));
        (*proposal)->umt = (MYREAL *) mycalloc (world->data->uepsites, sizeof (MYREAL));
        (*proposal)->umf = (MYREAL *) mycalloc (world->data->uepsites, sizeof (MYREAL));
    }
#endif
    
#ifdef BEAGLE
    (*proposal)->leftid   = 0;
    (*proposal)->rightid  = 0;
    (*proposal)->parentid = 0;
    reset_beagle(world->beagle);
#endif
}

///
/// sets new proposal structure using template gproposal, this function will be called most of the time instead of new_proposal()
void
set_new_proposal (proposal_fmt ** proposal, timelist_fmt * tl, world_fmt * world)
{
  (void) tl;
  (void) proposal;
  (void) world;
  //long mal = world->data->maxalleles[world->locus];
  //long oldsize=0;
  error("needs fixing");
#if 0
  long newsize =0;
  mutationmodel_fmt *s;
  long sumtips = world->sumtips;
#ifdef UEP    
    long j;
#endif
    long listsize = 2*(world->sumtips + 2);
    (*proposal)->likelihood = -HUGE;
    (*proposal)->sumtips = sumtips;
    // pointers and values to outside structures
    (*proposal)->world = world;
    (*proposal)->datatype = world->options->datatype;
    (*proposal)->numpop = world->numpop;
    (*proposal)->endsite = world->data->seq[0]->endsite;
    (*proposal)->fracchange = world->data->seq[0]->fracchange;
    (*proposal)->param0 = world->param0;
    (*proposal)->root = world->root;
    (*proposal)->migration_model = world->options->migration_model;
    // precalculated values
    (*proposal)->mig0list = world->mig0list;
    (*proposal)->design0list = world->design0list;
    newsize = 4 * listsize;
    if((*proposal)->nodedata == NULL)
        (*proposal)->nodedata = (node **) mycalloc (newsize, sizeof (node *));
    else
        (*proposal)->nodedata = (node **) myrealloc ((*proposal)->nodedata, newsize * sizeof (node *));
    //xcode  oldsize = (*proposal)->listsize;
    memset((*proposal)->nodedata, 0, (size_t) newsize * sizeof (node *));
    (*proposal)->aboveorigin = (*proposal)->nodedata;
    (*proposal)->bordernodes = (*proposal)->nodedata + listsize;
    (*proposal)->line_f =   (*proposal)->bordernodes + listsize;
    (*proposal)->line_t_allocsize = listsize;
    (*proposal)->line_t =  (*proposal)->line_f + listsize;
    // mf holds also mt array  

    (*proposal)->mf = (MYREAL **) mycalloc((world->numsubloci[world->locus]*2),sizeof(MYREAL*));
    (*proposal)->mt = (*proposal)->mf + world->numsubloci[world->locus];
    long sublocus;
    long endsite;
    //long rcategs;                                                                                                 
    long sublocistart = world->sublocistarts[world->locus];
    long sublociend = world->sublocistarts[world->locus+1];
    (*proposal)->xf = (xarray_fmt *) mycalloc(world->numsubloci[world->locus],sizeof(xarray_fmt));
    (*proposal)->xt = (xarray_fmt *) mycalloc(world->numsubloci[world->locus],sizeof(xarray_fmt));
    for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
      {
        s = &world->mutationmodels[sublocus];
        const long xs = sublocus - sublocistart;
        if (strchr (SEQUENCETYPES, s->datatype))
          {
            endsite = s->numpatterns + s->addon;
            (*proposal)->mf[xs] = (MYREAL *) mycalloc(endsite,sizeof(MYREAL));
            (*proposal)->mt[xs] = (MYREAL *) mycalloc(endsite,sizeof(MYREAL));
            allocate_xseq(&(*proposal)->xf[xs], endsite, s->numsiterates);
            allocate_xseq(&(*proposal)->xt[xs], endsite, s->numsiterates);
          }
        else
          {
            long mal = s->maxalleles + 1;
            (*proposal)->mf[xs] = (MYREAL *) mycalloc(mal,sizeof(MYREAL));
            (*proposal)->mt[xs] = (MYREAL *) mycalloc(mal,sizeof(MYREAL));
            (*proposal)->xf[xs].a = (MYREAL *) mycalloc (mal, sizeof (MYREAL));
            (*proposal)->xt[xs].a = (MYREAL *) mycalloc (mal, sizeof (MYREAL));
          }
      }
    (*proposal)->old_migr_table_counter = 4 * sumtips /* 100 */;
    (*proposal)->old_migr_table_counter2 = 4 * sumtips /* 100 */;
    (*proposal)->migr_table =
      (migr_table_fmt *) mycalloc ((*proposal)->old_migr_table_counter,
                                 sizeof (migr_table_fmt));
                                 
    (*proposal)->migr_table2 =
    (migr_table_fmt *) mycalloc ((*proposal)->old_migr_table_counter2,
                                 sizeof (migr_table_fmt));
    (*proposal)->migr_table_counter = 0;
    (*proposal)->migr_table_counter2 = 0;
#ifdef UEP
    if (world->options->uep)
    {
        (*proposal)->ueplike =
        (MYREAL **) mycalloc (world->data->uepsites, sizeof (MYREAL *));
        (*proposal)->ueplike[0] =
        (MYREAL *) mycalloc (world->numpop * world->data->uepsites,
                             sizeof (MYREAL));
        for (j = 1; j < world->data->uepsites; ++j)
            (*proposal)->ueplike[j] = (*proposal)->ueplike[0] + j * world->numpop;
        
        (*proposal)->ut.s = (pair *) mycalloc (world->data->uepsites, sizeof (pair));
        (*proposal)->uf.s = (pair *) mycalloc (world->data->uepsites, sizeof (pair));
        (*proposal)->umt = (MYREAL *) mycalloc (world->data->uepsites, sizeof (MYREAL));
        (*proposal)->umf = (MYREAL *) mycalloc (world->data->uepsites, sizeof (MYREAL));
    }
#endif
    
#ifdef BEAGLE
    (*proposal)->leftid   = 0;
    (*proposal)->rightid  = 0;
    (*proposal)->parentid = 0;
    reset_beagle(world->beagle);
#endif
    world->has_proposal_first=FALSE;
#endif /* if zero */
}


void jumblenodes (node ** s, long n)
{
    node **temp;
    
    long i, rr, tn = n;
    
    temp = (node **) mycalloc (n, sizeof (node *));
    memcpy (temp, s, sizeof (node *) * (size_t) n);
    for (i = 0; i < n && tn > 0; i++)
    {
        s[i] = temp[rr = RANDINT (0L, tn - 1)];
        temp[rr] = temp[tn - 1];
        tn--;
    }
    myfree(temp);
}

void traverse_check(node *theNode)
{
  node *right = theNode->next;
  if (theNode != NULL)
    {
      if (theNode->type != 't')
        {
	  while(!right->top)
	    {
	      if(right->back!=NULL)		
		//	  if (theNode->next->back != NULL)
		{
		  if(theNode->tyme < showtop(right->back)->tyme)
		    {
		      printf("Problem in traverse_check: id=%li time=%f type=%c time-up-next=%f type_up=%c\n",
			     theNode->id, theNode->tyme, theNode->type, right->back->tyme, right->back->type);
		      error("time conflict");
		    }
		  traverse_check (right->back);
		}
	      right = right->next;
	    }
	  //            if (theNode->type != 'm' && theNode->type != 'd'&& theNode->next->next->back != NULL)
	  //{
	  //    if(theNode->tyme < theNode->next->next->back->tyme)
	  //    {
	  //        printf("Problem in traverse_check: id=%li time=%f type=%c time-up-next-next=%f\n",
	  //               theNode->id, theNode->tyme, theNode->type, theNode->next->next->back->tyme);
	  //        error("time conflict");
	  //    }
	  //    
	  //    traverse_check (theNode->next->next->back);
	  //}
        }
    }
    else
    {
        error("Node is NULL????");
    }
}

///
/// tags the lineage that was added the last time the tree was changed.
void traverse_tagnew(node *theNode, node *origin)
{
  node *right=NULL;
    if (theNode != NULL)
    {
      if(!theNode->top)
	theNode = showtop(theNode);
      theNode->visited = FALSE;
      if (theNode->type != 't')
        {
	  right = theNode->next;
	  while(!right->top)
	    {		
	      traverse_tagnew (right->back, origin);       
	      //if (theNode->type != 'm' && theNode->type != 'd' && theNode->type && theNode->next->next->back != NULL)
	      //{
	      //traverse_tagnew (theNode->next->next->back, origin);
	      right = right->next;
	    }
        }
      if(theNode==origin)
        {
	  while(theNode != NULL && theNode->type != 'r')
            {
	      theNode = showtop(showtop(theNode)->back);
	      theNode->visited = TRUE;
            }
        }
    }
    else
    {
        error("a missing node pointer encountered, aborted");
    }
}

///
/// pick a node at random from the available list of nodes, this code will even 
/// work with datasets where there are only two nodes.
void
chooseOrigin (proposal_fmt * proposal)
{
    long elem = 0;
    // oldelem = (proposal->sumtips * 2.);
    node *tmp=NULL, **goal;
    //032110 goal = (node **) mycalloc (oldelem, sizeof (node *));
    //#ifdef DEBUG
    traverse_check(crawlback (proposal->root->next));
    //#endif
    //032110 traverseAllNodes (crawlback (proposal->root->next), &goal, &elem, &oldelem,
    //032110                  NO_MIGR_NODES);
    goal = proposal->world->nodep;
    elem = proposal->world->sumtips * 2;
    switch(elem)
    {
        case 0:
            error("problem with choosing a node for the tree change [choosOrigin()]");
            //break;
        case 1:
        case 2:
            tmp = goal[0];
            break;
        default:
            tmp = goal[RANDINT (0L, elem - 2)];
            while(tmp->back->type =='r')
                tmp = goal[RANDINT (0L, elem - 2)];
            break;
    }
    //032110 myfree(goal);
    proposal->origin = tmp;
    if (proposal->origin != showtop (crawlback (proposal->root->next)))
    {
        proposal->oback = showtop (crawlback (proposal->origin));
        proposal->osister = showsister (proposal->origin);
        if (proposal->oback != showtop (crawlback (proposal->root->next)))
        {
            proposal->ocousin = showsister (proposal->oback);
        }
        else
        {
            proposal->ocousin = NULL;
        }
    }
    if (proposal->origin == NULL)
        error ("Designation of origin for branch removal failed");
}

void add_proposal_to_tl(proposal_fmt *proposal, timelist_fmt ** timevector)
{
  size_t ii;
  long slot = (*timevector)->T;
  long lastfrom;
  size_t oldalloc;
  if ((*timevector)->allocT < ((*timevector)->allocT + proposal->migr_table_counter + proposal->migr_table_counter2))
    {
      oldalloc = (*timevector)->allocT;
      (*timevector)->allocT = (*timevector)->allocT + proposal->migr_table_counter + proposal->migr_table_counter2 + HUNDRED;
      (*timevector)->tl = (vtlist *) myrealloc((*timevector)->tl, sizeof(vtlist) * (*timevector)->allocT);
      allocate_lineages (timevector, 0, proposal->numpop);
    }
  else
      allocate_lineages (timevector, 0, proposal->numpop);
  
  vtlist * tls = &(*timevector)->tl[slot];
  tls->eventnode = proposal->origin;
  tls->age = proposal->origin->tyme;
  tls->from = proposal->origin->actualpop;
  tls->to = proposal->origin->pop;
  tls->type = proposal->origin->type;
  slot++;
  lastfrom = tls->from;
  for (ii=0;ii<proposal->migr_table_counter;ii++,slot++)
    {
      //printf("mig: %li -> %li %c %f\n", proposal->migr_table[ii].from, proposal->migr_table[ii].to, proposal->migr_table[ii].event, proposal->migr_table[ii].time);
      tls = &(*timevector)->tl[slot];
      tls->eventnode = NULL;
      tls->age = proposal->migr_table[ii].time;
      tls->from = proposal->migr_table[ii].from;
      tls->to = proposal->migr_table[ii].to;
      tls->type = proposal->migr_table[ii].event;
    }
  for (ii=0;ii<proposal->migr_table_counter2;ii++,slot++)
    {
      //printf("mig2: %li -> %li %c %f\n", proposal->migr_table2[ii].from, proposal->migr_table2[ii].to, proposal->migr_table2[ii].event, proposal->migr_table2[ii].time);
      tls = &(*timevector)->tl[slot];
      tls->eventnode = NULL;
      tls->age = proposal->migr_table2[ii].time;
      tls->from = proposal->migr_table2[ii].from;
      tls->to = proposal->migr_table2[ii].to;
      tls->type = proposal->migr_table2[ii].event;
    }
  if (proposal->migr_table_counter > 0)
    lastfrom = proposal->migr_table[proposal->migr_table_counter-1].from;
  if (proposal->migr_table_counter2 > 0)
    lastfrom = proposal->migr_table2[proposal->migr_table_counter2-1].from;
  tls = &(*timevector)->tl[slot]; //coalescence of orgin lineage into tree (or migrate2 table)
  tls->eventnode = NULL;
  tls->age = proposal->time;
  tls->from = lastfrom;
  tls->to = lastfrom;
  tls->type = 'i';
  slot++;
  (*timevector)->T = slot;
  if ((*timevector)->T > (*timevector)->allocT)
    error("timevector lineages failed to allocate mcmc1.c: 705");
}

///
/// fix a tymelist based on the lineage in proposal assuming proposal->origin is a tip
/// used for assignment probg_treetimes calculation
// this is the inverse of construct_localtimelist()
void fix_timelist(timelist_fmt *timevector, proposal_fmt *proposal)
{
  size_t ii;
  //printf("tips: %li, mig# %li, mig2# %li timelslice %li time: %f\n",proposal->sumtips,proposal->migr_table_counter, proposal->migr_table_counter2, proposal->timeslice,proposal->time);
  //printf("tyme: %f actualpop %li pop %li\n",proposal->origin->tyme,proposal->origin->actualpop, proposal->origin->pop);
  allocate_lineages (&timevector, 0, proposal->numpop);
  add_proposal_to_tl(proposal, &timevector);
  qsort ((void *) (*timevector).tl, (size_t) (*timevector).T, sizeof (vtlist), agecmp);
  add_partlineages(proposal->numpop, &timevector, proposal->world);
  for(ii=0;ii<timevector->T;ii++)
  {
    char type;
    long actualpop;
    long pop;
    if (timevector->tl[ii].eventnode != NULL)
      {
	type = timevector->tl[ii].eventnode->type;
	actualpop = timevector->tl[ii].eventnode->actualpop;
	pop = timevector->tl[ii].eventnode->pop;
      }
    else
      {
	type = timevector->tl[ii].type;
	actualpop = -1;
	pop = -1;
      }
    //printf("%i> ii=%li %0.5f  (%li->%li)[]%li, %li]  %c %3li\n",myID, ii, timevector->tl[ii].age,
    //	   timevector->tl[ii].from, timevector->tl[ii].to,actualpop,pop,
    //	   type,timevector->tl[ii].lineages[0]);
      //size_t jj;
      //for(jj = 1;jj < proposal->numpop; jj++)
      //printf(" %3li",timevector->tl[ii].lineages[jj]);
      //printf("\n");
      //printf(" backtyme:%f dif:%f\n",showtop(showtop(timevector->tl[ii].eventnode)->back)->tyme,dif=showtop(timevector->tl[ii].eventnode->back)->tyme-timevector->tl[ii].eventnode->tyme);
  }
  //error("check stop\n");
}

///
/// construct a list of ordered times with pointers to the tree
/// using the origin to prune the old time list
void
construct_localtimelist (timelist_fmt * timevector, proposal_fmt * proposal)
{
#ifdef TREEDEBUG1
    double dif;
    long ii;
#endif
    long z = 0;
    long oz = proposal->aboveorigin_allocsize;
    //    long numpop = timevector->numpop = proposal->numpop;
    traverseAllNodes (crawlback (proposal->origin)->back,
                      &proposal->aboveorigin, &z, &oz, WITH_MIGR_NODES,&proposal->divlist);
    proposal->aboveorigin_allocsize = oz;
    proposal->aboveorigin[z++] = proposal->oback;
    proposal->aboveorigin[z] = NULL;
    prune_timelist(&proposal->world->treetimes[0],timevector, proposal->aboveorigin, proposal, proposal->world->numpop);
#ifdef DEBUGXX
    printf("%i> in construct_local_timelist(): %li (%li)\n", myID, timevector->T, proposal->world->treetimes[0].T);
#endif
    add_partlineages(proposal->numpop, &timevector, proposal->world);
    //preparation for timeparam: add_skyline(world,&timevector)
#ifdef TREEDEBUG1
    for(ii=0;ii<timevector->T-2;ii++)
    {
        printf("%i> ii=%li %0.5f  (%li->%li) %c %3li",myID, ii, timevector->tl[ii].age, 
               timevector->tl[ii].from, timevector->tl[ii].to, 
               timevector->tl[ii].eventnode->type,timevector->tl[ii].lineages[0]);
        long jj;
        for(jj = 1;jj < proposal->numpop; jj++)
            printf(" %3li",timevector->tl[ii].lineages[jj]);
        printf(" backtyme:%f dif:%f\n",showtop(showtop(timevector->tl[ii].eventnode)->back)->tyme,dif=showtop(timevector->tl[ii].eventnode->back)->tyme-timevector->tl[ii].eventnode->tyme);
    }
#endif
}



void record_above_origin_divergence(div_fmt **divlist, long from, long to)
{ 
  if ( (*divlist)->div_elem >= (*divlist)->div_allocsize)
    {
      (*divlist)->div_allocsize += HUNDRED;
      (*divlist)->divlist = (longpair *) realloc((*divlist)->divlist, sizeof(longpair) * (size_t) (*divlist)->div_allocsize);
    }
  (*divlist)->divlist[(*divlist)->div_elem][0]=from;
  (*divlist)->divlist[(*divlist)->div_elem][1]=to;
  (*divlist)->div_elem += 1;
}

void allocate_nodelist(node ***nodelist, long *oldnode_elem, long elements)
{
  long olde = (*oldnode_elem);
  (*oldnode_elem) += elements;
  (*nodelist) = (node **) realloc(*nodelist, sizeof(node *) * (size_t) (*oldnode_elem));
  memset(*nodelist + olde, 0,sizeof(node*) * (size_t) elements);
}

void test_traverse_error(node * theNode)
{
  if(theNode->type != 'r' && theNode->tyme > showtop(theNode->back)->tyme)
    {
      printf("%i> Tree traverse failed\nwith node (%c) time %10.10f and node (%c) time %10.10f (diff=%20.20f)\n",
	     myID, theNode->type,  theNode->tyme, showtop(theNode->back)->type, showtop(theNode->back)->tyme,theNode->tyme - showtop(theNode->back)->tyme);
      error("traverseAllNodes() failed");
    }
}

/*----------------------------------------------------------------------------
 finds all nodes in a tree starting at the root node and crawling up 
 to the tips in a recursive fashion, writing nodeptrs in the nodelist vector
 the flag include_migration is 1 if we want to touch the migration nodes too,
 otherwise =0 -> jump over the migration nodes. for convenience we define the 
 the macros NO_MIGR_NODES=0 and WITH_MIGR_NODES=1 in the definitions.h file
 */
void
traverseAllNodes (node * theNode, node *** nodelist, long *node_elem,
                  long *oldnode_elem, int include_migration, div_fmt **divlist)
{
  node *right = NULL;
  test_traverse_error(theNode);
  if (include_migration == NO_MIGR_NODES)
    {
      if (theNode->type != 't')
        {
	  right = theNode->next;
	  while(!right->top)
	    {
	      node *testnode = crawlback(right);
	      if (testnode != NULL)
		traverseAllNodes (testnode, nodelist, node_elem, oldnode_elem, NO_MIGR_NODES, divlist);
	      right = right->next;	      
	      //node *testnode2 = crawlback (theNode->next->next);
	      //if (testnode2 != NULL && theNode->type != 'm' && theNode->type != 'd')
	      //traverseAllNodes (testnode2, nodelist, node_elem, oldnode_elem, NO_MIGR_NODES, divlist);
	    }
	  if (*node_elem >= *oldnode_elem-3)
	    allocate_nodelist(nodelist,oldnode_elem, HUNDRED);
	  (*nodelist)[(*node_elem)] = theNode;
	  (*node_elem) += 1;
        }
      else
        {
	  if (*node_elem >= *oldnode_elem-3)
	    allocate_nodelist(nodelist,oldnode_elem, HUNDRED);
	  (*nodelist)[(*node_elem)] = theNode;
	  (*node_elem) += 1;
        }
    }
  else
    {
      if (theNode->type != 't')
        {
	  right = theNode->next;
	  while(!right->top)
	    {
	      if (right->back != NULL)
		traverseAllNodes (right->back, nodelist, node_elem,
				  oldnode_elem, WITH_MIGR_NODES, divlist);
	      //if (divlist != NULL && theNode->type == 'd')
	      //record_above_origin_divergence(divlist, theNode->pop,theNode->actualpop);
	      //if (theNode->type == 'i' && theNode->next->next->back != NULL)
	      //traverseAllNodes (theNode->next->next->back, nodelist, node_elem,
	      //	      oldnode_elem, WITH_MIGR_NODES, divlist);
	      right = right->next;
	    }
	  if (*node_elem >= *oldnode_elem-3)
	    allocate_nodelist(nodelist,oldnode_elem, HUNDRED);
	  (*nodelist)[(*node_elem)++] = theNode;
        }
      else
        {
	  if (*node_elem >= *oldnode_elem-3)
	    allocate_nodelist(nodelist,oldnode_elem, HUNDRED);
	  (*nodelist)[(*node_elem)++] = theNode;
        }
    }
}

///
/// copies elements from the old timelist (oldtv)
/// into the new local timeslist (newtv) checking whether they
/// are used after the removal of the residual tree that is above the origin (ptr)
void
prune_timelist (timelist_fmt * oldtv, timelist_fmt * newtv, node ** ptr, proposal_fmt * proposal, long numpop)
{
  (void) numpop;
  (void) proposal;
    register long i    = 0;
    register long j    = 0;
    register node * thenode;
    long          slot = 0;
    vtlist        * tls;
    node *thisptr;
    slot=0;
    // go through all slices in old time list and
    // extract node pointers
    for (i = 0; i < (*oldtv).T; i++)
    {
        j = 0;
        thenode = (*oldtv).tl[i].eventnode;
	thisptr = ptr[0];
        while (thenode != thisptr) 
        {
	  if (thisptr == NULL)
	    {
	      // if the comparison reaches the end of the ptr list that holds an NULL element at the end
	      // then the node needs to be present in the new time list.
	      tls = &((*newtv).tl[slot]);
	      tls->eventnode = thenode;
	      tls->age       = (*oldtv).tl[i].age;
	      tls->interval  = (*oldtv).tl[i].interval;
	      if(thenode!=NULL)
		{

		  tls->from      = thenode->pop;
		  tls->to        = thenode->actualpop;
		}
	      tls->slice     = slot;
	      slot++;	    
	      break;
	    }
	  j++;
	  thisptr = ptr[j];
	}        
    }
    (*newtv).T = slot;
    qsort ((void *) (*newtv).tl, (size_t) (*newtv).T, sizeof (vtlist), agecmp);
}

/* replaces nodepointers in list 1 with NULL if they are present in list 2
 returns the first NULL slot in the array.
 */
long
xor (node ** ptrl1, node ** ptrl2)
{
    long i = 0, j = 0, slot = -1;
    /* assumes that there is an NULL element at the end */
    for (i = 0; ptrl1[i] != NULL; j = 0, i++)
    {
        while ((ptrl1[i] != ptrl2[j]) && (ptrl2[j] != NULL))
            j++;
        if (ptrl2[j] != NULL)
        {
            if (slot == -1)
                slot = i;
            ptrl1[i] = NULL;
        }
    }
    return slot;
}

/* migrate() fills the PROPOSAL->MIGRATION_TABLE in the tree or
 PROPOSAL->MIGRATION_TABLE2 when at the bottom of the tree
 
 PROPOSAL proposal-scratchpad
 UP       node above (younger)
 MIGR_TABLE_COUNTER migration array counter, 
 will increase by one during execution
 AIR      if true standard execution, if false updating the last 
 lineage in the residual tree.

 2013: addition of a second type of 'migration' which is equivalent to
 speciation array records now also the event
 */
int
migrate (proposal_fmt * proposal, node * up, char event, long from, long to)
{
    long tmp;
    long numpop = proposal->numpop;
    long i = proposal->migr_table_counter;
    migr_table_fmt *array = proposal->migr_table;
    if (i > MIGRATION_LIMIT * numpop || numpop < 2)
    {
        return 0;
    }
    if (i > 0)
        array[i].to = array[i - 1].from;
    else
        array[i].to = up->pop;

    if (array[i].to != to)
      {
	error("Aborted when adding a migration or divergence event");
      }
    if(event=='d')
      {
      	species_fmt * s = get_fixed_species_model(from, to, proposal->world->species_model, proposal->world->species_model_size);
	if (s!= NULL)
	  tmp = s->from;
	else
	  return 0;
	if (tmp == -1)
	  return 0;
      }
    else
      {
	tmp = from;
	event = 'm';
	if (tmp == -1)
	  return 0;
      }

    array[i].event = event;
    //#if DEBUG
    //if (i>0)
    //  {
    //	if (array[i-1].event == 'd' && event == 'm')
    //	  error("GOPF");
    //  }
    //#endif
    if(tmp < numpop)
        array[i].from = tmp;
    else
        return 0;
    
    //DEBUG
    // printf("i: %3li %li %li\n",i,proposal->migr_table[i].from,proposal->migr_table[i].to);
    array[i++].time = proposal->time;
    if(proposal->time < proposal->origin->tyme)
    {
        error("in migrate() wrong time found");
    }
    if (i >= proposal->old_migr_table_counter)
    {
        proposal->old_migr_table_counter += 10;
        proposal->migr_table =
        (migr_table_fmt *) myrealloc (proposal->migr_table,
                                      sizeof (migr_table_fmt) *
                                      (size_t) (proposal->old_migr_table_counter));
    }
    proposal->migr_table_counter = i;
    return 1;
}

// these function changed compared to 3.7 because we evaluate eventtime_single for each parameter
// and know to and from values
int
migrateb (proposal_fmt * proposal, node * up, char event, long from, long to)
{
  const long numpop = proposal->numpop;
  long tmp;
    long i = proposal->migr_table_counter2;
    migr_table_fmt *array = proposal->migr_table2;
    if (i > MIGRATION_LIMIT * numpop)
    {
        return 0;
    }
    if (i > 0)
        array[i].to = array[i - 1].from;
    else
        array[i].to = up->pop;
    if (array[i].to != to)
      {
	error("Aborted when adding a migration or divergence event");
      }
    if(event=='d')
      {
      	species_fmt * s = get_fixed_species_model(from, to, proposal->world->species_model, proposal->world->species_model_size);
	if (s!= NULL)
	  tmp = s->from;
	else
	  return 0;
	if (tmp == -1)
	  return 0;
      }
    else
      {
	tmp = from;
	event = 'm';
	if (tmp == -1)
	  return 0;
      }
    array[i].event = event;
    if(tmp < numpop)
        array[i].from = tmp;
    else
        return 0;
    
    array[i++].time = proposal->time;

    if (i >= proposal->old_migr_table_counter2)
    {
        proposal->old_migr_table_counter2 += 10;
        proposal->migr_table2 =
        (migr_table_fmt *) myrealloc (proposal->migr_table2,
                                      sizeof (migr_table_fmt) *
                                      (size_t) (proposal->old_migr_table_counter2));
    }
    proposal->migr_table_counter2 = i;
    return 1;
}
/* OLD CODE
int
migrateb (proposal_fmt * proposal, node * up, char event)
{
  const long numpop = proposal->numpop;
  long tmp;
    long i = proposal->migr_table_counter2;
    migr_table_fmt *array = proposal->migr_table2;
    if (i > MIGRATION_LIMIT * numpop)
    {
        return 0;
    }
    if (i > 0)
        array[i].to = array[i - 1].from;
    else
        array[i].to = up->pop;

    if(event=='d')
      {
      	tmp = speciation_from (array[i].to, proposal);
	if (tmp == -1)
	  return 0;
      }
    else
      {
	tmp = migration_from (array[i].to, proposal);
	event = 'm';
	if (tmp == -1)
	  return 0;
      }
    array[i].event = event;
    if(tmp < numpop)
        array[i].from = tmp;
    else
        return 0;
    
    array[i++].time = proposal->time;

    if (i >= proposal->old_migr_table_counter2)
    {
        proposal->old_migr_table_counter2 += 10;
        proposal->migr_table2 =
        (migr_table_fmt *) myrealloc (proposal->migr_table2,
                                      sizeof (migr_table_fmt) *
                                      (size_t) (proposal->old_migr_table_counter2));
    }
    proposal->migr_table_counter2 = i;
    return 1;
}
*/
/* migration_from() returns the FROM population when there was a migration
 TO        population to migrate to (or species offspring (from:ancestor --> to:ofsspringspecies)
 PROPOSAL  proposal-scratchpad
 */
long
migration_from (long to, proposal_fmt * proposal)
{
    long ii = 0;
    MYREAL *r = proposal->world->migproblist[to];
    MYREAL rr = UNIF_RANDUM ();
    while (rr >  r[ii] && ii < proposal->world->numpop-1)
    {
        ii++;
    }
    if (ii < to)
        return ii;
    else
        return ++ii;
}

void
chooseTarget (proposal_fmt * proposal, timelist_fmt * timevector,
              node ** bordernodes, long *bordernum)
{
    long actualpop = -99;
    node *rb = crawlback (proposal->root->next);
    *bordernum = 0;
    proposal->target = NULL;
    proposal->realtarget = NULL;
    if (proposal->migr_table_counter == 0)
        actualpop = proposal->origin->pop;
    else
        actualpop = proposal->migr_table[proposal->migr_table_counter - 1].from;
    if (rb->tyme < proposal->time)
    {
        error ("Wrong Time for action in chooseTarget()\n");
    }
    //printf("\n\n\n");
    findbordernodes (rb, proposal, actualpop, &bordernodes, 
		     &proposal->bordernodes_allocsize, bordernum,
                     &(*timevector).tl, (*timevector).T);
    if (*bordernum > 0)
    {
        // found elegible lineages to coalesce to
        proposal->target = bordernodes[RANDINT (0L, (*bordernum) - 1)];
        if (proposal->target != rb)
        {
            proposal->tsister = showsister (proposal->target);
            proposal->realtsister = crawlback (proposal->tsister)->back;
        }
        else
            proposal->tsister = NULL;
        proposal->realtarget = proposal->target;
        if (proposal->target->type == 'm' || proposal->target->type == 'd')
            proposal->target = crawlback (showtop (proposal->target)->next);
    }
    else
    {
        // no lineage to coalesce to was found.
        proposal->target = NULL;
        proposal->tsister = NULL;
        proposal->realtsister = NULL;
        proposal->realtarget = NULL;
    }
}

void findbordernodes_recursive (node * theNode, proposal_fmt * proposal, long pop,
                 node *** bordernodes, long *allocsize, long *bordernum, vtlist ** tyme,
                 long gte)
{
    node *tmp, *back;
    // we search on the old tree and reaching the oback node that was excised from
    // from the residual tree we jump over it and go up the sister branch and back is
    // going further down than the oback node
    if (theNode == proposal->oback)
    {
        tmp = showtop (crawlback (proposal->osister)->back);
        back = showtop (proposal->oback->back);
    }
    else
    {
        tmp = showtop (theNode);
        back = showtop (theNode->back);
    }
    //printf("%i> fb: %p %f <? %f <? %f (%c%c)[back%li tmp%li]",myID,proposal, back->tyme,proposal->time,tmp->tyme,back->type, tmp->type, back->actualpop,tmp->pop);  
    if (pop == tmp->pop && pop == back->actualpop && tmp->tyme < proposal->time
        && back->tyme > proposal->time)
    {
        // lineage end points bracket the the new proposed time,
        // therefore this branch needs to be included
      if (*bordernum >= *allocsize-3)
	{
	  allocate_nodelist(bordernodes, allocsize, HUNDRED);      
	}
      (*bordernodes)[*bordernum] = tmp;
      *bordernum += 1;
      //printf("****\n");
      return;
    }
    else
    {
      //printf("\n");
      // the lineage does not fit the populations (it may fit the time, so)
        // the if here checks this
        if (back->tyme < proposal->time)
        {
            // this makes sure that we looked at all lineages that
            // may be electable to receive a coalescent have been looked
            // at, if we reach here the time of the proposal is older than
            // all available lineages on this branch and we safely can stop
            // following this branch.
            return;
        }
        // if we are still working on interior nodes and the proposal time is still 
        // younger than then lineage bracketing nodes we need to continue our search
        // up in the tree following right and left branches, but need to stop when 
        // we encounter a tip node (this is important with dated tips because these can
        // be interspersed in the timelist.
        if (tmp->type != 't')
	  {
	    if (tmp->next->back != NULL)
	      {
		//node *tmpup = crawlback(tmp->next);
		//if (tmpup->tyme > proposal->time && back->tyme > proposal->time)
		//  {
		//  if(tmpup->type == 'i')
		//    {
		//findbordernodes (tmpup->next->back, proposal, pop, bordernodes, allocsize,
		//		 bordernum, tyme, gte);
		//findbordernodes (tmpup->next->next->back, proposal, pop,
		//		 bordernodes, allocsize, bordernum, tyme, gte);
		//    }
		//}
		//else
		//{
		    findbordernodes (tmp->next->back, proposal, pop, bordernodes, allocsize,
				     bordernum, tyme, gte);
		    //}
	      }
	    warning("needs fix for multi-branches\n");
	    if (tmp->type == 'i' && tmp->next->next->back != NULL)
	      findbordernodes (tmp->next->next->back, proposal, pop,
			       bordernodes, allocsize, bordernum, tyme, gte);
	  }
#ifdef DEBUGXXXX
	else
	  {
	    if (pop == tmp->pop && pop == back->actualpop && tmp->tyme < proposal->time
		&& back->tyme > proposal->time)
	      {
		// lineage end points bracket the the new proposed time,
		// therefore this branch needs to be included
		if (*bordernum >= *allocsize-3)
		  {
		    allocate_nodelist(bordernodes, allocsize, HUNDRED);      
		  }
		(*bordernodes)[*bordernum] = tmp;
		*bordernum += 1;
		//printf("***t\n");
		return;
	      }
	  }
#endif
    }
}

void add_node_border(node *tmp, node ***bordernodes,long *bordernum, long *allocsize)
{
  if (*bordernum >= *allocsize-3)
    {
      allocate_nodelist(bordernodes, allocsize, HUNDRED);      
    }
  (*bordernodes)[*bordernum] = tmp;
  *bordernum += 1;
}

boolean is_in_bracket(node *tmp, node*back, MYREAL thetime, long pop)
{
  if (pop == tmp->pop && pop == back->actualpop && tmp->tyme < thetime
      && back->tyme > thetime)
    return TRUE;
  else
    return FALSE;
}

void increase_stack(node ***stack, long *allocsize,  long newallocsize)
{
  if (newallocsize < *allocsize)
    return;
  *stack = (node **) myrealloc(*stack, sizeof(node **) * (size_t) newallocsize);
  memset(*stack + *allocsize, 0, sizeof(node **) * (size_t) (newallocsize - *allocsize));
  *allocsize = newallocsize;
}

//void
//findbordernodes_iterative (node * theNode, proposal_fmt * proposal, long pop,
//                 node *** bordernodes, long *allocsize, long *bordernum, vtlist ** tyme,
//                 long gte)
void
findbordernodes (node * theNode, proposal_fmt * proposal, long pop,
                 node *** bordernodes, long *allocsize, long *bordernum, vtlist ** tyme,
                 long gte)
{
  (void) gte;
  (void) tyme;
  //boolean backthreading = FALSE;
  //boolean done = FALSE;

  // non recursive version
  //http://codereview.stackexchange.com/questions/47932/recursion-vs-iteration-of-tree-structure


  node * tmp = theNode;
  node *back;
  if (theNode == NULL)
    return;
  
  // Create an empty stack and push root to it
  node ** stack = (node **) calloc(HUNDRED,sizeof(node*));
  long stackallocsize = HUNDRED;
  long stacksize = 0;
  stack[0] = tmp;
  
  /* Pop all items one by one. Do following for every popped item
     a) print it
     b) push its right child
     c) push its left child
     Note that right child is pushed first so that left is processed first */
  while (stacksize >= 0)
      {
        // Pop the top item from stack and print it
        tmp = stack[stacksize];
	stacksize--;
	if (tmp==NULL)
	  continue;
	if (tmp == proposal->oback)
          {
            //tmp = showtop (crawlback (proposal->osister)->back);
            tmp = showtop (proposal->osister);
	    //            back = showtop (proposal->oback->back);
            back = showtop (crawlback(proposal->oback));
          }
        else
          {
            tmp = showtop (tmp);
	    //            back = showtop (tmp->back);
            back = showtop (crawlback(tmp));
          }
	//if (tmp->type == 'm')
	//  error("gopf");

        //printf ("@@%c id=%li %f-%f %li [%f]", tmp->type, tmp->id, tmp->tyme, back->tyme, stacksize+1, proposal->time);
	node *newtmp = tmp;
	node *oldtmp = tmp;
	while (newtmp->tyme < proposal->time && newtmp != back)
	  {
	    oldtmp = newtmp;
	    newtmp = showtop(showtop(newtmp)->back);
	  }
	boolean yes = is_in_bracket(oldtmp,newtmp,proposal->time,pop);
	//printf("%c\n", yes ? '*' : ' ');
	if (yes)
	  add_node_border(oldtmp, bordernodes,bordernum, allocsize);
        if (tmp->next && !tmp->tip)
	  {
	    stacksize++;
	    increase_stack(&stack, &stackallocsize,  stacksize);
	    stack[stacksize] = showtop(crawlback(tmp->next));
	  }
        if (tmp->next!=NULL && tmp->next->next && !tmp->tip)
	  {
	    stacksize++;
	    increase_stack(&stack, &stackallocsize,  stacksize);
	    stack[stacksize] = showtop(crawlback(tmp->next->next));
	  }
      }
  myfree(stack);
}





#ifdef DEBUGXXXTEST

    node *tmp, *back;
    // we search on the old tree and reaching the oback node that was excised from
    // from the residual tree we jump over it and go up the sister branch and back is
    // going further down than the oback node
    boolean done = FALSE;
    node * tmp = theNode;
    while(!done)
      {	
	if (tmp == proposal->oback)
	  {
	    tmp = showtop (crawlback (proposal->osister)->back);
	    back = showtop (proposal->oback->back);
	  }
	else
	  {
	    tmp = showtop (tmp);
	    back = showtop (tmp->back);
	  }
	if (is_in_bracket(tmp,back,proposal->time,pop))
	  add_node_border(bordernodes,bordernum,allocsize);
	if (tmp->tyme < proposal->time)
	  {
	    tmp = showtop(crawlback(tmp)->next->next->back);
	    back = crawlback(tmp);
	    continue;
	  }
	else
	  {
	    
	  }

	if (pop == tmp->pop && pop == back->actualpop && tmp->tyme < proposal->time
	    && back->tyme > proposal->time)
	  {
	    // lineage end points bracket the the new proposed time,
	    // therefore this branch needs to be included
	    if (*bordernum >= *allocsize-3)
	      {
		allocate_nodelist(bordernodes, allocsize, HUNDRED);      
	      }
	    (*bordernodes)[*bordernum] = tmp;
	    *bordernum += 1;
	    tmp = crawlback(tmp)
	  }
	else
	  {
	    // the lineage does not fit the populations (it may fit the time, so)
	    // the if here checks this
	    if (back->tyme < proposal->time)
	      {
		// this makes sure that we looked at all lineages that
		// may be electable to receive a coalescent have been looked
		// at, if we reach here the time of the proposal is older than
		// all available lineages on this branch and we safely can stop
		// following this branch.
		return;
	      }
	    // if we are still working on interior nodes and the proposal time is still 
	    // younger than then lineage bracketing nodes we need to continue our search
	    // up in the tree following right and left branches, but need to stop when 
	    // we encounter a tip node (this is important with dated tips because these can
	    // be interspersed in the timelist.
	    if (tmp->type != 't')
	      {
		back = tmp;
		tmp = showtop(tmp->next->back);
		continue;
		tmp = tmp->next->next->back;
	      }
	  }
      }
}
#endif

#ifdef TESTING2
///
/// freeing the proposal structure, this is the replacement of the real free_proposal method
/// and does not allocate, but reuses old memory
=======
//void free_proposal(proposal_fmt *proposal)/
//{
    // do nothing
//}
void reset_simple_proposal_variables(proposal_fmt **proposal)
{
    (*proposal)->mig_removed = FALSE;
    (*proposal)->rr = 0.0;
    (*proposal)->origin = NULL;
    (*proposal)->target = NULL;
    (*proposal)->realtarget = NULL;
    (*proposal)->tsister = NULL;
    (*proposal)->realtsister = NULL;
    (*proposal)->osister = NULL;
    (*proposal)->realosister = NULL;
    (*proposal)->ocousin = NULL;
    (*proposal)->realocousin = NULL;
    (*proposal)->oback = NULL;
    (*proposal)->realoback = NULL;
    //
    (*proposal)->connect = NULL;
    (*proposal)->likelihood = 0.0;
    (*proposal)->time = 0.0;
    (*proposal)->v = 0.0;
    (*proposal)->vs = 0.0;
    //
#ifdef UEP
    (*proposal)->ueplikelihood = 0.0;
#endif
    (*proposal)->migr_table_counter=0;
    (*proposal)->migr_table_counter2=0;
    (*proposal)->timeslice = 0;
    //
    (*proposal)->treelen = 0.0;
#ifdef BEAGLE
    (*proposal)->parentid = 0;
    (*proposal)->leftid = 0;
    (*proposal)->rightid = 0;
#endif
}

void
reset_proposal (proposal_fmt ** proposal, world_fmt *world)
{
  mutationmodel_fmt *s;
  //world_fmt *world = proposal->world;
  const long listsize = 2 * (world->sumtips + 2);
  const long newsize = 4 * listsize;
  const long mal = (*proposal)->world->data->maxalleles[(*proposal)->world->locus]+1;
  (*proposal)->likelihood = -HUGE;
  // pointers and values to outside structures
  (*proposal)->world = world;
  (*proposal)->datatype = world->options->datatype;
  (*proposal)->sumtips = world->sumtips;
  (*proposal)->numpop = world->numpop;
  (*proposal)->endsite = world->data->seq[0]->endsite;
  (*proposal)->fracchange = world->data->seq[0]->fracchange;
  (*proposal)->param0 = world->param0;
  (*proposal)->root = world->root;
  (*proposal)->migration_model = world->options->migration_model;
  // precalculated values
  (*proposal)->mig0list = world->mig0list;
  (*proposal)->design0list = world->design0list;
  (*proposal)->listsize =  listsize;
  // ..line_t are also reset with this

//  memset(proposal->nodedata,0,sizeof(node *) * (2 * listsize + 2 * sumtips)); 
  long locus = world->locus;
  const long sublocistart = world->sublocistarts[locus];
  const long sublociend   = world->sublocistarts[locus+1];
  long sublocus;
  long endsite;
  long rcategs;
  for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
    {
      long i;
      for(i=0;i<newsize;i++)
        (*proposal)->nodedata[i]=NULL;
      //  memset((*proposal)->nodedata,0,sizeof(node *) * newsize); 
      memset((*proposal)->mf,0, sizeof(MYREAL)*(2 * (*proposal)->endsite)); // (*proposal)->mt is also freed with this
      // resetting pointers
      reset_simple_proposal_variables(proposal);
      
      if (strchr (SEQUENCETYPES, (*proposal)->datatype))
	{
	  s = &world->mutationmodels[sublocus];
	  endsite = s->numpatterns;
	  rcategs = s->numsiterates;
	  const long xs = sublocus - sublocistart;
	  memset((*proposal)->mf[xs],0, sizeof(MYREAL)*(endsite)); 
	  memset((*proposal)->mt[xs],0, sizeof(MYREAL)*(endsite)); 
	  zero_xseq(&(*proposal)->xf[xs], endsite, rcategs);
	  zero_xseq(&(*proposal)->xt[xs], endsite, rcategs);
	}
    }
  memset((*proposal)->migr_table, 0, sizeof(migr_table_fmt) * (*proposal)->old_migr_table_counter);
  memset((*proposal)->migr_table2, 0, sizeof(migr_table_fmt) * (*proposal)->old_migr_table_counter2);
#ifdef UEP
    
    if ((*proposal)->world->options->uep)
    {
        memset((*proposal)->ueplike, 0, sizeof(MYREAL) * world->numpop * world->data->uepsites);
        for (j = 1; j < world->data->uepsites; ++j)
            (*(*proposal))->ueplike[j] = (*(*proposal))->ueplike[0] + j * world->numpop;
        
        memset((*proposal)->uf.s, 0, world->data->uepsites * sizeof(pair));
        memset((*proposal)->ut.s, 0, world->data->uepsites * sizeof(pair));
        memset((*proposal)->umf, 0, world->data->uepsites * sizeof(MYREAL));
        memset((*proposal)->umt, 0, world->data->uepsites * sizeof(MYREAL));
    }
#endif
}
#endif
///
/// freeing the proposal structure, this will be replaced by a function that resets permanently
/// allocated structure [reset_proposal()]
#ifdef TESTING2
void
free_masterproposal (proposal_fmt * proposal)
{
  //implement  free_masterproposal(proposal);
}
#else
void
free_masterproposal (proposal_fmt * proposal)
{
  mutationmodel_fmt *s;
  long sublocus;
  long locus = proposal->world->locus;
  long sublocistart = proposal->world->sublocistarts[locus];
  long sublociend = proposal->world->sublocistarts[locus+1];
  
  //static long count=0;
  myfree(proposal->aboveorigin);
  myfree(proposal->bordernodes);
  myfree(proposal->line_f);
  myfree(proposal->line_t);
  myfree(proposal->divlist->divlist);
  myfree(proposal->divlist);
  myfree(proposal->param0save);
  //  printf("%li ",count++);fflush(stdout);
  //myfree(proposal->nodedata); // ..line_t are also freed with this
    for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
      {
	s = &proposal->world->mutationmodels[sublocus];
	const long xs = sublocus - sublocistart;
	myfree(proposal->mf[xs]); 
	myfree(proposal->mt[xs]); 
	if (strchr (SEQUENCETYPES, s->datatype))
	  {
	    myfree(proposal->xf[xs].s[0]);
	    myfree(proposal->xt[xs].s[0]);
	    myfree(proposal->xf[xs].s);
	    myfree(proposal->xt[xs].s);
	  }
	else
	  {
	    myfree(proposal->xf[xs].a);
	    myfree(proposal->xt[xs].a);
	  }
      }
    myfree(proposal->mf); 
    myfree(proposal->mt); 

    myfree(proposal->xt);
    myfree(proposal->xf);
    myfree(proposal->migr_table);
    myfree(proposal->migr_table2);
#ifdef UEP
    
    if (proposal->world->options->uep)
    {
        myfree(proposal->ueplike[0]);
        myfree(proposal->ueplike);
        myfree(proposal->uf.s);
        myfree(proposal->ut.s);
        myfree(proposal->umf);
        myfree(proposal->umt);
        
    }
#endif
    myfree(proposal);
}

#endif  /*TESTING2*/


void
free_timevector (timelist_fmt * timevector)
{
    myfree(timevector->lineages);
    myfree(timevector->tl);
    myfree(timevector);
    //printf("%i> free timevector\n",myID);
}

///
/// using a global timevector to save malloc/free pair calls and 
/// so to save time
void free_timevector_new (timelist_fmt * timevector)
{
  (void) timevector;
    //    memset(timevector->lineages, 0, sizeof(long) * );
    // (timevector->tl);
}

/*----------------------------------------------------------*
 * rejection/acceptance of the new tree according to the likelihood
 * and an acceptance ratio which is higher the better the
 * likelihood values are (-> Metropolis)
 */
boolean
acceptlike (world_fmt * world, proposal_fmt * proposal, long g,
            timelist_fmt * tymelist, boolean assign)
{
    //proposal changes when chain is stuck, to accept some new trees without
    //considering data to remove out of sticky area [hopefully]
    //static long not_accepted = 0;
    
  //with assign we use the tymelist and the proposal to recalculate the tymelist for the new tree probg.
  // we will change tymelist as a side-effect asuming that we do not need to use the tymelist further 
    const long limit = MIGRATION_LIMIT * world->numpop;
    
    long rm  = 0L;
    long rmc = rmigrcount (proposal);
    
    MYREAL rr;
    MYREAL expo;
    
#ifdef UEP
    node *first;
#endif
    
    rm =  proposal->migr_table_counter + proposal->migr_table_counter2 
    + world->migration_counts - rmc;
    
    if (rm > limit )//|| proposal->divlist->div_elem > 0)
    {
      //should we use this at all     if (world->cold)
      //warning ("migration limit (%li) exceeded: %li\n",
      //	 MIGRATION_LIMIT * world->numpop, rm);
#ifdef BEAGLE
      //printf("Reset to: %f\n",world->likelihood[g]);
      //set_beagle_dirty(proposal->origin,proposal->target,showtop(world->root->next->back));
      //calcLnL(world, world->beagle->instance);
#endif
      return FALSE;
    }
#ifdef UEP
    if (world->options->uep)
    {
        first = first_uep2 (proposal->world->root->next->back,
                            proposal->world->root->next->back,
                            proposal->world->data->uepsites);
        proposal->firstuep = first_uep (first, proposal->world->root,
                                        proposal->world->data->uepsites);
        proposal->ueplikelihood = pseudo_ueplikelihood (world, proposal);
        proposal->likelihood = pseudotreelikelihood (world, proposal);
    }
    else
    {
        proposal->likelihood = pseudotreelikelihood (world, proposal);
    }
#else
    //    printf("DEBUG: world->likellihood[%li]=%f heat=%f\n",g,world->likelihood[g],world->heat);
    proposal->likelihood = pseudotreelikelihood (world, proposal);
    //    if (world->cold)
    //  printf("DEBUG [%s]: old=%f proposal=%f heat=%f\n",world->in_burnin ? "burnin" : "sample",world->likelihood[g], proposal->likelihood, world->heat);
#endif
    
    //if(proposal->likelihood <= -HUGE &&  world->likelihood[g] <= -HUGE)
    //  {
	//	if(!world->in_burnin)
	//  {
    //long burnin = world->options->burn_in;
    //world->options->burn_in *= 10;
    //fprintf(stdout,"likelihood is -inf, insert a burn-in phase\n");
    //	    burnin_chain(world);
    //world->options->burn_in = burnin;
    //  }
    //	return TRUE;
    // }
    //@@@@@@@@@@@@@@@@@@@@@ <<<<<<
#ifdef DEBUGX
    if(world->cold)
      {
    	printf("DEBUG %i> like=%f proplike=%f diff=%f\n",myID, world->likelihood[g], proposal->likelihood, -world->likelihood[g]+ proposal->likelihood);
      }
#endif
    double oldlike = world->likelihood[g];
    double oldprobg = 0.0;
    double newprobg = 0.0;
    double newp = 0.0;
    double oldp = 0.0;
    if (assign)
      {
	fix_timelist(tymelist, proposal);
	//calc probg(theta|g,an)
	oldprobg =  probg_treetimes(world);
	newprobg =  probg_treetimes_local(world, tymelist);
	
      }
    newp = proposal->likelihood + newprobg;
    oldp = oldlike + oldprobg;
    if (oldp <= newp) //ratio is > 1
      {
        return TRUE;
    }
    rr = LOG (RANDUM ());
    if (!world->options->heating)
    {
	if (rr < (newp - oldp)) // r < ratio
        {
            return TRUE;
        }
    }
    else
    {
	if (rr < ((newp - oldp) * world->heat))
        {
            return TRUE;
        }
    }
    if(world->options->prioralone)
    {
        return TRUE;
    }
    //#ifdef DEBUG
    //if (world->heat < 0.1)
    //  {
    //warning("weird: %f ?< (%f-%f)*%f\n",rr, proposal->likelihood, world->likelihood[g], 
    //		world->heat);
    // }
    //#endif
    return FALSE;
}

long
rmigrcount (proposal_fmt * proposal)
{
    node *p;
    long count = 0;
    for (p = proposal->origin; p != proposal->oback; p = showtop (p->back))
    {
        if (p->type == 'm')
            count++;
    }
    return count;
}


///
/// generates the time interval to the next event (migration or coalescence)
/// and also decides what event happens.
/// This funcion is modified by the rate of the mutation rate that can be
/// estimated in a Bayesian context, and be fixed in a ML context
/// this rate is only influencing the time and not the eventtype
/// because the rate cancels out of the ratio that chooses the event.
MYREAL
eventtime (proposal_fmt * proposal, long pop, vtlist * tentry, char *event)
{
    //    static boolean mig_force = TRUE;
    MYREAL  interval;
    MYREAL  lines;
    MYREAL  denom;
    MYREAL  invdenom;
    MYREAL  inheritance = proposal->world->options->inheritance_scalars[proposal->world->locus];
    MYREAL  rate = proposal->world->options->mu_rates[proposal->world->locus];
    MYREAL  invrate = 1./ rate;
    long    timeslice = tentry->timeslice;
    long    addition=0;
    long    timepop = (proposal->world->numpop2+addition)*timeslice + pop;
    MYREAL  mm         = proposal->mig0list[pop] * invrate;
    MYREAL  timethetarate = proposal->world->timek[timepop] * proposal->param0[pop];
    MYREAL  *skyparam = proposal->world->timek+(proposal->world->numpop2+addition)*timeslice;
    lines   = 2.0 * tentry->lineages[pop];
    // param0 could be partitioned over time so that we we have param0[timeslice][pop]
    // makeing proposal for particular timeslices instead of all, one could have equal weights
    // for all timeslices or then weigh them using lineages? or an arbitrary prior.
    // we also need a vector for times to establish the size of the timeslices
    if(!proposal->world->options->skyline)
      denom    = mm + (lines * (1.0 / (inheritance*rate*timethetarate)));
    else
      {
	long msta = proposal->world->mstart[pop];
        long msto = proposal->world->mend[pop];
	mm = 0.0;
	long i;
        for (i = msta; i < msto; i++)
	  {
	    if(skyparam[i]<=0.0)
	      {
		warning("skyparam[%li]=%f rate=%f\n",i, skyparam[i],rate);
		error("mismatch with skyparams\n");
	      }
	    else
	      mm += proposal->world->data->geo[i] * proposal->world->param0[i] * skyparam[i]/rate;
	  }
	denom    = mm + (lines * (1.0 / (inheritance*rate*timethetarate)));
      }
    invdenom = 1.0 / denom;
    interval =  (-(LOG (/*nu = */UNIF_RANDUM ())) * invdenom);
    if(interval < 0.0 || isnan(interval))
      {
	error("abort in eventtime()");
      }
    //could be used to calculate p(gn|go): proposal->nu += nu; //adds the probabilities (=random numbers);
    //[overcautious]
    //if (interval < 0.0)
    //    error("interval time is negative");
    if (lines > 0.0)
    {
        if ((UNIF_RANDUM ()) < (mm * invdenom))
        {
            *event = 'm';
            return interval;
        }
        else
        {
            *event = 'c';
            return interval;
        }
    }
    else
    {
        //      printf("mcmc1.c 1653 pop = %li, lines = %li\n", pop, tentry->lineages[pop]);
        *event = 'm';
        return interval;
    }
}

/*--------------------------------------------------------*
 * showsister() 
 * find the sisternode, by going down the branch and up on 
 * the other side again, neglecting the migration nodes.
 */
node *
showsister (node * theNode)
{
    node *tmp = crawlback (theNode);
    
    if (tmp->next->top)
    {
        return crawlback (tmp->next->next);
    }
    else
    {
        if (tmp->next->next->top)
        {
            return crawlback (tmp->next);
        }
        else
        {
            error ("error in treestructure, cannot find sisternode\n");
        }
    }
    return NULL;
}

void
count_migrations (node * p, long *count)
{
  switch(p->type)
    {
    case 'm':
    case 'd':
      *count += 1;
      count_migrations (p->next->back, count);
      break;
    case 'i':
    case 'r':
      count_migrations (p->next->back, count);
      count_migrations (p->next->next->back, count);
    }
}

MYREAL
prob_tree (world_fmt * world, timelist_fmt * tyme)
{
    long j, pop;
    MYREAL mm, cc, ss = 0;
    long tlfrom;
    boolean usem = world->options->usem;
    MYREAL pk;
    const MYREAL mu_rate = world->options->mu_rates[world->locus];
    //    const MYREAL lmu_rate = world->options->lmu_rates[world->locus];
    tyme->tl[0].interval = tyme->tl[0].age;
    for (j = 1; j < tyme->T; j++)
    {
        tyme->tl[j].interval = tyme->tl[j].age - tyme->tl[j - 1].age;
    }
    for (j = 0; j < tyme->T - 1; j++)
    {
        mm = cc = 0.0;
        for (pop = 0; pop < world->numpop; pop++)
        {
            mm += tyme->tl[j].lineages[pop] * world->mig0list[pop];
            cc += tyme->tl[j].lineages[pop] * (tyme->tl[j].lineages[pop] -
					       1) / (world->param0[pop] * mu_rate);
        }
        ss += -(tyme->tl[j].interval) * (mm + cc);
	tlfrom = tyme->tl[j].from;
        if (tlfrom == tyme->tl[j].to)
            ss += LOG2 - LOG (world->param0[tlfrom]*mu_rate);
        else
	  { //@@
	    pk = usem ? world->param0[tlfrom]/mu_rate : world->param0[tlfrom] /(mu_rate * world->param0[pop]); 
            ss += LOG (pk);
	  }
    }
    return ss;
}
