// assignment of individuals to populations
// works only for a per locus base and then may be summarized
// for all loci: each locus is treated indpendently.
// this needs numpop assuming that we know the number of populations
// to simplify the issue with the migration rate matrix
// TODO:
// - admixture so that each allele is treated independently, currently
//   for msat an individual with a.b will be changing state for both alleles
//   (not include admixture) but sequence data will treat each allele independently
// - currently we use a flat uniform prior to pick between populations, thus 
//   assuming we know the number of populations AND also that each has the same prob
//   this should change and use a Dirichlet process prior
//  RESULTS:
//   seems to work for migrate-n debug, and migrate-n-mpi debug, both without heating
//   heating confuses the nodes on the tree I guess a link in the assignment database 
//   is not updated correctly. [check this 2022!]
/*
(c) Peter Beerli, Tallahassee 2013/2022

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
 */
#include "migration.h"
#include "random.h"
#include "tools.h"
#include "tree.h"
#include "sighandler.h"
#include "bayes.h"
#include "migrate_mpi.h"

#include <string.h>  /* for strdup()*/


extern node *showsister (node * theNode);

extern int myID;

//#ifndef HAS_INDIX
#define INDIX(a,b,c) ((a)*(b)+(c))
//#define HAS_INDIX
//#endif
// functions 
void fill_world_unassigned(world_fmt *world);
long find_in_unassignedDB(char *nayme, unassigned_fmt **db);
void empty_world_unassigned(world_fmt *world);
void remove_node_assigndb(world_fmt * world,node *p);
void reset_all_assigned_nodes(world_fmt *world);
void reassign_all_individuals(world_fmt *world);
void swap_unassigned_nodecollection(world_fmt *world1, world_fmt *world2);
void reassign_individual(node * node1, long newpop, long otherpop);
void record_assignment(long locus, world_fmt * world);
void copy_assignment(world_fmt *world1, world_fmt *world2);
void add_unassigned(node *p, long locus, world_fmt *world);
 long reset_unassigned(long index, node *p,  long locus, world_fmt *world);
void allocate_unassigned(node *p,  long locus, world_fmt *world);
void  set_unassigned(node * p, world_fmt * world);
void update_assignment(world_fmt *world);
void report_unassigned(FILE *file, world_fmt *world);
long chooseUnassigned (proposal_fmt * proposal);
#ifdef MPI
void get_assignments (world_fmt * world, option_fmt * options);
#endif
long assign_bypopfreq(world_fmt * world);
long assign_bypastfreq(node * thenode, long numpop);

// function implementations
void fill_world_unassigned(world_fmt *world)
{
  //assignment
  world->unassigned = (unassigned_fmt **) mycalloc(1,sizeof(unassigned_fmt*));
  world->unassigned[0] = (unassigned_fmt *) mycalloc(1,sizeof(unassigned_fmt));
  world->unassigned[0]->index=0;
  world->unassigned[0]->key=(char *) mycalloc(7,sizeof(char));
  world->unassigned[0]->accept=(long *) mycalloc(world->loci,sizeof(long));
  world->unassigned[0]->trials=(long *) mycalloc(world->loci,sizeof(long));
  strcpy(world->unassigned[0]->key,"%$^&!@");
  //  world->has_unassigned = options->has_unassigned;
  world->unassignednum=1;
  world->assigncount = 0;
}

void empty_world_unassigned(world_fmt *world)
{
  long no;
  if (!world->has_unassigned)
    return;
  long i = world->unassignednum;
  while (i>1)
    {
      i -= 1;
      free(world->unassigned[i]->key);
      free(world->unassigned[i]->probloc);
      free(world->unassigned[i]->accept);
      free(world->unassigned[i]->trials);
      for (no=0;no<world->unassigned[i]->ploidy;no++)
	if (world->unassigned[i]->thenodes[no] != NULL)
	    free(world->unassigned[i]->thenodes[no]->freqs);
      free(world->unassigned[i]->thenodes);
      free(world->unassigned[i]);
    }
  world->unassignednum=1;
  world->unassigned[0]->next = NULL;
}

// this may not work after all because the reporting needs all individuals
void remove_node_assigndb(world_fmt * world,node *p)
{
  long index;
  long i;
  index = find_in_unassignedDB(p->truename,world->unassigned);
  if (index != UNKNOWN)
    {
      unassigned_fmt *t = world->unassigned[index];
      for (i= (long) t->ploidy-1; i >= 0; i--)
	{
	  if(p == t->thenodes[i])
	    {
	      if (i!= (long) t->ploidy-1)
		{
		  t->thenodes[i] = t->thenodes[t->ploidy-1];
		  t->ploidy -= 1;
		}
	      else
		{
		  t->thenodes[i] = NULL;
		  t->ploidy -= 1;
		}
	      break;
	    }
	}
    }
}

void reset_all_assigned_nodes(world_fmt *world)
{
   long i;
   long j;
  for (i = 1;i < world->unassignednum; i++)
    {
      if (world->unassigned[i]->key != NULL)
	{
	  for(j=0;j<world->unassigned[i]->ploidy;j++)
	    world->unassigned[i]->thenodes[j] = NULL;
	}
    }
}

// used with heating
void reassign_all_individuals(world_fmt *world)
{
   long i, j;
  //  char *nayme;
  char *key;
   long ploidy;
   long testploidy;
  node **thenodes;
  node *thenode;
#ifdef DEBUG
  printf("%i> reassign_all_individuals()\n",myID);
#endif  
  for (i=1;i<world->unassignednum;i++)
    {
      if (world->unassigned[i]->key != NULL)
	key = world->unassigned[i]->key;
      else
	continue;
      thenodes = world->unassigned[i]->thenodes;
      if(thenodes==NULL)
	continue;
      ploidy = world->unassigned[i]->ploidy;
      testploidy = 0;
      for (j=0; j<world->sumtips;j++)
	{
	  // get the node from the tree (this may be a tree that was hotter or colder before)
	  thenode = world->nodep[j];
	  if (thenode->type != 't')
	    {
	      error("failed to show a top in ind assign for heating");
	    }
	  if (thenode->truename==NULL)
	    {
	      error("true name missing");
	    }
	  else
	    {
	      if (!strcmp(thenode->truename,key))
		{
		  thenodes[testploidy]=thenode;
		  testploidy++;
		}
	    }
	  // if the namye occures twice (e.g. ploidy=2) then stop
	  // because we already found the node(s) and do not need to search further 
	  if (testploidy == ploidy)
	    break;
	}
    }
}

void swap_unassigned_nodecollection(world_fmt *world1, world_fmt *world2)
{
#ifdef DEBUG
  printf("%i> swap unassigned: heat1=%f heat2=%f\n",myID, world1->heat,world2->heat);
#endif
  reassign_all_individuals(world1);
  reassign_all_individuals(world2);
}


//used for origin which is dangling so no migration events need to be considered
void reassign_individual(node * node1, long newpop, long otherpop)
{
  if (node1->type != 't')
    return;
  if (node1->truepop == UNKNOWN)
    {
      if (node1->actualpop != newpop)
	{
	  node1->actualpop = newpop;
	  node1->pop = newpop;
	  node1->freqs[newpop] += 1.0;
	}
      if (otherpop != -1)
	{
	  node1->freqs[otherpop] -= 1.0;
	}
    }
}


// records a node and its population assignment
void record_assignment( long locus, world_fmt * world)
{
   long j;
   long index;
  node *p;
  const  long len=world->unassignednum;
  for (index=1; index < len; index++)
    {
      for(j=0;j<world->unassigned[index]->ploidy;j++)
      	{
	  // this assumes all nodes talk about the same individuals
	  p = world->unassigned[index]->thenodes[j];
	  if (p!=NULL)
	    world->unassigned[index]->probloc[INDIX(world->numpop, locus,( long) p->actualpop)] += 1;
	}
    }
}
  
// for swapping temperatures
// this uses the swapped trees to match the unassigned
// up this will change the population assignment of the 
// tips 
void copy_assignment(world_fmt *world1, world_fmt *world2)
{
   long i,j;
  unassigned_fmt *db;
  node * node1, *node2;
#ifdef DEBUG
  printf("%i> WRONG?????? copy unassigned: heat1=%f heat2=%f\n",myID, world1->heat,world2->heat);
#endif
  for (i=0;i<world1->unassignednum;i++)
    {
      db = world1->unassigned[i];
      for(j=0;j<world1->unassigned[i]->ploidy;j++)
	{
	  node1 = db->thenodes[j];
	  node2 = world2->unassigned[i]->thenodes[j];
	  long temppop = node2->actualpop;
	  node2->actualpop = node2->pop = node1->actualpop;
	  node1->actualpop = node1->pop = temppop;
	}  
    }
}	      

long find_in_unassignedDB(char *nayme, unassigned_fmt **db)
{
  unassigned_fmt *d = NULL;

  if (db == NULL)
    error("unassigned DB not allocated");
  else
    d = *db; // static analyze disliked db[0]
  if (d==NULL)
    return UNKNOWN;

  while (strcmp(nayme,d->key))
    {
      if (d->next != NULL)
	d = d->next;
      else
	return UNKNOWN;// -1
    }
  return d->index; 
}

void add_unassigned(node *p,  long locus, world_fmt *world)
{
  const  long numpop = world->numpop;
  long newpop ;
  if (world->has_unassignedpoplist)
    newpop = ( long) world->unassignedpoplist[RANDINT(0,(long) world->unassignedpoplistnum-1)];
  else
    newpop = ( long) RANDINT(0,(long) numpop-1);
  unassigned_fmt *temp;
  if (p->type != 't')
    error("node is not tip");

  temp = (unassigned_fmt *) mycalloc(1,sizeof(unassigned_fmt));
  temp->key = strdup(p->truename);
  temp->accept=(long *) mycalloc(world->loci,sizeof(long));
  temp->trials=(long *) mycalloc(world->loci,sizeof(long));
  temp->probloc = (MYREAL *) mycalloc((size_t) (world->loci * world->numpop), sizeof(MYREAL));
  //temp->probloc[INDIX(world->numpop,locus,newpop)] = 0.0;
  temp->thenodes = (node **) mycalloc(1,sizeof(node *));

  if (p->freqs == NULL)
    {
      p->freqs = (double*) calloc((size_t) numpop,sizeof(double));
      int i;
      for (i=0;i<numpop;i++) p->freqs[i]=1.0;
    }
  else
    {
      p->freqs = (double*) realloc(p->freqs, ((size_t) numpop) * sizeof(double));
      int i;
      for (i=0;i<numpop;i++) p->freqs[i]=1.0;
    }
  temp->thenodes[0] = p;

  temp->index = (long) world->unassignednum;
  temp->lastlocus = locus;
  temp->lastpop = newpop;
  temp->ploidy=1;

  world->unassigned[temp->index-1]->next = temp;
  world->unassignednum++;
  world->unassigned = (unassigned_fmt **) myrealloc(world->unassigned,sizeof(unassigned_fmt*)* (size_t) world->unassignednum);
  world->unassigned[temp->index] = temp;
#ifdef DEBUG 
  printf("%i> added %s (%li) to unassigned db [heat=%f,replicate=%li]\n",myID, temp->key,temp->index, world->heat,
	 world->replicate);
#endif
  //no return value used ever: return world->unassignednum;
}

// resets the values for a new locus! to zero,
// but also checks whether this is a diploid etc that will need to have 
// symmetric updates of population labels for the reassignment
 long reset_unassigned(long index, node *p,  long locus, world_fmt *world)
{
  long  replicate = world->replicate;
  const  long numpop = world->numpop;
  //long pop;
  unassigned_fmt **db = world->unassigned;
  unassigned_fmt *temp = db[index];
  long newpop;
  if (world->has_unassignedpoplist)
    newpop = ( long) world->unassignedpoplist[RANDINT(0,(long) world->unassignedpoplistnum-1)];
  else
    newpop = ( long) RANDINT(0,(long) numpop-1);

#ifdef DEBUG
  printf("%i> reset_unassigned()\n", myID);
#endif
  if (p->type != 't')
    error("node is not tip");
  if(strcmp(p->truename,temp->key))
    error("assignment Database is incorrect -- check reset_unassigned()");
  if (temp->lastlocus == locus )
    {
      if (temp->lastreplicate == replicate)
	{
	  temp->thenodes = (node **) myrealloc(temp->thenodes, (temp->ploidy+1)*sizeof(node *));

	  if (p->freqs == NULL)
	    {
	      p->freqs = (double*) calloc((size_t) numpop,sizeof(double));
	      int i;
	      for (i=0;i<numpop;i++) p->freqs[i]=1.0;
	    }
	  else
	    {
	      p->freqs = (double*) realloc(p->freqs, ((size_t) numpop) * sizeof(double));
	      int i;
	      for (i=0;i<numpop;i++) p->freqs[i]=1.0;
	    }
	  temp->thenodes[temp->ploidy] = p;
	  temp->ploidy++;
	  //      printf("changed %s (%li) [ploidy=%li] in unassigned-db\n",temp->key,temp->index,temp->ploidy);
	}
      else
	{
	  // new replicate adding to the earlier replicate, but needs new nodes
	  temp->ploidy=1;
	  if (p->freqs == NULL)
            {
              p->freqs = (double*) calloc((size_t) numpop,sizeof(double));
              int i;
              for (i=0;i<numpop;i++) p->freqs[i]=1.0;
            }
          else
            {
              p->freqs = (double*) realloc(p->freqs, ((size_t) numpop) * sizeof(double));
              int i;
              for (i=0;i<numpop;i++) p->freqs[i]=1.0;
            }

	  temp->thenodes[0] = p;	  
	}
    }
  else
    {
      //for (pop=0;pop<numpop;pop++)
      //temp->probloc[INDIX(numpop,locus,pop)] = 0.0;
      temp->lastpop = newpop;
      temp->ploidy = 1;
      if (p->freqs == NULL)
	{
	  p->freqs = (double*) calloc((size_t) numpop,sizeof(double));
	  int i;
	  for (i=0;i<numpop;i++) p->freqs[i]=1.0;
	}
      else
	{
	  p->freqs = (double*) realloc(p->freqs, ((size_t) numpop) * sizeof(double));
	  int i;
	  for (i=0;i<numpop;i++) p->freqs[i]=1.0;
	}
      temp->thenodes[0] = p;
    }
  temp->lastlocus = locus;
  temp->lastreplicate = replicate;
  return temp->lastpop;
}


void allocate_unassigned(node *p,  long locus, world_fmt *world)
{
  long index;
   long newpop;
  index = find_in_unassignedDB(p->truename,world->unassigned);
  if (index == UNKNOWN)
    {
      /*newpop =*/ add_unassigned(p, locus, world); // static-analyze Value stored to 'newpop' is never read
    }
  else
    {
      newpop = reset_unassigned(index, p, locus, world);
      p->actualpop = p->pop = (long) newpop; //syncs the individual for diploid data
    }
}

void  set_unassigned(node * p, world_fmt * world)
{
   long locus = world->locus;
  if (world->has_unassigned)
    {
      if (p->truename[0] == '?')
	{
	  allocate_unassigned(p, locus,world);
	  p->truepop = UNKNOWN;
	} 
    }
}

long assign_bypopfreq(world_fmt * world)
{
  long i,j;
  double s=0.0;
  node **nodelist;
  long nl = 0;
  double *freqs;
  nodelist = (node **) mycalloc ((world->sumtips + 1), sizeof (node *));
  freqs = (double *) mycalloc (world->numpop, sizeof (double));
  find_tips (world->root, nodelist, &nl);
  for (j=0;j<world->numpop;j++)
    {
      for (i=0; i<nl; i++)
	{
	  if (nodelist[i]->actualpop == j)
	    freqs[j] += 1.0;
	}
    }
  for (j=0;j<world->numpop;j++)
    {
      if (freqs[j]<BIGEPSILON)
	{
	  freqs[j]=TENTH;
	}
      s += freqs[j];
    }
  for (j=0;j<world->numpop;j++)
    {
      freqs[j] /= s;
    }
  return random_from_freqlist(freqs, world->numpop);
}

long assign_bypastfreq(node * thenode, long numpop)
{
  double *freqs = thenode->freqs; 
  long pop = random_from_freqlist(freqs, numpop);
  //freqs[pop] += 1;
  return pop;
}

void update_assignment(world_fmt *world)
{
   long i;
  long r = RANDINT(1 , (long) world->unassignednum-1);
  node ** nodes = world->unassigned[r]->thenodes;
   long ploidy   = world->unassigned[r]->ploidy;
   long newpop;
   if (world->has_unassignedfreq)
     //newpop = assign_bypopfreq(world);
     newpop = assign_bypastfreq(nodes[r], world->numpop);
   else
     {
       if (world->has_unassignedpoplist)
	 newpop = ( long) world->unassignedpoplist[RANDINT(0,(long) world->unassignedpoplistnum-1)];
       else
	 newpop = ( long) RANDINT(0,(long) world->numpop-1);
       //nodes[r]->freqs[newpop] += 1.0;
     }   
  for (i=0; i<ploidy;i++)
    {
      reassign_individual(nodes[i], newpop, -1);
#ifdef DEBUG
      printf("%i> assignment: old:%li new:%li\n", myID, nodes[i]->actualpop, newpop);
#endif
    }
}

void report_unassigned(FILE *file, world_fmt *world)
{
   long i;
   long idi;
   long locus;
   long pop;
  MYREAL sum = 0.0;
  MYREAL lsum;
  MYREAL *total;
  MYREAL totalsum;
  char *key;
  MYREAL maxtotal;
  MYREAL val;
  total = (MYREAL *) mycalloc(world->numpop, sizeof(MYREAL));
  if (world->has_unassigned)
    {
      fprintf(file,"\n\n\n\nAssignment of Individuals to populations\n");
      fprintf(file,"========================================\n\n");
      fprintf(file,"Individual Locus    ");
      for (pop=0;pop<world->numpop;pop++)
	{
	  fprintf(file,"%5li ",pop+1);
	}
      fprintf(file,"\n");
      fprintf(file,"---------- -------- ");
      for (pop=0;pop<world->numpop;pop++)
	{
	  fprintf(file,"----- ");
	}
      fprintf(file,"\n");
      for (i=1;i<world->unassignednum;i++)
	{
	  key = world->unassigned[i]->key;
	  memset(total,0,sizeof(MYREAL)* (size_t) world->numpop);
	  for (locus=0;locus<world->loci;locus++)
	    {
	      fprintf(file, "%-10.10s ",key);
	      sum = 0.0;
	      for (pop=0;pop<world->numpop;pop++)
		{
		  val = world->unassigned[i]->probloc[INDIX(world->numpop,locus,pop)];
		  if(val==0.0)
		    {
		      val = world->unassigned[i]->probloc[INDIX(world->numpop,locus,pop)]=SMALL_VALUE;
		    }
		  sum += val;
		}
	      fprintf(file," %8li ",locus+1);
	      for (pop=0;pop<world->numpop;pop++)
		{
		  idi = INDIX(world->numpop,locus,pop);
		  fprintf(file, "%5.3f ", world->unassigned[i]->probloc[idi]/sum);
		  total[pop] += log(world->unassigned[i]->probloc[idi]/sum);
		  
		}
	      fprintf(file, "[%li/%li = %f]\n", world->unassigned[i]->accept[locus],world->unassigned[i]->trials[locus],
		      (double) world->unassigned[i]->accept[locus]/world->unassigned[i]->trials[locus]);	   
	    }
	  maxtotal = (MYREAL) -HUGE;
	  for (pop=0;pop<world->numpop;pop++)
	    {
	      if (maxtotal < total[pop])
		maxtotal = total[pop];
	    }
	  totalsum = 0.0;
	  for (pop=0;pop<world->numpop;pop++)
	    {
	      totalsum += exp(total[pop]-maxtotal);
	    }
	  lsum = log(totalsum);
	  fprintf(file, "%-10.10s       All ",key);
	  for (pop=0;pop<world->numpop;pop++)
	    {
	      fprintf(file, "%5.3f ", exp(total[pop]-maxtotal-lsum));
	    }
	  fprintf(file, "\n");	  
	}
    }
  myfree(total);
}

///
/// pick a node at random from the unassigned list
long chooseUnassigned (proposal_fmt * proposal)
{
  //long elem = 0;
  //node *tmp=NULL;//, **goal;
  long index= -1;
  long j = -1;
  long counter = 0;
  if(proposal->world->unassignednum>0)
    {
      while(index<0 && counter++ < 100)
	{
	  index = RANDINT(1,(long) proposal->world->unassignednum-1);
	  if (proposal->world->unassigned[index]->thenodes[0]==NULL)
	    index = -1;
	}
      if (counter > 99)
	return FAILURE;
      else
	counter = 0;
    }
  else
    return FAILURE;
  while(j== -1 && counter++ < 100)
    {
      if(proposal->world->unassigned[index]->ploidy>0)
	j = RANDINT(0, (long) proposal->world->unassigned[index]->ploidy-1);
      else
	{
	  j = -1;
	  index = RANDINT(1,(long) proposal->world->unassignednum-1);
	}
    }
  if (counter>99)
    return index; //FAILURE; is -1

  proposal->origin = proposal->world->unassigned[index]->thenodes[j];
  proposal->world->unassigned[index]->trials[proposal->world->locus] += 1;
  if (!proposal->origin->tip)
    error("unassigned individual is not a tip?????");
  if (proposal->origin->back == NULL)
    {
	warning("proposal->orig is not part of the tree");
	error("");
      }
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
    return index ; // SUCCESS;
}


///
/// save all population assignments from the worker nodes
#ifdef MPI
void get_assignments (world_fmt * world, option_fmt * options)
{
    long maxreplicate = (options->replicate
                         && options->replicatenum >
                         0) ? options->replicatenum : 1;
    
    if (myID == MASTER && world->has_unassigned)
    {
        mpi_results_master (MIGMPI_ASSIGN, world, maxreplicate,
                            unpack_assign_buffer);
    }
}
#endif
