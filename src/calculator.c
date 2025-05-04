/*
(c) Peter Beerli 2010-12 Tallhassee
beerli@fsu.edu

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
#ifdef BEAGLE
#include "beagle.h"
#include "migration.h"
#include "tree.h"
#include "tools.h"
#ifdef PRETTY
#include "pretty.h"
#endif
#include "sighandler.h"

//#define BEAGLE_PARTIALS 7

extern long myID;
extern int use_beagle_gpu;
extern int use_beagle_dynamicscale;
extern int use_beagle_autoscale;
extern int use_beagle_manualscale;

void debug_partials(int nodenum, mutationmodel_fmt *s, long instance,beagle_fmt *beagle, world_fmt *world);
void debug_beagle(beagle_fmt *beagle);
long new_id(long id, long sumtips);
void reset_beagle(beagle_fmt *beagle);
/*-----------------------------------------------------------------------------
|	This function sets up the beagle library and initializes all data members.
*/
void init_beagle(world_fmt *world, long locus)
{
  int numscalingbuffers;
  beagle_fmt *beagle = world->beagle;
  beagle->scaling = !world->options->fastlike;
  if(beagle->scaling)
    {
      numscalingbuffers = (4 *(world->sumtips));
      beagle->ievectrans = FALSE;
      beagle->logscalers = TRUE;
      beagle->eigencomplex = FALSE;
      beagle->dynamicscaling = FALSE;
      beagle->autoscaling = TRUE;
      beagle->requireDoublePrecision = FALSE;
      beagle->requireSSE = TRUE;
    }
  else
    {
      numscalingbuffers = BEAGLE_OP_NONE;
      beagle->ievectrans = FALSE;
      beagle->logscalers = FALSE;
      beagle->eigencomplex = FALSE;
      beagle->dynamicscaling = FALSE;
      beagle->autoscaling = FALSE;
      beagle->requireDoublePrecision = FALSE;
      beagle->requireSSE = TRUE;
    }

  beagle->instance_handle    = (int *) mycalloc(world->numsubloci[locus],sizeof(int));
  beagle->numallocoperations = 200;
  beagle->operations         = (BeagleOperation *) mycalloc(beagle->numallocoperations,sizeof(BeagleOperation));
  beagle->branch_indices     = (int *) mycalloc(beagle->numallocoperations * 2, sizeof(int));
  beagle->branch_lengths     = (double *) mycalloc(beagle->numallocoperations * 2, sizeof(double));
  if(numscalingbuffers!=BEAGLE_OP_NONE)
    {
      beagle->scalingfactorsindices   = (int *) mycalloc(numscalingbuffers, sizeof(int));
    }
  beagle->scalingfactorscount  = 0; //actual count of which scale factors are updated? this was world->sumtips-1;
  long sublocus;
  mutationmodel_fmt *s;
  const long sublocistart = world->sublocistarts[locus];
  const long sublociend   = world->sublocistarts[locus+1];
  for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
    {
      s = &world->mutationmodels[sublocus];
      s->numallocpartials    = s->numsiterates * world->sumtips * s->numstates * s->numpatterns;
      s->partials            = (double *) mycalloc(s->numallocpartials, sizeof(double));
      s->numpartials = 0;
    }      
  //  beagle->numallocallyweights = world->mutationmodels[0].numsites;
  //beagle->allyweights = (int *) mycalloc(beagle->numallocallyweights,sizeof(int));
  reset_beagle(beagle);
} 
/*-----------------------------------------------------------------------------
|	This function reinitializes the beagle library.
*/
void reinit_beagle(world_fmt *world, long locus)
{
  beagle_fmt *beagle = world->beagle;
  beagle->scaling = !world->options->fastlike;
  long sublocus;
  mutationmodel_fmt *s;
  const long sublocistart = world->sublocistarts[locus];
  const long sublociend   = world->sublocistarts[locus+1];
  for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
    {
      s = &world->mutationmodels[sublocus];
      s->numpartials = 0;
    }
  beagle->numoperations = 0;
  beagle->scalingfactorscount  = 0; //actual count of which scale factors are updated? this was world->sumtips-1;      
  reset_beagle(beagle);
} 

/*-----------------------------------------------------------------------------
|	This function deletes the beagle library and destroys all data members.
*/
void destroy_beagle(world_fmt *world, long locus)
{
  int numscalingbuffers = world->options->fastlike ? BEAGLE_OP_NONE : (4 * world->sumtips);
  beagle_fmt *beagle = world->beagle;
  myfree(beagle->instance_handle);
  long sublocus;
  mutationmodel_fmt *s;
  const long sublocistart = world->sublocistarts[locus];
  const long sublociend   = world->sublocistarts[locus+1];
  for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
    {
      s = &world->mutationmodels[sublocus];
      myfree(s->partials);
    }

  if(numscalingbuffers!=BEAGLE_OP_NONE)
    {
      myfree(beagle->scalingfactorsindices);
    }
  myfree(beagle->branch_lengths);
  myfree(beagle->operations);
  myfree(beagle->branch_indices);
  //  myfree(beagle->allyweights);
} 


void print_beagle_available_resources(world_fmt *world)
{
  int i;
  BeagleResourceList* rList;
  rList = beagleGetResourceList();
  fprintf(stdout, "Available resources for likelihood calculator:\n");
  for (i = 0; i < rList->length; i++) 
    {
      fprintf(stdout, "\tResource %i:\n\t\tName : %s\n", i, rList->list[i].name);
      fprintf(stdout, "\t\tDesc : %s\n\n", rList->list[i].description);
    }
}


void print_beagle_resources(BeagleInstanceDetails instDetails, world_fmt *world)
{
  if(myID==MASTER)
    {
      if(world->options->progress)
	{
	  fprintf(stdout, "Likelihood calculator using HMSBEAGLE\n");
	  fprintf(stdout, "-------------------------------------\n");
	  fprintf(stdout, "\tRsrc Name : %s\n", instDetails.resourceName);
	  fprintf(stdout, "\tImpl Name : %s\n", instDetails.implName);
	  fprintf(stdout, "\n");        
	}
      fprintf(world->outfile, "Likelihood calculator using HMSBEAGLE\n");
      fprintf(world->outfile, "-------------------------------------\n");
      fprintf(world->outfile, "\tRsrc Name : %s\n", instDetails.resourceName);
      fprintf(world->outfile, "\tImpl Name : %s\n", instDetails.implName);
      fprintf(world->outfile, "\n");        
#ifdef PRETTY
      print_beagle_resources_pretty(instDetails);
#endif
    }
}

void  set_beagle_instances(world_fmt *world, long locus)
{
  int code;
  beagle_fmt *beagle = world->beagle;
  int resource = use_beagle_gpu;//cpu=0, gpu=1
  boolean ievectrans=beagle->ievectrans;
  boolean logscalers=beagle->logscalers;
  boolean eigencomplex=beagle->eigencomplex;
  boolean dynamicscaling=beagle->dynamicscaling;
  boolean autoscaling=beagle->autoscaling;
  boolean requireDoublePrecision=beagle->requireDoublePrecision;
  boolean requireSSE=beagle->requireSSE;

  static boolean done=FALSE;


  int numscalingbuffers;
  
  long sublocus;
  const long sublocistart = world->sublocistarts[locus];
  const long sublociend   = world->sublocistarts[locus+1];
  BeagleInstanceDetails instDetails;

  numscalingbuffers = (beagle->scaling ? (4 *(world->sumtips)) : 0);

  for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
    {
      const long xs = sublocus - sublocistart;
      mutationmodel_fmt *s = &world->mutationmodels[sublocus];
      beagle->instance_handle[xs] = beagleCreateInstance(
	       world->sumtips,            // number of tips
	       2 * (2 * world->sumtips-1), // total buffers (= total nodes in tree) twice to accomodate rejections.
	       0, // number of compact state representation [funny tips]
	       s->numstates,  // number of states (nucleotides etc)
	       s->numpatterns + s->addon,  // number of site patterns
	       1, // number of rate matrices eigen-decomp buffers [eigencount, but no clue why this could be bigger than 1]
	       2*(2 * world->sumtips - 1),// (=number of branches)
	       s->numsiterates,           // categoryCount
	       numscalingbuffers,         // number of scaling buffers (times eigencount)
	       &resource,		   // List of resources 
	       1,			   // Number of resources
	       //	       BEAGLE_FLAG_VECTOR_SSE | (beagle->scaling ? BEAGLE_FLAG_SCALING_AUTO : 0),  // preferenceFlags
	       //(ievectrans ? BEAGLE_FLAG_INVEVEC_TRANSPOSED : BEAGLE_FLAG_INVEVEC_STANDARD) |
                (logscalers ? BEAGLE_FLAG_SCALERS_LOG : BEAGLE_FLAG_SCALERS_RAW) |
                (eigencomplex ? BEAGLE_FLAG_EIGEN_COMPLEX : BEAGLE_FLAG_EIGEN_REAL) |
                (dynamicscaling ? BEAGLE_FLAG_SCALING_DYNAMIC : 0) | 
                (autoscaling ? BEAGLE_FLAG_SCALING_AUTO : 0) |
                (requireDoublePrecision ? BEAGLE_FLAG_PRECISION_DOUBLE : BEAGLE_FLAG_PRECISION_SINGLE) |
                (requireSSE ? BEAGLE_FLAG_VECTOR_SSE : BEAGLE_FLAG_VECTOR_NONE),	  /**< Bit-flags indicating required implementation characteristics, see BeagleFlags (input) */
	       0L,			   // requirementFlags (see BeagleFlags) 
	       &instDetails
	     );
      if (beagle->instance_handle[xs] < 0)
	{
	  usererror("beagleCreateInstance returned a negative instance handle (and that's not good)");
      	}
      else
	{
	  code = beagleSetStateFrequencies(beagle->instance_handle[xs],
				    0, //state frequencies index
				    s->basefreqs);
	  if (code != 0)
	    usererror("setStateFrequencies encountered a problem");

	  code = beagleSetCategoryWeights(beagle->instance_handle[xs],
				    0, //category weights  index
				    s->siteprobs);
	  if (code != 0)
	    usererror("setCategoryWeights encountered a problem");

	  code = beagleSetCategoryRates(beagle->instance_handle[xs], 
					s->siterates);
	  if (code != 0)
	    usererror("setCategoryRates encountered a problem");

	  code = beagleSetPatternWeights(beagle->instance_handle[xs], s->aliasweight);
	  if (code != 0)
	    usererror("setPatternWeights encountered a problem");

	  if(!done)
	    {
	      done=TRUE;
	      print_beagle_resources(instDetails, world);
	    }
	  if (!(instDetails.flags & BEAGLE_FLAG_SCALING_AUTO))
	    autoscaling = FALSE;
	}
    }
}


long
set_branch_index (node * p,  long *bid)
{
  long bbid = *bid;
    if (p->type != 't')
    {
      bbid = set_branch_index (crawlback (p->next), &bbid);
      bbid = set_branch_index (crawlback (p->next->next), &bbid);
    }
    p->bid = bbid;
    //printf("%li -- %li\n",p->id, bbid);
    return bbid+1;
}    /* set_branch_index */


void adjust_beagle(beagle_fmt *beagle)
{
  beagle->numallocoperations = beagle->numoperations + 10;
  beagle->operations = (BeagleOperation *) myrealloc(beagle->operations, sizeof(BeagleOperation)  * beagle->numallocoperations);
  beagle->branch_lengths   = (double *) myrealloc(beagle->branch_lengths, sizeof(double) * 2 * beagle->numallocoperations);
  beagle->branch_indices = (int *) myrealloc(beagle->branch_indices, beagle->numallocoperations * 2 * sizeof(int));
  //beagle->numallocbranches = beagle->numallocoperations * 2;
}


void prepare_beagle_instances(node *theNode, node * left, node *right, beagle_fmt *beagle)
{
  const boolean scaling = beagle->scaling;
  const long parentid = theNode->id;
  const long rightid  = right->id;
  const long leftid   = left->id;
 
  const long rightbid  = right->bid;
  const long leftbid   = left->bid;

  const double leftbranch = left->v;
  const double rightbranch = right->v;
  long ii = beagle->numoperations;
  long i = ii;
  long j = ii * 2;

  int *branch_indices = beagle->branch_indices;
  double *branch_lengths = beagle->branch_lengths;
  BeagleOperation *operations     = beagle->operations;
  int *scalingfactorsindices = beagle->scalingfactorsindices;

  if(beagle->numoperations  >=  beagle->numallocoperations)
    {
      adjust_beagle(beagle);
    }
  branch_indices[j] = leftbid;
  branch_indices[j+1] = rightbid;

  operations[i].destinationPartials   = parentid;
  operations[i].destinationScaleWrite = scaling ? parentid : BEAGLE_OP_NONE; 
  operations[i].destinationScaleRead  = BEAGLE_OP_NONE;//scaling ? parentid : BEAGLE_OP_NONE; 
  operations[i].child1Partials  = leftid;
  operations[i].child1TransitionMatrix = leftbid;
  operations[i].child2Partials = rightid;
  operations[i].child2TransitionMatrix = rightbid;

  branch_lengths[j]     = leftbranch; 
  branch_lengths[j+1]   = rightbranch; 
  if(scaling)
    {
      scalingfactorsindices[ii] = parentid;
      beagle->scalingfactorscount += 1;
    }
  beagle->numbranches += 2;
  beagle->numoperations++;
#ifdef BEAGLEDEBUG
  printf("%li> TREE: %i {%i %i %i %i %i %i %i} {%i %i} {%f %f}\n",myID, beagle->numoperations,
	 beagle->operations[i],   
	 beagle->operations[i+1] ,
	 beagle->operations[i+2] ,
	 beagle->operations[i+3] ,
	 beagle->operations[i+4] ,
	 beagle->operations[i+5] ,
	 beagle->operations[i+6] ,
	 beagle->branch_indices[j],
	 beagle->branch_indices[j+1],                 
	 beagle->branch_lengths[j],
	 beagle->branch_lengths[j+1]);
#endif
}

void prepare_beagle_instances_proposal(proposal_fmt *proposal, long trueparentid, long leftid, long leftbid, double leftbranch, 
				       long rightid, long rightbid, double rightbranch, beagle_fmt *beagle)
{
  long ii= beagle->numoperations;
  long i = ii;
  long j = ii * 2;
  const long nodep_boundary = proposal->world->sumtips * 2 - 1;
  const boolean scaling = beagle->scaling;
  long parentid = trueparentid > nodep_boundary ? trueparentid - nodep_boundary : trueparentid + nodep_boundary; 

  int *branch_indices = beagle->branch_indices;
  double *branch_lengths = beagle->branch_lengths;
  BeagleOperation *operations     = beagle->operations;
  int *scalingfactorsindices = beagle->scalingfactorsindices;

  if(beagle->numoperations  >=  beagle->numallocoperations)
    {
      adjust_beagle(beagle);
    }
#ifdef BEAGLE_DEBUG
  printf("PINT: (%li, %li) l:%li,%f, r:%li,%f\n",trueparentid,parentid,leftid,leftbranch,rightid,rightbranch);
#endif
  operations[i].destinationPartials   = parentid;
  operations[i].destinationScaleWrite = scaling ? parentid : BEAGLE_OP_NONE; 
  operations[i].destinationScaleRead  = BEAGLE_OP_NONE;//scaling ? parentid : BEAGLE_OP_NONE; 
  operations[i].child1Partials  = leftid;
  operations[i].child1TransitionMatrix = leftbid;
  operations[i].child2Partials = rightid;
  operations[i].child2TransitionMatrix = rightbid;

  branch_indices[j] = leftbid;
  branch_indices[j+1] = rightbid;

  branch_lengths[j]     = leftbranch; 
  branch_lengths[j+1]   = rightbranch; 

  if(scaling)
    {
      scalingfactorsindices[ii] = parentid;
      beagle->scalingfactorscount += 1;
    }
  beagle->numbranches += 2;
  beagle->numoperations++;
#ifdef BEAGLEDEBUG
  printf("%li> PROP: %i {%i %i %i %i %i %i %i} {%i %i} {%f %f}\n",myID,beagle->numoperations,
	 beagle->operations[i],   
	 beagle->operations[i+1] ,
	 beagle->operations[i+2] ,
	 beagle->operations[i+3] ,
	 beagle->operations[i+4] ,
	 beagle->operations[i+5] ,
	 beagle->operations[i+6] ,
	 beagle->branch_indices[j],
	 beagle->branch_indices[j+1],                 
	 beagle->branch_lengths[j],
	 beagle->branch_lengths[j+1]);
#endif
}

///
/// prepares conditional likelihood down to the root, this assumes that the newtree and the oldtree
/// are the same once this routine is used
void evaluate_beagle_instances_proposal (proposal_fmt * proposal,
				    node * mother,  
				    node * newdaughter, long newdaughter_id, long newdaughter_bid, 
				    MYREAL v)
{
    node *nn = NULL, *d1 = NULL, *d2 = NULL, *oldnn = NULL;

    if (mother->type != 'r')
    {
        children (mother, &d1, &d2);
        if (d1 == newdaughter)
	  {
	    d1 = d2;
	  }
	long id0, id1, id2;
	long bid1, bid2;
	double v1, v2;
	id0  = mother->id; 
	id1  = new_id(newdaughter_id, proposal->world->sumtips);
	bid1 = newdaughter_bid;
	v1 = v;
	id2  = d1->id;
	bid2 = d1->bid;
	v2 = d1->v;
	prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);

        oldnn = mother;
        nn = showtop (crawlback (mother));
        while (nn->type != 'r')
        {
            children (nn, &d1, &d2);
            if (d1 == oldnn)
	      {
		d1 = d2;
		d2 = oldnn;
	      }
	    if(d2->id != oldnn->id)
	      {
		warning("One of the children is not what it should be\n");
	      }
	    id0  = nn->id; 
	    id1  = new_id(d2->id,proposal->world->sumtips);;
	    bid1 = d2->bid;
	    v1 = d2->v;
	    id2  = d1->id;
	    bid2 = d1->bid;
	    v2 = d1->v;
	    prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);

            oldnn = nn;
            nn = showtop (crawlback (nn));
        }
    }
}


void fill_beagle_instances(world_fmt *world, long locus)
{
  beagle_fmt *beagle = world->beagle;
  //unsigned long i;
  //unsigned long ii;
  unsigned long j;
  //unsigned long zz;
  unsigned long site;
  unsigned int numpatterns;
  unsigned int numstates;
  //unsigned int numweights=0;
  int code;
  // set partials for all tipnodes
  // use z to advance through the partial array
  unsigned long z = 0L;
  mutationmodel_fmt *s;
  long instance;
  long sublocus;
  const long sublocistart = world->sublocistarts[locus];
  const long sublociend   = world->sublocistarts[locus+1];
  for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
    {
      const long xs = sublocus - sublocistart;
      instance = beagle->instance_handle[xs];
      s = &world->mutationmodels[sublocus];
      numstates = s->numstates;
      numpatterns = s->numpatterns + s->addon;
      if(s->partials!=NULL)
	{
	  s->numallocpartials = world->sumtips * numpatterns * numstates;
	  s->partials = (double *) myrealloc(s->partials, sizeof(double) * s->numallocpartials);
	}
      else
	{
	  s->numallocpartials = world->sumtips * numpatterns * numstates;
	  s->partials = (double *) mycalloc( s->numallocpartials,sizeof(double));
	}
      s->numpartials=0;
      z = 0L;
      for (j = 0; j < world->sumtips; ++j)
	{
	  for(site=0;site<numpatterns;site++)
	    {
	      memcpy(&s->partials[z], 
		     &(world->nodep[j]->x[xs].s[site][0][0]),
		     sizeof(double) * numstates);
	      s->numpartials++;
	      z += numstates;
	    }
	  code = beagleSetTipPartials(instance,       // instance
				      j,				// indicator for tips
				      &s->partials[z - numpatterns*numstates]);// inPartials
	  if (code != 0)
	    usererror("setTipPartials encountered a problem");
	}
#ifdef BEAGLEDEBUG
      // debug_eigensystem(ii, beagle, world);
#endif
      code = beagleSetEigenDecomposition(instance,		  // instance
				   0,					  // eigenIndex,
				   (const double *)s->eigenvectormatrix,	// inEigenVectors,
				   (const double *)s->inverseeigenvectormatrix,// inInverseEigenVectors,
				   s->eigenvalues); // inEigenValues      
      if (code != 0)
	usererror("setEigenDecomposition encountered a problem");
    }
}
/*-----------------------------------------------------------------------------
|	Calculates the log likelihood by calling the beagle functions
|	updateTransitionMatrices, updatePartials and calculateEdgeLogLikelihoods.
*/
double calcLnL(world_fmt *world, int scalingFactorIndex)
{
  beagle_fmt *beagle = world->beagle;
  const boolean scaling = beagle->scaling;
  int cumulativeScalingFactorIndex = scaling ? scalingFactorIndex : BEAGLE_OP_NONE;
  boolean accumulate_on_the_fly = scaling ? TRUE : FALSE;
  double logL = 0.0;
  //unsigned long i;
  unsigned long j;
  //unsigned long z;
  //unsigned long z;
  long locus = world->locus;
  const long snp_addon = 4;
  //long ii;
  //int max = 0;
  double *patternloglike = (double *) calloc(snp_addon+world->maxnumpattern[locus],sizeof(double));
  int rootIndex = beagle->operations[beagle->numoperations-1].destinationPartials;
  //  printf("Scalingfactorindex=%i\n",cumulativeScalingFactorIndex);
  long sublocus;
  const long sublocistart = world->sublocistarts[locus];
  const long sublociend   = world->sublocistarts[locus+1];
  for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
    {
      const long xs = sublocus - sublocistart;
      long instance = beagle->instance_handle[xs];
      mutationmodel_fmt *s = &world->mutationmodels[sublocus];
      const long numpatterns = s->addon+s->numpatterns;
      patternloglike = myrealloc(patternloglike,  numpatterns * sizeof(double));
      memset(patternloglike, 0, numpatterns * sizeof(double));
      int code = beagleUpdateTransitionMatrices(instance,	                	// instance,
						0,					// eigenIndex,
						(const int *) beagle->branch_indices,	// indicators transitionrates for each branch,
						NULL, 			                // firstDerivativeIndices,
						NULL,					// secondDervativeIndices,
						beagle->branch_lengths,		        // edgeLengths,
						beagle->numbranches);			// number branches to update, count

      if (code != 0)
	usererror("updateTransitionMatrices encountered a problem");
      
      if(beagle->scaling)
	{
	  beagleResetScaleFactors(instance,cumulativeScalingFactorIndex);

	  code = beagleUpdatePartials(instance,	// instance
				      beagle->operations,		                // operations
				      beagle->numoperations,				// operationCount
				      (accumulate_on_the_fly ?
				       cumulativeScalingFactorIndex : BEAGLE_OP_NONE)); // adjust scalers on the fly

	  //	  beagleResetScaleFactors(instance, cumulativeScalingFactorIndex);	

	  beagleAccumulateScaleFactors(instance,
				       beagle->scalingfactorsindices,
				       beagle->scalingfactorscount,
				       cumulativeScalingFactorIndex);
	}
      else
	{
	  code = beagleUpdatePartials(instance,                  	// instance
				      beagle->operations,		                // operations
				      beagle->numoperations,				// operationCount
				      BEAGLE_OP_NONE);                                // no scaling
	}
      if (code != 0)
	usererror("updatePartials with scaling encountered a problem");
      //      debug_beagle(beagle);
      double logLbeagle;
      int weightindex = 0;
      int stateFrequencyindex = 0;
      code = beagleCalculateRootLogLikelihoods(instance,                      // instance
					       (const int *) &rootIndex,      // bufferIndices
					       &weightindex,
					       &stateFrequencyindex,
					       &cumulativeScalingFactorIndex, //scalingfactors index,
					       1,                             // count is this correct
					       &logLbeagle);               // outLogLikelihoods
      code = beagleGetSiteLogLikelihoods(instance,                      // instance
					 patternloglike);               // outLogLikelihoods
      if (code != 0)
	usererror("calculateRootLogLikelihoods encountered a problem");
#ifdef BEAGLEDEBUG
      printf("after ----------------");
      debug_partials((int) rootIndex, s, instance, beagle, world);
#endif
#undef BEAGLEDEBUG
      if(s->addon>0)
	{
	  double maxinvariants = -HUGE;
	  double invariants    =  0.0;
	  const long added = s->numpatterns + s->addon;
	  for (j = s->numpatterns; j < added; j++) 
	    {
	      double pat = patternloglike[j];
	      if(pat > maxinvariants)
		{
		  maxinvariants = pat;
		}
	    }
	  for (j = s->numpatterns; j < added; j++) 
	    {
	      invariants += exp(patternloglike[j]-maxinvariants);
	    }
	  invariants = log(1.0 - invariants * exp(maxinvariants));
	  for (j = 0; j < s->numpatterns; j++) 
	    {
	      logL += s->aliasweight[j] * (patternloglike[j]-invariants);
	    }
	}
      else
	{
	  for (j = 0; j < s->numpatterns; j++) 
	    {
	      logL += s->aliasweight[j] * patternloglike[j];
	    }
	}
    }
  //  printf("%f\n",logL);
#ifdef BEAGLEDEBUG
  printf("DEBUG: Log LnL=%f\n",logL);
  debug_beagle(beagle);
#endif
  myfree(patternloglike);
  //  myfree(outlike);
  return logL; 
}


double force_beagle_recalculate(world_fmt *world, long locus)
{
  reset_beagle(world->beagle);
  set_all_dirty(world->root->next, crawlback (world->root->next), world, locus);
  smooth (world->root->next, crawlback (world->root->next), world, locus);
  return treelikelihood(world);
}

///
/// switches the final scaler when scaling is enabled
void change_beagle_scalingIndex(beagle_fmt *beagle)
{
  //if(beagle->scalingfactorsindices == beagle->scalingfactorsindices1)
  // {
      //beagle->scalingfactorsindices = beagle->scalingfactorsindices2;
  // }
  //else
  //{
      //beagle->scalingfactorsindices = beagle->scalingfactorsindices1;
  //}
  int tmp = beagle->scalingIndex1;
  beagle->scalingIndex1 =  beagle->scalingIndex2;
  beagle->scalingIndex2 = tmp;
}

///
/// start out at the root node and traverse up and match the ids with the new ids
void change_beagle(node *theNode, beagle_fmt *beagle, long sumtips)
{
  // move node ids;
  // this works only just after the pseudolikelihood step + the acceptlike with TRUE outcome
  long i;
  //long tmp;
  long parentid;
  long oldparentid;
  long nodep_boundary = sumtips * 2 - 1;
  if (theNode->type != 't')
    {
      change_beagle(crawlback(theNode->next), beagle,sumtips);
      change_beagle(crawlback(theNode->next->next),beagle, sumtips);
    }
  for(i=0;i<beagle->numoperations;i++)
    {
      parentid = beagle->operations[i].destinationPartials;
      oldparentid = parentid < nodep_boundary ? parentid + nodep_boundary : parentid - nodep_boundary; 
      if(theNode->id == oldparentid)
	{
	  theNode->id = parentid;
	}
    }
}

long new_id(long id, long sumtips)
{
  const long nodep_boundary = 2 * sumtips - 1;
  return (id < nodep_boundary ? id + nodep_boundary : id - nodep_boundary); 
}

void reset_beagle(beagle_fmt *beagle)
{
  beagle->numoperations = 0;
  beagle->numbranches   = 0;
  beagle->scalingfactorscount = 0;
}




void beagle_stop(world_fmt **universe, long usize)
{
  long u;
  for(u=0;u< usize; u++)
    {
      world_fmt *world = universe[u];
      beagle_fmt *beagle = world->beagle;
      unsigned long i;
      long locus = world->locus;
      for(i=0; i<world->numsubloci[locus]; i++)
	{
	  beagleFinalizeInstance(beagle->instance_handle[i]);
	}
    }
}


void set_beagle_dirty(node *origin, node *target, node *mrca)
{
  node *nn = origin;
  while((nn=showtop(crawlback(nn)))!=mrca)
    {
      set_dirty(nn);
    }
  nn = target;
  while((nn=showtop(crawlback(nn)))!=mrca)
    {
      set_dirty(nn);
    }
}

void debug_beagle(beagle_fmt *beagle)
{
  long i;
  printf("----beagle content---------------------\n");
  printf("operations    : %p\n",beagle->operations);
  printf("     alloc    : %i\n",beagle->numallocoperations);
  printf("       num    : %i\n",beagle->numoperations);
  for(i=0;i<beagle->numoperations;i++)
    printf("%i ",beagle->operations[i].destinationPartials);
  printf("\nbranch indices: %p\n",beagle->branch_indices);
  printf("     alloc    : %i\n",beagle->numallocoperations * 2);
  printf("       num    : %i\n",beagle->numbranches);
  printf("----beagle content end-----------------\n\n");
}

//BEAGLE_DLLEXPORT int beagleGetPartials(int instance,
//int bufferIndex,
//int scaleIndex,
//double* outPartials);
#ifdef BEAGLEDEBUG
void debug_partials(int nodenum, mutationmodel_fmt *s, long instance,beagle_fmt *beagle, world_fmt *world)
 {
   int i;
   double *outlike = (double *) calloc((s->numpatterns+s->addon)*s->numsiterates*s->numstates, sizeof(double));
   beagleGetPartials(instance,nodenum,BEAGLE_OP_NONE, outlike);
   printf("-----------------------\n");
   long z=0;
   for(i=0;i < s->numpatterns; i++)
     {
       printf("%i: {%g, %g, %g, %g}\n",i, outlike[z],outlike[z+1],outlike[z+2],outlike[z+3]);
       if(outlike[z]>1.0 || outlike[z+1] > 1.0 || outlike[z+2] > 1.0 || outlike[z+3] > 1.0)
	 error("failed");
	    z += 4;
     }
   myfree(outlike);
 }

void debug_eigensystem(long i, beagle_fmt *beagle, world_fmt *world)
 {
   long nn = world->mutationmodels[i].numstates * world->mutationmodels[i].numstates;
   int j;
   printf("Eigenvalues:");
   for(j=0;j<nn; j++)
	  {
	    printf(" %f",world->mutationmodels[i].eigenvalues[j]);
	  }  
   printf("\n\n");
 }
#endif /*BEAGLEDEBUG*/
#else
// we do not use BEAGLE
typedef double t_is_empty;
#endif /*BEAGLE*/
