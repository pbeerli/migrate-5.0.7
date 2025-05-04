/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 M C M C 2   R O U T I N E S 
 
 Tree changing routines
 
 
 Peter Beerli 1996, Seattle
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2009 Peter Beerli, Tallahassee FL
 
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
 
$Id: mcmc2.c 2158 2013-04-29 01:56:20Z beerli $
-------------------------------------------------------*/
/* \file mcmc2.c
Contains low-level routines to manipulate the treechanger
*/

#include "migration.h"
#include "sighandler.h"
#include "tree.h"
#include "uep.h"
#include "tools.h"
#include "mcmc.h"
#include "migrate_mpi.h"
#ifdef BEAGLE
#include "calculator.h"
#endif
#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif

#define comment(a) FPRINTF(stderr,a)


/* prototypes ------------------------------------------- */
/* changes the tree to the proposed tree by rebuilding parts
   of the tree and also inserts new migration events. */
void coalesce1p (proposal_fmt * proposal);
/* and its pretend version */
void pretendcoalesce1p (proposal_fmt * proposal);

/* private functions */
/* family of tree manipulation routines: dependent on the position of
   target and origin in the tree different routines are used.
   oback: there is no change in the tree topology.
   ocousin: target is the cousin-node, the branch is inserted in its
   former sister branch.
   rbcoa: target is the bottommost node, the ripped branch will be 
   bring the new bottommost node.
   stancoa: all other cases. */
void coat_oback (proposal_fmt * proposal);
void coat_ocousin (proposal_fmt * proposal);
void coat_obbcoa (proposal_fmt * proposal);
void coat_rbcoa (proposal_fmt * proposal);
void coat_stancoa (proposal_fmt * proposal);
/* pretend versions of the above */
void target_oback (proposal_fmt * proposal);
void target_obbcoa (proposal_fmt * proposal);
void target_ocousin (proposal_fmt * proposal);
void target_rbcoa (proposal_fmt * proposal);
void target_stancoa (proposal_fmt * proposal);
/* subroutines in target_stancoa */
void t_s_upper (proposal_fmt * proposal, node * connect, node * obb);
node *t_s_obbncon (proposal_fmt * proposal, node * obb);
void t_s_tncon (proposal_fmt * proposal, node * obb);
void t_s_tcon (proposal_fmt * proposal, node * uppern);
long erase_migr_nodes (world_fmt *world, node * up);
node *findcrossing (node ** ptrl1, node ** ptrl2);
void connectnodes (world_fmt *world, node * mother, node * brother, node * sister);
void gotoroot (node * origin, node ** ptrlist, long *ptrallocsize);
void adjust (node * theNode, MYREAL tyme, long level);
void localevaluate (node * mother);
void copy_x (proposal_fmt * proposal, xarray_fmt *axx1, xarray_fmt *axx2);
void fix_root_pop (world_fmt *world, node * p);
void pseudoevaluate (proposal_fmt * proposal, xarray_fmt *x, MYREAL **lx,
                     node * mother, node * newdaughter, MYREAL v);
node *crawl_down (node * theNode, MYREAL tyme);
void scalecopy(MYREAL **copy, MYREAL **original, world_fmt *world);

void old_target_stancoa (proposal_fmt * proposal);
void erase_migr_nodes2 (world_fmt *world, node * up);

//##
/*=========================================================*/
void scalecopy(MYREAL **copy, MYREAL **original, world_fmt *world)
{
  const long locus        = world->locus;
  const long sublocistart = world->sublocistarts[locus];
  const long sublociend   = world->sublocistarts[locus+1];
  long sublocus;
  for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
    {
      mutationmodel_fmt *s = &world->mutationmodels[sublocus];
      const long xs = sublocus - sublocistart;
      memcpy(copy[xs],original[xs], sizeof(MYREAL) * (size_t) (s->numpatterns + s->addon));
    }
}


void
coalesce1p (proposal_fmt * proposal)
{
#ifdef MESS
  int skip=1;
#endif
  //static long count=0;
    node *obb = NULL;
    erase_migr_nodes (proposal->world, proposal->origin);
    if (proposal->migr_table_counter2 > 0) //|| proposal->mig_removed)
      erase_migr_nodes (proposal->world, proposal->realtarget);
    if (proposal->target == proposal->ocousin)
    {
        coat_ocousin (proposal);
        if (proposal->oback->tyme > showtop (proposal->oback->back)->tyme
                && crawlback (proposal->oback)->type != 'r')
            error ("problem with time encountered"); /* was >= ** */
    }
    else
    {
        if (((proposal->target == proposal->oback)
                || (proposal->target == proposal->osister)
                || (proposal->target == proposal->origin)))
        {
            //xcode obb = proposal->oback->back;
            coat_oback (proposal);
            if (proposal->oback->tyme > showtop (proposal->oback->back)->tyme)
            {   /* was >= ** */
                if (proposal->oback->back->type == 'r')
                {
                    adjust (proposal->root, proposal->oback->tyme + 10000, 0L);
                    comment ("Root time adjusted: was to small");
                }
                else
                    error ("problem with time encountered");
            }
        }
        else
        {
            if (crawlback (proposal->target)->type == 'r')
            {
                coat_rbcoa (proposal);
                if (proposal->oback->tyme >
                        showtop (proposal->oback->back)->tyme)
                {  /* was >= ** */
                    if (proposal->oback->back->type == 'r')
                    {
                        adjust (proposal->root, proposal->oback->tyme + 10000,
                                0L);
                        comment ("Root time adjusted: was to small");
                    }
                    else
                        error ("problem with time encountered");
                }
            }
            else
            {
                obb = showtop (crawlback (proposal->oback));
                if ((obb != NULL) && (obb == proposal->target))
                {
                    coat_obbcoa (proposal);
                    if (proposal->oback->tyme > showtop (proposal->oback->back)->tyme) /* was >= ** */
                        error ("problem with time encountered");
                }
                else
                {
                    coat_stancoa (proposal);
                    if (proposal->oback->tyme > showtop (proposal->oback->back)->tyme) /* was >= ** */
                        error ("problem with time encountered");
                }
            }
        }
    }
    set_pop (proposal->oback, proposal->origin->actualpop,
             proposal->origin->actualpop);
  //  free_migr_node(proposal->origin, crawlback(proposal->origin));
    //    if(++count % 10000 == 0)
    //  {
    //	printf("%li heat %f> %li %li %li %li\n", count, proposal->world->heat, proposal->migr_table_counter, proposal->migr_table_counter2, proposal->world->node_collection_count, proposal->world->node_collection_allocated);
    //
    //  }

    if (proposal->migr_table_counter > 0)
    {
        set_pop (proposal->oback,
                 proposal->migr_table[proposal->migr_table_counter - 1].from,
                 proposal->migr_table[proposal->migr_table_counter - 1].from);
        insert_migr_node (proposal->world, proposal->origin,
                          crawlback (proposal->origin), proposal->migr_table,
                          &(proposal->migr_table_counter));
    }
    proposal->migr_table_counter = 0;
    if (proposal->migr_table_counter2 > 0)
    {
        insert_migr_node (proposal->world, proposal->realtarget,
                          crawlback (proposal->target), proposal->migr_table2,
                          &(proposal->migr_table_counter2));
        proposal->migr_table_counter2 = 0;
    }
    /*    if (proposal->oback == crawlback(proposal->root->next)) */
    fix_root_pop (proposal->world, crawlback (proposal->root->next));
    //    if(proposal->world->in_last_chain)
    traverse_tagnew(crawlback (proposal->root->next),proposal->origin);
    traverse_check(crawlback (proposal->root->next));
#ifdef MESS
    //traverse_checker(proposal->root);
    if(skip==0)
      {
	FPRINTF (stdout,
		 "\n[& Comment: Locus %li, best log likelihood = %f (to file),heat=%f]\n",
		 proposal->world->locus + 1, proposal->world->likelihood[proposal->world->G], proposal->world->heat);
	
	debugtreeout (stdout, crawlback (proposal->world->root->next),
		      crawlback (proposal->world->root->next), 0);
      }
#endif
}



void
pretendcoalesce1p (proposal_fmt * proposal)
{
    if (proposal->target == proposal->ocousin)
    {
        target_ocousin (proposal);
    }
    else
    {
        if (((proposal->target == proposal->oback)
                || (proposal->target == proposal->osister)
                || (proposal->target == proposal->origin)))
        {
            target_oback (proposal);
        }
        else
        {
            if (crawlback (proposal->target)->type == 'r')
            {
                target_rbcoa (proposal);
            }
            else
            {
                if ((showtop (crawlback (proposal->oback)) != NULL)
                        && (showtop (crawlback (proposal->oback)) ==
                            proposal->target))
                {
                    target_obbcoa (proposal);
                }
                else
                {
                    target_stancoa (proposal);
                }
            }
        }
    }
}

/*================================================================*/

void
coat_oback (proposal_fmt * proposal)
{
    node *obb, *obm;
#ifdef TREEDEBUG
    printf("coat_oback\n");
#endif
    obb = showtop (crawlback (proposal->oback));
    connectnodes (proposal->world,showtop (proposal->oback->back), proposal->osister, NULL);
    adjust_time (proposal->oback, proposal->time);
    obm = showtop (crawl_down (proposal->osister, proposal->time)->back);
    connectnodes (proposal->world, proposal->oback, proposal->origin, proposal->osister);
    proposal->oback->back = NULL;
    if (obm->type == 'm' || obm->type == 'd')
      {
	connectnodes (proposal->world,obm, proposal->oback, NULL);
      }
    
    if (obb->type == 'r')
      {
        connectnodes (proposal->world,obb, proposal->oback, NULL);
	erase_migr_nodes(proposal->world,proposal->oback);
      }
    else
    {
      connectnodes (proposal->world,obb, proposal->oback, proposal->ocousin);
    }
    adjust (obb, obb->tyme, 2L);
    set_dirty (proposal->oback);
    localevaluate (proposal->oback);
}

void
coat_ocousin (proposal_fmt * proposal)
{
    node *obb, *rcback;
#ifdef TREEDEBUG
    printf("coat_ocousin\n");
#endif

    obb = showtop (crawlback (proposal->oback));
    connectnodes (proposal->world,showtop (proposal->oback->back), proposal->osister, NULL);
    rcback = showtop (crawl_down (proposal->ocousin, proposal->time)->back);
    adjust_time (proposal->oback, proposal->time);
    proposal->oback->back = NULL;
    connectnodes (proposal->world,rcback, proposal->oback, NULL);
    connectnodes (proposal->world,proposal->oback, proposal->ocousin, proposal->origin);
    connectnodes (proposal->world,obb, proposal->oback, proposal->osister);
    adjust (obb, obb->tyme, 2L);
    set_dirty (proposal->oback);
    set_dirty (obb);
    if (crawlback (obb)->type != 'r')
        localevaluate (showtop (crawlback (obb)));
}

void
coat_rbcoa (proposal_fmt * proposal)
{
    node *obb, *root;
#ifdef TREEDEBUG
    printf("coat_rbcoa\n");
#endif

    if (proposal->ocousin != NULL)
    {
        obb = showtop (crawlback (proposal->oback));
        root = proposal->root;
        connectnodes (proposal->world,showtop (proposal->oback->back), proposal->osister, NULL);
        connectnodes (proposal->world,obb, proposal->ocousin, proposal->osister);
        adjust (obb, obb->tyme, 2L);
        adjust_time (proposal->oback, proposal->time);
        proposal->oback->back = NULL;
        connectnodes (proposal->world,proposal->oback, proposal->origin, proposal->target);
        connectnodes (proposal->world,root, proposal->oback, NULL);
        adjust (proposal->oback, proposal->time, 2L);
        set_dirty (obb);
        set_dirty (proposal->oback);
        localevaluate (showtop (crawlback (obb)));
    }
    else
    {
        error ("error in coat_rbcoa\n");
    }
}

void
coat_obbcoa (proposal_fmt * proposal)
{
  char x;
  node /**obb,*/  * proposalack;
#ifdef TREEDEBUG
    printf("coat_obbcoa\n");
#endif

    if (proposal->ocousin != NULL)
    {
        /*   obb = showtop(crawlback(proposal->oback)); */
        connectnodes (proposal->world,showtop (proposal->oback->back), proposal->osister, NULL);
        // proposalack = showtop (crawl_down (proposal->target, proposal->time)->back);
	       proposalack = showtop (crawlback (proposal->target));
        adjust_time (proposal->oback, proposal->time);
        connectnodes (proposal->world,proposal->target, proposal->ocousin, proposal->osister);
        adjust (proposal->target, proposal->target->tyme, 1L);
        proposal->oback->back = NULL;
        adjust_time (proposal->oback, proposal->time);
	x = showtop (proposal->realtarget->back)->type;
        if (x  == 'm' || x == 'd')
        {
            connectnodes (proposal->world,showtop (proposal->realtarget->back), proposal->oback,
                          NULL);
        }
        connectnodes (proposal->world,proposal->oback, proposal->origin, proposal->target);
        connectnodes (proposal->world,proposalack, proposal->tsister, proposal->oback);
        adjust (proposalack, proposalack->tyme, 2L);
        set_dirty (proposal->target);
        set_dirty (proposal->oback);
        localevaluate (proposalack);
    }
    else
    {
        error ("error in coat_obbcoa\n");
    }
}

void
coat_stancoa (proposal_fmt * proposal)
{
    node *obb, *proposalack;
#ifdef TREEDEBUG
    printf("coat_stancoa\n");
#endif

    if (proposal->ocousin != NULL)
    {
        obb = showtop (crawlback (proposal->oback));
        //proposalack = showtop (crawl_down (proposal->target, proposal->time)->back);
	proposalack = showtop (crawlback (proposal->target));
        connectnodes (proposal->world,showtop (proposal->oback->back), proposal->osister, NULL);
        proposal->oback->back = NULL;
        connectnodes (proposal->world,obb, proposal->osister, proposal->ocousin);
        adjust_time (proposal->oback, proposal->time);
	char x = showtop (proposal->realtarget->back)->type;
        if (x  == 'm' || x == 'd')
	  {
            connectnodes (proposal->world,showtop (proposal->realtarget->back), proposal->oback,
                          NULL);
	  }
        connectnodes (proposal->world,proposal->oback, proposal->origin, proposal->target);
        connectnodes (proposal->world,proposalack, proposal->oback, proposal->tsister);
        adjust (proposalack, proposalack->tyme, 3L);
        adjust (obb, obb->tyme, 3L);
        set_dirty (obb);
        set_dirty (proposal->oback);
        localevaluate (proposalack);
        localevaluate (obb);
    }
    else
    {
        adjust_time (proposal->oback, proposal->time);
        obb = showtop (crawlback (proposal->oback));
        connectnodes (proposal->world,showtop (proposal->oback->back), proposal->osister, NULL);
        //proposalack = showtop (crawl_down (proposal->target, proposal->time)->back);
	       proposalack = showtop (crawlback (proposal->target));
        adjust_time (proposal->oback, proposal->time);
        proposal->oback->back = NULL;
	char x = showtop (proposal->realtarget->back)->type;
        if (x == 'm' || x == 'd')
            connectnodes (proposal->world,showtop (proposal->realtarget->back), proposal->oback,
                          NULL);
        connectnodes (proposal->world,proposal->oback, proposal->origin, proposal->target);
        connectnodes (proposal->world,proposalack, proposal->oback, proposal->tsister);
        /* either proposalack is osister or a descendent of her
           other case are disallowed by time constraints */
        connectnodes (proposal->world,obb, proposal->osister, NULL);
        adjust (proposalack, proposalack->tyme, 3L);
        adjust (obb, obb->tyme, 3L);
        set_dirty (proposal->oback);
        localevaluate (proposalack);

    }
}

///
/// pretends to connect the nodes so that the new tree looks like
/// insert the line on the same branch osister -- obb
/// therefore the target is oback
/// @verbatim
/// ori     osis
///  |       |
///  +-oback-+   <--- target
///      |
///      |
/// @endverbatim
void target_oback (proposal_fmt * proposal)
{
    node *obb = showtop (crawlback (proposal->oback));
#ifdef TREEDEBUG
    printf("target_oback\n");
#endif
#ifndef BEAGLE
    copy_x (proposal, proposal->xf, proposal->origin->x);

#ifdef UEP

    if(proposal->world->options->uep)
        copy_uepx(proposal,proposal->uf, proposal->origin->ux);
#endif

    scalecopy (proposal->mf, proposal->origin->scale, proposal->world);
#endif
    proposal->v = (proposal->time - proposal->origin->tyme); // time from origin to new oback
    proposal->vs = (proposal->time - proposal->osister->tyme); // time adjusted for time from osister to oback
#ifdef BEAGLE
    long id0, id1, id2;
    long bid1, bid2;
    double v1, v2;
    id0  = proposal->oback->id;
    id1  = proposal->origin->id;
    bid1 = proposal->origin->bid;
    v1   = proposal->v;
    id2  = proposal->osister->id;
    bid2 = proposal->osister->bid;
    v2   = proposal->vs;
    prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
    evaluate_beagle_instances_proposal (proposal,obb, proposal->oback, id0, proposal->oback->bid,
					fabs (proposal->time - obb->tyme));
#else
    pseudonuview (proposal, proposal->xf, proposal->mf, proposal->v,
                  proposal->osister->x, proposal->osister->scale, proposal->vs);
    pseudoevaluate (proposal, proposal->xf, proposal->mf, obb, proposal->oback,
                    fabs (proposal->time - obb->tyme));
#endif
}

///
/// prepares or calculates the conditional likelihoods when the new tree looks like this
/// @verbatim
///   ori   ocousin
///    |       |
///    +-oback-+     osis
///       |           |
///       +----obb----+
///             |
/// @endverbatim
void target_ocousin (proposal_fmt * proposal)
{
    node *obb;
#ifdef TREEDEBUG
    printf("target_ocousin\n");
#endif
    obb = showtop (crawlback (proposal->oback));
#ifndef BEAGLE
    copy_x (proposal, proposal->xf, proposal->origin->x);
#ifdef UEP

    if(proposal->world->options->uep)
        copy_uepx(proposal,proposal->uf, proposal->origin->ux);
#endif
    scalecopy (proposal->mf, proposal->origin->scale, proposal->world);
    copy_x (proposal, proposal->xt, proposal->ocousin->x);
#ifdef UEP
    if(proposal->world->options->uep)
        copy_uepx(proposal,proposal->ut, proposal->ocousin->ux);
#endif
    scalecopy (proposal->mt, proposal->ocousin->scale, proposal->world);
#endif /*if not define BEAGLE*/
    proposal->v = (proposal->time - proposal->origin->tyme);// time to oback
    proposal->vs = fabs (proposal->ocousin->tyme - proposal->time); //time from ocousin to new oback
#ifdef BEAGLE
    long id0, id1, id2;
    long bid1, bid2;
    double v1, v2;
    id0  = proposal->oback->id;
    id1  = proposal->origin->id;
    bid1 = proposal->origin->bid;
    v1   = proposal->v;
    id2  = proposal->ocousin->id;
    bid2 = proposal->ocousin->bid;
    v2   = proposal->vs;
    prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
#else
    pseudonuview (proposal, proposal->xf, proposal->mf, proposal->v,
                  proposal->xt, proposal->mt, proposal->vs);
#endif
    proposal->v = fabs (proposal->time - obb->tyme); //time from new oback to obb
    proposal->vs = obb->tyme - proposal->osister->tyme; // time from obb to osister (this has also changed 
#ifdef BEAGLE
    id0  = obb->id;
    id1  = new_id(proposal->oback->id,proposal->world->sumtips); // connect to the oback/ori/ocousin trio
    bid1 = proposal->oback->bid;
    v1   = proposal->v;
    id2  = proposal->osister->id;
    bid2 = proposal->osister->bid;
    v2   = proposal->vs;
    prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
    // connects to the root from obb
    evaluate_beagle_instances_proposal (proposal,
					showtop (crawlback (obb)), obb, obb->id, obb->bid,
					(showtop (crawlback (obb))->tyme - obb->tyme));

#else
    pseudonuview (proposal, proposal->xf, proposal->mf, proposal->v,
                  proposal->osister->x, proposal->osister->scale, proposal->vs);
    pseudoevaluate (proposal, proposal->xf, proposal->mf,
                    showtop (crawlback (obb)), obb,
                    (showtop (crawlback (obb))->tyme - obb->tyme));
#endif
}

///
/// calculates the conditional likelihoods when new tree looks like
/// @verbatim
/// ori      1 (target)    osis     ocous     
///  |       |              |       |       
///  +-oback-+              +--obb--+
///      |                      |
///                            to root (which is 1=target)
/// @endverbatim
void target_rbcoa (proposal_fmt * proposal)
{
    node *obb;
#ifdef TREEDEBUG
    printf("target_rbcoa\n");
#endif
    if (proposal->ocousin != NULL)
    {
      obb = showtop (crawlback (proposal->oback));
      proposal->vs = obb->tyme - proposal->osister->tyme;
      proposal->v = obb->tyme - proposal->ocousin->tyme;
#ifndef BEAGLE
      copy_x (proposal, proposal->xt, proposal->osister->x);
#ifdef UEP
      
      if(proposal->world->options->uep)
	copy_uepx(proposal,proposal->ut, proposal->osister->ux);
#endif
      
      scalecopy (proposal->mt, proposal->osister->scale, proposal->world);	     
#endif /*if not BEAGLE*/
#ifdef BEAGLE
      long id0, id1, id2;
      long bid1, bid2;
      double v1, v2;
      id0  = obb->id;
      id1  = proposal->osister->id;
      bid1 = proposal->osister->bid;
      v1   = proposal->vs;
      id2  = proposal->ocousin->id;
      bid2 = proposal->ocousin->bid;
      v2   = proposal->v;
      prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
      //  evaluate to the root
     evaluate_beagle_instances_proposal(proposal, 
					 showtop (crawlback (obb)), obb, obb->id, obb->bid,
					 showtop (crawlback (obb))->tyme - obb->tyme);
#else
      pseudonuview (proposal, proposal->xt, proposal->mt, proposal->vs,
		    proposal->ocousin->x, proposal->ocousin->scale,
		    proposal->v);
      pseudoevaluate (proposal, proposal->xt, proposal->mt,
		      showtop (crawlback (obb)), obb,
		      showtop (crawlback (obb))->tyme - obb->tyme);
#endif
      proposal->v = proposal->time - proposal->origin->tyme;
      proposal->vs = proposal->time - proposal->target->tyme;
#ifndef BEAGLE
      copy_x (proposal, proposal->xf, proposal->origin->x);
#ifdef UEP
      
      if(proposal->world->options->uep)
	copy_uepx(proposal,proposal->uf, proposal->origin->ux);
#endif
      
      scalecopy (proposal->mf, proposal->origin->scale, proposal->world);

#endif /*if not BEAGLE*/
#ifdef BEAGLE
      id0  = proposal->oback->id;
      id1  = proposal->origin->id;
      bid1 = proposal->origin->bid;
      v1   = proposal->v;
      id2  = new_id(proposal->target->id, proposal->world->sumtips);//target == root->next->back
      bid2 = proposal->target->bid;
      v2   = proposal->vs;
      prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
#else
      pseudonuview (proposal, proposal->xf, proposal->mf, proposal->v,
		    proposal->xt, proposal->mt, proposal->vs);
#endif
    }
    else
      {
        error ("error in target_rbcoa\n");
      }
}

///
/// calculates conditional likelihood for this
/// @verbatim
///    osis     ocous
///      |       |   
///      +--obb--+	
/// ori      |                            
///  |       |               
///  +-oback-+           
///      |               
///     obbb                           
/// @endverbatim
void target_obbcoa (proposal_fmt * proposal)
{
  node *obb, *obbb;
#ifdef TREEDEBUG
  printf("target_obbcoa\n");
#endif
  if (proposal->ocousin != NULL)
    {
      obb = showtop (crawlback (proposal->oback));
      proposal->vs = obb->tyme - proposal->osister->tyme;
      proposal->v = obb->tyme - proposal->ocousin->tyme;
#ifndef BEAGLE
      copy_x (proposal, proposal->xt, proposal->osister->x);
#ifdef UEP
      if(proposal->world->options->uep)
	copy_uepx(proposal,proposal->ut, proposal->osister->ux);
#endif
      scalecopy (proposal->mt, proposal->osister->scale, proposal->world);
#endif /*if not BEAGLE*/
#ifdef BEAGLE
    long id0, id1, id2;
    long bid1, bid2;
    double v1, v2;
    id0  = obb->id;
    id1  = proposal->osister->id;
    bid1 = proposal->osister->bid;
    v1   = proposal->vs;
    id2  = proposal->ocousin->id;
    bid2 = proposal->ocousin->bid;
    v2   = proposal->v;
    prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
#else
        pseudonuview (proposal, proposal->xt, proposal->mt, proposal->vs,
                      proposal->ocousin->x, proposal->ocousin->scale,
                      proposal->v);
#endif
        proposal->v = proposal->time - proposal->origin->tyme;
        proposal->vs = fabs (obb->tyme - proposal->time);
#ifndef BEAGLE
        copy_x (proposal, proposal->xf, proposal->origin->x);
#ifdef UEP
        if(proposal->world->options->uep)
            copy_uepx(proposal,proposal->uf, proposal->origin->ux);
#endif
	scalecopy (proposal->mf, proposal->origin->scale, proposal->world);
#endif /*if not BEAGLE*/
#ifdef BEAGLE
    id0  = proposal->oback->id;
    id1  = proposal->origin->id;
    bid1 = proposal->origin->bid;
    v1   = proposal->v;
    id2  = new_id(obb->id, proposal->world->sumtips);
    bid2 = obb->bid;
    v2   = proposal->vs;
    prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
    obbb = showtop (crawlback (obb));
    proposal->v = fabs (obbb->tyme - proposal->time);
    /*check this carefully: fakes obb and uses id of oback which is inserted between obb and obbb*/
    /*this should guarantee the proper evaluation but feeding the right ids into beagle [091509 this seems OK]*/
    evaluate_beagle_instances_proposal (proposal,obbb, obb, proposal->oback->id, proposal->oback->bid,
					proposal->v);
#else
        pseudonuview (proposal, proposal->xf, proposal->mf, proposal->v,
                      proposal->xt, proposal->mt, proposal->vs);
        obbb = showtop (crawlback (obb));
        proposal->v = fabs (obbb->tyme - proposal->time);
        pseudoevaluate (proposal, proposal->xf, proposal->mf, obbb, /* is this correct??????*/ obb,
                        proposal->v);
#endif
    }
    else
    {
        error ("error in target_obbcoa\n");
    }
}

///
/// calculates all standard cases
/// @verbatim
/// ori    target         osis     ocous    ori    target                       ori    target        		    
///  |       |              |       |        |       |              	         |       |              		    
///  +-oback-+              +--obb--+	     +-oback-+              	         +-oback-+              		    
///      |                      |	         |                  	             |                  		    
///     ...                    ...	        ...			            ...				    
///      |                      |	         |			             |				    
///      +--------connect-------+	        osis                  ocous         ocous                  osis	    	    
///                 |                            |                      |            |                      |	    	    
///					         +---------obb----------+            +---------obb----------+	    	    
///					                    |                                   |                             
///                                                     t_s_upper()                           t_s_upper()  
/// @endverbatim
void target_stancoa (proposal_fmt * proposal)
{
    node *obb, *nn, *oldnn, *d1, *d2, *proposalack, *uppern = NULL; /* o R */
    node **double_line;
#ifdef TREEDEBUG
    printf("target_stancoa\n");
#endif
    double_line = (node **) mycalloc (2, sizeof (node *));
    if (proposal->ocousin != NULL)
      {
        obb = showtop (crawlback (proposal->oback));
        gotoroot (proposal->target, proposal->line_t, &proposal->line_t_allocsize); // looks at target --> .... ---> root
        double_line[0] = proposal->osister;          
        /*findcrossing needs a last NULL element */
        if (findcrossing (double_line, proposal->line_t) != NULL)
	  {
	    //printf("osister ");
            t_s_upper (proposal, proposal->osister, obb); // osister is one of the ancestors of target
#ifdef BEAGLEDEBUG
	    debug_beagle(proposal->world->beagle);
#endif
	  }
        else
	  {
            double_line[0] = proposal->ocousin;
            /*findcrossing needs a last NULL element */
            if (findcrossing (double_line, proposal->line_t) != NULL)
	      {
		//		printf("ocousin ");
                t_s_upper (proposal, proposal->ocousin, obb); //ocousin is one of the ancestors of target
#ifdef BEAGLEDEBUG
		debug_beagle(proposal->world->beagle);
#endif
	      }
            else
	      {
                gotoroot (obb, proposal->line_f, &proposal->line_f_allocsize); //goto the root from obb, the node that is below origin and stable
                proposal->connect = findcrossing (proposal->line_f, proposal->line_t); // node that is ancestor of obb and target
                proposal->vs = obb->tyme - proposal->osister->tyme;
                proposal->v = obb->tyme - proposal->ocousin->tyme;
#ifndef BEAGLE
                copy_x (proposal, proposal->xt, proposal->osister->x);
#ifdef UEP
		
                if(proposal->world->options->uep)
		  copy_uepx(proposal,proposal->ut, proposal->osister->ux);
#endif
		
                scalecopy (proposal->mt, proposal->osister->scale, proposal->world);
#endif /*if not BEAGLE*/
		if (obb!=proposal->connect) // handles conditionals down to connect@@
		  { 
		    uppern = t_s_obbncon (proposal, obb);//returns the last node before connect, from the obb-side
		  }
		else
		  {
		    // this should not happen because target_oback() should cover this case
		    warning("this should not happen because this is covered by target_oback()");
		    uppern = proposal->connect;
		  }
                proposal->v = proposal->time - proposal->origin->tyme;
#ifndef BEAGLE
                copy_x (proposal, proposal->xf, proposal->origin->x);
#ifdef UEP
		
                if(proposal->world->options->uep)
		  copy_uepx(proposal,proposal->uf, proposal->origin->ux);
#endif
		
                scalecopy (proposal->mf, proposal->origin->scale, proposal->world);
#endif /*if not BEAGLE*/
                if (proposal->target != proposal->connect)
		  {
		    //                    t_s_tncon (proposal, obb);
		    t_s_tncon (proposal, uppern);//target is not connect
		  }
                else
		  {
		    //warning("t_s_tcon() should not be called, right?");
                    t_s_tcon (proposal, uppern); // target is connect, on the old tree obb-->target<--...
		    // the oback will connect below target and target needs also to be updated 
		  }
	      }
	  }
      }
    else
      { // ocousin is NULL, only for two tip trees?
        proposal->v = proposal->time - proposal->origin->tyme;
        proposal->vs = proposal->time - proposal->target->tyme;
#ifndef BEAGLE
        copy_x (proposal, proposal->xf, proposal->origin->x);
#ifdef UEP
	
        if(proposal->world->options->uep)
	  copy_uepx(proposal,proposal->uf, proposal->origin->ux);
#endif
	
        scalecopy (proposal->mf, proposal->origin->scale, proposal->world);
        copy_x (proposal, proposal->xt, proposal->target->x);
#ifdef UEP
	
        if(proposal->world->options->uep)
	  copy_uepx(proposal,proposal->ut, proposal->target->ux);
#endif	
        scalecopy (proposal->mt, proposal->target->scale,proposal->world);
        //CHECK HERE      memcpy(proposal->mt,proposal->osister->scale, sizeof(MYREAL)*proposal->endsite);
#endif /*is not BEAGLE*/
#ifdef BEAGLE
	long id0, id1, id2;
	long bid1, bid2;
	double v1, v2;
	id0  = proposal->oback->id;
	id1  = proposal->origin->id;
	bid1 = proposal->origin->bid;
	v1   = proposal->v;
	id2  = proposal->target->id;
	bid2 = proposal->target->bid;
	v2   = proposal->vs;
	prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
#else
        pseudonuview (proposal, proposal->xf, proposal->mf, proposal->v,
                      proposal->xt, proposal->mt, proposal->vs);
#endif
	proposalack = showtop (crawlback (proposal->target));
        proposal->v = fabs (proposalack->tyme - proposal->time);
        if (proposalack != proposal->oback)
	  {
            children (proposalack, &d1, &d2);
            if (d1 == proposal->target)
	      {
		d1 = d2;
		d2 = proposal->target;
	      }
#ifdef BEAGLE
	    id0  = proposalack->id;
	    id1  = new_id(proposal->oback->id,proposal->world->sumtips);
	    bid1 = proposal->oback->bid;
	    v1   = proposal->v;
	    id2  = d1->id;
	    bid2 = d1->bid;
	    v2   = d1->v;
	    prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
#else
            pseudonuview (proposal, proposal->xf, proposal->mf, proposal->v,
                          d1->x, d1->scale, d1->v);
#endif
            oldnn = proposalack;
            nn = showtop (crawlback (oldnn));
            while (nn != proposal->oback)
	      {
		//warning("why are we here?");// because oback is the root-back node
                children (nn, &d1, &d2);
                if (d1 == oldnn)
		  {
		    d1 = d2;
		    d2 = oldnn;
		  }
#ifdef BEAGLE
		id0  = nn->id;
		id1  = new_id(oldnn->id,proposal->world->sumtips);
		bid1 = oldnn->bid;
		v1   = oldnn->v;
		id2  = d1->id;
		bid2 = d1->bid;
		v2   = d1->v;
		prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
#else
                pseudonuview (proposal, proposal->xf, proposal->mf, oldnn->v,
                              d1->x, d1->scale, d1->v);
#endif
                oldnn = nn;
                nn = showtop (crawlback (nn));
	      }
	  }
      }
    myfree(double_line);
}

/// safe copy
void old_target_stancoa (proposal_fmt * proposal)
{
    node *obb, *nn, *oldnn, *d1, *d2, *proposalack, *uppern = NULL; /* o R */
    node **double_line;
#ifdef TREEDEBUG
    printf("target_stancoa\n");
#endif
    double_line = (node **) mycalloc (2, sizeof (node *));
    if (proposal->ocousin != NULL)
      {
        obb = showtop (crawlback (proposal->oback));
        gotoroot (proposal->target, proposal->line_t, &proposal->line_t_allocsize);
        double_line[0] = proposal->osister;
        /*findcrossing needs a last NULL element */
        if (findcrossing (double_line, proposal->line_t) != NULL)
	  {
            t_s_upper (proposal, proposal->osister, obb);
	  }
        else
	  {
            double_line[0] = proposal->ocousin;
            /*findcrossing needs a last NULL element */
            if (findcrossing (double_line, proposal->line_t) != NULL)
	      {
                t_s_upper (proposal, proposal->ocousin, obb);
	      }
            else
	      {
                gotoroot (obb, proposal->line_f, &proposal->line_f_allocsize);
                proposal->connect =
		  findcrossing (proposal->line_f, proposal->line_t);
                proposal->vs = obb->tyme - proposal->osister->tyme;
                proposal->v = obb->tyme - proposal->ocousin->tyme;
                copy_x (proposal, proposal->xt, proposal->osister->x);
#ifdef UEP
		
                if(proposal->world->options->uep)
		  copy_uepx(proposal,proposal->ut, proposal->osister->ux);
#endif
		
                scalecopy (proposal->mt, proposal->osister->scale, proposal->world);
		if (obb!=proposal->connect) { 
		  uppern = t_s_obbncon (proposal, obb);
                }
		else
		  {
		    uppern = proposal->connect;
		  }
                proposal->v = proposal->time - proposal->origin->tyme;
                copy_x (proposal, proposal->xf, proposal->origin->x);
#ifdef UEP
		
                if(proposal->world->options->uep)
		  copy_uepx(proposal,proposal->uf, proposal->origin->ux);
#endif
		
                scalecopy (proposal->mf, proposal->origin->scale, proposal->world);
                if (proposal->target != proposal->connect)
		  {
		    //                    t_s_tncon (proposal, obb);
		    t_s_tncon (proposal, uppern);
		  }
                else
		  {
                    t_s_tcon (proposal, uppern);
		  }
	      }
	  }
      }
    else
      {
        proposal->v = proposal->time - proposal->origin->tyme;
        proposal->vs = proposal->time - proposal->target->tyme;
        copy_x (proposal, proposal->xf, proposal->origin->x);
#ifdef UEP
	
        if(proposal->world->options->uep)
	  copy_uepx(proposal,proposal->uf, proposal->origin->ux);
#endif
	
        scalecopy (proposal->mf, proposal->origin->scale,proposal->world);
        copy_x (proposal, proposal->xt, proposal->target->x);
#ifdef UEP
	
        if(proposal->world->options->uep)
	  copy_uepx(proposal,proposal->ut, proposal->target->ux);
#endif
	
        scalecopy (proposal->mt, proposal->target->scale, proposal->world);

        //CHECK HERE      memcpy(proposal->mt,proposal->osister->scale, sizeof(MYREAL)*proposal->endsite);
#ifdef BEAGLE
	long id0, id1, id2;
	long bid1, bid2;
	double v1, v2;
	id0  = proposal->oback->id;
	id1  = proposal->origin->id;
	bid1 = proposal->origin->bid;
	v1   = proposal->v;
	id2  = proposal->target->id;
	bid2 = proposal->target->bid;
	v2   = proposal->vs;
	prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
#else
        pseudonuview (proposal, proposal->xf, proposal->mf, proposal->v,
                      proposal->xt, proposal->mt, proposal->vs);
#endif
        proposal->v =
	  fabs (((proposalack =
		  showtop (crawlback (proposal->target)))->tyme -
		 proposal->time));
        if (proposalack != proposal->oback)
	  {
            children (proposalack, &d1, &d2);
            if (d1 == proposal->target)
	      {
		d1 = d2;
		d2 = proposal->target;
	      }
#ifdef BEAGLE
	    id0  = proposalack->id;
	    id1  = new_id(proposal->oback->id,proposal->world->sumtips);
	    bid1 = proposal->oback->bid;
	    v1   = proposal->v;
	    id2  = d1->id;
	    bid2 = d1->bid;
	    v2   = d1->v;
	    prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
#else
            pseudonuview (proposal, proposal->xf, proposal->mf, proposal->v,
                          d1->x, d1->scale, d1->v);
#endif
            oldnn = proposalack;
            nn = showtop (crawlback (oldnn));
            while (nn != proposal->oback)
	      {
                children (nn, &d1, &d2);
                if (d1 == oldnn)
		  {
		    d1 = d2;
		    d2 = oldnn;
		  }
#ifdef BEAGLE
		id0  = nn->id;
		id1  = new_id(oldnn->id,proposal->world->sumtips);
		bid1 = oldnn->bid;
		v1   = oldnn->v;
		id2  = d1->id;
		bid2 = d1->bid;
		v2   = d1->v;
		prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
#else
                pseudonuview (proposal, proposal->xf, proposal->mf, oldnn->v,
                              d1->x, d1->scale, d1->v);
#endif
                oldnn = nn;
                nn = showtop (crawlback (nn));
	      }
	  }
      }
    myfree(double_line);
}
// end safe copy


void
t_s_upper (proposal_fmt * proposal, node * connect, node * obb)
{
    node *nn, *oldnn, *d1, *d2;
#ifdef TREEDEBUG
    printf("target_stancoa:t_s_upper\n");
#endif
    proposal->v = proposal->time - proposal->origin->tyme;
    proposal->vs = proposal->time - proposal->target->tyme;
#ifndef BEAGLE
    copy_x (proposal, proposal->xf, proposal->origin->x);
#ifdef UEP

    if(proposal->world->options->uep)
        copy_uepx(proposal,proposal->uf, proposal->origin->ux);
#endif

    scalecopy (proposal->mf, proposal->origin->scale, proposal->world);
#endif /*is not BEAGLE*/
#ifdef BEAGLE
    //  orig    target
    //    |       |
    //    +-oback-*
    long id0, id1, id2;
    long bid1, bid2;
    double v1, v2;
    id0  = proposal->oback->id;
    id1  = proposal->origin->id;
    bid1 = proposal->origin->bid;
    v1   = proposal->v;
    id2  = proposal->target->id;
    bid2 = proposal->target->bid;
    v2   = proposal->vs;
    prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
#else
    pseudonuview (proposal, proposal->xf, proposal->mf, proposal->v,
                  proposal->target->x, proposal->target->scale, proposal->vs);
#endif
    nn = showtop (crawlback (proposal->target));
    oldnn = proposal->target;
    proposal->vs = fabs (nn->tyme - proposal->time);
#ifdef BEAGLE
    //  orig    target=oldnn
    //    |       |
    //    +-oback-*                                d1
    //        |                                    |
    //        *-nn(target_orig leads to this node)-+
    if(nn!=connect)
      {
	children (nn, &d1, &d2);
	if (d1 == oldnn)
	  {
	    d1 = d2;
	    d2 = oldnn;
	  }
	id0  = nn->id;
	id1  = new_id(proposal->oback->id, proposal->world->sumtips);
	bid1 = proposal->oback->bid;
	v1   = proposal->vs;
	id2  = d1->id;
	bid2 = d1->bid;
	v2   = d1->v;
	prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
	proposal->vs = nn->v;
    //if (nn == connect)
    //  {
	oldnn = nn;
	nn = showtop (crawlback (nn));
	//  }
	//else
	//  {
      }
    //31 31 31 3 13 6 10 15 15 15 23 10 2 9 16 16 16 12 11 15 12 
    //29 29 29 9 7 0 1 27 27 27 29 11 8 14 28 28 28 27 15 13 4 31 31 31 11 3 3 13 16 16 16 31 8 32 12 

#endif  
    while (nn != connect && nn->type!='r')
      {
	children (nn, &d1, &d2);
	if (d1 == oldnn)
	  {
	    d1 = d2;
	    d2 = oldnn;
	  }
#ifdef BEAGLE
	id0  = nn->id;
	id1  = new_id(oldnn->id,proposal->sumtips);
	bid1 = oldnn->bid;
	v1   = proposal->vs;
	id2  = d1->id;
	bid2 = d1->bid;
	v2   = d1->v;
	prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
#else
	pseudonuview (proposal, proposal->xf, proposal->mf, proposal->vs, d1->x,
		      d1->scale, d1->v);
#endif
	proposal->vs = nn->v;
	oldnn = nn;
	nn = showtop (crawlback (nn));
      }
    if(nn->type!='r')
       {
	 children (nn, &d1, &d2);
	if (d1 == oldnn)
	    {
	      d1 = d2;
	      d2 = oldnn;
	    } 
#ifdef BEAGLE
	// getting the conditiionals of the connection node
	id0  = nn->id;
	if(nn==connect && oldnn==proposal->target)
	  {
	    id1  = new_id(proposal->oback->id, proposal->world->sumtips);
	    bid1 = proposal->oback->bid;
	  }
	else
	  {
	    id1  = new_id(oldnn->id, proposal->world->sumtips);
	    bid1 = oldnn->bid;
	  }
	v1   = proposal->vs;
	id2  = d1->id;
	bid2 = d1->bid;
	v2   = d1->v;
	prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
#else
	pseudonuview (proposal, proposal->xf, proposal->mf, proposal->vs, d1->x,
		      d1->scale, d1->v);
#endif
       }
    proposal->v = obb->tyme - proposal->osister->tyme;
    proposal->vs = obb->tyme - proposal->ocousin->tyme;
    if (connect == proposal->ocousin)
      {
#ifdef BEAGLE
	//12 12 12 0 1 4 0 32 32 32 12 11 5 5 33 33 33 32 12 1 6 31 31 31 33 16 6 10 11 11 11 31 8 13 4 9 9 9 11 3 27 15 

	id0  = obb->id;
	id1  = proposal->osister->id;
	bid1 = proposal->osister->bid;
	v1   = proposal->v;
	id2  = new_id(proposal->ocousin->id,proposal->sumtips);
	bid2 = proposal->ocousin->bid;
	v2   = proposal->vs;
	prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
#else
	pseudonuview (proposal, proposal->xf, proposal->mf, proposal->vs,
		      proposal->osister->x, proposal->osister->scale,
		      proposal->v);
#endif
	  }
    else
      {
#ifdef BEAGLE
	id0  = obb->id;
	if(connect==proposal->osister)
	  id1  = new_id(proposal->osister->id,proposal->sumtips);
	else
	  id1 = proposal->osister->id;
	bid1 = proposal->osister->bid;
	v1   = proposal->v;
	id2  = proposal->ocousin->id;
	bid2 = proposal->ocousin->bid;
	v2   = proposal->vs;
	prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
#else
	pseudonuview (proposal, proposal->xf, proposal->mf, proposal->v,
		      proposal->ocousin->x, proposal->ocousin->scale,
		      proposal->vs);
#endif
      }
#ifdef BEAGLE
	//}
    evaluate_beagle_instances_proposal (proposal, showtop (crawlback (obb)), obb, obb->id, obb->bid,  obb->v);
#else
    pseudoevaluate(proposal,proposal->xf,proposal->mf,showtop (crawlback (obb)), obb, obb->v );
#endif
}

///
/// calculates conditionals from obb to connection, used in target_stancoa()
node * t_s_obbncon (proposal_fmt * proposal, node * obb)
{

    node *nn, *oldnn, *d1, *d2;
    /*     FPRINTF(stdout,"t_s_obbncon\n"); */
#ifdef BEAGLE
    long id0, id1, id2;
    long bid1, bid2;
    double v1, v2;
#ifdef TREEDEBUG
    printf("target_stancoa:t_s_obbncon\n");
#endif
    id0  = obb->id;
    id1  = proposal->osister->id;
    bid1 = proposal->osister->bid;
    v1   = proposal->vs;
    id2  = proposal->ocousin->id;
    bid2 = proposal->ocousin->bid;
    v2   = proposal->v;
    prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
#else
    pseudonuview (proposal, proposal->xt, proposal->mt, proposal->vs,
                  proposal->ocousin->x, proposal->ocousin->scale, proposal->v);
#endif
    nn = showtop (crawlback (obb));
    oldnn = obb;
    proposal->vs = fabs (nn->tyme - obb->tyme);
    while (nn != proposal->connect)
    {
        children (nn, &d1, &d2);
        if (d1 == oldnn)
	  {
	    d1 = d2;
	    d2 = oldnn;
	  }
#ifdef BEAGLE
	id0  = nn->id; // will be offset in prepare_beagle_instances_proposal()
	id1  = new_id(oldnn->id, proposal->world->sumtips);
	bid1 = oldnn->bid;
	v1   = proposal->vs;
	id2  = d1->id;
	bid2 = d1->bid;
	v2   = d1->v;
	prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
#else
        pseudonuview (proposal, proposal->xt, proposal->mt, proposal->vs, d1->x,
                      d1->scale, d1->v);
#endif
        proposal->vs = nn->v;
        oldnn = nn;
        nn = showtop (crawlback (nn));
    }
    return oldnn;
}

///
/// calculates lines from target ---> connect <--- obb
///                                      |
void t_s_tncon (proposal_fmt * proposal, node * fromobb)
{
    node *nn, *oldnn, *d1, *d2;
#ifdef BEAGLE
    long id0, id1, id2;
    long bid1, bid2;
    double v1, v2;
#ifdef TREEDEBUG
    printf("target_stancoa:t_s_tncon\n");
#endif
    id0  = proposal->oback->id;//id0 will point to the new id+nodep_boundary
    id1  = proposal->origin->id;
    bid1 = proposal->origin->bid;
    v1   = proposal->v;
    id2  = proposal->target->id;
    bid2 = proposal->target->bid;
    v2   = proposal->time - proposal->target->tyme;
    // connect oback, origin, target
    prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
#else
    pseudonuview (proposal, proposal->xf, proposal->mf, proposal->v,
                  proposal->target->x, proposal->target->scale,
                  proposal->time - proposal->target->tyme);
#endif
    nn = showtop (crawlback (proposal->target));
    oldnn = proposal->target;    
    // calculates the time between the newly inserted oback and the back of target
    proposal->v = fabs (nn->tyme - proposal->time);
#ifdef BEAGLE
    children (nn, &d1, &d2);
    if (d1 == oldnn)
      {
	d1 = d2;
      }
    id0  = nn->id;
    id1  = new_id(proposal->oback->id,proposal->world->sumtips);
    bid1 = proposal->oback->bid;
    v1   = proposal->v;
    if(d1!=fromobb)
      id2  = d1->id;
    else
      id2  = new_id(d1->id,proposal->sumtips);
    bid2 = d1->bid;
    v2   = d1->v;
    prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
    proposal->v = nn->v;
    oldnn = nn;
    if(nn != proposal->connect)
      {
	nn = showtop (crawlback (nn));
#endif
	while (nn != proposal->connect)
	  {
	    children (nn, &d1, &d2);
	    if (d1 == oldnn)
	      {
		d1 = d2;
		d2 = oldnn;
	      }
#ifdef BEAGLE
	    id0  = nn->id;
	    id1  = new_id(oldnn->id,proposal->world->sumtips);
	    bid1 = oldnn->bid;
	    v1   = oldnn->v;
	    id2  = d1->id;
	    bid2 = d1->bid;
	    v2   = d1->v;
	    prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
#else
	    pseudonuview (proposal, proposal->xf, proposal->mf, proposal->v, d1->x, d1->scale, d1->v);
#endif
	    proposal->v = nn->v;
	    oldnn = nn;
	    nn = showtop (crawlback (nn));
	  }
	// why obb here this is target ---> root
	// nn is connect now
	if (fromobb != proposal->connect)
	  {
#ifdef BEAGLE
	    id0  = nn->id;
	    id1  = new_id(oldnn->id, proposal->world->sumtips);
	    bid1 = oldnn->bid;
	    v1   = proposal->v;
	    id2  = new_id(fromobb->id, proposal->world->sumtips);
	    bid2 = fromobb->bid;
	    v2   = fromobb->v;
	    prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
#else
	    pseudonuview (proposal, proposal->xf, proposal->mf, proposal->v,
			  proposal->xt, proposal->mt, proposal->vs);
#endif
	  }
#ifdef BEAGLE
      }
#endif
    proposal->v = proposal->connect->v;
#ifdef BEAGLE
    //this assume that we have finished all above the connect node
    //and the connect node is a daughter in the evaluation
    evaluate_beagle_instances_proposal (proposal,
					showtop (crawlback (proposal->connect)), proposal->connect, 
					proposal->connect->id, proposal->connect->bid,
					proposal->v);
#else
    pseudoevaluate (proposal, proposal->xf, proposal->mf,
                    showtop (crawlback (proposal->connect)), proposal->connect,
                    proposal->v);
#endif
}

void
t_s_tcon (proposal_fmt * proposal, node * uppern)
{
    node *obb, *d1, *d2;
#ifdef TREEDEBUG
    printf("target_stancoa:t_s_tcon\n");
#endif
    children (proposal->target, &d1, &d2);
    if (d1 == uppern)
      {
	d1 = d2;
	d2 = uppern;
      }
    proposal->vs = fabs (proposal->target->tyme - showtop (uppern)->tyme);
#ifdef BEAGLE
    long id0, id1, id2;
    long bid1, bid2;
    double v1, v2;
    id0  = proposal->target->id;
    id1  = new_id(uppern->id, proposal->world->sumtips);
    bid1 = uppern->bid;
    v1   = proposal->vs;
    id2  = d1->id;
    bid2 = d1->bid;
    v2   = d1->v;
    prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
#else
    pseudonuview (proposal, proposal->xt, proposal->mt, proposal->vs, d1->x, d1->scale, d1->v); /* eval target */
#endif
    proposal->vs = proposal->time - proposal->target->tyme;
#ifdef BEAGLE
    id0  = proposal->oback->id;
    id1  = proposal->origin->id;
    bid1 = proposal->origin->bid;
    v1   = proposal->v;
    id2  = new_id(proposal->target->id,proposal->sumtips);
    bid2 = proposal->target->bid;
    v2   = proposal->vs;
    prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
#else
    pseudonuview (proposal, proposal->xf, proposal->mf, proposal->v, proposal->xt, proposal->mt, proposal->vs); /*eval oback */
#endif
    obb = showtop (crawlback (proposal->target));
    proposal->v = fabs (obb->tyme - proposal->time);
    if (obb->type == 'r')
        return;
    else
    {
        children (obb, &d1, &d2);
        if (d1 == proposal->target)
	  {
	    d1 = d2;
	    d2 = proposal->target;
	  }
#ifdef BEAGLE
	id0  = obb->id;
	id1  = new_id(proposal->oback->id, proposal->world->sumtips);
	bid1 = proposal->oback->bid;
	v1   = proposal->v;
	id2  = d1->id;
	bid2 = d1->bid;
	v2   = d1->v;
	prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
	evaluate_beagle_instances_proposal (proposal,
					showtop (crawlback (obb)), obb, obb->id, obb->bid,
					showtop (crawlback (obb))->tyme - obb->tyme);
#else
        pseudonuview (proposal, proposal->xf, proposal->mf, proposal->v, d1->x,
                      d1->scale, d1->v);
        pseudoevaluate (proposal, proposal->xf, proposal->mf,
                        showtop (crawlback (obb)), obb, obb->v);
#endif
    }
}


/*-----------------------------------------------------------------
erase all migration nodes between proposal->origin and proposal->oback 
*/
long
erase_migr_nodes (world_fmt *world, node * up)
{
    long deleted = 0;
    node *theNode, *down, *oldNode;
    down = crawlback (up);
    theNode = showtop(up->back);
    //    printf("-------------\nup=%li\n",up->id);
    while (theNode->type == 'm' || theNode->type == 'd')
    {

        oldNode = theNode;
        theNode = showtop(theNode->back);
	free_mignodelet(oldNode,world);
	deleted++;
	  //printf("-");
    }
    //printf("down=%li deleted=%li\n-------------\n",down->id, 2*deleted);
    down->back = up;
    up->back = down;
    //if(world->heat < 0.001)
    //  printf("%i> heat=%f deleted=%li\n",myID,1./world->heat,deleted);
    return deleted;
}

void erase_migr_nodes2 (world_fmt *world, node * up)
{
  int skip=0;
    long deleted = 0;
    node *theNode, *oldNode;
    //node *down
    if(skip!=0)
      return;
    //xcode down = crawlback (up);
    theNode = showtop(up->back);
    //    printf("-------------\nup=%li\n",up->id);
    while (theNode->type == 'm' || theNode->type == 'd')
    {
        oldNode = theNode;
        theNode = showtop(theNode->back);
	free_mignodelet(oldNode,world);
	deleted++;
	  //printf("-");
    }
}
/*-------------------------------------------------------
Connects and adjusts three nodes, the first is the 
mother node and the next two are child nodes, if one
of the child nodes is NULL it just connects mother with
the not-NULL child.
*/
void
connectnodes (world_fmt *world, node * mother, node * brother, node * sister)
{
#ifndef MESS
  (void) world;
#endif
    node *tmp;
    if(mother->top!=1)
      error("mother is not on top");
    if(brother !=0 && brother->top!=1)
      error("brother is not on top");
    if(sister!=NULL && sister->top!=1)
      error("sister is not on top");
    if ((sister != NULL) && (brother != NULL))
    {
        if ((mother == brother) || (mother == sister)
                || (sister == brother))
        {
            error ("connectnodes() conflict");
        }
        tmp = crawl_down (brother, mother->tyme);
#ifdef MESS
	if(tmp->back != NULL && (tmp->back->type == 'm' || tmp->back->type == 'd') && showtop(tmp->back)->tyme > mother->tyme)
	  {
	    printf("new code 1 \n");
	    erase_migr_nodes2(world,showtop(tmp));
	  }
#endif 
        mother->next->back = tmp;
        tmp->back = mother->next;
	
        tmp = crawl_down (sister, mother->tyme);
#ifdef MESS
	if(tmp->back != NULL && (tmp->back->type == 'm' || tmp->back->type == 'd') && showtop(tmp->back)->tyme > mother->tyme)
	  {
	    printf("new code 2\n");
	    erase_migr_nodes2(world,showtop(tmp));
	  }
#endif

        mother->next->next->back = tmp;
        tmp->back = mother->next->next;
    }
    else
    {
        if (sister == NULL)
        {
            tmp = crawl_down (brother, mother->tyme);
#ifdef MESS
	    if(tmp->back != NULL && (tmp->back->type == 'm' || tmp->back->type == 'd') && showtop(tmp->back)->tyme > mother->tyme)
	  {
	    printf("new code 3\n");
	    erase_migr_nodes2(world,showtop(tmp));
	  }
#endif

            mother->next->back = tmp;
            tmp->back = mother->next;
        }
        else
        {

            if (brother == NULL)
            {
                if (mother->type == 'm' || mother->type == 'd' )
                {
                    tmp = crawl_down (sister, mother->tyme);
#ifdef MESS
		    if(tmp->back != NULL && (tmp->back->type == 'm' || tmp->back->type == 'd') && showtop(tmp->back)->tyme > mother->tyme)
		      {
			printf("new code 4 \n");
			erase_migr_nodes2(world,showtop(tmp));
		      }
#endif
                    mother->next->back = tmp;
                    tmp->back = mother->next;
                }
                else
                {
                    tmp = crawl_down (sister, mother->tyme);
#ifdef MESS
		    if(tmp->back != NULL && (tmp->back->type == 'm' || tmp->back->type == 'd') && showtop(tmp->back)->tyme > mother->tyme)
	  {
	    printf("new code 5\n");
	    erase_migr_nodes2(world,showtop(tmp));
	  }
#endif

                    mother->next->next->back = tmp;
                    tmp->back = mother->next->next;
                }
            }
            else
            {
                error ("Single, lonely rootnode detected in connectenodes()\n");
            }
        }
    }
}

void
gotoroot (node * origin, node ** ptrlist, long *ptrallocsize)
{
    node *theNode;
    long i = 0;

    for (theNode = origin;
            (theNode != NULL) && (crawlback (theNode)->type != 'r');
            theNode = showtop (crawlback (theNode)))
    {
      if (i >= *ptrallocsize-3)
	allocate_nodelist(&ptrlist,ptrallocsize, HUNDRED);
      ptrlist[i++] = theNode;
    }
    ptrlist[i] = theNode;  /* adds root->back to the list */
    ptrlist[i+1] = NULL;   /* guarantees an end for the array)*/
}

void
adjust (node * theNode, MYREAL tyme, long level)
{
    if (level < 0 || theNode == NULL)
        return;

    theNode->tyme = tyme;
    if (theNode->type == 'r')
    {
        theNode->v = 1;
        theNode->length = MYREAL_MAX;
        if (theNode->next->back != NULL)
            adjust (crawlback (theNode->next), crawlback (theNode->next)->tyme,
                    level);
        if (theNode->next->next->back != NULL)
            adjust (crawlback (theNode->next->next),
                    crawlback (theNode->next->next)->tyme, level);
    }
    else if (theNode->type == 't')
    {
      //#ifndef TESTINGDATE
      //theNode->tyme = 0.0;
      //#endif
        theNode->length = lengthof (theNode);
        ltov (theNode);
        return;
    }
    else if ((theNode->type != 't'))
    {
        if (theNode->type != 'm' && theNode->type != 'd')
        {
            theNode->length = lengthof (theNode);
            ltov (theNode);
            if (theNode->next->back != NULL)
            {
                theNode->next->tyme = tyme;
                theNode->next->v = theNode->v;
                theNode->next->length = theNode->length;
                adjust (crawlback (theNode->next),
                        crawlback (theNode->next)->tyme, level - 1);
            }
            if (theNode->next->next->back != NULL)
            {
                theNode->next->next->tyme = tyme;
                theNode->next->next->v = theNode->v;
                theNode->next->next->length = theNode->length;
                adjust (crawlback (theNode->next->next),
                        crawlback (theNode->next->next)->tyme, level - 1);
            }
        }
        else
            adjust (crawlback (theNode->next), crawlback (theNode->next)->tyme,
                    level);
    }
}


/* calculates x-array down the tree assuming that only one line
   is affected by the change of the sub-likelihoods above that line
   BUT does NOT calculate the tree-likelihood as evaluate() does.
   THIS CHANGES THE TREE
 */
void
localevaluate (node * mother)
{
    node *nn = NULL;

    if (mother->type != 'r')
    {
        set_dirty (mother);
        for (nn = mother; crawlback (nn)->type != 'r';
                nn = showtop (crawlback (nn)))
        {
            set_dirty (nn);
        }
    }
}


///
/// copy the sequence data
void
copy_x (proposal_fmt * proposal, xarray_fmt *axx1, xarray_fmt *axx2)
{
  world_fmt *world  = proposal->world;
  long locus        = world->locus;
  mutationmodel_fmt *s;
  long sublocus;
  long sublocistart = world->sublocistarts[locus];
  long sublociend   = world->sublocistarts[locus+1];
  for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
    {
      s = &world->mutationmodels[sublocus];
      long endsite      = s->numpatterns + s->addon;
      long rcategs      = s->numsiterates;
      //long *maxalleles  = s->maxalleles;
      const long xs = sublocus - sublocistart;
      xarray_fmt  xx1   = axx1[xs];
      xarray_fmt  xx2   = axx2[xs];
      const size_t sitelikesize = sizeof (sitelike) * (size_t) rcategs;
      long i;
      
      switch (s->datatype)
	{
	case 'a':
	case 'b':
	case 'm':
	  memcpy (xx1.a, xx2.a, sizeof (MYREAL) * (size_t) s->maxalleles);
	  break;
	case 's':
	case 'n':
	case 'h':
	case 'u':
	case 'f':
	  for (i = 0; i < endsite; i++)
	    {
	      memcpy (xx1.s[i], xx2.s[i], sitelikesize);
	      long ii;
	      long sum=0;
	      for(ii=0;ii<4;ii++)
		{
		  if(xx2.s[i][0][ii] == 0.0)
		    sum += 1;
		}
	      if(sum == 4)
		error("copy_x copied only zeroes");
	    }
	  break;
	}
    }
}

void
fix_root_pop (world_fmt *world, node * p)
{
    if (crawlback (p) != p->back)
    {
      erase_migr_nodes (world, p);
    }
    if (p->back->actualpop != p->actualpop)
    {
        p->back->pop = p->back->next->pop = p->back->next->next->pop = p->pop;
        p->back->actualpop = p->back->next->actualpop = p->actualpop;
        p->back->next->next->actualpop = p->actualpop;
    }
}


/* transfers an x-array down the tree assuming that only one line
   is affected by the change of the sub-likelihoods above that line
   BUT does NOT calculate the tree-likelihood as evaluate() does.
   DOES NOT CHANGE THE TREE
 */
void
pseudoevaluate (proposal_fmt * proposal, xarray_fmt *ax, MYREAL **lx,
                node * mother, node * newdaughter, MYREAL v)
{
    node *nn = NULL, *d1 = NULL, *d2 = NULL, *oldnn = NULL;

    if (mother->type != 'r')
    {
        children (mother, &d1, &d2);
        if (d1 == newdaughter)
	  {
	    d1 = d2;
	    d2 = newdaughter;
	  }
#ifdef BEAGLE
	long id0, id1, id2;
	long bid1, bid2;
	double v1, v2;
	id0  = mother->id; 
	id1  = new_id(proposal->oback->id,proposal->world->sumtips);
	bid1 = proposal->oback->bid;
	v1 = v;
	id2  = d1->id;
	bid2 = d1->bid;
	v2 = d1->v;
	prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
#else
        pseudonuview (proposal, ax, lx, v, d1->x, d1->scale, d1->v);
#endif
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
#ifdef BEAGLE
	    id0  = nn->id; 
	    id1  = new_id(d2->id,proposal->world->sumtips);
	    bid1 = d2->bid;
	    v1 = d2->v;
	    id2  = d1->id;
	    bid2 = d1->bid;
	    v2 = d1->v;
	    prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);
#else
            pseudonuview (proposal, ax, lx, oldnn->v, d1->x, d1->scale, d1->v);
#endif
            oldnn = nn;
            nn = showtop (crawlback (nn));
        }
    }
}

node *
findcrossing (node ** ptrl1, node ** ptrl2)
{
    long i = 0, j = 0;

    /* assumes that there is an NULL element at the end */
    for (i = 0; ptrl1[i] != NULL; j = 0, i++)
    {
        while ((ptrl2[j] != NULL) && (ptrl1[i] != ptrl2[j]))
            j++;
        if (ptrl2[j] != NULL)
        {
            break;
        }
    }
    return ptrl1[i];
}

node *
crawl_down (node * theNode, MYREAL tyme)
{
  node *tN = showtop(theNode);
  node *otmp, *tmp = tN->back;

    otmp = tN;
    if (tmp == NULL)
        return otmp;
    while ((tmp->type == 'm' || tmp->type == 'd') && showtop (tmp)->tyme < tyme)
    {
#ifdef TREEDEBUG1
      printf("Crawldown: start=%x<%x> (%li<%li>)  %li %f %f\n",theNode, tN, theNode->id, tN->id,  tmp->id, showtop (tmp)->tyme, tyme);
#endif
        otmp = tmp->next;
        tmp = tmp->next->back;
        if (tmp == NULL)
            return otmp;
    }
    return otmp;
}
