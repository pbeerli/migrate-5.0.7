/*------------------------------------------------------
 Maximum likelihood estimation
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm
 -------------------------------------------------------
 Unique event polymorphism   R O U T I N E S
 
 Peter Beerli 2000, Seattle
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
 
$Id: uep.c 2158 2013-04-29 01:56:20Z beerli $
*/
/* \file uep.c 
Unique event polymorphism manipulation
*/
#ifdef UEP
#include "migration.h"
#include "tools.h"
#include "random.h"
#include "uep.h"
#include "sighandler.h"
#include "migrate_mpi.h"

#ifdef DMALLOC_DEBUG
#include <dmalloc.h>
#endif

#define AVERAGE    0  /* used for time of uep */
#define POPULATION 1  /* used for time of uep given that it raised in */
/*   pop x */

node *first_uep (node * p, node * nop, long uepsites);
boolean uep_check (node * root, int *uep, long uepsites);
node *first_uep2 (node * p, node * root, long uepsites);
MYREAL sum_ueplike (proposal_fmt * proposal);

void print_uep (FILE * out, world_fmt * world, long copies);
MYREAL *calc_uep_popprob (world_fmt * world, long copies, long rootstate);
MYREAL *calc_uep_timeprob (world_fmt * world, long copies, int type,
                           long rootstate);
long normalize_uep (world_fmt * world);
void analyze_uep (world_fmt * world);
void print_ancestor (FILE * out, world_fmt * world, long copies);
void update_uepanc (world_fmt * world);
void check_uep_root (node * p, world_fmt * world);

void
setup_uep (world_fmt * world, option_fmt * options)
{
    if (options->uep)
    {
        alloc_ueplike (&(world->ueplike), world->data->uepsites,
                       &(world->ueplikestore), world->options->lsteps,
                       &(world->ueprootstore), world->numpop, world->loci);
        alloc_ueptime (&(world->ueptime), world->data->uepsites,
                       &(world->ueptimestore), world->options->lsteps,
                       world->numpop, world->loci);
        world->oldrootuep =
            (long *) mycalloc (world->data->uepsites, sizeof (long));

        alloc_uepanc (&(world->uepanc), world->data->uepsites, world->numpop,
                      world->loci);
    }
}

void
destroy_uep (world_fmt * world, long storesize)
{
    long loci = world->loci;
    long locus;
    long i;
    long j;
  if (world->options->uep)
    {
      myfree(world->uepanc[0]);
      myfree(world->uepanc);
      myfree(world->ueplike[0]);
      myfree(world->ueplike);
      for (locus = 0; locus < loci; locus++)
	{
	  for (i = 0; i < storesize; ++i)
	    {
	      for (j = 0; j < world->options->lsteps; ++j)
		{
		  myfree(world->ueplikestore[locus][i][j]);
		  myfree(world->ueptime[j].populations);
		  myfree(world->ueptime[j].ueptime);
		}
	      myfree(world->ueplikestore[locus][i]);
	      myfree(world->ueptimestore[locus][i]);
	      myfree(world->ueprootstore[locus][i]);
	    }
	  myfree(world->ueplikestore[locus]);
	  myfree(world->ueptimestore[locus]);
	  myfree(world->ueprootstore[locus]);
	}
      myfree(world->ueplikestore);
      myfree(world->ueptimestore);
      myfree(world->ueprootstore);
      for (j = 0; j < world->options->lsteps; ++j)
	{
	  myfree(world->ueptime[j].populations);
	  myfree(world->ueptime[j].ueptime);
	}
      myfree(world->ueptime);
      myfree(world->oldrootuep);    
    }
}

void
pseudonu_twostate (proposal_fmt * proposal,
                   ueparray_fmt * xx1, MYREAL *lx1, MYREAL v1,
                   ueparray_fmt * xx2, MYREAL *lx2, MYREAL v2)
{
    //  static long count=0;
    long j;
    MYREAL pija11, pija12, pija21, pija22;
    //  MYREAL x3m = -MYREAL_MAX;
    MYREAL r = proposal->world->options->uepmu;
    MYREAL s = proposal->world->options->uepnu;
    MYREAL rs = r + s;
    pair *x1;
    pair *x2;
    v1 = EXP (-v1 * rs);
    v2 = EXP (-v2 * rs);
    r = r / rs;
    s = s / rs;
    for (j = 0; j < proposal->world->data->uepsites; j++)
    {
        x1 = &(xx1->s[j]);
        x2 = &(xx2->s[j]);
        pija11 = (v1 * r + s) * (*x1)[0] + (s - v1 * s) * (*x1)[1];
        pija12 = (v2 * r + s) * (*x2)[0] + (s - v2 * s) * (*x2)[1];
        pija21 = (s - v1 * s) * (*x1)[0] + (v1 * s + r) * (*x1)[1];
        pija22 = (s - v2 * s) * (*x2)[0] + (v2 * s + r) * (*x2)[1];
        (*x1)[0] = pija11 * pija12;
        (*x1)[1] = pija21 * pija22;
    }
    //*lx1 += lx2;
    //count++;
    //if (count == SCALEINTERVAL)
    // {
    // count = 0;
    //if ((*xx1)[0] > x3m)
    //      x3m = (*xx1)[0];
    //      if ((*xx1)[1] > x3m)
    //      x3m = (*xx1)[1];
    //      (*xx1)[0] /= x3m;
    //     (*xx1)[1] /= x3m;
    //      *lx1 += log (x3m);
    //    }
}


void
twostate_nuview (node * mother, world_fmt * world, const long locus)
{
    //  static long count=0;
    long j;
    node *d1 = NULL, *d2 = NULL;
    MYREAL v1, v2;
    MYREAL pija11, pija12, pija21, pija22;
    pair *xx1, *xx2;
    //  MYREAL x3m = -MYREAL_MAX;
    pair *xx3;
    MYREAL r = world->options->uepmu;
    MYREAL s = world->options->uepnu;
    MYREAL rs = r + s;
    children (mother, &d1, &d2);
    v1 = EXP (-d1->v * rs);
    v2 = EXP (-d2->v * rs);
    r = r / rs;
    s = s / rs;
    for (j = 0; j < world->data->uepsites; j++)
    {
        xx1 = &d1->ux.s[j];
        xx2 = &d2->ux.s[j];
        xx3 = &mother->ux.s[j];

        pija11 = (v1 * r + s) * (*xx1)[0] + (s - v1 * s) * (*xx1)[1];
        pija12 = (v2 * r + s) * (*xx2)[0] + (s - v2 * s) * (*xx2)[1];
        (*xx3)[0] = pija11 * pija12;
        pija21 = (s - v1 * s) * (*xx1)[0] + (v1 * s + r) * (*xx1)[1];
        pija22 = (s - v2 * s) * (*xx2)[0] + (v2 * s + r) * (*xx2)[1];
        (*xx3)[0] = pija21 * pija22;
    }
    //  count++;
    //  mother->scale[0] = d1->scale[0] + d2->scale[0];
    //  if (count == SCALEINTERVAL)
    //    {
    //      count = 0;
    //      if (xx3[0] > x3m)
    //      x3m = xx3[0];
    //      if (xx3[1] > x3m)
    //      x3m = xx3[1];
    //      xx3[0] /= x3m;
    //      xx3[1] /= x3m;
    //      mother->scale[0] += log (x3m);
    //    }
}

void
alloc_uepanc (long ***uepanc, long size, long numpop, long loci)
{
    long locus;
    (*uepanc) = (long **) mycalloc (loci, sizeof (long *));
    (*uepanc)[0] = (long *) mycalloc (loci * 2 * size * numpop, sizeof (long));
    for (locus = 1; locus < loci; locus++)
      (*uepanc)[locus] = (*uepanc)[0] + 2 * size * numpop;
}

void
alloc_ueplike (MYREAL ***ueplike, long size,
               MYREAL *****ueplikestore, long storesize,
               long ****ueprootstore, long numpop, long loci)
{
    long i, j, locus;
    (*ueplike) = (MYREAL **) mycalloc (size, sizeof (MYREAL *));
    (*ueplike)[0] = (MYREAL *) mycalloc (size * numpop, sizeof (MYREAL));
    for (i = 1; i < size; ++i)
        (*ueplike)[i] = (*ueplike)[0] + i * numpop;

    (*ueprootstore) = (long ***) mycalloc (loci, sizeof (long **));
    (*ueplikestore) = (MYREAL ****) mycalloc (loci, sizeof (MYREAL ***));

    for (locus = 0; locus < loci; locus++)
    {
        (*ueprootstore)[locus] = (long **) mycalloc (storesize, sizeof (long *));
        (*ueplikestore)[locus] =
            (MYREAL ***) mycalloc (storesize, sizeof (MYREAL **));
        for (i = 0; i < storesize; ++i)
        {
            (*ueprootstore)[locus][i] = (long *) mycalloc (size, sizeof (long));
            (*ueplikestore)[locus][i] =
                (MYREAL **) mycalloc (size, sizeof (MYREAL *));
            for (j = 0; j < size; ++j)
                (*ueplikestore)[locus][i][j] =
                    (MYREAL *) mycalloc (numpop, sizeof (MYREAL));
        }
    }
}

void
alloc_ueptime (ueptime_fmt ** ueptime, long size,
               ueptime_fmt **** ueptimestore,
               long storesize, long numpop, long loci)
{
    long i, j, locus;
    (*ueptime) = (ueptime_fmt *) mycalloc (size, sizeof (ueptime_fmt));
    for (j = 0; j < size; ++j)
    {
        (*ueptime)[j].size = 3;
        (*ueptime)[j].populations = (long *) mycalloc (3, sizeof (long));
        (*ueptime)[j].ueptime = (MYREAL *) mycalloc (3, sizeof (MYREAL));
    }
    (*ueptimestore) = (ueptime_fmt ***) mycalloc (loci, sizeof (ueptime_fmt **));
    for (locus = 0; locus < loci; locus++)
    {
        (*ueptimestore)[locus] = (ueptime_fmt **) mycalloc (storesize,
                                 sizeof (ueptime_fmt
                                         *));
        (*ueptimestore)[locus][0] =
            (ueptime_fmt *) mycalloc (size * storesize, sizeof (ueptime_fmt));
        for (i = 1; i < storesize; ++i)
        {
            (*ueptimestore)[locus][i] = (*ueptimestore)[locus][0] + i * size;
        }
        for (i = 0; i < storesize; ++i)
        {
            for (j = 0; j < size; ++j)
            {
                (*ueptimestore)[locus][i][j].size = 3;
                (*ueptimestore)[locus][i][j].populations =
                    (long *) mycalloc (3, sizeof (long));
                (*ueptimestore)[locus][i][j].ueptime =
                    (MYREAL *) mycalloc (3, sizeof (MYREAL));
            }
        }
    }
}

void
allocate_uep (node * p, world_fmt * world, char datatype, boolean withtips)
{
    if (p->type != 't')
    {
        if (p->next->back != NULL)
            allocate_uep (crawlback (p->next), world, datatype, withtips);
        if (p->next->next->back != NULL)
            allocate_uep (crawlback (p->next->next), world, datatype, withtips);
        //      p->ux.s = (pair *) mycalloc (world->data->uepsites, sizeof (pair));
    }
    else
    {
        if (withtips)
        {
            //   p->ux.s = (pair *) mycalloc (world->data->uepsites, sizeof (pair));
        }
    }
}

void
constrain_distance_uep (int **uep, long uepsites, MYREAL **distm, long tips)
{
    long i, j;
    //  long ueppatt = howmany_uep_pattern(uep,uepsites, tips);
    //  long ueppatt = 2;
    for (i = 0; i < tips; ++i)
    {
        for (j = 0; j < i; ++j)
        {
            if (uep[i][0] == uep[j][0])
            {
                distm[i][j] = 0.;
                distm[j][i] = 0.;
            }
            else
            {
                distm[i][j] = 1.;
                distm[j][i] = 1.;
            }
        }
    }
}

void
set_uepvalues (data_fmt * data, node ** treenode, long tii, long ii)
{
    long j;
    for (j = 0; j < data->uepsites; j++)
    {
        switch (data->uep[ii][j])
        {
        case '0':
            treenode[tii]->ux.s[j][0] = 1.;
            treenode[tii]->ux.s[j][1] = 0.;
            break;
        case '1':
            treenode[tii]->ux.s[j][0] = 0.;
            treenode[tii]->ux.s[j][1] = 1.;
            break;
        case '?':
        default:
            treenode[tii]->ux.s[j][0] = treenode[tii]->ux.s[j][1] = 1.;
            break;
        }
    }
}

void
make_uep_values (world_fmt * world, data_fmt * data, long locus)
{
    long pop, ind, ii = 0, poptips;
    long halfuep = world->sumtips / 2;
    long tmp = 0;
    node **treenode = world->nodep;
    if (strchr (SEQUENCETYPES, world->options->datatype))
    {
        for (ii = 0; ii < world->sumtips; ii++)
        {
            set_uepvalues (data, treenode, ii, ii);
        }
    }
    else
    {
        for (pop = 0; pop < world->numpop; pop++)
        {
            poptips = data->numind[pop][0] / 2;
            for (ind = 0; ind < poptips; ind++)
            {
                set_uepvalues (data, treenode, tmp + ind, ii);
                set_uepvalues (data, treenode, tmp + ind + poptips,
                               halfuep + ii);
                ii++;
            }
            tmp += poptips;
        }
    }
}

MYREAL
pseudo_tl_uep (ueparray_fmt * xx1, ueparray_fmt * xx2, MYREAL v1, MYREAL v2,
               proposal_fmt * proposal, world_fmt * world)
{
    MYREAL summ;
    long i;
    MYREAL tterm;
    pair *x1;
    summ = 0.0;
    for (i = 0; i < world->data->uepsites; i++)
    {
        x1 = &(xx1->s[i]);
        tterm =
            world->options->uepfreq0 * (*x1)[0] + world->options->uepfreq1 * (*x1)[1];
        summ += (LOG (tterm) + proposal->mf[i]);
    }
    return summ;
}


MYREAL
treelike_uep (world_fmt * world, long locus)
{
    MYREAL tterm;
    MYREAL summ;
    long i;
    //  MYREAL scale;
    node *p;
    pair *x1;
    p = crawlback (world->root->next);
    summ = 0.0;
    for (i = 0; i < world->data->uepsites; i++)
    {
        x1 = &(p->ux.s[i]);
        //        scale = p->uepscale[i];
        tterm =
            world->options->uepfreq0 * (*x1)[0] + world->options->uepfreq1 * (*x1)[1];
        summ += (LOG (tterm)); //+ scale);
    }
    return summ;
}    /* treelikelihood */


// uses the conditional likelihoods to set
// the actual state at a node
void
update_uep (node * p, world_fmt * world)
{
    long i;
    if (p->type != 't')
    {
        if (p->next->back != NULL)
            update_uep (crawlback (p->next), world);
        if (p->next->next->back != NULL)
            update_uep (crawlback (p->next->next), world);

        if (p->type != 'r')
        {
            for (i = 0; i < world->data->uepsites; ++i)
                p->uep[i] = (p->ux.s[i][0] < p->ux.s[i][1]) ? '1' : '0';
        }
    }
}

void
check_uep_root (node * p, world_fmt * world)
{
    long i;
    for (i = 0; i < world->data->uepsites; ++i)
    {
        if (p->uep[i] > 1)
        {
            p->uep[i] = (char) RANDINT (0, 1);
        }
        world->oldrootuep[i] = p->uep[i];
    }
    //  memcpy(world->oldrootuep,p->uep,sizeof(long)*world->data->uepsites);
}

boolean
is_success_pseudo_uep (proposal_fmt * proposal)
{
    boolean o, ob, t, tb;
    long uepsites = proposal->world->data->uepsites;
    //if(!proposal->world->in_last_chain && proposal->world->options->uep_last))
    //  return TRUE;
    node *root = proposal->world->root->next->back;

    if (proposal->world->options->uep)
    {
        o = uep_check (root, proposal->origin->uep, uepsites);
        ob = uep_check (root, proposal->oback->uep, uepsites);
        t = uep_check (root, proposal->target->uep, uepsites);
        tb =
            uep_check (root, showtop (crawlback (proposal->target))->uep,
                       uepsites);
        if ((o && t && ob && tb) || (!o && !t && !ob && !tb))
            return TRUE;
        if (o && ob && t && !tb)
            return TRUE;
        if (o && !ob && t && tb)
            return TRUE;
        if (o && !ob && !t && !tb)
            return TRUE;
        if (!o && !ob && t && !tb)
            return TRUE;
    }
    return FALSE;
}

boolean
is_success_pseudo_uepOLD (proposal_fmt * proposal)
{
    boolean o, ob, t;
    //if(!proposal->world->in_last_chain && proposal->world->options->uep_last))
    //  return TRUE;
    node *root = proposal->world->root->next->back;
    if (proposal->world->options->uep)
    {
        o =
            uep_check (root, proposal->origin->uep,
                       proposal->world->data->uepsites);
        ob =
            uep_check (root, proposal->oback->uep,
                       proposal->world->data->uepsites);
        t =
            uep_check (root, proposal->target->uep,
                       proposal->world->data->uepsites);
        if (o && t)
            return TRUE;
        else
        {
            if (o && !t && !ob)
            {
                return TRUE;
            }
            else
            {
                if (!o && !t && !ob)
                    return TRUE;
            }
        }
    }
    return FALSE;
}

boolean
uep_check (node * root, int *uep, long uepsites)
{
    long i;
    for (i = 0; i < uepsites; ++i)
        if (uep[i] == root->uep[i])
            return FALSE;
    return TRUE;
}

boolean
uep_compare (int *target, int *origin, long uepsites)
{
    long i;

    for (i = 0; i < uepsites; ++i)
        if (target[i] != origin[i])
            return FALSE;
    return TRUE;
}


void
store_uep (world_fmt * world)
{
    long i;
    if (world->options->uep)
    {
        if (world->in_last_chain)
        {
            for (i = 0; i < world->data->uepsites; ++i)
            {
                memcpy (world->ueplikestore[world->locus][world->G][i],
                        world->ueplike[i], sizeof (MYREAL) * world->numpop);
                memcpy (world->ueptimestore[world->locus][world->G][i].
                        populations, world->ueptime[i].populations,
                        sizeof (long) * world->ueptime[i].size);
                memcpy (world->ueptimestore[world->locus][world->G][i].ueptime,
                        world->ueptime[i].ueptime,
                        sizeof (MYREAL) * world->ueptime[i].size);
                world->ueptimestore[world->locus][world->G][i].size =
                    world->ueptime[i].size;
            }
        }
    }
}

void
print_uep (FILE * out, world_fmt * world, long copies)
{
    MYREAL *prob, *tyme, *timep;
    long pop;
    FPRINTF (out, "\n\n\nUnique event polymorphism\n");
    FPRINTF (out, "=========================\n\n");
    FPRINTF (out, "Probabilities and time of UEP\n");
    FPRINTF (out, "Ancestral State of UEP is either 0 or 1\n");
    FPRINTF (out, "-------------------------------------------------\n");
    FPRINTF (out, "Populations  Probability   Time(*)       Time(**)\n");
    FPRINTF (out, "-------------------------------------------------\n");
    prob = calc_uep_popprob (world, copies, -1); // don't care if anc uep=1/0
    tyme = calc_uep_timeprob (world, copies, AVERAGE, -1);
    timep = calc_uep_timeprob (world, copies, POPULATION, -1);
    for (pop = 0; pop < world->numpop; ++pop)
    {
        FPRINTF (out, "  %3li        %2.6f      %2.6f      %2.6f\n",
                 pop + 1, prob[pop], tyme[pop], timep[pop]);
    }
    FPRINTF (out, "-------------------------------------------------\n");
    FPRINTF (out, "(*)  Assuming that mutation arose on branch to UEP\n");
    FPRINTF (out, "(**) Mutation arose in specific population.\n");

    FPRINTF (out, "\nAncestral State of UEP is 0\n");
    FPRINTF (out, "-------------------------------------------------\n");
    FPRINTF (out, "Populations  Probability   Time(*)       Time(**)\n");
    FPRINTF (out, "-------------------------------------------------\n");
    prob = calc_uep_popprob (world, copies, 1); // don't care if anc uep=1
    tyme = calc_uep_timeprob (world, copies, AVERAGE, 1);
    timep = calc_uep_timeprob (world, copies, POPULATION, 1);
    for (pop = 0; pop < world->numpop; ++pop)
    {
        FPRINTF (out, "  %3li        %2.6f      %2.6f      %2.6f\n",
                 pop + 1, prob[pop], tyme[pop], timep[pop]);
    }
    FPRINTF (out, "-------------------------------------------------\n");
    FPRINTF (out, "Probability=Prob(is in population | uep_anc=0)\n");
    FPRINTF (out, "(*), (**) see above\n");

    FPRINTF (out, "\nAncestral State of UEP is 1\n");
    FPRINTF (out, "-------------------------------------------------\n");
    FPRINTF (out, "Populations  Probability   Time(*)       Time(**)\n");
    FPRINTF (out, "-------------------------------------------------\n");
    prob = calc_uep_popprob (world, copies, 0); // don't care if anc uep=0
    tyme = calc_uep_timeprob (world, copies, AVERAGE, 0);
    timep = calc_uep_timeprob (world, copies, POPULATION, 0);
    for (pop = 0; pop < world->numpop; ++pop)
    {
        FPRINTF (out, "  %3li        %2.6f      %2.6f      %2.6f\n",
                 pop + 1, prob[pop], tyme[pop], timep[pop]);
    }
    FPRINTF (out, "-------------------------------------------------\n");
    FPRINTF (out, "Probability=Prob(is in population | uep_anc=1)\n");
    FPRINTF (out, "(*), (**) see above\n");
    myfree(prob);
    myfree(tyme);
    myfree(timep);
}

MYREAL *
calc_uep_popprob (world_fmt * world, long copies, long rootstate)
{
    long tree, pop, j, locus;
    MYREAL *result;
    result = (MYREAL *) mycalloc (world->numpop, sizeof (MYREAL));
    for (locus = 0; locus < world->loci; locus++)
    {
        for (tree = 0; tree < world->atl[world->repstop - 1][locus].T - 1;
                tree++)
        {
            for (pop = 0; pop < world->numpop; pop++)
            {
                for (j = 0; j < world->data->uepsites; ++j)
                {
                    if (rootstate != world->ueprootstore[locus][tree][j])
                        result[pop] += world->ueplikestore[locus][tree][j][pop];
                }
            }
        }
    }
    for (pop = 0; pop < world->numpop; pop++)
    {
        result[pop] /= copies;
    }
    return result;
}

MYREAL
ueptime_report (ueptime_fmt * ueptime, int type, long pop)
{
    long i;
    MYREAL summ = 0.0;
    long elem = 0;

    switch (type)
    {
    case AVERAGE:
        return (ueptime->ueptime[ueptime->size - 1] + ueptime->ueptime[0]) / 2.;
    case POPULATION:
        summ = 0.0;
        elem = 0;
        for (i = 0; i < ueptime->size - 1; ++i)
        {
            if (ueptime->populations[i] == pop)
            {
                summ += ueptime->ueptime[i] + ueptime->ueptime[i + 1];
                elem += 2;
            }
        }
        if (elem > 0)
            return summ / (MYREAL) elem;
        else
            return 0.0;
        break;
    }
    return 0.0;
}

MYREAL *
calc_uep_timeprob (world_fmt * world, long copies, int type, long rootstate)
{
    long pop, tree, j, locus;
    MYREAL *result;
    result = (MYREAL *) mycalloc (world->numpop, sizeof (MYREAL));
    for (locus = 0; locus < world->loci; locus++)
    {
        for (tree = 0; tree < world->atl[world->repstop - 1][locus].T - 1;
                tree++)
        {
            for (pop = 0; pop < world->numpop; pop++)
            {
                for (j = 0; j < world->data->uepsites; ++j)
                {
                    if (rootstate != world->ueprootstore[locus][tree][j])
                        result[pop] += world->ueplikestore[locus][tree][j][pop] *
                                       ueptime_report (&world->ueptimestore[locus][tree][j],
                                                       type, pop);
                }
            }
        }
    }
    for (pop = 0; pop < world->numpop; pop++)
    {
        result[pop] /= copies;
    }
    return result;
}

long
normalize_uep (world_fmt * world)
{
    long tree, pop, j, locus;
    MYREAL summ;
    long copies = 0;
    for (locus = 0; locus < world->loci; locus++)
    {
        for (tree = 0; tree < world->atl[world->repstop - 1][locus].T - 1;
                tree++)
        {
            summ = 0.0;
            copies += world->atl[world->repstop - 1][locus].tl[tree].copies;
            for (pop = 0; pop < world->numpop; pop++)
            {
                for (j = 0; j < world->data->uepsites; ++j)
                    summ += world->ueplikestore[locus][tree][j][pop];
            }
            for (pop = 0; pop < world->numpop; pop++)
            {
                for (j = 0; j < world->data->uepsites; ++j)
                    world->ueplikestore[locus][tree][j][pop] *=
                        ((MYREAL) world->atl[world->repstop - 1][locus].tl[tree].
                         copies) / summ;
            }
        }
    }
    return copies;
}

void
analyze_uep (world_fmt * world)
{
    long copies;
    copies = normalize_uep (world);
    print_uep (world->outfile, world, copies);
    print_ancestor (world->outfile, world, copies);
}

void
show_uep_store (world_fmt * world)
{
    long i, j, z, steps, locus;
    FILE *file;
    char filename[20] = "uepout";
    openfile (&file, filename, "w+",  NULL);
    for (locus = 0; locus < world->loci; locus++)
    {
        for (steps = 0; steps < world->G; ++steps)
        {
            for (i = 0; i < world->numpop; ++i)
            {
                for (j = 0; j < world->data->uepsites; ++j)
                    FPRINTF (file, "%f ",
                             world->ueplikestore[locus][steps][j][i]);
            }
            FPRINTF (file, "\n");
        }
        for (steps = 0; steps < world->G; ++steps)
        {
            for (j = 0; j < world->data->uepsites; ++j)
            {
                FPRINTF (file, "%li ",
                         world->ueptimestore[locus][steps][j].size);
                for (z = 0; z < world->ueptimestore[locus][steps][j].size; ++z)
                {
                    FPRINTF (file, "%li ",
                             world->ueptimestore[locus][steps][j].
                             populations[z]);
                }
                FPRINTF (file, "\n");
                FPRINTF (file, "%li ",
                         world->ueptimestore[locus][steps][j].size);
                for (z = 0; z < world->ueptimestore[locus][steps][j].size; ++z)
                {
                    FPRINTF (file, "%f ",
                             world->ueptimestore[locus][steps][j].ueptime[z]);
                }
                FPRINTF (file, "\n");
            }
        }
    }
}

void
print_ancestor (FILE * out, world_fmt * world, long copies)
{
    char uepsite[10];
    long us, pop, half;
    MYREAL subtotal = 0.0;
    MYREAL total = 0;
    long locus;
    MYREAL ancsum = 0.;
    for (locus = 0; locus < world->loci; locus++)
    {
        for (us = 0; us < 2 * world->numpop * world->data->uepsites; us++)
            total += world->uepanc[locus][us];
    }
    FPRINTF (out, "\n\nUEP probabilities at the MRCA\n");
    FPRINTF (out, "--------------------------------------\n");
    FPRINTF (out, "UEP allele    Population   Probability\n");
    FPRINTF (out, "--------------------------------------\n");
    half = world->numpop * world->data->uepsites;
    for (us = 0; us < world->data->uepsites; us++)
    {
        if (world->data->uepsites != 1)
            sprintf (uepsite, "%3li:", us + 1);
        else
            sprintf (uepsite, "%4s", "     ");
        for (pop = 0; pop < world->numpop; pop++)
        {
            ancsum = 0;
            for (locus = 0; locus < world->loci; locus++)
            {
                subtotal += world->uepanc[locus][pop * us + pop];
                ancsum += world->uepanc[locus][pop * us + pop];
            }
            FPRINTF (out, "%4s 0            %3li     %f\n", uepsite, pop + 1,
                     ancsum / total);
        }
        FPRINTF (out, "%4s 0            All     %f\n", uepsite,
                 subtotal / total);
        FPRINTF (out, "--------------------------------------\n");
        subtotal = 0.;
        for (pop = 0; pop < world->numpop; pop++)
        {
            ancsum = 0;
            for (locus = 0; locus < world->loci; locus++)
            {
                subtotal += world->uepanc[locus][half + pop * us + pop];
                ancsum += world->uepanc[locus][half + pop * us + pop];
            }
            FPRINTF (out, "%4s 1            %3li     %f\n", uepsite, pop + 1,
                     ancsum / total);
        }
        FPRINTF (out, "%4s 1            All     %f\n", uepsite,
                 subtotal / total);
        FPRINTF (out, "--------------------------------------\n");
    }
}

void
calc_ueplike (node * p, world_fmt * world, MYREAL **ueplike)
{
    boolean done;
    int temp;
    long i, save_i = 0;
    node *d1, *d2;
    if (p->type != 't')
    {
        if (p->next->back != NULL)
            calc_ueplike (crawlback (p->next), world, ueplike);
        if (p->next->next->back != NULL)
            calc_ueplike (crawlback (p->next->next), world, ueplike);

        if (p->type != 'r')
        {
            d1 = crawlback (p->next);
            d2 = crawlback (p->next->next);
            done = FALSE;
            for (i = 0; i < world->data->uepsites; ++i)
            {
                temp = (d1->uep[i]=='0'  ? 0 : 1 )+ (d2->uep[i]=='0' ? 0 : 1);
                if (temp == 1)
                {
                    if (!done) //shortcut as long we have only one UEP
                    {
		      if(world->options->verbose)
                        printf ("CALCUEP: at <%li> with time %f (%f,%f)\n",p->id,p->tyme,d1->tyme,d2->tyme);
                        collect_ueplike (p, d1, d2, i, world, ueplike[i]);
                        save_i = i;
                    }
                    else
                        memcpy (ueplike[i], ueplike[save_i],
                                sizeof (MYREAL) * world->numpop);
                    done = TRUE;
                }
            }
        }
    }
}

void
fill_ueptime (ueptime_fmt * ueptime, node * p, node * last)
{
    long i = 1;
    long count = 0;
    node *nod;

    ueptime->populations[0] = last->actualpop;
    ueptime->ueptime[0] = last->tyme;
    while ((nod = showtop (last->back)) != p)
    {
        last = nod;
        count++;
    }
    count += 2;
    if (count > ueptime->size)
    {
        ueptime->populations = (long *) myrealloc (ueptime->populations,
                               sizeof (long) * count);
        ueptime->ueptime = (MYREAL *) myrealloc (ueptime->ueptime,
                                               sizeof (MYREAL) * count);
        ueptime->size = count;
    }
    while ((nod = showtop (last->back)) != p)
    {
        ueptime->populations[i] = nod->actualpop;
        ueptime->ueptime[i] = nod->tyme;
        last = nod;
        i++;
    }
    ueptime->populations[i] = nod->actualpop;
    ueptime->ueptime[i] = nod->tyme;
    ueptime->size = i + 1;
}

void
collect_ueplike (node * p, node * d1, node * d2, long uepsite,
                 world_fmt * world, MYREAL *ueplike)
{
    node *nod, *last;
    //  long pop;
    //  MYREAL interval = world->root->next->back->tyme - p->tyme;
    //  printf("collect_ueplike: %f %f %f\n",world->root->next->back->tyme,p->tyme,interval);
    memset (ueplike, 0, sizeof (MYREAL) * world->numpop);
    //  if (d1->uep[uepsite] == 1)
    if (d1->uep[uepsite] != world->root->next->back->uep[uepsite])
        last = d1;
    else
        last = d2;
    fill_ueptime (&world->ueptime[uepsite], p, last);
    // = (p->tyme + last->tyme) / 2.;
    while ((nod = showtop (last->back)) != p)
    {
        //      printf("real: %li> %f %f\n",nod->id,nod->tyme, nod->tyme - last->tyme);
        ueplike[nod->actualpop] += nod->tyme - last->tyme;
        last = nod;
    }
    ueplike[nod->actualpop] += nod->tyme - last->tyme;
    //for(pop=0;pop<world->numpop;++pop)
    //{
    //  if(ueplike[pop]!=0.)
    //    ueplike[pop] *= EXP (world->treelen);
    //}
    //      FPRINTF(stdout,"%f %f %f %li - ",ueplike[pop], ueplike[pop]-interval,
    //      nod->tyme, nod->id);
    //  }
    // FPRINTF(stdout,"real\n");
}


node *
first_uep (node * p, node * nop, long uepsites)
{
    node *bt = showtop (crawlback (p));
    while (uep_compare (p->uep, bt->uep, uepsites) && bt != nop)
    {
        p = bt;
        bt = showtop (crawlback (p));
    }
    return p;
}

node *
first_uep2 (node * p, node * root, long uepsites)
{
    node *pn, *pnn;
    if (p->type != 'r' && p->type != 't')
    {
        pn = p->next;
        pnn = p->next->next;
        p = first_uep2 (crawlback (pn), root, uepsites);
        if (uep_check (root, p->uep, uepsites))
            return p;
        p = first_uep2 (crawlback (pnn), root, uepsites);
        if (uep_check (root, p->uep, uepsites))
            return p;
    }
    return p;
}

MYREAL
sum_ueplike (proposal_fmt * proposal)
{
    long u, i;
    MYREAL like = 1.;
    MYREAL poplike;

    for (u = 0; u < proposal->world->data->uepsites; ++u)
    {
        poplike = 0.;
        for (i = 0; i < proposal->world->numpop; ++i)
        {
            if (proposal->ueplike[u][i] != 0.0)
                poplike += proposal->ueplike[u][i];
        }
        like *= poplike;
    }
    if (proposal->world->options->ueprate > 0.0)
    {
        return LOG (like) -
               (proposal->world->options->ueprate * proposal->treelen);
    }
    else
        return LOG (like);
}


MYREAL
adjust_uep_base (proposal_fmt * proposal, MYREAL interval, MYREAL oldinterval)
{
    long u, pop;
    // MYREAL expT = EXP (-(proposal->treelen-proposal->world->treelen));
    for (u = 0; u < proposal->world->data->uepsites; ++u)
    {
        memcpy (proposal->ueplike[u], proposal->world->ueplike[u],
                sizeof (MYREAL) * proposal->world->numpop);
        //                  printf("adjust_base   ");
        for (pop = 0; pop < proposal->world->numpop; ++pop)
        {
            if (proposal->world->ueplike[u][pop] != 0.0)
                proposal->ueplike[u][pop] -= oldinterval - interval;
            //      proposal->ueplike[u][pop] *= expT;
            //                      printf("%f ", proposal->ueplike[u][pop]);
        }
        //            printf("- %f %f\n",interval, oldinterval);
    }
    return sum_ueplike (proposal);
}

MYREAL
adjust_uep_target (node * first, node * firstb, proposal_fmt * proposal,
                   MYREAL interval, MYREAL oldinterval)
{
    long u, pop;
    node *last, *p;
    MYREAL lasttime, endtime;
    //MYREAL expT = EXP (-(proposal->treelen-proposal->world->treelen));
    boolean o =
        uep_check (proposal->world->root->next->back, proposal->origin->uep,
                   proposal->world->data->uepsites);
    for (u = 0; u < proposal->world->data->uepsites; ++u)
    {
        memset (proposal->ueplike[u], 0,
                sizeof (MYREAL) * proposal->world->numpop);
        last = first;
        p = showtop (first->back);
        if (o)
        {
            lasttime = proposal->time;
            endtime = firstb->tyme;
        }
        else
        {
            lasttime = last->tyme;
            endtime = proposal->time;
        }
        while (p->tyme < lasttime)
        {
            last = p;
            p = showtop (last->back);
        }
        while (p->tyme <= endtime)
        {
            proposal->ueplike[u][p->actualpop] += p->tyme - lasttime;
            last = p;
            lasttime = p->tyme;
            p = showtop (last->back);
        }
        //                  printf("adjust_target ");
        for (pop = 0; pop < proposal->world->numpop; ++pop)
        {
            if (proposal->world->ueplike[u][pop] != 0.0)
                proposal->ueplike[u][pop] -= oldinterval - interval;
            //  proposal->ueplike[u][pop] *= expT;
            //               printf("%f ", proposal->ueplike[u][pop]);
        }
        //                 printf("- %f %f\n",interval, oldinterval);
    }
    return sum_ueplike (proposal);
}

MYREAL
adjust_uep_origin (node * first, node * firstb, proposal_fmt * proposal,
                   MYREAL interval, MYREAL oldinterval)
{
    long u, i, pop;
    MYREAL lasttime;
    //  node *last, *p;
    migr_table_fmt *mt = proposal->migr_table;
    long mtc = proposal->migr_table_counter;
    //MYREAL expT = EXP (-(proposal->treelen-proposal->world->treelen));
    //  printf("@@@@@@@ adjust branch below origin @@@@@@@@@\n");
    for (u = 0; u < proposal->world->data->uepsites; ++u)
    {
        // memcpy(proposal->ueplike[u], proposal->world->ueplike[u],
        //             sizeof(MYREAL)*proposal->world->numpop);
        memset (proposal->ueplike[u], 0,
                sizeof (MYREAL) * proposal->world->numpop);
        lasttime = first->tyme;
        for (i = 0; i < mtc; ++i)
        {
            proposal->ueplike[u][mt[i].to] += mt[i].time - lasttime;
            lasttime = mt[i].time;
        }
        proposal->ueplike[u][proposal->target->actualpop] +=
            proposal->time - lasttime;
        //                  printf("adjust_target ");
        for (pop = 0; pop < proposal->world->numpop; ++pop)
        {
            if (proposal->world->ueplike[u][pop] != 0.0)
                proposal->ueplike[u][pop] -= oldinterval - interval;
            //  proposal->ueplike[u][pop] *= expT;
            //               printf("%f ", proposal->ueplike[u][pop]);
        }
        //        printf("- %f %f\n",interval, oldinterval);
    }
    return sum_ueplike (proposal);
}

MYREAL
adjust_uep_oback (node * first, node * firstb, proposal_fmt * proposal,
                  MYREAL interval, MYREAL oldinterval)
{
    long u, pop;
    MYREAL lasttime = 0.;
    node *last, *p;
    //MYREAL expT = EXP (-(proposal->treelen-proposal->world->treelen));
    if (first->tyme < proposal->time) // target == osister
        first = proposal->oback;
    for (u = 0; u < proposal->world->data->uepsites; ++u)
    {
        //memcpy(proposal->ueplike[u], proposal->world->ueplike[u],
        //             sizeof(MYREAL)*proposal->world->numpop);
        memset (proposal->ueplike[u], 0,
                sizeof (MYREAL) * proposal->world->numpop);
        last = first;
        p = showtop (first->back);
        if (last != proposal->oback)
            lasttime = last->tyme;
        else
            lasttime = proposal->time;
        while (p != firstb)
        {
            proposal->ueplike[u][p->actualpop] += p->tyme - lasttime;
            last = p;
            p = showtop (last->back);
            if (p == proposal->oback)
                p = showtop (p->back);
            lasttime = last->tyme;
        }
        proposal->ueplike[u][p->actualpop] += p->tyme - lasttime;
        //       printf(" adjust_oback ");
        for (pop = 0; pop < proposal->world->numpop; ++pop)
        {
            if (proposal->world->ueplike[u][pop] != 0.0)
                proposal->ueplike[u][pop] -= oldinterval - interval;
            //  proposal->ueplike[u][pop] *= expT;
            //                printf("%f ", proposal->ueplike[u][pop]);
        }
        //                  printf("- %f %f\n",interval, oldinterval);
    }
    return sum_ueplike (proposal);
}

MYREAL
adjust_uep_osister (node * first, node * firstb, proposal_fmt * proposal,
                    MYREAL interval, MYREAL oldinterval)
{
    long u, pop;
    MYREAL lasttime = 0.;
    node *last, *p;
    //MYREAL expT = EXP (-(proposal->treelen-proposal->world->treelen));
    for (u = 0; u < proposal->world->data->uepsites; ++u)
    {
        memset (proposal->ueplike[u], 0,
                sizeof (MYREAL) * proposal->world->numpop);
        last = first;
        p = showtop (first->back);
        lasttime = last->tyme;
        while (p != firstb)
        {
            proposal->ueplike[u][p->actualpop] += p->tyme - lasttime;
            last = p;
            p = showtop (last->back);
            lasttime = last->tyme;
        }
        proposal->ueplike[u][p->actualpop] += p->tyme - lasttime;
        //            printf("adjust_osister ");
        for (pop = 0; pop < proposal->world->numpop; ++pop)
        {
            if (proposal->world->ueplike[u][pop] != 0.0)
                proposal->ueplike[u][pop] -= oldinterval - interval;
            //  proposal->ueplike[u][pop] *= expT;
            //        printf("%f ", proposal->ueplike[u][pop]);
        }
        // printf("- %f %f\n",interval, oldinterval);
    }
    return sum_ueplike (proposal);
}

MYREAL
adjust_uep_osistertarget (node * first, node * firstb,
                          proposal_fmt * proposal, MYREAL interval,
                          MYREAL oldinterval)
{
    long u, pop;
    MYREAL lasttime = 0.;
    node *last, *p;
    //MYREAL expT = EXP (-(proposal->treelen-proposal->world->treelen));
    for (u = 0; u < proposal->world->data->uepsites; ++u)
    {
        memset (proposal->ueplike[u], 0,
                sizeof (MYREAL) * proposal->world->numpop);
        last = first;
        p = showtop (first->back);
        lasttime = last->tyme;
        while (p != firstb)
        {
            proposal->ueplike[u][p->actualpop] += p->tyme - lasttime;
            last = p;
            p = showtop (last->back);
            if (last != proposal->oback)
                lasttime = last->tyme;
            else
                lasttime = proposal->time;
        }
        proposal->ueplike[u][p->actualpop] += p->tyme - lasttime;
        printf ("adjust_osistertarget ");
        for (pop = 0; pop < proposal->world->numpop; ++pop)
        {
            if (proposal->world->ueplike[u][pop] != 0.0)
                proposal->ueplike[u][pop] -= oldinterval - interval;
            //  proposal->ueplike[u][pop] *= expT;
            printf ("%f ", proposal->ueplike[u][pop]);
        }
        printf ("- %f %f\n", interval, oldinterval);
    }
    return sum_ueplike (proposal);
}

MYREAL
adjust_uep_osistertarget2 (node * first, node * firstb,
                           proposal_fmt * proposal, MYREAL interval,
                           MYREAL oldinterval)
{
    long u, pop, actualpop;
    MYREAL lasttime = 0., endtime;
    node *last, *p;
    // MYREAL expT = EXP (-(proposal->treelen-proposal->world->treelen));
    for (u = 0; u < proposal->world->data->uepsites; ++u)
    {
        memset (proposal->ueplike[u], 0,
                sizeof (MYREAL) * proposal->world->numpop);
        last = first;
        p = showtop (first->back);
        actualpop = p->actualpop;
        lasttime = last->tyme;
        endtime = proposal->time;
        while (p != proposal->realtarget)
        {
            proposal->ueplike[u][p->actualpop] += p->tyme - last->tyme;
            last = p;
            p = showtop (last->back);
        }
        proposal->ueplike[u][p->actualpop] += p->tyme - last->tyme;
        proposal->ueplike[u][p->pop] += endtime - p->tyme;
        //                  printf("adjust_osistertarget ");
        for (pop = 0; pop < proposal->world->numpop; ++pop)
        {
            if (proposal->world->ueplike[u][pop] != 0.0)
                proposal->ueplike[u][pop] -= oldinterval - interval;
            //  proposal->ueplike[u][pop] *= expT;
            //                printf("%f ", proposal->ueplike[u][pop]);
        }
        //      printf("- %f %f\n",interval, oldinterval);
    }
    return sum_ueplike (proposal);
}

MYREAL
pseudo_ueplikelihood (world_fmt * world, proposal_fmt * proposal)
{
    MYREAL interval = 0., oldinterval = 0.; //  node *first, *firstb;
    node *root = proposal->root->next->back;
    long uepsites = proposal->world->data->uepsites;

    boolean o, t, ob, tb;

    node *firstb;
    node *first;
    //memcpy(root->uep, world->oldrootuep, sizeof(long)*uepsites);
    o = uep_check (root, proposal->origin->uep, uepsites);
    t = uep_check (root, proposal->target->uep, uepsites);
    ob = uep_check (root, proposal->oback->uep, uepsites);
    tb = uep_check (root, showtop (crawlback (proposal->target))->uep,
                    uepsites);
    first = first_uep2 (root, root, uepsites);
    firstb =
        showtop (crawlback
                 (first = first_uep (first, proposal->world->root, uepsites)));
    proposal->firstuep = first;
    if (proposal->world->options->ueprate > 0.0)
        proposal->treelen =
            calc_pseudotreelength (proposal, proposal->world->treelen);
    /*  if(root==proposal->target)
       {
       oldinterval = root->tyme - firstb->tyme;
       interval = proposal->time - firstb->tyme;
       }
       else 
       {
       if(root==proposal->oback)
       {
       newroot = crawlback(proposal->oback->next);
       if(newroot==proposal->origin)
       newroot = crawlback(proposal->oback->next->next);
       oldinterval = root->tyme - firstb->tyme;
       if(proposal->time > newroot->tyme)
       interval = proposal->time - firstb->tyme;
       else
       interval = newroot->tyme - firstb->tyme;
       }
       else
       {
       oldinterval= root->tyme - firstb->tyme;
       interval = root->tyme - firstb->tyme;
       }
       } */
    if (o && t && ob && tb)
    {
        if (proposal->oback == first)
        {

            return adjust_uep_oback (proposal->osister, firstb, proposal,
                                     interval, oldinterval);
        }
        else
        {
            //      FPRINTF(stdout,"+");
            return adjust_uep_base (proposal, interval, oldinterval);
        }
    }
    else
    {
        if (!o && !t && !ob && !tb)
        {
            if (proposal->osister == first)
            {
                //      FPRINTF(stdout,"s");
                if (proposal->target == first)
                    return adjust_uep_osistertarget (proposal->osister, firstb,
                                                     proposal, interval,
                                                     oldinterval);
                else
                {
                    if (proposal->target == proposal->oback)
                        return adjust_uep_osistertarget2 (proposal->osister,
                                                          firstb, proposal,
                                                          interval, oldinterval);
                    else
                        return adjust_uep_osister (proposal->osister,
                                                   showtop (crawlback (firstb)),
                                                   proposal, interval,
                                                   oldinterval);
                }
            }
            else
            {
                //      FPRINTF(stdout,"-");
                return adjust_uep_base (proposal, interval, oldinterval);
            }
        }
        else
        {
            if (proposal->target == first)
            {
                //    if(proposal->oback == firstb)
                // adjust_uep_target_oback(first,firstb,proposal,
                //                      interval,oldinterval);
                //else
                return adjust_uep_target (first, firstb, proposal, interval,
                                          oldinterval);
            }
            else
            {
                if (proposal->oback == firstb && proposal->origin != first)
                    return adjust_uep_target (first, showtop (crawlback (firstb)),
                                              proposal, interval, oldinterval);
                else
                {
                    if (proposal->origin == first)
                        return adjust_uep_origin (proposal->origin,
                                                  proposal->oback, proposal,
                                                  interval, oldinterval);
                    //  else
                    //    error ("do not know what to do with oback!=firstb");
                }
            }
        }
    }
    return -MYREAL_MAX;  //never come here
}

MYREAL
ueplikelihood (world_fmt * world)
{
    long u, i;
    MYREAL like = 1.;
    MYREAL poplike = 0;
    calc_ueplike (world->root->next->back, world, world->ueplike);
    //  printf("update_likeli ");
    for (u = 0; u < world->data->uepsites; ++u)
    {
        poplike = 0.;
        for (i = 0; i < world->numpop; ++i)
        {
            poplike += world->ueplike[u][i];
            //      printf("%f ", world->ueplike[u][i]);
        }
        like *= poplike;
        //      printf("\n");
    }
    if (world->options->ueprate > 0.0)
    {
        return LOG (like) - (world->options->ueprate * world->treelen);
    }
    else
        return LOG (like);
}

void
show_uep_time (node * p, world_fmt * world)
{
    if (p->type != 'r')
    {
        if (p->type == 'm' || p->type == 'd')
        {
            printf ("%li> - %f %f\n", p->id, p->tyme,
                    world->root->next->back->tyme - p->tyme);
            show_uep_time (showtop (p->back), world);
        }
        else
        {
            printf ("%li> %c %f %f\n", p->id, p->uep[0], p->tyme,
                    world->root->next->back->tyme - p->tyme);
            show_uep_time (showtop (p->back), world);
        }
    }
}
void
show_uep_time2 (node * p, world_fmt * world)
{
    if (p->type == 't')
        return;
    if (p->type != 'r')
    {
        if (p->type == 'm' || p->type == 'd')
        {
            show_uep_time2 (showtop (p->next->back), world);
            printf ("%li> - %f %f\n", p->id, p->tyme,
                    world->root->next->back->tyme - p->tyme);
        }
        else
        {
            show_uep_time2 (showtop (p->next->back), world);
            show_uep_time2 (showtop (p->next->next->back), world);
            printf ("%li> %c %f %f\n", p->id, p->uep[0], p->tyme,
                    world->root->next->back->tyme - p->tyme);

        }
    }
}

void
update_uepanc (world_fmt * world)
{
    long i, half, pop, popi;
    half = world->data->uepsites * world->numpop;
    for (i = 0; i < world->data->uepsites; ++i)
    {
        pop = world->root->next->back->actualpop;
        popi = pop * i + pop;
        world->ueprootstore[world->locus][world->G][i] =
            world->root->next->back->uep[i];
        if (world->root->next->back->uep[i] == 0)
            world->uepanc[world->locus][popi] += 1;
        else
            world->uepanc[world->locus][half + popi] += 1;
    }
}

void
copy_uepx (proposal_fmt * proposal, ueparray_fmt xx1, ueparray_fmt xx2)
{
    memcpy (xx1.s, xx2.s, sizeof (pair)*proposal->world->data->uepsites);
}

#endif
typedef double uep_empty;





