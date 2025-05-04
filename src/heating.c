/* heating scheme with threads if available
 
 
part of migrate 
Peter Beerli 2000
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2008 Peter Beerli, Tallahassee FL
 
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
 
$Id: heating.c 2067 2012-07-27 20:59:32Z beerli $
 
*/
/*! \file heating.c 

Heating handler for threaded heating

*/

#include "migration.h"
#include "mcmc.h"
#include "sighandler.h"
#include "random.h"
#ifdef mPI
#include "migrate_mpi.h"
#endif

extern int myID; //identifier for nodes in MPI context

// functions 
void adjust_temperatures(world_fmt ** universe, long hchains, long step, long steps);
void adjust_temperatures_bounded(world_fmt ** universe, long hchains, long step, long steps);
void reset_weight(world_fmt *world);

#ifdef PTHREADS   /* ==================threaded version======== */

#include <pthread.h>
#include "heating.h"
#include "bayes.h"

void fill_tpool (tpool_t tpool, world_fmt ** universe, int universe_size);
void tpool_init (tpool_t * tpoolp, int num_worker_threads, int max_queue_size,
                 int do_not_block_when_full);
void thread_tree_update (void *thisworld);
int tpool_destroy (tpool_t tpool, int finish);
void tpool_thread (tpool_t tpool);
int tpool_synchronize (tpool_t tpool, int finish);
void allocate_thread_workp(long n, tpool_work_t ** workpointers);

extern void run_one_update(world_fmt *world); //defined in main.c

tpool_work_t **thread_work_pointers;

void
tpool_init (tpool_t * tpoolp, int num_worker_threads, int max_queue_size,
            int do_not_block_when_full)
{
    int i, rtn;
    tpool_t tpool;

    tpool = (tpool_t) mymalloc (sizeof (struct _tpool_t));

    tpool->num_threads = num_worker_threads;
    tpool->max_queue_size = max_queue_size;
    tpool->do_not_block_when_full = do_not_block_when_full;
    tpool->threads =
        (pthread_t *) mymalloc (sizeof (pthread_t) * num_worker_threads);

    tpool->cur_queue_size = 0;
    tpool->queue_head = NULL;
    tpool->queue_tail = NULL;
    tpool->queue_closed = 0;
    tpool->shutdown = 0;
    tpool->done = 0;
    if ((rtn = pthread_mutex_init (&(tpool->random_lock), NULL)) != 0)
        error ("pthread_mutex_init for tpool->random_lock failed");
    if ((rtn = pthread_mutex_init (&(tpool->queue_lock), NULL)) != 0)
        error ("pthread_mutex_init for tpool->queue_lock failed");
    if ((rtn = pthread_cond_init (&(tpool->queue_done), NULL)) != 0)
        error ("pthread_cond_init for tpool->queue_done failed");
    if ((rtn = pthread_cond_init (&(tpool->queue_not_empty), NULL)) != 0)
        error ("pthread_cond_init for tpool->queue_not_empty failed");
    if ((rtn = pthread_cond_init (&(tpool->queue_not_full), NULL)) != 0)
        error ("pthread_cond_init for tpool->queue_not_full failed");
    if ((rtn = pthread_cond_init (&(tpool->queue_empty), NULL)) != 0)
        error ("pthread_cond_init for tpool->queue_empty failed");

    for (i = 0; i < tpool->num_threads; i++)
    {
        if ((rtn =
                    pthread_create (&(tpool->threads[i]), NULL, (void *) tpool_thread,
                                    (void *) tpool)) != 0)
            error ("Thread pool creation failed with pthread_create");
    }
    *tpoolp = tpool;
    thread_work_pointers = (tpool_work_t **) mycalloc(num_worker_threads, sizeof(tpool_work_t *));
    allocate_thread_workp(num_worker_threads, thread_work_pointers);
}


void allocate_thread_workp(long n, tpool_work_t ** workpointers)
{
    // memory savings?
  long i;
  for(i=0;i<n;i++)
    (workpointers)[i] = (tpool_work_t *) mymalloc (sizeof (tpool_work_t));
}

int
tpool_add_work (tpool_t tpool, void *routine, void *arg, long z)
{
    tpool_work_t *workp;
    pthread_mutex_lock (&tpool->queue_lock);

    if ((tpool->cur_queue_size == tpool->max_queue_size)
            && tpool->do_not_block_when_full)
    {
        pthread_mutex_unlock (&tpool->queue_lock);
        return -1;
    }
    while ((tpool->cur_queue_size == tpool->max_queue_size)
            && (!tpool->shutdown || tpool->queue_closed))
    {
        pthread_cond_wait (&tpool->queue_not_full, &tpool->queue_lock);
    }

    if (tpool->shutdown || tpool->queue_closed)
    {
        pthread_mutex_unlock (&tpool->queue_lock);
        return -1;
    }

    //    workp = (tpool_work_t *) mymalloc (sizeof (tpool_work_t));
    workp = thread_work_pointers[z];
    workp->routine = routine;
    workp->arg = arg;
    workp->next = NULL;
    if (tpool->cur_queue_size == 0)
    {
        tpool->queue_tail = tpool->queue_head = workp;
        pthread_cond_broadcast (&tpool->queue_not_empty);
    }
    else
    {
        (tpool->queue_tail)->next = workp;
        tpool->queue_tail = workp;
    }
    tpool->cur_queue_size++;
    pthread_mutex_unlock (&tpool->queue_lock);
    return 1;
}

/// \brief put work for heated chains on to the heated-chain-stack
///
/// queue the work needed to be done into the stack of heated chains
/// if there is an error the routien will abort with an 
/// error "filling heating queue failed"
void fill_tpool (tpool_t tpool, world_fmt ** universe, int universe_size)
{
    int i;
    for (i = 0; i < universe_size; i++)
    {
      if (!tpool_add_work (tpool, thread_tree_update, (void *) universe[i], i))
            error ("Filling heating queue failed");
    }
}

void
tpool_thread (tpool_t tpool)
{
    tpool_work_t *myworkp;

    for (;;)
    {
        pthread_mutex_lock (&(tpool->queue_lock));
        while ((tpool->cur_queue_size == 0) && (!tpool->shutdown))
        {
            pthread_cond_wait (&(tpool->queue_not_empty), &(tpool->queue_lock));
        }
        if (tpool->shutdown)
        {
            pthread_mutex_unlock (&tpool->queue_lock);
            pthread_exit (NULL);
        }
        myworkp = tpool->queue_head;
        tpool->cur_queue_size--;
        if (tpool->cur_queue_size == 0)
            tpool->queue_head = tpool->queue_tail = NULL;
        else
            tpool->queue_head = myworkp->next;

        if ((!tpool->do_not_block_when_full)
                && (tpool->cur_queue_size == (tpool->max_queue_size - 1)))
            pthread_cond_broadcast (&(tpool->queue_not_full));

        if (tpool->cur_queue_size == 0)
            pthread_cond_signal (&(tpool->queue_empty));

        pthread_mutex_unlock (&tpool->queue_lock);
        (*(myworkp->routine)) (myworkp->arg);
        pthread_mutex_lock (&tpool->queue_lock);
        tpool->done++;
        //      printf("-%i",tpool->done);
        pthread_cond_signal (&tpool->queue_done);
        pthread_mutex_unlock (&tpool->queue_lock);
        //myfree(myworkp);
    }
}

/// \brief updates the heated tree and records acceptance
///
/// updates the tree during threaded heating scheme 
void thread_tree_update (void *thisworld)
{
    long i;
    world_fmt * world = (world_fmt *) thisworld;
    // the version commented out would increase the runtime by the interval
    // using the heating interval in heated_swap() is not as efficient in thread but would
    // do the right thing, but just not checking that often.
    //    for (i = 0; i < world->options->heating_interval; i++)
    //{
		run_one_update(thisworld);
		//  }
}

///
/// synchronizes the heated chains and controls the locking of chains
int
tpool_synchronize (tpool_t tpool, int fullnumber)
{
    int rtn;
    //  printf("\nSync ");
    if ((rtn = pthread_mutex_lock (&(tpool->queue_lock))) != 0)
        error ("pthread_mutex_lock failed in tpool_synchronize()");
    //  printf("{%i} ",tpool->done);
    while (tpool->done < fullnumber)
    {
        //      printf("(%i [%i]) ",tpool->done,   tpool->cur_queue_size);
        if ((rtn =
                    pthread_cond_wait (&(tpool->queue_done),
                                       &(tpool->queue_lock))) != 0)
            error ("pthread_cond_wait failed in tpool_synchronize()");
    }
    tpool->done = 0;

    if ((rtn = pthread_mutex_unlock (&(tpool->queue_lock))) != 0)
        error ("pthread_mutex_unlock failed in tpool_destroy()");
    return 0;
}

int
tpool_destroy (tpool_t tpool, int finish)
{
    int i, rtn;
    tpool_work_t *cur_nodep;

    if ((rtn = pthread_mutex_lock (&(tpool->queue_lock))) != 0)
        error ("pthread_mutex_lock failed in tpool_destroy()");

    if (tpool->queue_closed || tpool->shutdown)
    {
        if ((rtn = pthread_mutex_unlock (&(tpool->queue_lock))) != 0)
            error ("pthread_mutex_unlock failed in tpool_destroy()");
        return 0;
    }

    tpool->queue_closed = 1;

    if (finish == 1)
    {
        while (tpool->cur_queue_size != 0)
        {
            if ((rtn =
                        pthread_cond_wait (&(tpool->queue_empty),
                                           &(tpool->queue_lock))) != 0)
                error ("pthread_cond_wait failed in tpool_destroy()");
        }
    }
    tpool->shutdown = 1;

    if ((rtn = pthread_mutex_unlock (&(tpool->queue_lock))) != 0)
        error ("pthread_mutex_unlock failed in tpool_destroy()");

    if ((rtn = pthread_cond_broadcast (&(tpool->queue_not_empty))) != 0)
        error ("pthread_cond_broadcast failed for queue_not_empty()");
    if ((rtn = pthread_cond_broadcast (&(tpool->queue_not_full))) != 0)
        error ("pthread_cond_broadcast failed for queue_not_full");

    for (i = 0; i < tpool->num_threads; i++)
    {
        if ((rtn = pthread_join (tpool->threads[i], NULL)) != 0)
            error ("pthread_join failed in tpool_destroy()");
    }

    myfree(tpool->threads);
    while (tpool->queue_head != NULL)
    {
        cur_nodep = tpool->queue_head->next;
        tpool->queue_head = tpool->queue_head->next;
        myfree(cur_nodep);
    }
    myfree(tpool);
    return 0;
}

#else /* ==================NOT threaded version======== */



#endif

#define HEATCHECKINTERVAL 1000
#define HEATSWAPLOW 0
#define HEATSWAPHIGH 10
void adjust_temperatures(world_fmt ** universe, long hchains, long step, long steps)
{
    long i;
    const MYREAL bigger = 1.1;
    const MYREAL smaller = 0.9;
    const MYREAL corrsum = steps / HEATCHECKINTERVAL;
    MYREAL *delta;
    if(step == 0)
    {
        for(i=0; i< hchains; i++)
        {
            universe[i]->treeswapcount=0;
            if(steps < HEATCHECKINTERVAL)
                universe[i]->averageheat  = 1./universe[i]->heat;
            else
                universe[i]->averageheat  = 0.0;
        }
    }
    else
    {
        if ( (step % HEATCHECKINTERVAL ) == 0)
        {
            universe[0]->averageheat= 1.0;
            delta = (MYREAL *) mycalloc(hchains,sizeof(MYREAL));
	    // FPRINTF(stdout,"\n%f %li\n",1./universe[0]->heat,universe[0]->treeswapcount);
            for(i=1; i< hchains; i++)
            {
                universe[i]->averageheat += (1./(universe[i]->heat * corrsum ));
		//FPRINTF(stdout,"%f %li\n",1./universe[i]->heat,universe[i]->treeswapcount);
                delta[i-1] = 1./universe[i]->heat - 1./universe[i-1]->heat;
		if(delta[i-1] < EPSILON)
		  delta[i-1] = EPSILON;
            }
            for(i=1; i< hchains; i++)
            {
                if(universe[i-1]->treeswapcount <= HEATSWAPLOW)
                    delta[i-1] *= smaller;
                else
                {
                    if(delta[i-1]<1000 && universe[i-1]->treeswapcount >= HEATSWAPHIGH)
			  delta[i-1] *= bigger;
                }
            }
            universe[0]->heat = 1.;
            universe[0]->treeswapcount=0;
	    for(i=1; i< hchains; i++)
	      {
		universe[i]->heat = 1./(1./universe[i-1]->heat + delta[i-1]);
		universe[i]->treeswapcount=0;
	      }
            myfree(delta);
        }
    }
}
///
/// adjust temperatures using a lower and upper bound (temperature=1 and temperature=highest)
void adjust_temperatures_bounded(world_fmt ** universe, long hchains, long step, long steps)
{
  //const MYREAL corrsum = steps / HEATCHECKINTERVAL;
    long i;
    long deltasum;
    MYREAL negheat = 0.0;
    MYREAL *delta;
    if(step == 0)
    {
        for(i=0; i< hchains; i++)
        {
            universe[i]->treeswapcount=0;
            if(steps < HEATCHECKINTERVAL)
                universe[i]->averageheat  = 1./universe[i]->heat;
            else
                universe[i]->averageheat  = 0.0;
        }
    }
    else
    {
        if ( (step % HEATCHECKINTERVAL ) == 0)
        {
            universe[0]->averageheat= 1.0;
            delta = (MYREAL *) mycalloc(hchains,sizeof(MYREAL));
#ifdef DEBUG
	    fprintf(stdout,"\n%i> chain 1: %f %f %li\n",myID, 1./universe[0]->heat,universe[0]->averageheat, universe[0]->treeswapcount);
#endif
	    deltasum = 0;
            for(i=1; i< hchains; i++)
            {
	      universe[i]->averageheat += HEATCHECKINTERVAL * (1./universe[i]->heat - universe[i]->averageheat) / step;
#ifdef DEBUG
	      fprintf(stdout,"%i> chain %li: %f %f %li (step=%li (%f))\n", myID, i+1, 1./universe[i]->heat,universe[i]->averageheat, universe[i]->treeswapcount, step, (MYREAL) HEATCHECKINTERVAL/step);
#endif
                delta[i-1] = (MYREAL) universe[i]->treeswapcount;
		if(delta[i-1] < 1.0)
		  delta[i-1] = 1.0;
		deltasum += (long) delta[i-1];
            }
	    
            universe[0]->heat = 1.;
	    universe[0]->treeswapcount=0;
            for(i=1; i < hchains-1; i++)
            {
	      negheat += delta[i-1]/deltasum;
	      universe[i]->heat = 1. - negheat;
	      universe[i]->treeswapcount=0;
	    }
            myfree(delta);
        }
    }
}


void reset_weight(world_fmt *world)
{
  long i;
  seqmodel_fmt *seq = world->data->seq[0]; 
  long endsite = seq->endsite;
  long heated_end = (long) (endsite - endsite * 1./world->heat);
  //MYREAL v;
  memcpy(seq->aliasweight, seq->savealiasweight, sizeof(long)* (size_t) endsite);
  for(i=0;i<heated_end; i++)
    {
    //xcode   v = RANDINT(0,endsite); 
      seq->aliasweight[0] = 0;
    }
}
