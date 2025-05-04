/*! \File migrate_mpi.c */
/* MPI parts for migrate
   started November 2000, Seattle
   Peter Beerli beerli@fsu.edu
 
   
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2014 Peter Beerli, Tallahassee FL
 
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
 
 
$Id: migrate_mpi.c 2170 2013-09-19 12:08:27Z beerli $
*/
#include "migration.h"
#ifdef MPI
#include "tools.h"
#include "assignment.h"
#include "sighandler.h"
#include "migrate_mpi.h"
#include "pretty.h"
#include "options.h"
#include "tree.h"
#include "world.h"
#include "data.h"
#include "laguerre.h"
#include "random.h"
#ifdef UEP
#include "uep.h"
#endif
#include "bayes.h"

#include "haplotype.h"
#include "mutationmodel.h"

#ifndef WINDOWS
#include <unistd.h>
#endif

/*should go into profile.h*/
#define GRIDSIZE 9

#ifdef PRETTY
extern double page_height;
#endif

extern int numcpu;

extern const MPI_Datatype mpisizeof;

extern void run_replicate(long locus,
                          long replicate,
                          world_fmt **universe,
                          option_fmt *options,
                          data_fmt *data, 
                          tpool_t * heating_pool,
                          int usize,
                          long *treefilepos,
                          long *Gmax);
extern void run_locus (world_fmt ** universe, int usize,
                       option_fmt * options, data_fmt * data,
                       tpool_t * heating_pool, long maxreplicate,
                       long locus, long *treefilepos, long *Gmax);

#ifndef HAS_INDIX
#define INDIX(a,b,c) ((a)*(b)+(c))
#define HAS_INDIX
#endif

void mpi_run_locus(world_fmt ** universe, int usize, option_fmt * options,
                   data_fmt * data, tpool_t * heating_pool, long maxreplicate,
                   long locus, long *treefilepos, long *Gmax);
void mpi_runreplicates_worker (world_fmt ** universe, int usize,
                               option_fmt * options, data_fmt * data,
                               tpool_t * heating_pool,
                               long *treefilepos, long *Gmax);
long pack_databuffer (data_fmt * data, option_fmt * options);
void unpack_databuffer (data_fmt * data, option_fmt * options, world_fmt *world);
void pack_allele_data (char **buffer, long *bufsize, long *allocbufsize, data_fmt * data,
                       long pop, long ind);
void pack_sequence_data (char **buffer, long *bufsize, long *allocbufsize, data_fmt * data,
                         long pop, long ind, long locus);

void mpi_resultsmaster (long sendtype, world_fmt * world,
                         long maxrep,
                         void (*unpack) (char *buffer, world_fmt * world,
                                         long locus, long maxrep,
                                         long numpop));

void mpi_results_worker (long bufs, world_fmt * world,
                         long maxrep,
                         long (*pack) (MYREAL **buffer, world_fmt * world,
                                       long locus, long maxrep, long numpop));
void assignloci_worker (world_fmt * world, option_fmt *options, long *Gmax);
void assign_worker_cleanup (void);
void swap_atl (long from, long to, world_fmt * world);
long pack_quantile (char **buffer, quantile_fmt quant, long n);
void unpack_quantile (char *buffer, quantile_fmt quant, long n);
long pack_failed_percentiles (char **buffer, boolean *failed, long n);
void unpack_failed_percentiles (char *buffer, boolean *failed, long n);

void handle_message(char *rawmessage,int sender, world_fmt *world);
void handle_mdim(float *values,long n, int sender, world_fmt * world);
void handle_burnin_message(char *rawmessage,int sender, world_fmt * world);
void handle_dataondemand(int sender,int tag,char *tempstr, world_fmt *world, option_fmt * options, data_fmt *data);
void set_filehandle(char *message, world_fmt *world,void **file, long *msgstart);

void mpi_receive_replicate( int sender, int tag, long locus, long replicate, world_fmt * world);

long unpack_single_bayes_buffer(MYREAL *buffer, bayes_fmt * bayes, world_fmt * world,long locus);
long pack_single_bayes_buffer(MYREAL **buffer, bayes_fmt *bayes, world_fmt *world,long locus);
long pack_single_bayes_buffer_part(MYREAL **buffer, bayes_fmt *bayes, world_fmt *world,long locus);

void unpack_hist_bayes_buffer(MYREAL *buffer, bayes_fmt *bayes, world_fmt *world, long locus);
long pack_hist_bayes_buffer(MYREAL **buffer, bayeshistogram_fmt *hist, world_fmt * world, long startposition);

long unpack_BF_buffer(MYREAL *buffer, long start, long locus, world_fmt * world);
long unpack_ess_buffer(MYREAL *buffer, long start, world_fmt *world);
long pack_BF_buffer(MYREAL **buffer, long start, long locus, world_fmt * world);
long pack_ess_buffer(MYREAL **buffer, long start, world_fmt *world);

long unpack_hyper_buffer(MYREAL *buffer, long start, world_fmt *world);
long pack_hyper_buffer(MYREAL **buffer, long start, world_fmt *world);

void unpack_sumfile_buffer (MYREAL *buffer, world_fmt * world,
                            long locus, long maxrep, long numpop);
void unpack_single_sumfile_buffer (MYREAL *buffer, timearchive_fmt **ta, world_fmt *world,
                                   long locus, long replicate, long numpop, long *startz);
long pack_sumfile_buffer (MYREAL **buffer, world_fmt * world,
                     long locus, long maxrep, long numpop);

long pack_single_sumfile_buffer(MYREAL **buffer, long z, long *allocbufsize, world_fmt * world,
                                long locus, long replicate, long numpop);

/// pack the assign buffer
long pack_seqerror_buffer(MYREAL **buffer, world_fmt * world,
			  long locus, long maxrep, long numpop);
void unpack_seqerror_buffer(MYREAL *buffer, world_fmt * world,
			  long locus, long maxrep, long numpop);


void mpi_send_replicate(int sender, long locus, long replicate, world_fmt * world);
long  mpi_send_stop_mcmc_lociworker(long numcpu, long loci);
long  mpi_send_stop_mcmc_replicateworker(long numcpu, long loci);
void mpi_send_stop_tag (int worker, world_fmt * world);
long  mpi_send_stop_mcmc_worker(long numcpu, long loci, MPI_Comm *comm, MPI_Request *irequests, MPI_Status *istatus, long id);
void mpi_send_stop_tag (int worker, world_fmt * world);
void send_receive_bayes_params(world_fmt *world, long locus);
void handle_replication(int sender,int tag,char *tempstr, world_fmt *world);
boolean in_mpistack(int sender, world_fmt *world);

void mpi_fprintf(FILE *file, const char *fmt, ...);
void mpi_fprintf2(FILE *file, long filesize, const char *fmt, ...);


#ifdef MPICHECK
#include <sys/resource.h>
void set_memory_limit(rlim_t softsize,rlim_t maxsize);
void check_memory_limit();
#endif

#ifdef PARALIO
boolean my_write_error;
#endif


////////////////////////////////////////////////////////////////////////////////
///
/// Controls all loci in the MPI implementation, uses a load balancing scheme and 
/// distributes the work on the waiting nodes
///
void
mpi_runloci_master (long loci, int *who, world_fmt *world, option_fmt * options, 
		    data_fmt *data, boolean options_readsum, boolean menu)
{
    int tag;
    int ll;
    int sender       = 0;
    boolean done     = FALSE;
    long alldone     = 0;
    long locusdone   = -1;
    long numsent     = 0;
    long tempstrsize = LINESIZE;
    long floattempstrsize = LINESIZE;
    long nbase       = loci + 1;
    long minnodes    = MIN((long) numcpu-1, (long) nbase-1);
    long locus;
    long newsize;
    long *twolongs;
    char *savetempstr;
    char *tempstr; 
    char *leadstr;
    float *temp;   
    MPI_Status status;
    MPI_Status *istatus;
    MPI_Request *irequests;

    irequests  = (MPI_Request *) mycalloc(minnodes, sizeof(MPI_Request));
    istatus    = (MPI_Status *) mycalloc(minnodes, sizeof(MPI_Status));
    twolongs   = (long *) mycalloc(TWO, sizeof(long));
    savetempstr= (char *) mycalloc(tempstrsize+1, sizeof(char));
    tempstr    = savetempstr;
    leadstr    = (char *) mycalloc(SMALLBUFSIZE, sizeof(char));
    temp       = (float *) mycalloc(floattempstrsize, sizeof(float));
    twolongs[1]= 0;

    for (locus = 0; locus < minnodes; locus++)
      {
	twolongs[0] = locus;
	ll = (int) (locus + 1);
	MYMPIISEND(twolongs, TWO, MPI_LONG, ll, ll, comm_world, &irequests[numsent]);
	numsent++;
      }
    // waits until all minodes nodes received their message
    MYMPIWAITALL(minnodes, irequests, istatus);
    locus = 0;   
    world->mpistacknum = 0;
    while((alldone < world->loci) || (world->mpistacknum < numcpu-1))
      {
        done=FALSE;
        while(!done)
	  {
	    memset(savetempstr,0,sizeof(char)*(size_t)(tempstrsize+1));
	    MYMPIRECV (leadstr, SMALLBUFSIZE, MPI_CHAR, (MYINT) MPI_ANY_SOURCE, 
		       (MYINT) MPI_ANY_TAG, comm_world, &status);
	    sender = status.MPI_SOURCE;
	    tag = status.MPI_TAG;
	    switch(leadstr[0])
	      {
	      case 'D': //request data [handle data on demand]
		handle_dataondemand(sender, tag, leadstr, world, options, data);
		break;

	      case 'M': //get message to print or file
		newsize = atol(leadstr+1);
		if (newsize > tempstrsize)
		  {
		    tempstrsize = newsize;
		    savetempstr = (char*) myrealloc(savetempstr,sizeof(char)*(size_t)(tempstrsize+1));
		  }
		tempstr = savetempstr;
		memset(tempstr,0,sizeof(char)*(size_t) (tempstrsize+1));	       
		MYMPIRECV (tempstr, (size_t) newsize, MPI_CHAR, sender, tag,
			   comm_world, &status);
		handle_message(tempstr, sender, world);
		break;
		
	      case 'Z': // reading/writing of the raw bayes posterior data , also used to guide
		// multiple replicates on when to start sampling
		newsize = atol(leadstr+1);
		if (newsize > floattempstrsize)
		  {
		    floattempstrsize = newsize;		  
		    temp = (float*) myrealloc(temp,sizeof(float)* (size_t) (floattempstrsize+1));
		  }
		MYMPIRECV (temp, (size_t) newsize, MPI_FLOAT, sender, tag, comm_world, &status);
		handle_mdim(temp,newsize,sender, world);
		break;
		
	      case 'B': //burnin stopping rule
		newsize = atol(leadstr+1);
		if (newsize > tempstrsize)
		  {
		    tempstrsize = newsize;
		    savetempstr = (char*) myrealloc(savetempstr,sizeof(char)*(size_t) (tempstrsize+1));
		  }
		tempstr = savetempstr;
		memset(tempstr,0,sizeof(char)*(size_t)(tempstrsize+1));
		MYMPIRECV (tempstr, (size_t) tempstrsize, MPI_CHAR, sender, tag,
			   comm_world, &status);
		handle_burnin_message(tempstr,sender-BURNTAG, world);
		break;

	      case 'R': //ignore first character and translate into locusnumber
                locusdone = atol(leadstr+1);
		//if negative this means locus had no data at all -> ignored 
		if(locusdone<0)
		  {
		    locusdone = -locusdone;
		    world->data->skiploci[locusdone]=TRUE;
		  }
		done=TRUE;
		++alldone;
		//printff("%i> MASTER: locus %li, %li of %li finished \n",myID, locusdone,alldone,world->loci);
                break;

	      case 'N': // need a replicator, this message is sent from a
		// locus-node to the master for distribution among nodes that
		// are waiting for work, the master send will delegate a 
		// replicate and tell also the replicator 
		// where it needs to send the final result.
#ifdef DEBUG_MPI
	        printf("%i> MASTER: received replicator request from n%i using string %s\n",myID, sender, tempstr);
#endif
		// add to the mpistack_request list
		handle_replication(sender,tag,leadstr,world);
		break;

	      case 'G': // tell master that worker can do work
		// node has free time and could do replicator work
#ifdef DEBUG_MPI
	        printf("%i> MASTER: accepts replicator n%i at stack position %li\n",myID, sender,world->mpistacknum);
#endif
		world->mpistack[world->mpistacknum] = sender;
		world->mpistacknum += 1;
		if(alldone>=world->loci && (world->mpistacknum >= (numcpu-1)) && world->mpistack_requestnum==0)
		  {
#ifdef DEBUG_MPI
		    printf("%i> MASTER: reached mpistacknum=%li\n",myID,world->mpistacknum);
#endif
		    done=TRUE;
		    continue;
		  }
		while(world->mpistacknum > 0 && world->mpistack_requestnum > 0)
		  {
		    world->mpistack_requestnum -= 1;
		    sender = world->mpistack_request[world->mpistack_requestnum].sender;
		    tag = world->mpistack_request[world->mpistack_requestnum].tag;
		    handle_replication(sender,tag,
		      world->mpistack_request[world->mpistack_requestnum].tempstr,world);
		  }
		break;

	      default: /* never go here under normal run condition */
#ifdef DEBUG_MPI
                fprintf(stderr,"%i> message=@%s@\n@%i@> sender=@%i@ tag=@%i@\n",
			myID,leadstr, myID,status.MPI_SOURCE,status.MPI_TAG);
		done=TRUE;
#else
                MPI_Finalize();
                error("DIED because of wrong message from worker");
#endif
                break;
	      }
	  }
        who[locusdone] = sender;
	// if not done with loci send another locus-work to sender (a node 1..numcpu)
        if (numsent < loci)
	  {
            twolongs[0]=numsent;
            MYMPISEND (twolongs, TWO, MPI_LONG, (MYINT) sender, (MYINT) numsent + 1, comm_world);
            numsent++;
	  }
        else
	  { 
            twolongs[0] = 0;
	    //tell workers to stop waiting for new loci
	    if(!in_mpistack(sender, world))
	      {
		MYMPISEND (twolongs, TWO, MPI_LONG, (MYINT) sender, (MYINT) 0, comm_world); 
	      }
	  }
	locus++;
      }
    // stop loci and/or replicate worker that had never the chance to work on a locus 
    // or replicate, but are still listening
    for(locus=0;locus < world->mpistacknum;locus++)
      {
#ifdef DEBUG_MPI        
	fprintf(stdout,"%i> sent kill to node %i\n",myID, world->mpistack[locus]);
#endif
	mpi_send_stop_tag(world->mpistack[locus], world);
      }
    myfree(twolongs);
    myfree(istatus);
    myfree(irequests);
    myfree(savetempstr);
    myfree(leadstr);
    myfree(temp);
}


///
/// worker nodes execute this function and serve the master
void mpi_runloci_worker (world_fmt ** universe, int usize,
                    option_fmt * options, data_fmt * data,
                    tpool_t * heating_pool, long maxreplicate,
                    long *treefilepos, long *Gmax)
{
    boolean done = FALSE;

    char *rawmessage;
    //char tempstr[255];
    long locus;
    long *twolongs;
    long rawmsgsize = 0;
    long nbase      = data->loci+1;
    MPI_Status status;   
    long i;
    twolongs        = (long *) mycalloc(TWO, sizeof(long));
    rawmessage      = (char *) mycalloc(SMALLBUFSIZE,sizeof(char));
    if(myID < nbase)
      {
        while (!done)
	  {
            MYMPIRECV (twolongs, TWO, MPI_LONG, MASTER, MPI_ANY_TAG,
		       comm_world, &status);
            locus = twolongs[0];
            if (status.MPI_TAG != 0) //stop condition
	      {
                mpi_run_locus(universe, usize, options, data, 
                              heating_pool, maxreplicate, locus, treefilepos, Gmax);  
		if(universe[0]->data->skiploci[locus]==FALSE)
		  rawmsgsize = 1 + sprintf(rawmessage,"R%li",locus);
		else
		  rawmsgsize = 1 + sprintf(rawmessage,"R-%li",locus);
                MYMPISEND (rawmessage, SMALLBUFSIZE, MPI_CHAR, (MYINT) MASTER, 
			   (MYINT) (locus + ONE), comm_world);
                /* we want to know what locus we worked for
                   - to control the work sent by master*/
                universe[0]->who[locidone++] = (int) locus;
	      }
            else
	      {
                done = TRUE;
		mpi_runreplicates_worker (universe, usize, options,  data, heating_pool, treefilepos, Gmax);
	      }
	  }
      }
    else
      {
        mpi_runreplicates_worker (universe, usize, options,  data, heating_pool, treefilepos, Gmax);
      }        
    myfree(twolongs);
    myfree(rawmessage);
}

///
/// each locus is responsible for replication farming and and reporting back to master
/// master <- locus-master <- locus-replicate-worker
/// generate genealogies for a single locus
/// \callgraph
void
mpi_run_locus(world_fmt ** universe, int usize, option_fmt * options,
          data_fmt * data, tpool_t * heating_pool, long maxreplicate,
          long locus, long *treefilepos, long *Gmax)
{
    boolean done=FALSE;    
    char   *tempstr;
    int     tag;
    int    *who;
    int     minnodes;                         // number of replicate-worker nodes
    int     sender          = 0;
    int     nbase           = (int) data->loci;    //number of loci-worker nodes
    int     senderlocus     = -1;
    long    replicate;
    long    i;
    long   *temp;
    long    numsent         = 1;//
    long    senderreplicate = 0;

    MPI_Request *irequests = NULL;

    MPI_Status status;
    MPI_Status *istatus    = NULL;
    
    temp    = (long *) mycalloc(TWO, sizeof(long));
    tempstr = (char *) mycalloc(LONGLINESIZE,sizeof(char));
    who     = (int *) mycalloc(maxreplicate,sizeof(int));
    
    if(maxreplicate>1)
      {
	// number of nodes available for replicate workers,
	// numcpu - nbase - 1 [for masternode] and maxreplicate-1 
	// [locus-worker is doing one replicate itself]  are limiting
	//minnodes =  maxreplicate; //MAX(0,MIN(maxreplicate-1,numcpu-nbase-1))+1;
	minnodes =  MAX(0,MIN(numcpu-nbase-1,maxreplicate-1))+1;
	irequests = (MPI_Request *) mycalloc(minnodes+1,sizeof(MPI_Request));
	istatus   = (MPI_Status *) mycalloc(minnodes+1,sizeof(MPI_Status));
	temp[0] = locus;
	// replicate 1 to maxreplicate should be worked on by other nodes
	// minnodes is the number of nodes free for work on replicates alone
	// so with 3 replicates to work on and 3 total nodes and a single locus we need
	// 1 master, 1 locus-worker and have 1 replicate-worker, the locus-worker and 
	// the replicate worker will both work on replicates, the locus-worker sends 
	// off 1 request (for the replicate worker), because we start the locus worker 
	// with replicate 0, the loop for the replicate workers starts at one
	// and to address this shift the MIN(..,minnodes+1) addresses that
	for (replicate = 1; replicate < minnodes; replicate++)//Cesky Krumlov 2013 replicate=1 replaced
	  {
            temp[1] = replicate;
	    sprintf(tempstr,"N%i %li %li", myID, locus, replicate);
            MYMPIISEND (tempstr, SMALLBUFSIZE, MPI_CHAR, 
			(MYINT) MASTER, (MYINT) (locus + 1 + REPTAG), comm_world, &irequests[numsent]);
            numsent++;   // counter of how many replicates are sent off-node
	  }
	run_replicate(locus, 0, universe, options, data,
		      heating_pool, usize,
		      treefilepos, Gmax); 
	who[0] = myID;
        MYMPIWAITALL(numsent,irequests, istatus); // wait for all replicators to finish
        // set replicate counter to 1 because locus-worker itself finished first replicate
        replicate=1;        
        done = FALSE;
        while(!done)
	  {
            // done=TRUE means that
            // no replicator worker is available yet
            // the loci-worker has to do all the work if numsent==1
            if(numsent==1)
	      {
                run_replicate(locus, replicate, universe, options, data, 
                              heating_pool, usize, treefilepos, Gmax);
                who[replicate] = myID;
                replicate++;
                if(replicate >= maxreplicate)
		  done=TRUE;
	      }            
            else
	      {
                memset(irequests,0,sizeof(int)*(size_t) minnodes);
                memset(istatus,0,sizeof(int)*minnodes);
                MYMPIRECV (tempstr, SMALLBUFSIZE, MPI_CHAR, MPI_ANY_SOURCE, 
			   (MYINT)(locus+1+ REPTAG), comm_world, &status);
                sender = status.MPI_SOURCE;  // repID of the replicator that did the work
                tag = status.MPI_TAG;        // tag is working locus + replicator tag 
                senderlocus = tag-REPTAG;    // locus that was worked on by sender  
                // test so that we can be sure that we got things from a valid replicator worker
                if(sender == myID)
		  {
                    fprintf(stdout,"%i, %i> DIE\n###########################################\n",myID, myRepID);
                    error("tried to send a replicate to myself using MPI -- this is not allowed\n");
		  }
                // test whether the locus is the same between locus-worker and replicator
                if(senderlocus-1 != locus)
		  warning("%i> !!!!!!!!!!!!! got wrong locus from worker myRepID=%i (my locus %i != its locus %i )\n",
                          myID, sender, locus, senderlocus-1);
                // receive only messages that are prefixed with 'R' and exit on all others
                if(tempstr[0]=='R')
		  {
                    //ignore first character and translate into repnumber
                    senderreplicate = atol(tempstr+1);
		  }
                else
		  {
		    fprintf(stderr,"%i> message=%s\n%i> sender=%i tag=%i\n",myID,tempstr, myID,
			    status.MPI_SOURCE,status.MPI_TAG);
                    error("DIED because of wrong message from worker");
		  }
		// record sender , this record should be filled at the end of this function
                who[senderreplicate] = sender;
                mpi_receive_replicate(sender, tag, locus, senderreplicate, universe[0]); 
                replicate++;   
                if(replicate >= maxreplicate)
		  {
                    done=TRUE;
		  }
                temp[0] = locus;
                if (numsent < maxreplicate) //at least one set was worked by the locus-worker
		  {
                    temp[1] = numsent;
		    sprintf(tempstr,"N%i %li %li", myID, locus, numsent);
		    MYMPIISEND (tempstr, SMALLBUFSIZE, MPI_CHAR, 
				(MYINT) MASTER, (MYINT) (locus + 1 + REPTAG), comm_world, &irequests);
                    numsent++;
                    MYMPIWAITALL(ONE, &irequests, &status);
		  }
	      }
	  }
	myfree(irequests);
	myfree(istatus);
      }
    else
      { /* no replicates */
        run_replicate(locus, 0, universe, options, data, 
                      heating_pool, usize,
                      treefilepos, Gmax);
      }    
    myfree(temp);
    myfree(tempstr);
    myfree(who);
    if(maxreplicate>1 && (numcpu-1) > universe[0]->loci)
      {
	const long hc = universe[0]->options->heated_chains;
	long t;
	if(universe[0]->options->heating)
	  {
	    for(t=0; t < hc; t++)
	      {
		long z = locus * hc + t;
		universe[0]->bf[z] /= maxreplicate;
		universe[0]->steppingstones[z] /= maxreplicate;
	      }
	  }
      }
#ifdef UEP
    if (options->uep)
      show_uep_store(universe[0]);
#endif
    if (options->replicate && options->replicatenum > 0)
      {
	(universe[0])->repkind = MULTIPLERUN;
      }
    if(!universe[0]->options->has_bayesmdimfile)
      calculate_credibility_interval(universe[0], locus);
    
    // cleanup
    if (options->heating)
    {
        for (i = 0; i < options->heated_chains; i++)
        {
	  free_tree(universe[i]->root, universe[i]);
        }
    }
    else
    {
      free_tree(universe[0]->root, universe[0]);        
    }
    bayes_reset(universe[0]);
}


///
/// run replicates on replicate-worker nodes
/// this function is called in main() and will work on any preset locus and any preset replicate
/// the calling function is responsible for the correct assignment to locus and replicate.
void
mpi_runreplicates_worker (world_fmt ** universe, int usize,
                    option_fmt * options, data_fmt * data,
                    tpool_t * heating_pool, 
                    long *treefilepos, long *Gmax)
{
    boolean done = FALSE;
    char *rawmessage;
    int sender;
    long nng;
    long *temp;
    long replicate;
    long locus;
    long rawmsgsize = 0;
    char *ready;
    MPI_Status status;
#ifdef IPROBE
    int notwaiting=0;
    char tempstr[255];
#endif
    ready = mycalloc(SMALLBUFSIZE,sizeof(char));
    temp = (long *) mycalloc(3, sizeof(long));
    rawmessage = (char *) mycalloc(STRSIZE,sizeof(char));
    sprintf(ready,"G%i",myID);//redy message to be sent to the master
    while (!done)
      {
	MYMPISEND (ready, SMALLBUFSIZE, MPI_CHAR, (MYINT) MASTER, myID, comm_world);
#ifdef IPROBE
	// this forces nodes that do not have work to sleep for 10 seconds and then check again
	// this reduces the load on the machine when the run is near finishing and only few nodes
	// do work, on multicore machine this will give more cycles to the working nodes.
	while(1)
	  {
	    MPI_Iprobe(MASTER, MPI_ANY_TAG, comm_world, &notwaiting, &status);
	    if(!notwaiting)
	      {
		//get_time (tempstr, "%H:%M:%S");
		//fprintf(stdout,"%i> replicate worker waits for master -- %s\n",myID, tempstr);
		sleep(10);
		continue;
	      }
	    else
	      break;
	  }
#endif
        MYMPIRECV (temp, 3, MPI_LONG, (MYINT) MASTER, MPI_ANY_TAG, comm_world, &status);
        sender = (int) temp[0];
        locus = temp[1];
        replicate = temp[2];
        if (status.MPI_TAG != 0) //stop condition
          {
	    nng=universe[0]->numpop2 + universe[0]->bayes->mu + 1 + universe[0]->species_model_size * 2 + universe[0]->grownum;	    
	    memset(universe[0]->accept_archive,0,sizeof(long)*2*nng);//resets also trials_archive
            run_replicate(locus, replicate, universe, options, data, heating_pool, usize,treefilepos, Gmax);
            rawmsgsize = 1 + sprintf(rawmessage,"R%li ",replicate);
            MYMPISEND (rawmessage, rawmsgsize, MPI_CHAR, (MYINT) sender, 
		       (MYINT) (locus+1+ REPTAG), comm_world);
            mpi_send_replicate(sender, locus, replicate, universe[0]);
	    universe[0]->bayes->numparams = 0;
          }
        else
          {
            done = TRUE;
          }
      }
    myfree(ready);
    myfree(temp);
    myfree(rawmessage);
}


void
assign_worker_cleanup (void)
{
    boolean done = FALSE;
    char *rawmessage;
    int sender;
    long *temp;
    char *ready;
    MPI_Status status;
    ready      = mycalloc(SMALLBUFSIZE,sizeof(char));
    temp       = (long *) mycalloc(3, sizeof(long));
    rawmessage = (char *) mycalloc(STRSIZE,sizeof(char));
    sprintf(ready,"G%i",myID);
    while (!done)
      {
	MYMPISEND (ready, SMALLBUFSIZE, MPI_CHAR, (MYINT) MASTER, myID, comm_world);
        MYMPIRECV (temp, 3, MPI_LONG, (MYINT) MASTER, MPI_ANY_TAG, comm_world, &status);
        sender = (int) temp[0];
        if (status.MPI_TAG != 0) //stop condition
          {
	    printf("%i> received real message but do not know what to do with this\n",myID);
          }
        else
          {
            done = TRUE;
          }
      }
    myfree(ready);
    myfree(temp);
    myfree(rawmessage);
}

//---------end replication in MPI


///
/// setup start parameters 
void
mpi_startparam_master(world_fmt * world)
{
    int tag;
    int numreceived = 0;

    long i;
    long sender;
    long workerloci = 0;

    MYREAL  *tmp;

    MPI_Status status;

    long nn  = world->numpop2+ world->bayes->mu + 1 + world->species_model_size * 2 + world->grownum;
    long nng = nn + 1;

    tmp = (MYREAL*) mycalloc(nng,sizeof(MYREAL));

    while (numreceived < world->loci)
      {
        MYMPIRECV (tmp, nng, mpisizeof, MPI_ANY_SOURCE,
                   MPI_ANY_TAG, comm_world, &status);
        sender = status.MPI_SOURCE;
        tag = status.MPI_TAG;
        workerloci=tmp[0];
        for(i=0; i< nn; i++)
	  {
	    world->param0[i] += tmp[i+1];
	  }
        numreceived+=workerloci;
      }
    for(i=0; i<nn; i++)
      world->param0[i] /= world->loci;
    
    myfree(tmp);
}

///
/// set start parameters for workers
void
mpi_startparam_worker (world_fmt * world)
{
    long ww;
    long repstart;
    long repstop;
    long r;
    long i;
    long locus;

    MYREAL *tmp;
    long nn  = world->numpop2+ world->bayes->mu + 1 + world->species_model_size * 2 + world->grownum;
    long nng = nn + 1;

    if(locidone>0)
      {
	tmp = (MYREAL*) mycalloc(nng,sizeof(MYREAL));
        set_replicates (world, world->repkind, world->options->replicatenum,
                        &repstart, &repstop);
        tmp[0]=(MYREAL)locidone;
        for (ww = 0; ww < locidone; ww++)
	  {
            locus = world->who[ww];
            for (r = repstart; r < repstop; r++)
	      {
                for(i=0; i < nn; i++)
		  tmp[i+1] += world->atl[r][locus].param[i];
	      }
	  }
        for(i=1; i < nng; i++)
	  {
	    tmp[i] /= locidone * (repstop-repstart);
	  }
	MYMPISEND (tmp, nng, mpisizeof, MASTER, myID, comm_world);
	myfree(tmp);
      }
}


///
/// orchestrates max(gmax) over all nodes
void
mpi_gmax_master (world_fmt * world, long *Gmax)
{
    int tag;
    int sender;
    int numreceived = 0;

    long tmp;

    MPI_Status status;

    *Gmax = 0.;

    MYMPIBCAST (Gmax, ONE, MPI_LONG, MASTER, comm_world);

    while (numreceived < MIN(world->loci, numcpu - 1))
      {
	MYMPIRECV (&tmp, ONE, MPI_LONG, MPI_ANY_SOURCE,
		   MPI_ANY_TAG, comm_world, &status);
	sender = status.MPI_SOURCE;
	tag = status.MPI_TAG;
        if (*Gmax < tmp)
	  {
            *Gmax = tmp;
	  }
        numreceived++;
      }
    //  do we need this barrier really?
    MYMPIBARRIER(comm_world);
}

///
/// returns the gmax values to the master to evaluate max(gmax)
void
mpi_gmax_worker (world_fmt * world)
{
    long ww;
    long repstart;
    long repstop;
    long r;
    long locus;
    long Gmax = 1;

    MYMPIBCAST (&Gmax, ONE, MPI_LONG, MASTER, comm_world);
#ifdef DEBUG_MPI
    printf("%i> locidone=%i\n",myID,locidone);
#endif
    if(locidone>0)
      {
	set_replicates (world, world->repkind, world->options->replicatenum,
			&repstart, &repstop);
	
	for (ww = 0; ww < locidone; ww++)
	  {
	    locus = world->who[ww];
	    for (r = repstart; r < repstop; r++)
	      {
		if (Gmax < world->atl[r][locus].T)
		  Gmax = world->atl[r][locus].T;
	      }
	  }
	MYMPISEND (&Gmax, ONE, MPI_LONG, MASTER, myID, comm_world);
      }
    //  do we need this barrier really?
    MYMPIBARRIER(comm_world);
}

///
/// first worker (myID=1, myRepId=0) will send stop-message to replication nodes
/// the comm_worker group's master is the first worker who has ID=0 in this group
/// as results we send messages to id=1..x in the comm_worker group, do not mix this 
/// with the id in comm_world that include the master (id=0 there).
long  mpi_send_stop_mcmc_worker(long numcpu, long loci, MPI_Comm *comm, 
				     MPI_Request *irequests, MPI_Status *istatus, long id)
{
    long twolongs[2];
    long *temp;
    long receiver;
    long sent      = 0;
    long xx        = (id==0) ? 0 : 1;
 
    temp = (long *) mycalloc(TWO, sizeof(long));
    twolongs[0]=0;
    twolongs[1]=0;

    for(receiver=loci+1-xx; receiver< numcpu-xx; receiver++)
      {
	MYMPIISEND (temp, TWO, MPI_LONG, (MYINT) receiver, 0, *comm, &irequests[sent]);
        sent++;
      }
    if(sent>0)
      {
	MYMPIWAITALL(sent,irequests, istatus); // wait for all replicators to finish
      }
    myfree(temp);
    return sent;
}

///
/// first worker (myID=1, myRepId=0) will send stop-message to replication nodes
/// the comm_worker group's master is the first worker who has ID=0 in this group
/// as results we send messages to id=1..x in the comm_worker group, do not mix this 
/// with the id in comm_world that include the master (id=0 there).
long  mpi_send_stop_mcmc_lociworker(long numcpu, long loci)
{

  long receiver;
  long *temp;
  long xx       = (myID==0) ? 0 : 1;
  long sent     = 0;
  long minnodes = MIN(numcpu,loci -1);
  
  
  MPI_Request *irequests;
  MPI_Status *istatus;
  irequests = (MPI_Request *) mycalloc(minnodes+1,sizeof(MPI_Request));
  istatus = (MPI_Status *) mycalloc(minnodes+1,sizeof(MPI_Status));
  
  temp = (long *) mycalloc(TWO, sizeof(long));
  
  for(receiver=loci+1-xx; receiver< numcpu-xx; receiver++)
    {
           error("disable because testing why openmpi breaks");
      //      MYMPIISEND (temp, TWO, MPI_LONG, (MYINT) receiver, 0, comm_workers, &irequests[sent]);
      sent++;
    }
  if(sent>0)
    MYMPIWAITALL(sent,irequests, istatus); // wait for all replicators to finish
  myfree(temp);
  myfree(irequests);
  myfree(istatus);
  return sent;
}

///
/// stops all replicate workers
long  mpi_send_stop_mcmc_replicateworker(long numcpu, long loci)
{

  long *temp;
  long receiver;
  long sent      = 0;
  long xx        = (myID==0) ? 0 : 1;
  long minnodes  = labs(numcpu - loci -1);

  MPI_Request *irequests;
  MPI_Status *istatus;
  
  irequests = (MPI_Request *) mycalloc(minnodes+1,sizeof(MPI_Request));
  istatus   = (MPI_Status *) mycalloc(minnodes+1,sizeof(MPI_Status));
  temp      = (long *) mycalloc(TWO, sizeof(long));

  for(receiver=loci+1-xx; receiver< numcpu-xx; receiver++)
    {
      error("disable because testing why openmpi breaks");
      //     MYMPIISEND (temp, 2, MPI_LONG, (MYINT) receiver, 0, comm_workers, &irequests[sent]);
      sent++;
    }
  if(sent>0)
    {
      MYMPIWAITALL(sent,irequests, istatus); // wait for all replicators to finish
    }
  myfree(temp);
  myfree(irequests);
  myfree(istatus);
  return sent;
}

///
/// sends a stop signal to all loci-worker
void
mpi_send_stop (world_fmt * world)
{
  long worker;
  long numelem  = world->numpop2 + (world->options->gamma ? 1 : 0);
  long numelem2 = 2 * numelem;
  
  MYREAL *temp;

  temp = (MYREAL *) mycalloc (numelem2+2, sizeof (MYREAL));
  temp[0] = MIGMPI_END;

  for (worker = 1; worker < numcpu; worker++)
    {
      MYMPISEND (temp, (MYINT) numelem2+2, mpisizeof, (MYINT) worker, (MYINT) 0, comm_world); //end of loci
    }
  myfree(temp);
}

///
/// sends a stop signal to a specific worker used in for assignloci_worker
void mpi_send_stop_tag (int worker, world_fmt * world)
{
  long *temp;

  temp = (long *) mycalloc (TWO, sizeof (MYREAL));

  temp[0] = 0;
  temp[1] = 0;
  MYMPISEND (temp, (MYINT) TWO, MPI_LONG, (MYINT) worker, (MYINT) 0, comm_world);

  myfree(temp);

}

///
/// sends a continue with further analysis to the workers
void
mpi_results_stop (void)
{
    long worker;
    long dummy = 0;

    for (worker = 1; worker < numcpu; worker++)
      {
        MYMPISEND (&dummy, ONE, MPI_LONG, (MYINT) worker, (MYINT) 0, comm_world);
      }
}

///
/// provides calculations for the master. calculates likelihoods, gradients, but also
/// delivers results, migration histories, bayesian results and more.
/// workers stay a long time in this function and answer requests from the master
void
mpi_maximize_worker (world_fmt * world, option_fmt *options, long kind, long rep)
{
  boolean done = FALSE;
  long locus;
  long repstart;
  long repstop;
  long Gmax;
  long numelem  =  world->numpop2 + (world->options->gamma ?  1 : 0);
  long numelem2 = numelem * 2;
  MYREAL *temp;
  MPI_Status status;
  nr_fmt *nr;
  helper_fmt helper;
#ifdef PARALIO
  char tempfilename[STRSIZE];
#endif
  temp = (MYREAL *) mycalloc (numelem2 + 2, sizeof (MYREAL));
  helper.xv = (MYREAL *) mycalloc (numelem2, sizeof (MYREAL));
  helper.expxv = (MYREAL *) mycalloc (numelem2, sizeof (MYREAL));
  helper.analystype = SINGLELOCUS; //this may be problematic byt seems only
  // to be used here for decision on gamma or not 
  nr = (nr_fmt *) mycalloc (1, sizeof (nr_fmt));
  set_replicates (world, world->repkind, rep, &repstart, &repstop);
  //which_calc_like (world->repkind);
  Gmax=1;
  //create_nr (nr, world, Gmax, 0, world->loci, world->repkind, repstart);

  while (!done)
    {
      MYMPIRECV (temp, numelem2+2, mpisizeof, MASTER, MPI_ANY_TAG,
		 comm_world, &status);
      locus = world->locus = status.MPI_TAG - 1;
      switch ((long) temp[0])
        {
        case MIGMPI_RESULT: // returns the results (stats for the loci worked on
	  mpi_results_worker ((long) temp[0], world, repstop, pack_result_buffer);
	  break;
	  //        case MIGMPI_PLOTPLANE: // returns the results (stats for the loci worked on
	  //mpi_results_worker ((long) temp[0], world, repstop, pack_plotplane_buffer);
	  //break;
        case MIGMPI_TREESPACE: // returns the treefile to the master
	  //	  printf("%i> print treespace to master\n",myID); 
	  mpi_results_worker ((long) temp[0], world, repstop, pack_treespace_buffer);
	  break;
        case MIGMPI_MIGHIST: // returns the histogram of events (top parts + obsolete parts)
	  mpi_results_worker ((long) temp[0], world, repstop, pack_mighist_buffer);
	  break;
        case MIGMPI_SKYLINE: // returns the histogram of events
	  //fprintf(stdout,"%i-------------------\n",myID);
	  //	  debug_skyline(world,"mpiresultsworker-skyline");
	  mpi_results_worker ((long) temp[0], world, repstop, pack_skyline_buffer);
	  //fprintf(stdout,"%i-------------------\n",myID);
	  break;
        case MIGMPI_BAYESHIST: // returns the histogram of the parameters
	  //printf("%i> send bayeshist\n",myID);
	  mpi_results_worker ((long) temp[0], world, repstop, pack_bayes_buffer);
	  break;
	case MIGMPI_SEQERROR:
	  mpi_results_worker((long) temp[0], world, repstop, pack_seqerror_buffer);
	  break;
	case MIGMPI_ASSIGN:
	  mpi_results_worker((long) temp[0], world, repstop, pack_assign_buffer);
	  break;
        case MIGMPI_HAPLOTYPING: // returns the histogram of the parameters
	  //	  printf("%i> send haplotypes\n",myID);
	  mpi_results_worker ((long) temp[0], world, repstop, pack_haplotypes_buffer);
	  break;
#ifdef PARALIO
	case MIGMPI_PARALIO:
	  fprintf(stdout,"%i> bayesmdimfile %s closing\n", myID, options->bayesmdimfilename);
	  MPI_File_close(&world->mpi_bayesmdimfile);
	  sprintf(tempfilename,"%s",options->bayesmdimfilename);
	  MYMPISEND (tempfilename, STRSIZE, MPI_CHAR, (MYINT) MASTER, (MYINT) myID, comm_world);
	  break;
#endif
        case MIGMPI_END:
	  done = TRUE;
	  break;
        default:
	  fprintf (stdout, "%i> Do not understand task --> exit \n", myID);
	  MPI_Finalize();
	  exit (0);
        }
    }
  myfree(temp);
  myfree(helper.xv);
  myfree(helper.expxv);
  //destroy_nr (nr, world);
}

///
/// sends out the options to all workers
void
broadcast_options_master (option_fmt * options, data_fmt *data)
{
  long allocbufsize = LONGLINESIZE;
  long bufsize;
  char *buffer;
#ifdef DEBUG
  char nowstr[LINESIZE];
  get_time (nowstr, "%c");
  printf("Master started option broadcast at %s\n", nowstr);
#endif
  //  char *buffermu;
  buffer = (char *) mycalloc (allocbufsize, sizeof (char));
  bufsize = save_options_buffer (&buffer, &allocbufsize, options, data);
  MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
  MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);
  myfree(buffer);
#ifdef DEBUG
  get_time (nowstr, "%c");
  printf("master finished broadcast of options %s\n",nowstr);
#endif
}

///
/// receives the options 
void
broadcast_options_worker (option_fmt * options)
{
  //  long i;
  //  char *temp;
  long bufsize;
  char *buffer, *sbuffer;
  MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
  buffer = (char *) mycalloc (bufsize + 1, sizeof (char));
  sbuffer = buffer;
  MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);
  if (strncmp(buffer,"MENUSTOP",8)==0)
    {
      MPI_Finalize();
      exit(0);
    }
  read_options_worker (&buffer, options);
  myfree(sbuffer); 
#ifdef DEBUG
  char nowstr[LINESIZE];
  get_time (nowstr, "%c");
  printf("%i > worker finished option broadcast receive at %s\n",myID, nowstr);
#endif
 
}


///
/// broadcasts the data to all nodes
void
broadcast_data_master (data_fmt * data, option_fmt * options)
{
  long bufsize;
  long allocbuffermusize = 1;
  char *buffermu;
#ifdef DEBUG
  char nowstr[LINESIZE];
  get_time (nowstr, "%c");
  printf("Master broadcasts data at %s\n",nowstr);
#endif
  bufsize = pack_databuffer (data, options);
  if(options->murates)
    {
      buffermu = (char *) mycalloc (1, sizeof (char));
      bufsize = save_mu_rates_buffer(&buffermu,&allocbuffermusize, options);
      MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
      MYMPIBCAST (buffermu, bufsize, MPI_CHAR, MASTER, comm_world);
      myfree(buffermu);
    }
#ifdef DEBUG
  get_time (nowstr, "%c");
  printf("Master finished broadcasts of data at %s\n",nowstr);
#endif
}

///
/// receives the data from the master
void
broadcast_data_worker (data_fmt * data, option_fmt * options, world_fmt *world)
{
  long bufsize=1;
  char *buffermu, *sbuffermu;
  char *temp;
  long i;
  unpack_databuffer (data, options, world);
  if(options->murates)
    {
      MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
      buffermu = (char *) mycalloc (bufsize + 1, sizeof (char));
      sbuffermu = buffermu;
      MYMPIBCAST (buffermu, bufsize, MPI_CHAR, MASTER, comm_world);
      //fprintf(stdout,"%i> received muratebuffer:>%s<\n",myID, buffermu);
      temp=strsep(&buffermu,"\n");
      options->muloci = atoi(temp);
      if(options->mu_rates==NULL)
	options->mu_rates = (MYREAL *) mycalloc(1+options->muloci, sizeof(MYREAL));
      for(i=0; i< options->muloci; i++)
	{
	  temp = strsep(&buffermu,"\n");
	  options->mu_rates[i] = atof(temp);
	}
      myfree(sbuffermu);
    }
#ifdef DEBUG
  char nowstr[LINESIZE];
  get_time (nowstr, "%c");
  printf("%i > worker finished data broadcast receive at %s\n",myID, nowstr);
#endif

}


boolean store_to_buffer_inloop(char **buffer, long *bufsize, char *smallbuffer, long *smallbufsize, long smallbufferallocsize, long threshold)
{
  if(*smallbufsize >= smallbufferallocsize - threshold)
  {
    *buffer = (char *) myrealloc (*buffer, sizeof (char) * (*bufsize + *smallbufsize));
    memcpy(*buffer + *bufsize, smallbuffer, *smallbufsize * sizeof(char));
    *bufsize += *smallbufsize;
    *smallbufsize = 0;
    return TRUE;
  }
  return FALSE;
}

boolean store_to_buffer_afterloop(char **buffer, long *bufsize, char *smallbuffer, long *smallbufsize, boolean done)
{
  if(!done)
  {
    *buffer = (char *) myrealloc (*buffer, sizeof (char) * (*bufsize + *smallbufsize+1));
    memset(*buffer+ *bufsize, '\0', (*smallbufsize+1) * sizeof(char));
    memcpy(*buffer + *bufsize, smallbuffer, *smallbufsize * sizeof(char));
    *bufsize += *smallbufsize;
    //printf("|%s|%li|%li\n",*buffer,*bufsize,*smallbufsize);
    //*buffer[*bufsize-1]='\0';
    *smallbufsize = 0;
    return TRUE;
  }
  return FALSE;
}


///
/// pack the data for shipment to the workers
long
pack_databuffer (data_fmt * data, option_fmt * options)
{
  //boolean done=FALSE;
  long locus, pop, ind;
  //long lo;
  long bufsize = 0;
  //long fpsize = 0;
  long biggest;
  long i;
#ifdef UEP
  long sumtips;
#endif
  char *buffer;
  long allocfpsize=LONGLINESIZE;
  
  buffer = (char *) mycalloc (allocfpsize,sizeof (char));
  bufsize = sprintf (buffer, "%c %li %li %li %li\n", options->datatype, (long) data->hasghost,
			  data->numpop, data->loci, data->allsubloci);
  // send header info
  bufsize += 1;
  MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
  MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);
  // send sites and loci information
  bufsize=0;
  for (locus = 0; locus < data->loci; locus++)
    {
      bufsize += sprintf (buffer+bufsize, "%li %li\n", data->subloci[locus], data->sublocistarts[locus]);
      if (bufsize >= allocfpsize-100)
	{
	  allocfpsize += LONGLINESIZE;
	  buffer = (char *) myrealloc (buffer, sizeof (char) * allocfpsize);
	}
    }
  bufsize += sprintf(buffer+bufsize,"%li\n", data->sublocistarts[locus]);
  for (locus = 0; locus < data->allsubloci; locus++)
    {
      int dataclass = 0;
      int scaling = 0;
      mutationmodel_fmt *s = &data->mutationmodels[locus];
      if(s->dataclass==SITEWORD)
	dataclass = 1;
      if(s->scaling==TRUE)
	scaling = 1;
      if (bufsize >= allocfpsize-500)
	{	  
      	  allocfpsize += LONGLINESIZE;
      	  buffer = (char *) myrealloc (buffer, sizeof (char) * allocfpsize);
	}
      if (s->datatype == '\0')
	s->datatype = options->datatype;
      bufsize += sprintf (buffer + bufsize, "%c %i %i %li %li\n",
			  s->datatype, dataclass, s->model, s->numpatterns, s->numsites);
      bufsize += sprintf (buffer + bufsize, "%li %li %li %f\n",
			  s->startsite, s->numstates, s->numsiterates,s->lambda); 
      bufsize += sprintf (buffer + bufsize, "%f %f %f %i %i\n",
			  s->parameters[0],s->parameters[1],s->parameters[2],(int) s->from_infile, (int) s->finished);
      bufsize += sprintf (buffer + bufsize, "%li %f %f %f %f\n",
			  s->weightsum, s->xi, s->xv, s->ttratio, s->fracchange);
      bufsize += sprintf (buffer + bufsize, "%f %f %li %f\n",
			  s->freq, s->freqlast,s->maxalleles,s->browniandefault);
      bufsize += sprintf (buffer + bufsize, "%li %li %li\n",
			  s->micro_threshold, s->microrange,s->microstart);
      bufsize += sprintf (buffer + bufsize, "%li %li %i %i %i\n",
			  s->addon, s->numcategs,(int) s->estimateseqerror, (int) s->seqerrorcombined, (int) s->scaling);
      bufsize += sprintf (buffer + bufsize, "%f %f %f %f\n", 
			  s->scoring_error[0], s->scoring_error[1], s->scoring_error[2], 
			  s->scoring_error[3]);
      for(i=0;i<s->numsites;i++)
	{
	  if (bufsize >= allocfpsize-100)
	    {
	      allocfpsize += LONGLINESIZE;
	      buffer = (char *) myrealloc (buffer, sizeof (char) * allocfpsize);
	    }
	  bufsize += sprintf (buffer + bufsize,"%f %li\n",s->aliasweight[i],s->alias[i]);
	}
    }
  // population data
  bufsize += 1;
  MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
  MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);


  // send population data, send locus data, then one individual at a time
  for (pop = 0; pop < data->numpop; pop++)
    {
      //send each population name 
      if (allocfpsize < LINESIZE)
	{
	  allocfpsize += LINESIZE;
	  buffer = (char *) myrealloc (buffer, sizeof (char) * allocfpsize);
	}
      bufsize = sprintf (buffer, "%s\n", data->popnames[pop]);//set to start
      bufsize += 1;
      MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
      MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);
      biggest = 0;
      bufsize = 0;
      // send number of ind for each locus per population
      for(locus=0; locus<data->loci; locus++)
        {
	  if (bufsize >= allocfpsize-100)
	    {
	      allocfpsize += LONGLINESIZE;
	      buffer = (char *) myrealloc (buffer, sizeof (char) * allocfpsize);
	    }
#ifdef DEBUG
	  printf("%i>numind send: %li %li %li %li\n", myID, pop, locus, data->numind[pop][locus], data->numalleles[pop][locus]);
#endif
	  bufsize += sprintf (buffer + bufsize, "%li %li\n", data->numind[pop][locus],data->numalleles[pop][locus]);
	  if(biggest < data->numind[pop][locus])
	    biggest = data->numind[pop][locus];
	}
      bufsize += 1;
      MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
      MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);
      bufsize = 0;
      // send individuals
#ifndef MPIDATAONDEMAND    
      if (!strchr (SEQUENCETYPES, options->datatype) && options->datatype != '@')
        {
	  for (ind = 0; ind < biggest; ind++)
            {
	      // assume that buffer is always larger than nmlength
	      // resetting bufsize to start!
	      bufsize = sprintf (buffer, "%*s\n",  (int) options->nmlength, data->indnames[pop][ind][0]);
	      pack_allele_data (&buffer, &bufsize, &allocfpsize, data, pop, ind);
	      bufsize += 1;
	      //printf("%i> pack_allele_data: bufsize=%li\n",myID,bufsize);
	      MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
	      MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);
	    }
        }
      else
        {
	  bufsize = 0;
	  for(locus=0;locus<data->loci; ++locus)
	    {
	      for (ind = 0; ind < data->numind[pop][locus]; ind++)
                {
		  // assume that buffer is always larger than nmlength
		  // is resetting bufsize to start!
		  bufsize = sprintf (buffer, "%*s\n", (int) options->nmlength,
				     data->indnames[pop][ind][locus]);
		  pack_sequence_data (&buffer, &bufsize, &allocfpsize, data, pop, ind, locus);
		  bufsize += 1;
		  MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
		  MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);
                }
            }
        }
#endif /*MPIDATAONDEMAND*/
    }
  // send geofile
  if (options->geo)
    {
      bufsize = 0;
      for (pop = 0; pop < data->numpop * data->numpop; pop++)
        {
	  if (bufsize >= allocfpsize-100)
	    {
	      allocfpsize += LONGLINESIZE;
	      buffer = (char *) myrealloc (buffer, sizeof (char) * allocfpsize);
	    }
	  bufsize += sprintf (buffer+ bufsize, "%f %f\n", data->geo[pop], data->lgeo[pop]);
        }
      bufsize += 1;
      MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
      MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);
    }
  // send randomsubset data
  bufsize = 0;
  for (pop = 0; pop < data->numpop; pop++)
    {
      for(locus=0;locus<data->loci;locus++)
	{
	  for(ind=0; ind < data->numind[pop][locus]; ind++)
	    {
	      if (bufsize >= allocfpsize-100)
		{
		  allocfpsize += LONGLINESIZE;
		  buffer = (char *) myrealloc (buffer, sizeof (char) * allocfpsize);
		}
	      bufsize += sprintf (buffer+ bufsize, "%li\n", data->shuffled[pop][locus][ind]);
	    }
        }
    }
  bufsize += 1;
  MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
  MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);
  bufsize = 0;
  // tipdate file
  if (options->has_datefile)
    {
      locus = 0;
      if(data->oneliner) // assumes all data for an individual is on one line
      {
	for (pop = 0; pop < data->numpop; pop++)
	  {
	    for(ind=0; ind < data->numind[pop][0]; ind++)
	      {
		  if (bufsize >= allocfpsize-200)
		    {
		      allocfpsize += LONGLINESIZE;
		      buffer = (char *) myrealloc (buffer, sizeof (char) * allocfpsize);
		    }
		  bufsize += sprintf (buffer+bufsize, "%s %f\n", data->sampledates[pop][locus][ind].name, data->sampledates[pop][locus][ind].date);
	      }
	  }
      }
      else
	{
	  // assumes locus and sublocus are the same!!!
	  for(locus=0;locus<data->loci;locus++)
	    {
	      for (pop = 0; pop < data->numpop; pop++)
		{
		  for(ind=0; ind < data->numind[pop][locus]; ind++)
		    {
		      if (bufsize >= allocfpsize-200)
			{
			  allocfpsize += LONGLINESIZE;
			  buffer = (char *) myrealloc (buffer, sizeof (char) * allocfpsize);
			}
		      bufsize += sprintf (buffer+bufsize, "%s %f\n", data->sampledates[pop][locus][ind].name, data->sampledates[pop][locus][ind].date);
		    }
		}
	    }
	}
      // assumes buffer is big enough
      bufsize += sprintf (buffer + bufsize, "%20.20f\n", data->maxsampledate);
      //      printf("%i> pack_data maxsampledate=%f %s\n",myID, data->maxsampledate, fp);
      //printf("%i> pack_data:%s\n""""""""""""""""""""""""\n\n\n",myID, buffer);
      bufsize += 1;
      MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
      MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);
    }
  // uepfile
#ifdef UEP
  if (options->uep)
    {
      sumtips = 0;
      bufsize = 0;
      for (pop = 0; pop < data->numpop; ++pop)
	sumtips += data->numind[pop][0];//Assumes UEP is matched by locus 1
      bufsize += sprintf (buffer+ bufsize, "%li %li\n", sumtips, data->uepsites);
      if (strchr (SEQUENCETYPES, options->datatype))
        {
	  for (pop = 0; sumtips; pop++)
            {
	      for (i = 0; i < data->uepsites; i++)
                {
		  if (bufsize >= allocfpsize-100)
		    {
		      allocfpsize += LONGLINESIZE;
		      buffer = (char *) myrealloc (buffer, sizeof (char) * allocfpsize);
		    }
		  bufsize += sprintf (buffer + bufsize, "%i\n", data->uep[pop][i]);
                }
	      done = store_to_buffer_afterloop(buffer, &bufsize, fp, &fpsize, done);
            }
        }
      else
        {
	  for (pop = 0; sumtips; pop++)
            {
	      for (i = 0; i < data->uepsites; i++)
                {
		  if (bufsize >= allocfpsize-100)
		    {
		      allocfpsize += LONGLINESIZE;
		      buffer = (char *) myrealloc (buffer, sizeof (char) * allocfpsize);
		    }
		  bufsize += sprintf (buffer + bufsize, "%i %i\n", data->uep[pop][i],
					  data->uep[pop + sumtips][i]);
                }
            }
        }
      bufsize += 1;
      MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
      MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);
    }
  
#endif
  // haplotyping distributor
  // data->individuals
  // data->numindividuals
  //
  long id;
  bufsize = 0;
  if(options->haplotyping)
    {
      for(locus=0;locus<data->loci;locus++)
	{
	  if (bufsize >= allocfpsize-100)
	    {
	      allocfpsize += LONGLINESIZE;
	      buffer = (char *) myrealloc (buffer, sizeof (char) * allocfpsize);
	    }
	  bufsize += sprintf (buffer + bufsize, "%li\n", data->numindividuals[locus]);
	}
      for(locus=0;locus<data->loci;locus++)
	{
	  for (id=0; id < data->numindividuals[locus]; id++)
	    {
	      if (bufsize >= allocfpsize-100)
		{
		  allocfpsize += LONGLINESIZE;
		  buffer = (char *) myrealloc (buffer, sizeof (char) * allocfpsize);
		}
	      bufsize += sprintf (buffer + bufsize, "%li\n%s\n%li %li %li\n",
				  data->individuals[locus][id].id,
				  data->individuals[locus][id].name,
				  data->individuals[locus][id].ind[0],
				  data->individuals[locus][id].ind[1],
				  data->individuals[locus][id].numstates);
	    }
	}
      bufsize += 1;
      printf("%i> %s\n",myID,buffer);
      MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
      MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);
    }
  myfree(buffer);
  return 0;
}


void
pack_allele_data (char **buffer, long *bufsize, long *allocbufsize, data_fmt * data, long pop,
                  long ind)
{
  //boolean done=FALSE;
  long locus;
  for (locus = 0; locus < data->allsubloci; locus++)
    {
      if (*bufsize > *allocbufsize-100)
	{
	  *allocbufsize += *bufsize;
	  *buffer = (char*) myrealloc(*buffer, *allocbufsize * sizeof(char));
	}
        *bufsize += sprintf (*buffer + *bufsize, "%s %s\n", data->yy[pop][ind][locus][0][0],
                                 data->yy[pop][ind][locus][1][0]);
    }
}

void
pack_sequence_data (char **buffer, long *bufsize, long *allocbufsize, data_fmt * data, long pop,
                    long ind, long locus)
{
    const long sublocistart = data->sublocistarts[locus];
    const long sublociend   = data->sublocistarts[locus+1];
    long sublocus;
    //boolean done = FALSE;
    for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
      {
	mutationmodel_fmt *s = &data->mutationmodels[sublocus];
	long site;
	//printf("(%li:",sublocus);
	for(site=0; site < s->numsites; site++)
	  {
	    if (*bufsize > *allocbufsize-s->numsites)
	      {
		*allocbufsize += s->numsites+2;
		*buffer = (char*) myrealloc(*buffer, *allocbufsize * sizeof(char));
	      }
	    *bufsize += sprintf (*buffer + *bufsize, "%s", data->yy[pop][ind][sublocus][0][site]);
	  }
	*bufsize += sprintf (*buffer + *bufsize, "\n"); // separates subloci
      }
}


// this function and get_data() do not mix well!
void
unpack_databuffer (data_fmt * data, option_fmt * options, world_fmt *world)
{
  //long len;
    long locus, pop, ind, i=0;
    long biggest;
#ifdef UEP
    long sumtips;
#endif
    long genomes;
    char *buffer;
    char *buf;
    char *input;
    long inputsize = LONGLINESIZE;
    char *name;
    long hasghost;
    long bufsize;
    long allocbufsize=LONGLINESIZE;
    buffer = (char *)   mycalloc(allocbufsize,sizeof(char));
    input = (char *) mycalloc (LONGLINESIZE, sizeof (char));
    name = (char *)  mycalloc (LONGLINESIZE, sizeof (char));
    // receive header info
    MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
    if (bufsize > allocbufsize)
      {
	allocbufsize = bufsize;
	buffer = (char*) myrealloc(buffer, allocbufsize * sizeof(char));
      }
    MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);
    buf = buffer;
    sgets_safe (&input, &inputsize, &buf);
    sscanf (input, "%c%li%li%li%li", &options->datatype, &hasghost, &data->numpop,
            &data->loci, &data->allsubloci);
    data->hasghost = (boolean) hasghost;
#ifdef DEBUG
    //printf("datatype=%c numpop=%li loci=%li\n",options->datatype,data->numpop,data->loci);
    printf("\nreading datatype etc: %c %li %li %li %li (input=%s)\n", options->datatype, (long) data->hasghost,
	   data->numpop, data->loci, data->allsubloci, input);
#endif
    if(options->datatype == '@')
      {
	data->oneliner = TRUE;
      }
    genomes = number_genomes(options->datatype);
 
    init_data_structure1 (&data);

    mutationmodel_fmt *s;

    if(data->subloci==NULL)
      {
	data->subloci = (long *) mycalloc(data->allsubloci,sizeof(long));
      }
    if(world->mutationmodels == NULL)
      {
	init_mutationmodel_first(world, data);
      }
    // receive site and loci information
    MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
    if (bufsize > allocbufsize)
      {
	allocbufsize = bufsize;
	buffer = (char*) myrealloc(buffer, allocbufsize * sizeof(char));
      }
    MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);    
    buf = buffer;
    for (locus = 0; locus < data->loci; locus++)
      {
	sgets_safe (&input, &inputsize, &buf);
	sscanf (input, "%li%li", &data->subloci[locus], &world->sublocistarts[locus]);
#ifdef DEBUG
	printf("%i> reading %li %li (input=%s)\n", myID, data->subloci[locus], world->sublocistarts[locus],input);
#endif
      }
    sgets_safe (&input, &inputsize, &buf);
    sscanf (input, "%li", &world->sublocistarts[locus]);
#ifdef DEBUG
    printf("%i> reading %li ----  last sublocistarts \n", myID, world->sublocistarts[locus]);
#endif
    for (locus = 0; locus < data->loci; locus++)
      {
	world->numsubloci[locus] = world->sublocistarts[locus+1] - world->sublocistarts[locus];
      }
    double fracchange;
    int dataclass;
    int scaling;
    long estimateseqerror;
    for (locus = 0; locus < data->allsubloci; locus++)
    {
      s = &(data->mutationmodels[locus]);
      sgets_safe (&input, &inputsize, &buf);
      sscanf (input,  "%c %i %i %li %li\n",
	      &s->datatype, &dataclass, &s->model, &s->numpatterns, &s->numsites);
      init_mutationmodel_readsites3(s, s->datatype, s->numsites);

      sgets_safe (&input, &inputsize, &buf);
      sscanf (input,  "%li %li %li %lf\n",
	      &s->startsite, &s->numstates, &s->numsiterates,&s->lambda); 
	
      sgets_safe (&input, &inputsize, &buf);
      int frominfile;
      int finished;
      sscanf (input,  "%lf %lf %lf %i %i\n",
	      &s->parameters[0],&s->parameters[1],&s->parameters[2],&frominfile, &finished);
      s->from_infile = (boolean) frominfile;
      s->finished = (boolean) finished;
#ifdef DEBUG
      //      printf("@@@@@@@@@@@@@@@@@@@@@@@@ %i> %p %i: 0:%f 1:%f 2:%f \n",myID, s, s->model, s->parameters[0],s->parameters[1],s->parameters[2]);
#endif	      
      sgets_safe (&input, &inputsize, &buf);
      sscanf (input,  "%li %lf %lf %lf %lf\n",
	      &s->weightsum, &s->xi, &s->xv, &s->ttratio, &fracchange);
      
      sgets_safe (&input, &inputsize, &buf);
      sscanf (input,  "%lf %lf %li %lf\n",
	      &s->freq, &s->freqlast,&s->maxalleles,&s->browniandefault);

      //if(s->datatype == 'b')
      //s->maxalleles = XBROWN_SIZE;
#ifdef DEBUG
      printf("@123456@@@@@@@@@@@@@@@@@ %i> datatype=%c maxalleles=%li \n",myID, s->datatype, s->maxalleles);
#endif	      

      sgets_safe (&input, &inputsize, &buf);
      sscanf (input,  "%li %li %li\n",
	      &s->micro_threshold, &s->microrange,&s->microstart);
	
      sgets_safe (&input, &inputsize, &buf);
      sscanf (input,  "%li %li %li %i %i\n",
		&s->addon, &s->numcategs,&estimateseqerror,&s->seqerrorcombined,&scaling);
	
      sgets_safe (&input, &inputsize, &buf);
      sscanf (input,  "%lf %lf %lf %lf\n", 
	      &s->scoring_error[0], &s->scoring_error[1], &s->scoring_error[2], 
	      &s->scoring_error[3]);
      
      
      s->dataclass = (dataclass==0) ? SITECHARACTER : SITEWORD;
      s->scaling = (scaling==0) ? FALSE : TRUE;
      s->estimateseqerror = (estimateseqerror==0) ? FALSE : TRUE;
      init_mutationmodel_readsites3(s, s->datatype, s->numsites);
      for(i=0;i<s->numsites;i++)
	{
	  sgets_safe (&input, &inputsize, &buf);
	  sscanf(input,"%lf%li", &s->aliasweight[i], &s->alias[i]);
	}
      s->numsiterates = 0; //reread of site rates from options and init
      set_siterates(locus,world,options);
      
#ifdef DEBUG
      printf("%i> locus=%li reading %c %li %li %li (input=%s)\n", myID, locus, s->datatype, s->numsites, s->startsite, s->addon,input);
#endif
      s->fracchange = (MYREAL) fracchange;
    }
    // population data
    for (pop = 0; pop < data->numpop; pop++)
    {
      // receive each population name
      MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
      if (bufsize > allocbufsize)
	{
	  allocbufsize = bufsize;
	  buffer = (char*) myrealloc(buffer, allocbufsize * sizeof(char));
	}
      MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);    
      buf = buffer;
      sgets_safe (&input, &inputsize, &buf);
      sscanf (input, "%s", data->popnames[pop]);
      biggest=0;
      // receive ind numbers for each locus
      MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
      if (bufsize > allocbufsize)
	{
	  allocbufsize = bufsize;
	  buffer = (char*) myrealloc(buffer, allocbufsize * sizeof(char));
	}
      MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);    
      buf = buffer;
#ifdef DEBUG
      //printf("%i>@@@%s@@@\n",myID,buffer);
#endif 
      for(locus=0; locus<data->loci; locus++)
        { 
	  sgets_safe (&input, &inputsize, &buf);
	  sscanf (input, "%li %li", &data->numind[pop][locus], &data->numalleles[pop][locus]);
#ifdef DEBUG
	  printf("%i>###%s@@@\n%i>numind receive: %li %li %li %li\n", myID,input,myID, pop, locus, data->numind[pop][locus],data->numalleles[pop][locus]);
#endif
	  if(biggest<data->numind[pop][locus])
	    biggest = data->numind[pop][locus];
	  //	    data->numalleles[pop][locus] = data->numind[pop][locus] * genomes;
        }
      init_data_structure2 (&data, options, world, pop);
      // receive individual data per pop
#ifndef MPIDATAONDEMAND    
        if (!strchr (SEQUENCETYPES, options->datatype) && options->datatype!='@')
        {
            for (ind = 0; ind < biggest; ind++)
            {
	      MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
	      if (bufsize > allocbufsize)
		{
		  allocbufsize = bufsize;
		  buffer = (char*) myrealloc(buffer, allocbufsize * sizeof(char));
		}
	      MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);    
	      buf = buffer;
	      sgets_safe (&input, &inputsize, &buf);
	      sscanf (input, "%s", data->indnames[pop][ind][0]);
	      for (locus = 0; locus < data->loci; locus++)
                {
		  sgets_safe (&input, &inputsize, &buf);
		  sscanf (input, "%s %s", data->yy[pop][ind][locus][0][0],
			  data->yy[pop][ind][locus][1][0]);
                }
            }
        }
        else
	  {
            for (locus = 0; locus < data->loci; locus++)
	      {
#ifdef DEBUG
		printf("reading over all individuals in locus %li: %li\n", locus, data->numind[pop][locus]);
#endif
                for (ind = 0; ind < data->numind[pop][locus]; ind++)
		  {
		    MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
		    if (bufsize > allocbufsize)
		      {
			allocbufsize = bufsize;
			buffer = (char*) myrealloc(buffer, allocbufsize * sizeof(char));
		      }
		    MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);    
		    buf = buffer;
		    sgets_safe (&input, &inputsize, &buf);
		    sprintf(data->indnames[pop][ind][locus],"%-*s", (int) options->nmlength, input);	
		    const long sublocistart = data->sublocistarts[locus];
		    const long sublociend   = data->sublocistarts[locus+1];
		    long sublocus;
		    for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
		      {
			mutationmodel_fmt *s = &data->mutationmodels[sublocus];
			// inputsize = 100 + s->numsites;
			// input =(char *) myrealloc (input, sizeof (char) * inputsize);
			// sgets_safe (&input,&inputsize, &buf);
			long site;
#ifdef DEBUG
			printf("%i> reading data %li %li %li: ",myID, ind, sublocus, s->numpatterns);
#endif
			//printf("reading sites %li: ",s->numsites);
			if(s->dataclass==SITECHARACTER)
			  {
			    for(site=0; site < s->numsites;site++)
			      {
				char c = *buf++;
				if (!(isspace(c) || isdigit(c)))
				  {
				    data->yy[pop][ind][sublocus][0][site][0]=c; //input[site];
				    //printf("%i> -%c- %li\n",myID,c, site);
				  }
				else
				  site--;
				//if(site < MIN(22,s->numsites))
				//  printf("%c",data->yy[pop][ind][sublocus][0][site][0]);
			      }
			    //printf(" ----> %22.22s\n",input);
			  }
			else
			  {
			    sgets_safe (&input, &inputsize, &buf);
			    strcpy( data->yy[pop][ind][sublocus][0][0],input);
			    printf("%s <--(input:%s) ---> Word data",data->yy[pop][ind][sublocus][0][0],input);
			  }
		      }
		  }
	      }
	  }
#endif /*MPIDATAONDEMAND*/    
    }
    // geofile
    data->geo =
        (MYREAL *) mycalloc (1, sizeof (MYREAL) * data->numpop * data->numpop);
    data->lgeo =
        (MYREAL *) mycalloc (1, sizeof (MYREAL) * data->numpop * data->numpop);
    if (!options->geo)
    {
        for (i = 0; i < data->numpop * data->numpop; i++)
            data->geo[i] = 1.0;
    }
    else
    {
      // receive geofile
      MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
      if (bufsize > allocbufsize)
	{
	  allocbufsize = bufsize;
	  buffer = (char*) myrealloc(buffer, allocbufsize * sizeof(char));
	}
      MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);    
      buf = buffer;
      for (pop = 0; pop < data->numpop * data->numpop; pop++)
        {
	  sgets_safe (&input, &inputsize, &buf);
#ifdef USE_MYREAL_FLOAT
	  sscanf (input, "%f%f", &data->geo[pop], &data->lgeo[pop]);
#else
	  sscanf (input, "%lf%lf", &data->geo[pop], &data->lgeo[pop]);
#endif
        }
    }
    
  // receive randomsubset
    data->shuffled = (long ***) mycalloc(data->numpop,sizeof(long **));
    MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
    if (bufsize > allocbufsize)
      {
	allocbufsize = bufsize;
	buffer = (char*) myrealloc(buffer, allocbufsize * sizeof(char));
      }
    MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);    
    buf = buffer;
    for (pop = 0; pop < data->numpop; pop++)
    {
      data->shuffled[pop] = (long **) mycalloc(data->loci,sizeof(long *));
      for(locus=0;locus<data->loci;locus++)
	{
	  data->shuffled[pop][locus] = (long *) mycalloc(data->numind[pop][locus],sizeof(long));
	  for(ind=0; ind < data->numind[pop][locus]; ind++)
	    {
	      sgets_safe (&input, &inputsize, &buf);
	      sscanf (input, "%li", &data->shuffled[pop][locus][ind]);
#ifdef DEBUG
	      // printf("%li-", data->shuffled[pop][locus][ind]);
#endif
	    }
#ifdef DEBUG
	  // printf("received shuffled\n");
#endif
	}
    }
  
    //receive date data
    data->sampledates = (tipdate_fmt ***) mycalloc (data->numpop, sizeof (tipdate_fmt **));
    data->sampledates[0] = (tipdate_fmt **) mycalloc (data->numpop * data->loci, sizeof (tipdate_fmt *));
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
    data->maxsampledate = 0.0;
    if (options->has_datefile)
      {
	MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
	if (bufsize > allocbufsize)
	  {
	    allocbufsize = bufsize;
	    buffer = (char*) myrealloc(buffer, allocbufsize * sizeof(char));
	  }
	MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);    
	buf = buffer;
	if(data->oneliner)
	  {
	    for (pop = 0; pop < data->numpop; pop++)
	      {
		for(ind=0; ind < data->numind[pop][0]; ind++)
		  {
		    //sgets_safe (&input, &inputsize, &buf);
		    sgets(input, LONGLINESIZE, &buf);
		    //printf("%i> %s\n",myID, input);
		    double date;
		    sscanf (input, "%s %lf", name, &date);
		    //printf("%i> new: %li %li %s %f\n",myID, ind, pop, name, date);
		    data->sampledates[pop][0][ind].date = (MYREAL) date;
		    //if (data->sampledates[pop][locus][ind].name == NULL)
		      data->sampledates[pop][0][ind].name = (char *) mycalloc(strlen(name)+1,sizeof(char));
		    strcpy(data->sampledates[pop][0][ind].name, name);
		    for(locus=1;locus<data->loci;locus++)
		      {
			data->sampledates[pop][locus][ind].date = (MYREAL) date;
			//if (data->sampledates[pop][locus][ind].name == NULL)
			  data->sampledates[pop][locus][ind].name = (char *) mycalloc(strlen(name)+1,sizeof(char));
			strcpy(data->sampledates[pop][locus][ind].name, name);
		      }
		  }
	      }
	  }
	else
	  {
	    for(locus=0;locus<data->loci;locus++)
	      {
		for (pop = 0; pop < data->numpop; pop++)
		  {
		    for(ind=0; ind < data->numind[pop][locus]; ind++)
		      {
			sgets (input, LONGLINESIZE, &buf);
			//printf("%i> %s\n",myID, input);
			double date;
			sscanf (input, "%s %lf", name, &date);
			//printf("%i> old: %li %li %s %f\n",myID, ind, pop, name, date);
			data->sampledates[pop][locus][ind].date = (MYREAL) date;
			data->sampledates[pop][locus][ind].name = (char *) mycalloc(strlen(name)+1,sizeof(char));
			strcpy(data->sampledates[pop][locus][ind].name, name);
		      }
		  }
	      }
	  }
	sgets_safe (&input, &inputsize, &buf);
	double dd;
	sscanf(input,"%lf",&dd);
	data->maxsampledate = (MYREAL) dd;
	
	//printf("%i> maxsampledate=%f %s\n''''''\n%s\n''''''''''''\n",myID, data->maxsampledate, input, buffer);
      }
    // uepfile
#ifdef UEP
    if (options->uep)
    {
	MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
	if (bufsize > allocbufsize)
	  {
	    allocbufsize = bufsize;
	    buffer = (char*) myrealloc(buffer, allocbufsize * sizeof(char));
	  }
	MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);    
	buf = buffer;
        sgets_safe (&input, &inputsize, &buf);
        sscanf (input, "%li%li", &sumtips, &data->uepsites);
        data->uep =
            (int **) mycalloc (number_genomes (options->datatype) * sumtips,
                             sizeof (int *));
        if (strchr (SEQUENCETYPES, options->datatype))
        {
            for (pop = 0; sumtips; pop++)
            {
                data->uep[i] = (int *) mycalloc (data->uepsites, sizeof (int));
                for (i = 0; i < data->uepsites; i++)
                {
                    sgets_safe (&input, &inputsize, &buf);
                    sscanf (input, "%i", &data->uep[pop][i]);
                }
            }
        }
        else
        {
            for (pop = 0; sumtips; pop++)
            {
                data->uep[i] = (int *) mycalloc (data->uepsites, sizeof (int));
                data->uep[i + sumtips] =
                    (int *) mycalloc (data->uepsites, sizeof (int));
                for (i = 0; i < data->uepsites; i++)
                {
                    sgets_safe (&input, &inputsize, &buf);
                    sscanf (input, "%i%i", &data->uep[pop][i],
                            &data->uep[pop + sumtips][i]);
                }
            }
        }

    }
#endif
    init_data_structure3 (data, options, world);

    // haplotyping
    if(options->haplotyping)
      {
	MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
	if (bufsize > allocbufsize)
	  {
	    allocbufsize = bufsize;
	    buffer = (char*) myrealloc(buffer, allocbufsize * sizeof(char));
	  }
	MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);    
	printf("hap %i> %s\n",myID,buffer);
	buf = buffer;
	if(data->numindividuals==NULL)
	  data->numindividuals = (long *) calloc(data->loci,sizeof(long));
	for(locus=0;locus<data->loci;locus++)
	  {
	    sgets (input, LONGLINESIZE, &buf);
	    sscanf (input, "%li", &data->numindividuals[locus]);
	  }
	if(data->individuals==NULL)
	  data->individuals = (individualDB_fmt **) calloc(data->loci,sizeof(individualDB_fmt *));
	for(locus=0;locus<data->loci;locus++)
	  {
	    data->individuals[locus] = (individualDB_fmt *) calloc(data->numindividuals[locus],sizeof(individualDB_fmt));
	    long id;
	    for (id=0; id < data->numindividuals[locus]; id++)
	      {
		sgets (input, LONGLINESIZE, &buf);
		sscanf (input, "%li", &data->individuals[locus][id].id);
		sgets (input, LONGLINESIZE, &buf);
		sscanf (input, "%s", data->individuals[locus][id].name);
		sgets (input, LONGLINESIZE, &buf);
		if (data->individuals[locus][id].ind == NULL)
		  data->individuals[locus][id].ind = mycalloc(TWO,sizeof(long));
		sscanf (input, "%li%li%li", &data->individuals[locus][id].ind[0],
			&data->individuals[locus][id].ind[1],
			&data->individuals[locus][id].numstates);
		data->individuals[locus][id].region = locus;
	      }
	  }
	// hook into data->indnames for world for haplotyping  
	world->indnames = data->indnames;
      }
    //
    switch (options->datatype)
    {
    case 'a':
      create_alleles (world, data, options);
        break;
    case 'b':
        for (pop = 0; pop < data->loci; pop++)
	  {
            data->maxalleles[pop] = XBROWN_SIZE;
	    world->mutationmodels[pop].maxalleles = XBROWN_SIZE;
	  }
        break;
    case 'm':
      create_alleles (world, data, options);
        for (pop = 0; pop < data->loci; pop++)
	  {
            data->maxalleles[pop] = options->micro_stepnum;
	  }
        break;
    }
#ifdef DEBUG
    printf("finished data transfer\n");
#endif
    myfree(input);
    myfree(name);
    myfree(buffer);
      }


void request_data(long pop,long ind, long locus, long sublocus,long allelenum, world_fmt *world, data_fmt *data, option_fmt * options, site_fmt ***datapart)
{
  //datapart is a pointer to an array of strings, it is assumed that it has NULL bytes allocated, it will receive
  //its memory from here
  long bufsize=0;
  mutationmodel_fmt *s = &world->mutationmodels[sublocus];
  long site;
  //long patterns = s->numpatterns;
  long numsites = s->numsites;
  MPI_Status  status;
  long tag=myID+ONDEMANDTAG;
  long len;
  char *buffer;
  char *buf;
  char *  p1 = (char *) mycalloc(LINESIZE,sizeof(char));
  long inputsize = LONGLINESIZE;
  char *  input = (char *) mycalloc(inputsize,sizeof(char));
  sprintf(p1,"D %li %li %li %li",pop,ind,sublocus,allelenum);
  //fprintf(stdout,"%i> requested data from master %s (pop,ind,subloc,allelnum)\n", myID, p1);
  //fflush(stdout);
  MYMPISEND (p1, SMALLBUFSIZE, MPI_CHAR, (MYINT) MASTER, (MYINT) tag, comm_world);
  MYMPIRECV (&bufsize, ONE, MPI_LONG, MASTER, (MYINT) tag, comm_world, &status);
  buffer = (char *) mycalloc (bufsize, sizeof (char));
  MYMPIRECV (buffer, bufsize, mpisizeof, MASTER, (MYINT) tag, comm_world, &status);
  buf = buffer;

  if (*datapart == NULL)
    {
#ifdef DEBUG
      printf("%i> allocation datapart: numsites=%li [%li]\n",myID, numsites, ind);
#endif
      *datapart = (site_fmt **) mycalloc (2,sizeof (site_fmt *));
      (*datapart)[0] = (site_fmt *) mycalloc (1+numsites, sizeof (site_fmt));
      if(!strchr(SEQUENCETYPES,options->datatype) && options->datatype!='@')
	(*datapart)[1] = (site_fmt *) mycalloc (1+numsites, sizeof (site_fmt));	  
      if(s->dataclass==SITEWORD)
	{
	  len = (options->allelenmlength+1);
	}
      else
	{
	  len = 2;
	}
      for(site=0;site<numsites;site++)
	{
	  (*datapart)[0][site] = (site_fmt) mycalloc (len, sizeof (char));
	  if(!strchr(SEQUENCETYPES,options->datatype) && options->datatype!='@')
	    (*datapart)[1][site] = (site_fmt) mycalloc (len, sizeof (char));
	}
    }
  else
    {
      error("datapart was not NULL, aborted in request_data()");
    }
  
  if (allelenum!=0)
    {
      sgets_safe (&input, &inputsize, &buf);
      sscanf (input, "%s", data->indnames[pop][ind][0]);
      sgets_safe (&input, &inputsize, &buf);
      sscanf (input, "%s %s", (*datapart)[0][0],(*datapart)[1][0]);
    }
  else
    {
      sgets_safe (&input, &inputsize, &buf);
      sprintf(data->indnames[pop][ind][locus],"%-*s", (int) options->nmlength, input);	

      sgets_safe (&input, &inputsize, &buf);
      if(s->basefreqs==NULL)
	{
	  s->basefreqs = (double*) mycalloc((s->numstates+BASEFREQLENGTH-4), sizeof(double));
	}
      else
	{
	  s->basefreqs = (double*) myrealloc(s->basefreqs,(s->numstates+BASEFREQLENGTH-4) * sizeof(double));
	}
      int tmp=0;
      sscanf(input,"%i %lf %lf %lf %lf", &tmp, &s->basefreqs[0],&s->basefreqs[1],&s->basefreqs[2],&s->basefreqs[3]);
      if (tmp==0)
	s->dataclass = SITECHARACTER;
      else
	s->dataclass = SITEWORD;
      
      if(s->dataclass==SITECHARACTER)
	{
	  for(site=0; site < s->numsites;site++)
	    {
	      char c = *buf++;
	      if (!(isspace(c) || isdigit(c)))
		{
		  (*datapart)[0][site][0]=c; //input[site];
		}
	      else
		site--;
	    }
	}
      else
	{
	  sgets_safe (&input, &inputsize, &buf);
	  strcpy((*datapart)[0][0],input);
	}
    }
  myfree(p1);
  myfree(buffer);
  myfree(input);
}


#if 0
///
/// unpacks results from buffer to fill data structures in master for final printout
void
unpack_plotplane_buffer (MYREAL *buffer, world_fmt * world,
                      long locus, long maxrep, long numpop)
{

}
///
/// pack plotplane into buffer to ship to master
long
pack_plotplane_buffer (MYREAL **buffer, world_fmt * world,
                    long locus, long maxrep, long numpop)
{
}
#endif


///
/// unpacks results from buffer to fill data structures in master for final printout
void
unpack_result_buffer (MYREAL *buffer, world_fmt * world,
                      long locus, long maxrep, long numpop)
{
  long z=0;
  long rep;
  long pop;
  long addon=0;
  long nn = world->numpop2+ world->bayes->mu + world->species_model_size * 2 + world->grownum;
  //  long numpop2 = world->numpop2;
  timearchive_fmt **atl = world->atl;
  MYREAL ***apg0 = world->apg0;

  if (maxrep > 1)
    addon = 1;

  for (rep = 0; rep < maxrep + addon; rep++)
    {
        atl[rep][locus].param_like = buffer[z++];
        for (pop = 0; pop < 4 * nn/*numpop2*/; pop++)//DEBUG PROBLEM why 4?
        {
            atl[rep][locus].parameters[pop] = buffer[z++];
        }
	//fprintf(stderr,"%i> unpacked result locus=%li replicate %li\n",myID,locus, rep);
    }
    // apg0
    for (rep = 0; rep < maxrep; rep++)
    {
        for (pop = 0; pop < world->options->lsteps; pop++)
        {
            apg0[rep][locus][pop] = buffer[z++];
        }
    }
    // BF material
    printf("unpack_result_buffer\n");
    z = unpack_BF_buffer(buffer, z, locus, world);
    // ESS material
    z = unpack_ess_buffer(buffer, z, world);
}

///
/// pack results into buffer to ship to master
long
pack_result_buffer (MYREAL **buffer, world_fmt * world,
                    long locus, long maxrep, long numpop)
{
  long rep;
  long pop;
  long bufsize;
  long z = 0;
  long addon = 0;
  long numpop2 = world->numpop2;
  long nn = numpop2 + world->bayes->mu + world->species_model_size * 2 + world->grownum;
  timearchive_fmt **atl = world->atl;
  MYREAL ***apg0 = world->apg0;
  
  if (maxrep > 1)
    addon = 1;
  
  bufsize = (maxrep+addon) + (maxrep+addon) * 4 * numpop2 + maxrep * world->options->lsteps + \
    world->options->heated_chains * world->loci + 5 * world->loci + 1 + nn * 2 + nn*6 + 2 * world->options->heated_chains;
  (*buffer) = (MYREAL *) myrealloc (*buffer, sizeof (MYREAL) * bufsize);
  memset (*buffer, 0, sizeof (MYREAL) * bufsize);
  
  for (rep = 0; rep < maxrep + addon; rep++)
    {
      (*buffer)[z++] = atl[rep][locus].param_like;
      for (pop = 0; pop < 4 * nn; pop++)
	{
	  (*buffer)[z++] =  atl[rep][locus].parameters[pop];
	}
    }
  // apg0
  for (rep = 0; rep < maxrep; rep++)
    {
      for (pop = 0; pop < world->options->lsteps; pop++)
        {
	  (*buffer)[z++] =  apg0[rep][locus][pop];
        }
    }
  // BF material
  if(!world->data->skiploci[locus])
    {
      z = pack_BF_buffer(buffer, z, locus, world);
    }
  // ESS material
  z = pack_ess_buffer(buffer, z, world);
  // hyper material
  z = pack_hyper_buffer(buffer, z, world);
#ifdef DEBUG_MPI
  fprintf(stdout,"DEBUG: %i> z=%li, bufsize=%li\n", myID, z, bufsize);
#endif
  if(bufsize >= z)
    {
      bufsize = z;
    }
  else
    {
      fprintf(stderr,"%i> bufsize=%li < z=%li\n",myID,bufsize,z);
      error("pack_results tried to stuff to much into buffer in pack_result_buffer()\n");
    }
  return bufsize;
}




///
/// unpacks replicate samples of migration events, adds numbers and part-vectors
/// to the final array per locus, this function is only used with replicates over
/// multiple loci.
void
unpack_mighist_replicate_buffer_old (MYREAL *buffer, world_fmt * world,
                       long locus, long numpop)
{
  long i; 
  long j;
  long z = 0;
  long nummighist;
  long nummighistold;
  mighistloci_fmt *aa;
  aa = &world->mighistloci[locus];
  nummighist = (long) buffer[z++];
  nummighistold = aa->mighistnum;
  aa->mighistnum += nummighist;
  if(aa->allocsize <= aa->mighistnum)
    {
      aa->mighist = (mighist_fmt *) myrealloc (aa->mighist, sizeof (mighist_fmt) *(aa->mighistnum+1));
      for(j=aa->allocsize; j<=aa->mighistnum; j++)
        {
	  aa->mighist[j].allocsize=1;
	  aa->mighist[j].migeventsize=0;
	  aa->mighist[j].migevents =
	    (migevent_fmt *) mycalloc (1,  sizeof (migevent_fmt) );
        }
      aa->allocsize = aa->mighistnum+1;
    }
  for (j = nummighistold; j < aa->mighistnum; j++)
    {
      aa->mighist[j].copies = (long) buffer[z++];
      aa->mighist[j].weight = (long) buffer[z++];
      aa->mighist[j].migeventsize = (long) buffer[z++];
      aa->mighist[j].allocsize = aa->mighist[j].migeventsize + 1;
      aa->mighist[j].migevents = (migevent_fmt *) myrealloc (aa->mighist[j].migevents,
							     sizeof (migevent_fmt) *
							     aa->mighist[j].allocsize);
      for (i = 0; i < aa->mighist[j].migeventsize; i++)
        {
	  aa->mighist[j].migevents[i].age = buffer[z++];
	  aa->mighist[j].migevents[i].from = (long) buffer[z++];
	  aa->mighist[j].migevents[i].to = (long) buffer[z++];
	  aa->mighist[j].migevents[i].sumlines = (long) buffer[z++];
        }
    }
}
///
/// unpacks replicate samples of migration events, adds numbers and part-vectors
/// to the final array per locus, this function is only used with replicates over
/// multiple loci.
void
unpack_mighist_replicate_buffer(MYREAL *buffer, world_fmt * world,
                       long locus, long numpop)
{
  long i;
  long pop;
  long z         = 0;
  long          oldeventbinnum;
  long         * eventbinnum = NULL;
  duo         ** eventbins;
  mighistloci_fmt *aa;
  long npall = world->numpop2 + world->bayes->mu + world->species_model_size * 2 + world->grownum;

  aa = &world->mighistloci[locus];
  eventbins = aa->migeventbins;
  eventbinnum = aa->migeventbinnum;
  

  for (pop = (world->options->mighist_all ? 0 : world->numpop); 
       pop <  npall; pop++)
    {
      oldeventbinnum = eventbinnum[pop];
      eventbinnum[pop] = (long) buffer[z++];
      if(oldeventbinnum < eventbinnum[pop])
	{
	  aa->migeventbins[pop] = (duo *) myrealloc(aa->migeventbins[pop], sizeof(duo) * eventbinnum[pop]);
	  memset(aa->migeventbins[pop][oldeventbinnum], 0, sizeof(duo) * (eventbinnum[pop] - oldeventbinnum));
	}
      for (i = 0; i < eventbinnum[pop]; i++)
	{
	  eventbins[pop][i][0] += buffer[z++];
	  eventbins[pop][i][1] += buffer[z++];
	}
    }
}


///
/// receive buffer with treespace content and unpack it.
/// Messing with format assignment to fit into the result_worker/master scheme
void
unpack_treespace_buffer (MYREAL *buffer, world_fmt * world,
                       long locus, long maxrep, long numpop)
{
  MYREAL like = -HUGE;
    char *input;
    char *sbuf;
    char *buf;
    long size  = strlen((char *) buffer)+1;
    buf = NULL;
    if(size<=1)
      {
	return;
      } 
    buf = (char *) mycalloc(size+1,sizeof(char));
    sbuf = buf;
    memcpy(buf, buffer,sizeof(char)*(size));
    if(world->options->treeprint==BEST)
      {
	input = strsep(&buf,"@");
	//printf("%i> in unpack_treespace_buffer() after strsep %s\n",myID,input);fflush(stdout);
	like = atof(input);
	if(like >= world->besttreelike[locus])
	  {
	    world->besttreelike[locus] = like;
	    //	if(world->treespacenum[locus] < size)
	    //  {
	    //printf("%i> in unpack_treespace_buffer() before treespace realloc  %s\n%s\n",myID,input,buf);
	    size = strlen(buf) + 1;
	    world->treespacealloc[locus] = size;
	    world->treespacenum[locus] = size;
	    world->treespace[locus] = (char *) myrealloc(world->treespace[locus],
							 sizeof(char)*(size+1));
	    strcpy(world->treespace[locus],buf);
	  }
      }
    else
      {
	size = strlen(buf) + 1;
	world->treespacealloc[locus] = size;
	world->treespacenum[locus] = size;
	world->treespace[locus] = (char *) myrealloc(world->treespace[locus],
						     sizeof(char)*(size+1));
	strcpy(world->treespace[locus],buf);
      }
    myfree(sbuf);
}

///
/// pack treespace into buffer
/// ship to master; messing with format assignment because the standard transport buffer
/// is a double 
long pack_treespace_buffer (MYREAL **buffer, world_fmt * world,
			    long locus, long maxrep, long numpop)
{
  long thissize = strlen(world->treespace[locus])+1;
  long thisrealsize;
  char *input;
  char *ptr2, *ptr3;
  long pos = 0;
  MYREAL like = -HUGE;
  if(world->treespace[locus]==NULL)
    return 0;
  ptr2 = strstr(world->treespace[locus],"=");
  ptr3 = strstr(world->treespace[locus]," ]\n");
  pos  = ptr3-ptr2;
  
  input = mycalloc(LONGLINESIZE,sizeof(char));
  if(world->options->treeprint == BEST)
    {
      sprintf(input,"%-*s", (int) pos, ptr2+1);
      like = atof(input);
      pos   = sprintf(input,"%f @", like);
      // fprintf(stdout,"%i> in pack_treespace_buffer() filling %li size\n",myID, thissize);
      // the master routine for this expects a MYREAL *buffer, but
      // the tree is a string
      thissize += pos;
    }
  thisrealsize = (long) (1. + (thissize * sizeof(char) / sizeof(MYREAL)));
  (*buffer) = (MYREAL *) myrealloc(*buffer, (1+thisrealsize) * sizeof(MYREAL));
  memset(*buffer, 0, sizeof(MYREAL) * (1+thisrealsize));
  // the sizeof(char) is NO MISTAKE!
  sprintf((char*)(*buffer), "%s%s", input, world->treespace[locus]);
  return thisrealsize;
}


void
unpack_mighist_buffer (MYREAL *buffer, world_fmt * world,
                       long locus, long maxrep, long numpop)
{
  long i;
  long pop;
  long z         = 0;
  long           oldeventbinnum;
  long         * eventbinnum = NULL;
  duo         ** eventbins;
  mighistloci_fmt *aa;
  long npall = world->numpop2 + world->bayes->mu + world->species_model_size * 2 + world->grownum;
  aa = &world->mighistloci[locus];
  eventbins = aa->migeventbins;
  eventbinnum = aa->migeventbinnum;
  

  for (pop = (world->options->mighist_all ? 0 : world->numpop); 
       pop <  npall; pop++)
    {
      oldeventbinnum = eventbinnum[pop];
      eventbinnum[pop] = (long) buffer[z++];
      if(oldeventbinnum < eventbinnum[pop])
	{
	  aa->migeventbins[pop] = (duo *) myrealloc(aa->migeventbins[pop], sizeof(duo) * eventbinnum[pop]);
	  memset(aa->migeventbins[pop][oldeventbinnum], 0, sizeof(duo) * (eventbinnum[pop] - oldeventbinnum));
	}

      for (i = 0; i < eventbinnum[pop]; i++)
	{
	  eventbins[pop][i][0] = buffer[z++];
	  eventbins[pop][i][1] = buffer[z++];
	}
    }
}


///
/// pack migration events and coalescence events into  buffer to 
/// ship to master
long pack_mighist_buffer_old (MYREAL **buffer, world_fmt * world,
			  long locus, long maxrep, long numpop)
{
  long i;
  long j;
  long bufsize = 1;
  long z = 0;
  mighistloci_fmt *aa;
  long thin = 1;
  MYREAL *tmp;
  aa = &world->mighistloci[locus];
  for (j = 0; j < aa->mighistnum; j++)
    {
      //aa->mighist[j].migeventsize = 1;
      bufsize += aa->mighist[j].migeventsize;
    }
  // trial code to avoid a crash due to not enough memory available for buffer
  if((tmp = (MYREAL *) myrealloc ((*buffer), sizeof (MYREAL) * (bufsize+1))) == NULL)
    {
      do
	{
	  thin *= 10;
	  warning("Recovery mode for the migration events: thinning the events list by a factor of %li",thin);
	  bufsize = 0;
	  for (j = 0; j < aa->mighistnum; j+=thin)
	    {
	      if(j<aa->mighistnum)
		bufsize += aa->mighist[j].migeventsize;
	    }
	  tmp = (MYREAL *) myrealloc ((*buffer), sizeof (MYREAL) * (bufsize+1));
	} 
      while(tmp == NULL || thin < 10000);
    }
  (*buffer) = tmp;
  memset (*buffer, 0, sizeof (MYREAL) * (bufsize+1));
  
  (*buffer)[z++] = (MYREAL)(aa->mighistnum/thin);
  for (j = 0; j < aa->mighistnum; j+=thin)
    {
      if(j < aa->mighistnum)
	{
#ifdef DEBUG_MPI
	  printf ("%i> packmighistbuffer: %li %li %li\n", myID, aa->mighist[j].copies,
		  aa->mighist[j].weight, aa->mighist[j].migeventsize);
#endif
	  (*buffer)[z++] = (MYREAL) aa->mighist[j].copies;
	  (*buffer)[z++] = (MYREAL) aa->mighist[j].weight;
	  (*buffer)[z++] = (MYREAL) aa->mighist[j].migeventsize;
	  
	  for (i = 0; i < aa->mighist[j].migeventsize; i++)
	    {
	      (*buffer)[z++] = aa->mighist[j].migevents[i].age;
	      (*buffer)[z++] = (MYREAL) aa->mighist[j].migevents[i].from;
	      (*buffer)[z++] = (MYREAL) aa->mighist[j].migevents[i].to;
	      (*buffer)[z++] = (MYREAL) aa->mighist[j].migevents[i].sumlines;
	    }
	}
    }
  if(bufsize < z)
    {
      fprintf(stderr,"%i> bufsize=%li < z=%li\n",myID,bufsize,z);
      error("pack_mighist_buffer() failed\n");
    }
  return z;
}
///
/// pack migration events and coalescence events into  buffer to 
/// ship to master
long pack_mighist_buffer (MYREAL **buffer, world_fmt * world,
			  long locus, long maxrep, long numpop)
{
  long i;
  long pop;
  long z         = 0;
  long bufsize   = 0;
  long         * eventbinnum = NULL;
  duo         ** eventbins;
  mighistloci_fmt *aa;
  long npall = world->numpop2 + world->bayes->mu + world->species_model_size * 2 + world->grownum;

  aa = &world->mighistloci[locus];
  eventbins = aa->migeventbins;
  eventbinnum = aa->migeventbinnum;
  for (pop = (world->options->mighist_all ? 0 : world->numpop); 
       pop <  npall; pop++)
    {
      bufsize += 1 + 2 * eventbinnum[pop];
    }
  (*buffer) = (MYREAL *) myrealloc ((*buffer), sizeof (MYREAL) * (bufsize+1));

  for (pop = (world->options->mighist_all ? 0 : world->numpop); 
       pop <  npall; pop++)
    {
      (*buffer)[z++] = (MYREAL)(eventbinnum[pop]);
      for (i = 0; i < eventbinnum[pop]; i++)
	{
	  (*buffer)[z++] = (MYREAL) eventbins[pop][i][0];
	  (*buffer)[z++] = (MYREAL) eventbins[pop][i][1];
	}
    }
  if(bufsize < z)
    error("pack_mighist_buffer() has a memory problem");
  return z;
}

///
/// receive buffer with skyline content and unpack it
void
unpack_skyline_buffer (MYREAL *buffer, world_fmt * world,
                       long locus, long maxrep, long numpop)
{
  //const   MYREAL invmaxrep = 1./((world->options->replicate
  //                       && world->options->replicatenum >
  //				  0) ? world->options->replicatenum : 1);
  //const long np2 = numpop * numpop;
  long i;
  long j;
  long templ;
  long allocsize;
  long z   = 0;
  long *receive_eventbinnum;
  MYREAL temp;
  mighistloci_fmt *aa;
  long npall = world->numpop2 + world->bayes->mu + world->species_model_size * 2 + world->grownum;
  receive_eventbinnum = (long *) mycalloc(npall, sizeof(long));

  aa = &world->mighistloci[locus];
  // read the number of bins for each parameter
  for(j=0; j< npall; j++)
    {
      templ = (long) buffer[z++];
      receive_eventbinnum[j] = templ;
      if(templ >= aa->eventbinnum[j])
	{
	  allocsize = templ + 1;
	  aa->eventbins[j] = (tetra *) myrealloc (aa->eventbins[j], sizeof (tetra) * allocsize);
	  //memset(aa->eventbins[j] + aa->eventbinnum[j],0,sizeof(tetra)*(templ-aa->eventbinnum[j]));
	  for(i=aa->eventbinnum[j];i < allocsize; i++)
	    {
	      aa->eventbins[j][i][0] = 0.F;
	      aa->eventbins[j][i][1] = 0.F;
	      aa->eventbins[j][i][2] = 0.F;
	      aa->eventbins[j][i][3] = 0.F;
	      aa->eventbins[j][i][4] = 0.F;
	      aa->eventbins[j][i][5] = 0.F;
	      aa->eventbins[j][i][6] = 0.F;
	      aa->eventbins[j][i][7] = 0.F;
	      aa->eventbins[j][i][8] = 0.F;
	    }
	  aa->eventbinnum[j] = allocsize;
	  //debug_skyline(world,"increased eventbins in unpack-skyline");
	}
    }
  // read time width of bins
  temp = (MYREAL) buffer[z++];
  if(temp - world->options->eventbinsize > EPSILON)
    error("problem with bins size transmission in unpack_skyline...");
  for(j=0; j< npall; j++)
    {
      for (i = 0; i < receive_eventbinnum[j]; i++)
	{
	  aa->eventbins[j][i][0] += (float) buffer[z++];// * invmaxrep;
	  aa->eventbins[j][i][1] += (float) buffer[z++];// * invmaxrep;
	  aa->eventbins[j][i][2] += (float) buffer[z++];// * invmaxrep;
	  aa->eventbins[j][i][3] += (float) buffer[z++];// * invmaxrep;
	  aa->eventbins[j][i][4] += (float) buffer[z++];// * invmaxrep;
	  aa->eventbins[j][i][5] += (float) buffer[z++];// * invmaxrep;
	}
    }
    //debug_skyline(world,"after unpack_skyline_buffer");
    myfree(receive_eventbinnum);
}

///
/// pack skyline into  buffer to 
/// ship to master
long pack_skyline_buffer (MYREAL **buffer, world_fmt * world,
                     long locus, long maxrep, long numpop)
{
  long j, i;
  long z = 0L;
  mighistloci_fmt *aa;
  long npall = world->numpop2 + world->bayes->mu + world->species_model_size * 2 + world->grownum;
  long bufsize = npall + 1;
  aa = &world->mighistloci[locus];
  //printf("Buffer in pack_skyline()=%f (%p)(%p)\n",(*buffer)[0], (*buffer), buffer);
  for (j = 0; j < npall; j++)
    {
      //      printf("%i> bufsize=%li eventbinnum=%li\n",myID, bufsize, aa->eventbinnum[j]);
      bufsize += 6 * aa->eventbinnum[j];
    }
  //  fprintf(stdout,"%i> bufsize in pack-skyline-buffer %li for locus %li\n", myID, bufsize, locus);
  //myfree(*buffer);
  (*buffer) = (MYREAL *) myrealloc ((*buffer), sizeof (MYREAL) * (bufsize));
  //(*buffer) = (MYREAL *) mycalloc(bufsize,sizeof (MYREAL));
  //  memset (*buffer, 0, sizeof (MYREAL) * (bufsize));
  // record how many bins in for each parameter
  for (j = 0; j < npall; j++)
    {
      (*buffer)[z++] = (MYREAL) aa->eventbinnum[j];
    }
  // time width of bins 
  (*buffer)[z++] = (MYREAL) world->options->eventbinsize;
  for (j = 0; j < npall; j++)
    {
      for (i = 0; i < aa->eventbinnum[j]; i++)
	{
	  (*buffer)[z++] =  (MYREAL) aa->eventbins[j][i][0];
	  (*buffer)[z++] =  (MYREAL) aa->eventbins[j][i][1];
	  (*buffer)[z++] =  (MYREAL) aa->eventbins[j][i][2];
	  (*buffer)[z++] =  (MYREAL) aa->eventbins[j][i][3];
	  (*buffer)[z++] =  (MYREAL) aa->eventbins[j][i][4];
	  (*buffer)[z++] =  (MYREAL) aa->eventbins[j][i][5];
	}
    }
  if(bufsize < z)
    {
      fprintf(stderr,"error: bufsize=%li, z=%li\n", bufsize, z);
      error("pack_skyline_buffer: bufsize is too small");
    }
    return z;
}


///
/// unpack bayes parameters to fit into mpi_results_master()
void
unpack_bayes_buffer (MYREAL *buffer, world_fmt * world,
                       long locus, long maxrep, long numpop)
{
  // this combines the single_bayes and hist_bayes buffer unpacking
  unpack_hist_bayes_buffer(buffer, world->bayes, world, locus); 
}

///
/// pack bayes parameters to fit into mpi_results_worker()
long 
pack_bayes_buffer (MYREAL **buffer, world_fmt * world,
		   long locus, long maxrep, long numpop)
{
  long i;
  long bufsize;
  long z; 
  long sizec       = sizeof(MYREAL);
  //  long numbins     = 0;
  bayes_fmt *bayes = world->bayes;
  long np2         = world->numpop2;
  long npp         = np2 + bayes->mu + world->species_model_size * 2 + world->grownum;
  bayeshistogram_fmt  *hist;
  hist = &(world->bayes->histogram[locus]);
  // memory for pack_single_bayes_buffer_part()
  // locus and numparams + accept and trial of genealogy
  bufsize = npp + 11*npp;
  //O bufsize =  4; 
  for(i=0; i < npp; i++)
    {
      bufsize += (3 * npp * hist->bins[i]);
    }
  
  //  // hist_bayes_buffer:
  //// max buffer memory needed is (npp + 11*npp + (3 * npp * hist->bins[i])
  //for(i=0; i < npp; i++)
  //  {
  //    if (!world->options->has_bayesmdimfile)
  //	{
  //	  //                   datastore           +11
  //	  bufsize += 11; 
  //	  
  //	  if(bayes->map[i][1] != INVALID)
  //	    {
  //	      //for each parameter
  //	      //                   bins number          +1
  //	      //                   bins*3              +3*bins
  //	      bufsize += 1 + 3 * hist->bins[i]; //set50, set95, result per bin 
  //	    }
  //	}
  //  }
  //bufsize += npp*npp;
  // pack_BF_buffer
  //   hmscale, hm    +2
  //   bf: number of heated chains +world->options->heated_chains
  //   steppingstone+scalar    +2*world->options->heated_chains
  bufsize += 4 + 3 * world->options->heated_chains + 1;
  // Autoarchive, ESS buffer:     parameters           +2*(numpop2+loci)
  //                              genealogy            +2
  bufsize += 2*(2*npp + 2);
  //test printf("%i> bufsize=%li (npp=%li, numbins=%li)\n",myID, bufsize, npp, numbins);

  if(bayes->hyperprior)
    {
      bufsize += npp*6;
    }
  (*buffer) = (MYREAL *) myrealloc(*buffer, sizeof(MYREAL) * (bufsize));
  memset (*buffer, 0,sizec * bufsize);
  
  z = pack_single_bayes_buffer_part(buffer,world->bayes,world,locus);

  //printf("%i> z=%li (%li)\n",myID,z, 2 + 2 * (np2+1));
  if(!world->options->has_bayesmdimfile)
    z = pack_hist_bayes_buffer(buffer, hist, world, z);
  // BF material
  if(!world->data->skiploci[locus])
    {
      // fprintf(stderr,"%i> pack_result_buffer(): packed BF result locus=%li replicate %li\n",myID,locus, -1);
      z = pack_BF_buffer(buffer, z, locus, world);
      if(z > bufsize)
	{
	  fprintf(stderr,"%i> ERROR: allocated bufsize=%li is smaller than used bufsize=%li\n",myID, bufsize, z);
	  error("buffer allocation overflowed");
	}
    }
  // ESS material
  z = pack_ess_buffer(buffer, z, world);
  if(z > bufsize)
    {
      fprintf(stderr,"%i> ERROR: allocated bufsize=%li is smaller than used bufsize=%li\n",myID, bufsize, z);
      error("buffer allocation overflowed");
    }
  z = pack_hyper_buffer(buffer,z,world);
  if(z > bufsize)
    {
      fprintf(stderr,"%i> ERROR: allocated bufsize=%li is smaller than used bufsize=%li\n",myID, bufsize, z);
      error("buffer allocation overflowed");
    }
  return z;
}

///
/// unpack the bayes histogram
void unpack_hist_bayes_buffer(MYREAL *buffer, bayes_fmt *bayes, world_fmt *world, long locus)
{
    long                j, i;
    long                pa;
    long                z = 0;
    long                numbins = 0;
    long                pnum;
    long                tmp1, tmp2;
    long                tmplocus;
    bayeshistogram_fmt  *hist;
    long                total = 0;
    long                np2 = world->numpop2; 
    long                npp = np2 + bayes->mu + world->species_model_size * 2 + world->grownum; 
    long                npp11 = 11 * npp;
    long                loci = world->loci > 1 ?  world->loci -1 : 0;
    //
    // begin unpack_single_bayes_buffer_part
    tmplocus = (long) buffer[z++]; 
    pnum = (long) buffer[z++];
#ifdef DEBUG
    const MYREAL        rat = numcpu / (world->maxreplicate * world->loci);
    fprintf (stdout, "%i> unpack_hist_bayes_buffer() received pnum=%li, rat=%f\n", myID, pnum, rat);fflush(stdout);
#endif
    if(tmplocus!=locus)
    {
        bayes->numparams=0;
        locus = tmplocus;
    }

// Cesky Krumlov 2013    printf("%i> original accepted/trials:  [%li,%li,%li,%li,%li]/[%li,%li,%li,%li,%li]\n",
// Cesky Krumlov 2013	   myID, world->accept_archive[0],world->accept_archive[1],world->accept_archive[2],world->accept_archive[3],world->accept_archive[4],world->trials_archive[0],world->trials_archive[1],world->trials_archive[2],world->trials_archive[3],world->trials_archive[4]);

    for (j = 0; j < npp; ++j)
    {
      if(bayes->map[j][1] != INVALID)
	  {
	    world->accept_archive[j] += (tmp1 = (long) buffer[z++]);//*rat);
	    world->trials_archive[j] += (tmp2 = (long) buffer[z++]);//*rat);
#ifdef DEBUG
	    fprintf (stdout, "%i> received (acc %li) (trial %li) => (sumacc %li) (sumtrial %li)\n", myID, tmp1, tmp2, world->accept_archive[j], world->trials_archive[j]);
#endif        
	  }
    }
    // genealogy
    world->accept_archive[j] += (tmp1 = (long) buffer[z++]);//*rat);
    world->trials_archive[j] += (tmp2 = (long) buffer[z++]);//*rat);
#ifdef DEBUG
	    fprintf (stdout, "%i> genealogy received (acc %li) (trial %li) => (sumacc %li) (sumtrial %li)\n", myID, tmp1, tmp2, world->accept_archive[j], world->trials_archive[j]);
#endif        
    // end unpack_single_bayes_buffer_part
    //
    if(!world->options->has_bayesmdimfile)
      {
	bayes->histogram[locus].datastore = (MYREAL *) myrealloc(bayes->histogram[locus].datastore, sizeof(MYREAL) * (npp11));    
	memset(bayes->histogram[locus].datastore,0,sizeof(MYREAL)*(npp11));
	hist = &(bayes->histogram[locus]);
	hist->minima = hist->datastore;    // contains minimal values for each parameter
        hist->maxima = hist->datastore + npp;    // contains maximal values for each parameter
        hist->adjmaxima = hist->datastore + 2*npp;// holds maxima values used in histogram [are smaller than maxima]
	hist->cred50l  = hist->datastore + 3*npp;    // holds 50%-credibility margins (<all lower values>, 
	hist->cred50u = hist->datastore + 4*npp;   //<all high values>)
	hist->cred95l = hist->datastore + 5*npp;    // holds 95%-credibility margins (<all lower values>)
	hist->cred95u = hist->datastore + 6*npp;   //<all high values>)
	hist->modes = hist->datastore + 7*npp;    // holds 95%-credibility margins (<all lower values>, <all high values>)
	hist->medians = hist->datastore + 8*npp;
	hist->means = hist->datastore + 9*npp;            
	hist->stds = hist->datastore + 10*npp;            


	for(i = 0; i < npp; ++i)
	  {
	    if(bayes->map[i][1] != INVALID)
	      {
		hist->bins[i] = (long) buffer[z++];
		total += hist->bins[i];
	      }
	  }
#ifdef DEBUG
	printf("%i> @@@@ locus=%li histbin[0]=%li (hist[loci].bins[0]=%li)\n", myID, locus, hist->bins[0],world->bayes->histogram[loci].bins[0]);
#endif
	hist->binsum = total; 
	// this steps kills poor memory machines [setting results to floats may help a little]
	hist->results = (double *) mycalloc(total * npp, sizeof(double));
	hist->set95 = (char *) mycalloc(total * npp * 2 + 2, sizeof(char));
	hist->set50 = hist->set95 + (total * npp + 1);
	long valids=0;
	for(i = 0; i < npp; ++i)
	  {
	    if(bayes->map[i][1] != INVALID)
	      {
		valids++;
	      }
	  }
	memcpy(hist->datastore,buffer+z,sizeof(MYREAL)*11*npp);
	z += 11*npp;
		//for(j=0;j<11;j++)
		//  hist->datastore[11*i+j] = buffer[z++];
	numbins = 0;
	for(pa=0; pa < npp; pa++)
	  {
	    for(i=0;i<hist->bins[pa];i++)
	      {
		hist->set50[numbins + i] = (buffer[z++] < 1.0 ? '0' : '1'); 
		hist->set95[numbins + i] = (buffer[z++] < 1.0 ? '0' : '1'); 
		hist->results[numbins + i] = (double) buffer[z++];
	      }
	    numbins += hist->bins[pa];
	    //
	    // CHECK
	    world->bayes->histtotal[locus * npp + pa] = hist->bins[pa];
	    //
	  }
	if(hist->covariance==NULL)
	  {
	    doublevec2d(&hist->covariance,npp,npp);
	  }
	for(i=0;i<npp;i++)
	  {
	    for(j=0;j<npp;j++)
	      {
		hist->covariance[i][j] = buffer[z++];
	      }
	  }
      }
    // BF material
#ifdef DEBUG
      printf("unpack_hist_bayes\n");
#endif
    z = unpack_BF_buffer(buffer, z, locus, world);
    // ESS material
    z = unpack_ess_buffer(buffer, z, world);
    // hyperprior
    z = unpack_hyper_buffer(buffer,z, world);
    // non-smoothed data transfer
    if(!world->options->has_bayesmdimfile)
      {
	if (bayes->histogram[locus].results2 == NULL)
	  {
	    bayes->histogram[locus].results2 = (double *) mymalloc(sizeof(double)*bayes->histogram[locus].binsum);
	    memcpy(bayes->histogram[locus].results2,bayes->histogram[locus].results,sizeof(double)*bayes->histogram[locus].binsum);//@@@@@@was float --> is double correct here daetwil
	    if(bayes->histogram[locus].smoothed == NULL)
	      bayes->histogram[locus].smoothed = (boolean *) mycalloc(npp, sizeof(boolean));
	    else
	      memset(bayes->histogram[locus].smoothed,0,sizeof(boolean)*npp);
	    calc_hpd_credibility(world, locus, world->numpop2, world->numpop2 + world->bayes->mu+world->species_model_size*2 + world->grownum);
	  }
      }
}

///
/// Bayes factor material buffer unpacker
long unpack_BF_buffer(MYREAL *buffer, long start, long locus, world_fmt * world)
{
  long i;
  long z = start;
  long hc = world->options->heated_chains;
  MYREAL              temp;
  MYREAL              *ttemp;
  MYREAL              *htemp;
  ttemp = calloc(3 * world->loci + 3 * world->options->heated_chains + 2, sizeof(MYREAL));
  htemp = ttemp + world->loci;
  //atemp = htemp + world->loci;
  //
  // harmonic mean calculation
  temp = buffer[z++]; //hmscale
  if(temp < world->hmscale[locus])
    {      
      if(world->hm[locus]>0.0)
	world->hm[locus] *= EXP(temp - world->hmscale[locus]);
      world->hmscale[locus] = temp;
      //      printf("%i> locus=%li hmscale=%f hm=%f temp=%f\n",myID,locus, world->hmscale[locus], world->hm[locus], temp);
    }
  htemp[locus] = temp;//hmscale store
  world->hm[locus] += EXP(-htemp[locus]+world->hmscale[locus]) * buffer[z++];

  // thermodynamic integration
  for(i=0;i < hc; i++)
    { 
      world->bf[locus * hc + i] +=  buffer[z++];
#ifdef DEBUG 
      printf("%i> ****bf[%li + hc + %li,%li,%li]=%f\n",myID, locus, i,z-start,z, world->bf[locus*hc+i]);
#endif
    }
  // stepping stone calculation
  for(i=0;i < hc; i++)
    {
      double stone = buffer[z++];
      double scalar = buffer[z++];
      long ii = locus * hc + i;
      if(world->steppingstones[ii] !=0.0)
	{
	  world->steppingstones[ii] *= exp(world->steppingstone_scalars[ii]-scalar);
	  world->steppingstone_scalars[ii] = scalar;
	  world->steppingstones[ii] += stone;
	}
      else
	{
	  world->steppingstones[ii] = stone;
	  world->steppingstone_scalars[ii] = scalar;
	}
    }
  myfree(ttemp);
  return z;
}

long unpack_ess_buffer(MYREAL *buffer, long start, world_fmt *world)
{
  long pa;
  long z = start;
  static long n=1;
  //long np2 = world->numpop2;
  long npp = world->numpop2 + ((long) world->bayes->mu)  + world->species_model_size * 2 + world->grownum;

#ifdef DEBUG
  const MYREAL        rat = (MYREAL) numcpu / (world->maxreplicate * world->loci);
  fprintf (stdout, "%i> unpack_ess_buffer() rat=%i / (%li * %li) = %f\n", myID, numcpu, world->maxreplicate, world->loci, rat);fflush(stdout);
#endif
  // unpacking autocorrelation and ess buffer
  for(pa=0; pa < npp; pa++)
    {
      if(world->bayes->map[pa][1] != INVALID)
	world->auto_archive[pa] += (buffer[z++] - world->auto_archive[pa])/n;
    }
  //genealogy
  world->auto_archive[pa] += (buffer[z++] - world->auto_archive[pa])/n;
  n++;
  //    printf("%i>>>>>> autoarchive %f\n",myID, world->auto_archive[0]);
  for(pa=0; pa < npp; pa++)
    {
      if(world->bayes->map[pa][1] != INVALID)
	world->ess_archive[pa] += buffer[z++];// * rat;
    }
  //genealogy
  world->ess_archive[pa] += buffer[z++];// * rat;
  return z;
}

 /// unpacking hyper buffer
long unpack_hyper_buffer(MYREAL *buffer, long start, world_fmt *world)
{
  long i;
  long z = start;
  hyper_fmt *hyper = world->bayes->hyperp;
  const long np2 = world->numpop2;
  const long npp = np2 + ((long) world->bayes->mu) + world->species_model_size * 2 + world->grownum; 
  if (world->bayes->hyperprior)
    {
      for(i=0;i<npp;i++)
	{      
	  MYREAL mean = buffer[z++]; 
	  MYREAL meanstd = buffer[z++]; 
	  MYREAL alpha = buffer[z++]; 
	  MYREAL alphastd = buffer[z++]; 
	  long meann = (long) buffer[z++]; 
	  long alphan = (long) buffer[z++];
	  if(meann > 1 && alphan>1)
	    {
	      onepass_mean_std_end(&mean,&meanstd,&meann);
	      onepass_mean_std_end(&alpha,&alphastd,&alphan);
	      //printf("mean: %f %f %li\n",mean,meanstd,meann);
	      //printf("alpha: %f %f %li\n",alpha,alphastd,alphan);
	      //printf("oldmean: %f %f %li\n",hyper[i].mean,hyper[i].meanstd,hyper[i].meann);
	      //printf("oldalpha: %f %f %li\n",hyper[i].alpha,hyper[i].alphastd,hyper[i].alphan);
	      combine_meanstd(&hyper[i].mean, &hyper[i].meanstd,&hyper[i].meann, mean, meanstd, meann);         
	      combine_meanstd(&hyper[i].alpha, &hyper[i].alphastd,&hyper[i].alphan, alpha, alphastd, alphan);         
	      //printf("newmean: %f %f %li\n",hyper[i].mean,hyper[i].meanstd,hyper[i].meann);
	      //printf("newalpha: %f %f %li\n",hyper[i].alpha,hyper[i].alphastd,hyper[i].alphan);
	    }
	}
    }
  return z;
}
   
///
/// Bayes factor material buffer packer
long pack_BF_buffer(MYREAL **buffer, long start, long locus, world_fmt * world)
{
  // buffer memory needs are (2 + #heatedchains + 2 * #heatedchains)
  long i;
  long ii;
  long z  = start;
  const long hc = world->options->heated_chains;
  (*buffer)[z++] = world->hmscale[locus];
  (*buffer)[z++] = world->hm[locus];
  for(i=0; i < hc; i++)
    {
      (*buffer)[z++] = world->bf[locus * hc + i];
    }
  for(i=0; i < hc; i++)
    {
      (*buffer)[z++] = world->steppingstones[locus * hc + i];
      (*buffer)[z++] = world->steppingstone_scalars[locus * hc + i];
    }
  //printf("%i> locus=%li send hmscale=%f hm=%f bf=%f %f %f %f\n",myID, locus,world->hmscale[locus],world->hm[locus], world->bf[locus*hc],world->bf[locus*hc+1],world->bf[locus*hc+2],world->bf[locus*hc+3]);
#ifdef DEBUG
  printf("%i> packbuffer: z=%li - %li\n",myID,z-hc-2,z);
  for(ii=z-hc-2;ii<z;ii++)
    printf("%f ",(*buffer)[ii]);
  printf("\n");
#endif
  return z;
}

 /// packing autocorrelation and ess buffer
long pack_ess_buffer(MYREAL **buffer, long start, world_fmt *world)
{
  // buffer memory needs are (npp + 1) + (npp+1) (numpop2 + mu + speciesmodel*2)
  long i;
  long z = start;
  const long np2 = world->numpop2;
  const long npp = np2 + ((long) world->bayes->mu) + world->species_model_size * 2 + world->grownum; 
  // (*buffer) = (MYREAL *) myrealloc(*buffer, (start + ((npp + 1) + (npp+1) * npp))*sizeof(MYREAL));
  for(i=0;i<npp;i++)
    {
      if(world->bayes->map[i][1] != INVALID)
	(*buffer)[z++] = world->auto_archive[i];
    }
  // genealogy
  (*buffer)[z++] = world->auto_archive[i];
  for(i=0;i<npp;i++)
    {
      if(world->bayes->map[i][1] != INVALID)
	(*buffer)[z++] = world->ess_archive[i];
    }
  //genealogy
  (*buffer)[z++] = world->ess_archive[i];
  return z;
}


 /// packing hyper buffer
long pack_hyper_buffer(MYREAL **buffer, long start, world_fmt *world)
{
  long i;
  long z = start;
  const long np2 = world->numpop2;
  const long npp = np2 + ((long) world->bayes->mu) + world->species_model_size * 2 + world->grownum;
  if (world->bayes->hyperprior)
    {
      for(i=0;i<npp;i++)
	{      
	  (*buffer)[z++] = world->bayes->hyperp[i].mean;
	  (*buffer)[z++] = world->bayes->hyperp[i].meanstd;
	  (*buffer)[z++] = world->bayes->hyperp[i].alpha;
	  (*buffer)[z++] = world->bayes->hyperp[i].alphastd;
	  (*buffer)[z++] = (MYREAL) world->bayes->hyperp[i].meann;
	  (*buffer)[z++] = (MYREAL) world->bayes->hyperp[i].alphan;
	}
    }
  return z;
}

///
/// pack the bayes histogram
long pack_hist_bayes_buffer(MYREAL **buffer, bayeshistogram_fmt *hist, world_fmt * world, long startposition)
{
  // buffer memory needed is (npp + 11*npp + (3 * npp * hist->bins[i])
  long  j;
  long  i;
  long  npp     = world->numpop2 + world->bayes->mu + world->species_model_size * 2 + world->grownum; 
  long  numbins = 0;
  long  z       = startposition;
  bayes_fmt *bayes = world->bayes;
#ifdef DEBUG_MPI
    printf("%i> pack_hist_bayes_buffer: position=%li last value = %f numparams=%li npp=%li\n",myID, startposition, (z > 0) ? (*buffer)[startposition] : -9999., hist->numparam, npp);
#endif
#ifdef DEBUG
	    printf("%i>",myID);
#endif    
    for(i = 0; i < npp; ++i)
      {
	if(bayes->map[i][1] != INVALID)
          {
	    (*buffer)[z++] = (MYREAL) hist->bins[i];
#ifdef DEBUG
	    printf("%f ",(*buffer)[z-1]);
#endif
	  }
      }
    //    printf("%i> npp=%li, z=%li\n",myID,npp,z);
    // parameter and mu datastore
    memcpy((*buffer)+z,hist->datastore, sizeof(MYREAL)*11*npp);
    z += 11*npp;
    numbins = 0;
    for(i=0; i < npp; i++)
      {
	// bins are zero for "c" and "0" parameters
	//	if(bayes->map[i][1] != INVALID)
        //  {
	    for(j=0;j<hist->bins[i];j++)
	      {
		(*buffer)[z++] = (MYREAL) (hist->set50[numbins + j]=='1' ? 1.0 : 0.0);
		(*buffer)[z++] = (MYREAL) (hist->set95[numbins + j]=='1' ? 1.0 : 0.0);
		(*buffer)[z++] = (MYREAL) hist->results[numbins + j];//@#@#@# was results2
	      }
	    //  }
#ifdef DEBUG
	    printf("<%li elements> ",3*hist->bins[i]);
#endif
	numbins += hist->bins[i];
      }
    // pack covariance matrix
    for(i=0;i<npp;i++)
      {
	for(j=0;j<npp;j++)
	  {
	    (*buffer)[z++] = hist->covariance[i][j];
#ifdef DEBUG
	    printf("%f ",(*buffer)[z-1]);
#endif
	  }
      }
#ifdef DEBUG
	    printf("\n");
#endif

#ifdef DEBUG_MPI
    printf("%i> pack_hist_bayes_buffer: position=%li numbins=%li, last value = %f\n",myID, z, numbins, (*buffer)[z-1]);
#endif
    return z;
}

///
/// pack the halotypebuffer
long pack_haplotypes_buffer(MYREAL **buffer, world_fmt * world,
			    long locus, long maxrep, long numpop)
{
  //send back per locus or per replicate: distinction?
  //loci/regions: new entry
  //send key and count as list: 
  //replicates: addition of counts of existing states, and addition of nonexisiting states
  long z=0;
  long ind;
  long pop;
  long bufsize=0;
  double ratio=1.0;
  if(world->data->haplotyping_report)
    {
      long allocbufsize = LINESIZE;
      char *charbuf = (char *) calloc(allocbufsize, sizeof(char));
      printf("PACK HAPLOTYPES===================================================\n");
      char * newname;
      char * oldname;
      oldname = (char *) calloc(LINESIZE,sizeof(char));
      newname = (char *) calloc(LINESIZE,sizeof(char));
      if(world->haplotypes[locus] == NULL)
	{
	  save_haplotypes(world,locus);
	  reset_haplotypes(world,locus);
	}
      else
	{
	  warning("%i> found a filled haplotypes%li]\n",locus);
	}
      for(pop=0;pop<numpop;pop++)
	{
  for(ind=0; ind<world->data->numind[pop][locus]; ind++)	   
	    {
	      // hooks into data->indnames
	      //printf("%s\n",world->indnames[pop][ind][locus]);
	      long numchar = strcspn(world->indnames[pop][ind][locus],":");
	      strncpy(newname,world->indnames[pop][ind][locus], numchar);
	      newname[numchar]='\0';
	      long ll = strlen(newname);
	      if(!strncmp(oldname,newname,ll))
		continue;
	      else
		strncpy(oldname,newname,ll);
	      //printf("%s   %s\n",oldname,newname);	      
	      long id = find_inDB(world->indnames[pop][ind][locus], locus, world->haplotypes[locus], world->numhaplotypes[locus]);
	      //printf("id=%li\n",id);
	      if( id != -1)
		{
		  print_haplotypes2(&charbuf, &allocbufsize, &bufsize, world, locus, pop, id, FALSE);
		  //printf("%i>PACK:%s\n", myID, charbuf);
		  z++;
		}
	    }
	}
      ratio = (float) sizeof(char)/(float) sizeof(MYREAL);
      *buffer = (MYREAL *) realloc(*buffer, bufsize*ratio*sizeof(MYREAL)); 
      memcpy(*buffer,charbuf,bufsize*sizeof(char));
      myfree(charbuf);
    }
  return bufsize * ratio;
}

///
/// unpack the haplotype buffer
void unpack_haplotypes_buffer(MYREAL *buffer, world_fmt * world,
                                 long locus, long maxrep, long numpop)
{
  long readpos=0;
  char *word;
  char *longword;
  char *charbuf = (char*) buffer;
  //get buflength
  //get number of individuals
  if(world->data->haplotyping_report)
    {
      word = (char *) calloc(LINESIZE,sizeof(char));
      longword = (char *) calloc(LONGLINESIZE,sizeof(char));
      while(charbuf[readpos]!='\0')
	{
	  // find the id of the sample
	  readpos+=read_word_delim(charbuf+readpos,word,":",FALSE);//population label       
	  readpos+=read_word_delim(charbuf+readpos,word,"\t",FALSE); //read individual
	  printf("hap %i> |%s|\n",myID, word);
	  //winner
	  long id = find_inDB(word, locus, world->haplotypes[locus], world->numhaplotypes[locus]);
	  if( id != -1)
	    {                 
	      printf("hap %i> charbuf+%li:|%30.30s|\n",myID, readpos,charbuf+readpos);
	      readpos+=read_word_delim(charbuf+readpos,longword,"\n",FALSE); //read individual
	      individualDB_fmt *winner = &(world->haplotypes[locus][id]);
	      long pos = nondestruct_trim(&longword);
	      while(*(longword+pos)!='\0')
		{
		  pos += read_word_delim(longword+pos,word,":",FALSE);
		  winner->numstates = set_haplotype_states(winner,word);
		  pos+=read_word_delim(longword+pos,word," ",TRUE);		  // read count
		  long count = atol(word);
		  add_state_counter(&winner->hash, &winner->numhash, &winner->total1, 
				    winner->states1, winner->numstates, count);
		}
	    }
	  else
	    {
	      readpos+=read_word_delim(charbuf+readpos,longword,"\n",FALSE);
	    }
	}
      myfree(word);
      myfree(longword);
    }
}

///
/// pack the assign buffer
long pack_assign_buffer(MYREAL **buffer, world_fmt * world,
			long locus, long maxrep, long numpop)
{
  // this is unconventional and may break on some architectures but 
  // simplifies interaction among nodes, assignment and haplotyping
  // need to transfer character strings but all others simply transfer 
  // numbers, char * are translated into reals and converted back
  // alternative would be use structures but they are cunbersome because
  // the number of elements is not fixed in the structures.
  long i,  p;
  unassigned_fmt * v;
  long z=0;
  char zchar[9];
  long dim = world->loci * numpop;
  float ratio = (float) sizeof(char)/(float) sizeof(MYREAL);
  long bufsize = 0;
  char * charbuf = (char *) mycalloc(world->unassignednum*STRSIZE,sizeof(char));
  bufsize = sprintf(charbuf,"         ");
  for (i=1; i<world->unassignednum; i++)
    {
      v = world->unassigned[i];
      bufsize += sprintf(charbuf+bufsize,"\t%s",v->key);
    }

  z = (long)(bufsize * ratio + 5 + (ratio*8));
#ifdef DEBUG
  long startz = z;
#endif
  sprintf(zchar,"%07li ",z);
  memcpy(charbuf,zchar,sizeof(char)*8);
  *buffer = (MYREAL *) realloc(*buffer, (dim*world->unassignednum + z)*sizeof(MYREAL)); 
  memset(*buffer,0, (dim*world->unassignednum + z)*sizeof(MYREAL)); 
  memcpy(*buffer,charbuf,bufsize*sizeof(char));
#ifdef DEBUG
  printf("%i> charbuf: %s\n",myID, charbuf);
#endif
  // we push first all numbers into the buffer, then we add
  // the characters of the names memcopying them into the buffer towards the 
  // and send then a fake number of MYREALS  
  for (i=1; i<world->unassignednum; i++)
    {
	  for (p=0;p<numpop;p++)
	    {
	      (*buffer)[z++] = world->unassigned[i]->probloc[INDIX(numpop,locus,p)];
	    }
	  (*buffer)[z++]  = world->unassigned[i]->accept[locus];
	  (*buffer)[z++]  = world->unassigned[i]->trials[locus];
    }
  //printf("%i> z=%li bufsize*ratio=%f dim*world->assigednum=%li\n",myID, z,bufsize*ratio,dim*world->unassignednum);
  myfree(charbuf);  
  return bufsize*ratio + z;
}

void unpack_assign_buffer(MYREAL *buffer, world_fmt * world,
			  long locus, long maxrep, long numpop)
{
  // read buffer into matrix 
  long i;
  long dim = world->loci * numpop;
  char *charbuf = (char*) (buffer);
  char *key=NULL;
  //boolean done=FALSE;
  long index;
  unassigned_fmt *temp;
  long z=0;
  //long unassignednum;
  // carefully check names etc before assigning values to the database
#ifdef DEBUG
  printf("%i> %s\n",myID, charbuf);
#endif
  key = strsep(&charbuf,"\t");
  z = atol(key);

  while (charbuf!=NULL)
    {
      key = strsep(&charbuf,"\t");
      if (key==NULL)
	break;
      if (key[0]=='\0')
	break;
      index = find_in_unassignedDB(key, world->unassigned);
      if (index==UNKNOWN)
	{
	  temp = (unassigned_fmt *) mycalloc(1,sizeof(unassigned_fmt));
	  temp->key = strdup(key);
	  temp->accept=(long *) mycalloc(world->loci,sizeof(long));
	  temp->trials=(long *) mycalloc(world->loci,sizeof(long));
	  temp->probloc = (MYREAL*) mycalloc(dim,sizeof(MYREAL));
	  temp->index = world->unassignednum;
	  world->unassignednum++;
	  world->unassigned[temp->index-1]->next = temp;
	  world->unassigned = (unassigned_fmt **) realloc(world->unassigned,sizeof(unassigned_fmt*)*world->unassignednum);
	  world->unassigned[temp->index] = temp;
	}
      else
	{
	  temp = world->unassigned[index];
	}
      for(i=0;i<numpop; i++)
	{
	  temp->probloc[INDIX(numpop,locus,i)] += buffer[z++];
	}
      temp->accept[locus] += buffer[z++];
      temp->trials[locus] += buffer[z++];
    }
}

///
/// pack the assign buffer
long pack_seqerror_buffer(MYREAL **buffer, world_fmt * world,
			long locus, long maxrep, long numpop)
{
  boolean is_combined = world->seqerrorcombined;
  long mult = is_combined ? 1 : 4;
  //long i, l, p;
  //unassigned_fmt * v;
  //long z=0;
  long bufsize = world->seqerrorratesnum[locus];
  (*buffer) = (MYREAL *) realloc(*buffer, sizeof(MYREAL)*(1+bufsize*mult));
  (*buffer)[0] = (MYREAL) bufsize;
  memcpy((*buffer)+1, world->seqerrorrates[locus], sizeof(MYREAL)*(bufsize*mult));
  return bufsize;
}

void unpack_seqerror_buffer(MYREAL *buffer, world_fmt * world,
			  long locus, long maxrep, long numpop)
{
  boolean is_combined = world->seqerrorcombined;
  long mult = is_combined ? 1 : 4;
  long i,j;
  long bufsize = (long) buffer[0];
  long oldsize = world->seqerrorratesnum[locus];;
  world->seqerrorratesnum[locus] += bufsize;
  if(world->seqerrorallocnum[locus] < world->seqerrorratesnum[locus])
    world->seqerrorallocnum[locus] += world->seqerrorratesnum[locus];
  world->seqerrorrates[locus] = (MYREAL *) realloc(world->seqerrorrates[locus],sizeof(MYREAL)*world->seqerrorallocnum[locus]*mult);
  printf("%i> unpack seqerror [%li, %li] \n",myID, bufsize, world->seqerrorratesnum[locus]);
  for (i=1,j=oldsize;i<bufsize*mult+1;i++,j++)
    {
      world->seqerrorrates[locus][j] = buffer[i];
    }
}

///
/// unpack bayes parameter buffer, sent from replicant nodes AND lociworker to the master
/// the values will be simply added to the bayes->params, no records of replicates will be done.
long unpack_single_bayes_buffer(MYREAL *buffer,bayes_fmt * bayes, world_fmt * world,long locus)
{
  long i, j;
  long z = 0;
  long pnum;
  long tmp1, tmp2;
  long tmplocus;
  long allocparams = world->bayes->allocparams;
  //    long oldallocparams = world->bayes->allocparams;
  long repstart;
  long repstop;
  long npp = world->numpop2 + (world->bayes->mu)  + world->species_model_size * 2 + world->grownum; 
  long nn = npp + 2;
  set_replicates (world, world->repkind, world->options->replicatenum,
		  &repstart, &repstop);
  tmplocus = (long) buffer[z++];
  pnum = (long) buffer[z++];
  //fprintf (stdout, "%i> received locus=%li (oldlocus=%li) pnum=%li\n", myID, tmplocus, locus, pnum);
  if(tmplocus!=locus)
    world->bayes->numparams=0;
  
  pnum += world->bayes->numparams;
  if(pnum >=world->bayes->allocparams)
    {
      allocparams = pnum + 1;
      world->bayes->params = (MYREAL *) myrealloc(world->bayes->params,sizeof(MYREAL)*allocparams*nn);
    }
  world->bayes->allocparams = allocparams;
  for(i = world->bayes->numparams; i < pnum; ++i)
    {
      // the first element is log(p(d|g)p(g|param))
      (world->bayes->params+(nn*i))[0] = buffer[z++];
      // the second element is log(p(d|g))
      (world->bayes->params+(nn*i))[1] = buffer[z++];
      //fprintf (stdout, "%i> receive params line %li ", myID, i);
      for (j = 2; j < nn; ++j) 
	  {
	    if(bayes->map[j-2][1] != INVALID)
	      (world->bayes->params+(nn*i))[j] = buffer[z++];
	    //fprintf (stdout, "%f ", (world->bayes->params+(nn*i + 1))[j]);
	  }
      //fprintf (stdout, "\n");
    }
  world->bayes->numparams = pnum;
  // acceptance ratios are added to the ones we have already
  // parameter acceptances
  for (j = 0; j < npp; ++j)
    {
      if(bayes->map[j][1] != INVALID)
	{
	  tmp1 = (long) buffer[z++];
	  tmp2 = (long) buffer[z++];
	  world->accept_archive[j] += tmp1;
	  world->trials_archive[j] += tmp2;
	}
    }
  // the last acceptance is the one for the genealogies
  tmp1 = (long) buffer[z++];
  tmp2 = (long) buffer[z++];
  world->accept_archive[j] += tmp1;
  world->trials_archive[j] += tmp2;
  return z;
}

///
/// Pack bayes parameter buffer, sent from replicant nodes AND lociworker to the master
/// the values will be simply added to the bayes->params, no records of specific replicates are kept.
long pack_single_bayes_buffer(MYREAL **buffer, bayes_fmt *bayes, world_fmt *world,long locus)
{
    long i, j;
    long bufsize;
    long z = 0;
    long nn = 2 + world->numpop2 + (world->bayes->mu)  + world->species_model_size * 2 + world->grownum; 
    const long nng= nn-1;
    bufsize = 2 * (nng); //acceptance ratio: params + tree
    bufsize += 2; // loci + numparams
    bufsize += world->bayes->numparams * nn;
    bufsize += 3 + world->options->heated_chains + 1;
    //printf("%i> bufsize in pack_single_bayes_buffer()=%li\n",myID,bufsize);
    (*buffer) = (MYREAL *) myrealloc(*buffer, bufsize * sizeof(MYREAL));
    memset (*buffer, 0, bufsize * sizeof(MYREAL));
#ifdef DEBUG_MPI
    fprintf(stdout, "%i>>>>>\n  buffersize=%li, numparams=%li\n>>>>\n", 
    	    myID,
	    bufsize,
	    world->bayes->numparams
	    );
    fflush(stdout);
#endif

    (*buffer)[z++] = (MYREAL) locus;
    (*buffer)[z++] = (MYREAL) bayes->numparams;

    for(i = 0; i < world->bayes->numparams; ++i)
      {
	//the first and second elements are logprob                                                                                   
	(*buffer)[z++] = (bayes->params+(i*nn))[0];
	(*buffer)[z++] = (bayes->params+(i*nn))[1];
	for (j = 2; j < nn; ++j) //the first and second elements are logprob                                              
	  {
	    if(bayes->map[j-2][1] != INVALID)
	      (*buffer)[z++] = (bayes->params+(i*nn))[j];
	  }
      }
    // for the parameters                                                                                                           
    for (j = 0; j < nng-1; ++j)
      {
	if(bayes->map[j][1] != INVALID)
	  {
	    (*buffer)[z++] =  (MYREAL) world->accept_archive[j];
	    (*buffer)[z++] =  (MYREAL) world->trials_archive[j];
	  }
      }
    // for the genealogy                                                                                                            
    (*buffer)[z++] =  (MYREAL) world->accept_archive[j];
    (*buffer)[z++] =  (MYREAL) world->trials_archive[j];

    if(bufsize < z)
      error("buffer is too small in pack_single_bayes_buffer()\n");
    memset(world->accept_archive,0,sizeof(long)*2*nng);//Cesky Krumlov 2013
    return z;
}
///
/// Pack bayes parameter buffer, sent from replicant nodes AND lociworker to the master
/// the acceptance and trial values will be simply added.

long pack_single_bayes_buffer_part(MYREAL **buffer, bayes_fmt *bayes, world_fmt *world,long locus)
{
    long j;
    long z = 0;
    const long nng= world->numpop2 + world->bayes->mu + 1 + world->species_model_size * 2 + world->grownum; 

    (*buffer)[z++] = (MYREAL) locus;
    (*buffer)[z++] = (MYREAL) bayes->numparams;
#ifdef DEBUG
    printf("%i> %f %f ", myID, (*buffer)[z-2],(*buffer)[z-1]);
#endif
    //
    for (j = 0; j < nng-1; ++j)
    {
      if(bayes->map[j][1] != INVALID)
	{
	  (*buffer)[z++] = (MYREAL) world->accept_archive[j]; 
	  (*buffer)[z++] = (MYREAL)  world->trials_archive[j];
#ifdef DEBUG
	  printf("%i> %f %f ",myID, (*buffer)[z-2],(*buffer)[z-1]);
#endif
	}
    }
    // genealogy
    (*buffer)[z++] = (MYREAL) world->accept_archive[j]; 
    (*buffer)[z++] = (MYREAL)  world->trials_archive[j];
#ifdef DEBUG
    printf("%i> %f %f [pack_single_...part]\n", myID, (*buffer)[z-2],(*buffer)[z-1]);
#endif
    memset(world->accept_archive,0,sizeof(long)* 2 * nng); //removes accept and trials archive
    return z;
}



///
/// gather results (sumfiles, results, migrate-histogram, ..) from workers
void
mpi_results_master (long sendtype, world_fmt * world, long maxreplicate,
                    void (*unpack) (MYREAL *buffer, world_fmt * world,
                                    long locus, long maxrep, long numpop))
{
#ifdef DEBUG_MPI
  long ii;
#endif
    long numpop = world->numpop;
    long bufsize = 1;
    // maxreplicate > 1 ---> add 1 [this all needs careful checking]
    // MIGMPI_SUMFILE -----> 0 
    // MIGMPI_HIST    -----> 0
    // still strange? long addon = (maxreplicate>1) ? 1 : ((sendtype == MIGMPI_SUMFILE) ||  (sendtype == MIGMPI_MIGHIST) )? 0 : ((world->loci == 1) ? 0 : 1);
    //long addon = (maxreplicate>1) ? 1 : ((world->loci > 1) ? 1 : 0);
    long addon = 1;
    //    boolean done = FALSE;
    MYREAL *buffer;
    MYREAL *temp;
    int worker;
    long z, tag, sender;
    MPI_Status status;
    long numelem = world->numpop2 + (world->options->gamma ? 1 : 0);
    long numelem2 = 2 * numelem;

    temp = (MYREAL *) mycalloc (numelem2 + 2, sizeof (MYREAL));
    buffer = (MYREAL *) mycalloc (bufsize+1, sizeof (MYREAL));
    temp[0] = (MYREAL)sendtype;
    temp[1] = (MYREAL) bufsize;
    for (worker = 1; worker < MIN (world->loci + addon, numcpu); worker++)
    {
      //printf("%i> MASTER requests information from n%i for locus %i\n",myID, worker, worker-1);
        MYMPISEND (temp, numelem2 + 2, mpisizeof, worker, worker, comm_world);
    }
    z = 0;
    while (z < world->loci)
    {
        MYMPIRECV (&bufsize, ONE, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG,
                  comm_world, &status);
        buffer = (MYREAL *) myrealloc (buffer, sizeof (MYREAL) * (bufsize + 1));
        memset (buffer, 0, sizeof (MYREAL) * (bufsize + 1));
        sender = status.MPI_SOURCE;
        tag = status.MPI_TAG;
#ifdef DEBUG_MPI
	fprintf(stdout, "%i> z=%li worker=%li bufsize=%li -------------------------------------\n",myID, z, sender, bufsize);
#endif
        MYMPIRECV (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) tag, comm_world, &status);
#ifdef DEBUG_MPI
	fprintf(stdout,"%i>------------------------------------\nbuffer=",myID);
	for(ii=0;ii<bufsize;ii++)
	  fprintf(stdout," %f", buffer[ii]);
	fprintf(stdout,"\n%i>-------------------------------------\n",myID);
#endif
        (*unpack) (buffer, world, tag - 1, maxreplicate, numpop);
	//fprintf(stdout,"%i> unpacked bufsize=%li from node %li\n",myID,bufsize,sender); fflush(stdout);
        z++;
    }
    myfree(buffer);
    myfree(temp);
}

void
mpi_results_worker (long bufs, world_fmt * world, long maxrep,
                    long (*pack) (MYREAL **buffer, world_fmt * world,
                                  long locus, long maxrep, long numpop))
{
    long numpop = world->numpop;
    long ww, locus;
    MYREAL *allbuffer;
    long bufsize = 1;
    allbuffer = (MYREAL *) mycalloc(1, sizeof(MYREAL));
    //fprintf(stdout,"%i> locidone=%i\n",myID, locidone); fflush(stdout);
    for (ww = 0; ww < locidone; ww++)
    {
      locus = world->who[ww];
      bufsize = (*pack) (&allbuffer, world, locus, maxrep, numpop);
#ifdef DEBUG_MPI
      fprintf(stdout,"%i> locus=%li after pack bufsize=%li\n",myID, locus, bufsize); fflush(stdout);
#endif
      MYMPISEND (&bufsize, ONE, MPI_LONG, MASTER, (int) (locus + 1), comm_world);
      //fprintf(stdout,"%i> sending results from locus=%li using bufsize=%li to master \n",myID, locus, bufsize); fflush(stdout);
      MYMPISEND (allbuffer, bufsize, mpisizeof, MASTER, (int) (locus + 1), comm_world);
    }
    myfree(allbuffer);
}

void
mpi_broadcast_results (world_fmt * world, long loci,
                       long (*pack) (MYREAL **buffer, world_fmt * world,
                                     long locus, long maxrep, long numpop),
                       void (*unpack) (MYREAL *buffer, world_fmt * world,
                                       long locus, long maxrep, long numpop))
{
    long locus;
    // long addon = (world->loci == 1) 0 : 1;
    long bufsize=1;
#ifdef DEBUG_MPI
    char nowstr[STRSIZE];
#endif
    MYREAL *allbuffer = NULL;// = &world->buffer;

    long maxreplicate = (world->options->replicate
                         && world->options->replicatenum >
                         0) ? world->options->replicatenum : 1;
    allbuffer = (MYREAL *) mycalloc (1, sizeof (MYREAL));
#ifdef DEBUG_MPI
    get_time (nowstr, "%H:%M:%S");
    if(world->options->progress)
      fprintf(stdout, "%i> Redistributing the data\nResult parts [Time is %s]\n",myID, nowstr);
#endif
    for (locus = 0; locus < loci; locus++)
    {
        if (myID == MASTER)
        {
            bufsize =(*pack) (&allbuffer, world, locus, maxreplicate,
                              world->numpop);
            MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
            MYMPIBCAST (allbuffer, bufsize, mpisizeof, MASTER, comm_world);
#ifdef DEBUG_MPI
            printf("%i> Locus %li results sent\n",myID, locus);
#endif
        }
        else
        {
            MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
            allbuffer = (MYREAL *) myrealloc (allbuffer, sizeof (MYREAL) * bufsize + 1);
            MYMPIBCAST (allbuffer, bufsize, mpisizeof, MASTER, comm_world);
            (*unpack)(allbuffer, world, locus, maxreplicate,
                      world->numpop);
#ifdef DEBUG_MPI
            printf("%i> Locus %li results received\n",myID, locus);
#endif
        }
	//myfree(allbuffer);
	//allbuffer=NULL;
	//memset (allbuffer, 0, sizeof (char) * bufsize);
    }
    MYMPIBARRIER(comm_world);
    myfree(allbuffer);
}


/*
// send the data over all loci/replicates to all nodes
// including the master node, so that all nodes can then 
// start calculating profiles [see calc_profiles()]
//
void distribute_locidata(world_fmt *world)
{
  char *buffer;
  pack_loci_data(world, &buffer);
  MPI_allgather(buffer);
  unpack_loci_data(buffer, world);
  myfree(buffer);
}
 
void pack_loci_data(world_fmt *world, char **buffer)
{
  long replicates = world->options->repl
  *buffer = myrealloc(*buffer,LONGLINESIZE);
  hits = sscanf (input, "%li %li %li %li %li", &world->loci, &world->numpop, &world->numpop2, &tmp, &replicates);  
}
*/

// necessary for analyzing old sumfiles using MPI
//
// master is reusing  mpi_runloci_master()
void
assignloci_worker (world_fmt * world, option_fmt *options, long * Gmax)
{
    boolean done = FALSE;
    long locus;
    MPI_Status status;
    long * twolongs;
    char *locusstring;
    twolongs = (long *) mycalloc(TWO,sizeof(long));
    locusstring = (char *) mycalloc(SMALLBUFSIZE,sizeof(char));

    world->options->progress=FALSE;
    options->progress=FALSE;
    
    while (!done)
    {
        MYMPIRECV (twolongs, TWO, MPI_LONG, (MYINT) MASTER, (MYINT) MPI_ANY_TAG,
		   comm_world, &status); //from mpi_runloci_master() around line migrate_mpi.c:163
        if (status.MPI_TAG != 0) //stop condition
        {
	  locus = twolongs[0];
	  sprintf(locusstring,"R%li",locus);
#ifdef DEBUG_MPI
	  printf("%i>>>>>> received locus %li in assignloci_worker{}\n",myID,locus);
	  swap_atl (locus, locidone, world);
	  printf("%i>>>>>> will send locus %li (%s) in assignloci_worker{}\n",myID,locus,locusstring);
#else
	  swap_atl (locus, locidone, world);
#endif
	  MYMPISEND (locusstring, SMALLBUFSIZE, MPI_CHAR, (MYINT) MASTER, (MYINT) locus + 1, comm_world);
	  /* we want to know what locus we worked for
	     - to control the work sent by master
	     - to use in setup_parameter0() [combroyden2.c] */
	  world->who[locidone++] = locus;
	  //
	  if (options->replicate && options->replicatenum > 0)
            {
	      world->locus = locus;
	      world->repkind = MULTIPLERUN;
#ifdef LONGSUM
	      change_longsum_times (world);   //multi run replicates
#endif         /*LONGSUM*/
            }

	  //
        }
        else
	  {
            done = TRUE;
	    assign_worker_cleanup();
#ifdef DEBUG_MPI
	    fprintf(stdout,"%i> STOP: received stop from %i\n",myID, status.MPI_SOURCE);
#endif
	  }
    }
    myfree(locusstring);
}

void
swap_atl (long from, long to, world_fmt * world)
{
    long r;
    timearchive_fmt *tmp;
    for (r = 0; r < world->options->replicatenum; r++)
    {
        tmp = &world->atl[r][to];
        world->atl[r][to] = world->atl[r][from];
        world->atl[r][from] = *tmp;
    }
}


//#ifdef SLOWNET

///
/// checks whether a node is already in the mpistack list
boolean in_mpistack(int sender, world_fmt *world)
{
  long i;
  for (i=0;i < world->mpistacknum; i++)
    {
      if((world->mpistack[i]) == sender)
	return TRUE;
    }
  return FALSE;
}

void handle_dataondemand(int sender,int tag,char *tempstr, world_fmt *world, option_fmt * options, data_fmt *data)
{
  //long  pos=0;
  long pop;
  long ind;
  long sublocus;
  long allelenum;
  char *buffer;
  long allocbufsize=LONGLINESIZE;
  long bufsize;
  //sleep(100);
  buffer = (char *) mycalloc (allocbufsize,sizeof (char));
  sscanf(tempstr+1,"%li%li%li%li", &pop, &ind, &sublocus, &allelenum);
  bufsize = sprintf (buffer, "%*s\n",  (int) options->nmlength, data->indnames[pop][ind][0]);
  if(allelenum!=0)
    {
      if (bufsize > allocbufsize-100)
	{
	  allocbufsize += bufsize;
	  buffer = (char*) myrealloc(buffer, allocbufsize * sizeof(char));
	}
      bufsize += sprintf (buffer + bufsize, "%s %s\n", data->yy[pop][ind][sublocus][0][0],data->yy[pop][ind][sublocus][1][0]);
    }
  else
    {
	mutationmodel_fmt *s = &data->mutationmodels[sublocus];
	int tmp = (s->dataclass == SITECHARACTER ? 0 : 1);
	bufsize += sprintf (buffer + bufsize, "%i %f %f %f %f\n", tmp, s->basefreqs[0],s->basefreqs[1],s->basefreqs[2],s->basefreqs[3]);
	long site;
	for(site=0; site < s->numsites; site++)
	  {
	    if (bufsize > allocbufsize-s->numsites)
	      {
		allocbufsize += s->numsites+2;
		buffer = (char*) myrealloc(buffer, allocbufsize * sizeof(char));
	      }
	    bufsize += sprintf (buffer + bufsize, "%s", data->yy[pop][ind][sublocus][0][site]);
	    //printf("%s",data->yy[pop][ind][sublocus][0][site]);
	  }
	//printf("###%i####\n",myID);
    }
  MYMPISEND (&bufsize, 1, MPI_LONG, (MYINT) sender, (MYINT) tag, comm_world);
  MYMPISEND (buffer, bufsize, MPI_CHAR, (MYINT) sender, (MYINT) tag, comm_world);
  myfree(buffer);
}


void handle_replication(int sender,int tag,char *tempstr, world_fmt *world)
{
  //long  pos=0;
  long i;
  long * temp;
  int realsender;
  long locus;
  long replicate;
  int replicator;
  //  boolean from_locus_sender = (sender <= world->loci);
  //  boolean from_replicator = (sender > world-> loci);
  sscanf(tempstr+1,"%i%li%li", &realsender, &locus, &replicate);
  //if(from_replicator)
  //  warning("these guys should not send to here");
  temp = (long *) mycalloc(3, sizeof(long));
  temp[0] = sender;
  temp[1] = locus;
  temp[2] = replicate;
  if(world->mpistacknum > 0)
    {
      world->mpistacknum -= 1;
      replicator = world->mpistack[world->mpistacknum];
      //printf("%i> checking out replicator:[%4li] %i\n",myID, world->mpistacknum,replicator);
      MYMPISEND(temp, 3, MPI_LONG, replicator, (MYINT) tag, comm_world);
    }
  else
    {
      if(world->mpistack_requestnum>=world->mpistack_request_numalloc)
	{

	  world->mpistack_request = (mpirequest_fmt *) myrealloc(world->mpistack_request,
					      sizeof(mpirequest_fmt)*(world->mpistack_requestnum+10));
	  for(i=world->mpistack_request_numalloc; i < world->mpistack_requestnum + 10; i++)
	    {
	      	world->mpistack_request[i].tempstr = (char *) mycalloc(SMALLBUFSIZE,sizeof(char));
	    }
	  world->mpistack_request_numalloc = world->mpistack_requestnum + 10;
	}
      strcpy(world->mpistack_request[world->mpistack_requestnum].tempstr, tempstr);
      world->mpistack_request[world->mpistack_requestnum].sender = sender;
      world->mpistack_request[world->mpistack_requestnum].tag = tag;
      world->mpistack_requestnum += 1;
#ifdef DEBUG
      printf("%i> MASTER: added request from n%i to request-stack position %li\n",myID, sender, world->mpistack_requestnum);
      fflush(stdout);
#endif
    }
  myfree(temp);
}

void handle_mdim(float *values,long n, int sender, world_fmt * world)
{
    register long i;
    register long z;

    int digits;
#ifdef ZNZ
    znzFile file = world->bayesmdimfile;
#else
    FILE *file = world->bayesmdimfile;
#endif
    register char *temp;
    long addition=0;
    long skypnum = world->timeelements*(world->numpop2+addition);

    temp = (char *) mycalloc(n*20+20*skypnum+LINESIZE,sizeof(char)); 
    z = sprintf(temp,"%li\t%li\t%li\t%f\t%f\t%f\t%f\t%li\t%f",
	   (long) values[0], (long) values[1]+1,(long) values[2]+1,
	    values[3], values[4], values[5], values[6],(long) values[7], values[8]);
    for(i=9;i<n; i++)
      {
	//fprintf (stdout, "%f (%li) ",values[i], i);
	if(i > (world->numpop2 +  world->bayes->mu))
	  {
	    z += sprintf(temp+z,"\t%f",values[i]);
	    continue;
	  }
	if(fabs(values[i]) < SMALLEPSILON)
	  {
	    z += sprintf(temp+z,"\t0");
	  }
	else
	  {
	    digits = (long) log10(values[i]);
	    switch(digits)
	      {
	      case -8:
	      case -6:
	      case -5:
		//		fmt = 10;
		z += sprintf(temp + z,"\t%.10f", values[i]);
		break;
	      case -4:
	      case -3:
		//fmt = 8;
		z += sprintf(temp + z,"\t%.8f", values[i]);
		break;
	      case -2:
	      case -1:
		//		fmt= 5;
		z += sprintf(temp + z,"\t%.5f", values[i]);
		break;
	      case 0:
	      case 1:
		//		fmt = 4;
		z += sprintf(temp + z,"\t%.4f", values[i]);
		break;
	      case 2:
		//		fmt = 2;
		z += sprintf(temp + z,"\t%.2f", values[i]);
		break;
	      case 3:
		//		fmt = 1;
		z += sprintf(temp + z,"\t%.1f", values[i]);
		break;
	      case 4:
	      case 5:
	      case 6:
	      case 7:
	      case 8:
		//		fmt = 0;
		z += sprintf(temp + z,"\t%.0f", values[i]);
		break;
	      default:
		if(digits<-8)
		  {
		    //		    fmt=20;
		    z += sprintf(temp + z,"\t%.20f", values[i]);
		  }
		else
		  {		   
		    //		    fmt = 5;
		    z += sprintf(temp + z,"\t%.5f", values[i]);
		  }
		break;
	      }
#ifdef DEBUG	    
	    fprintf(stdout,"\t%f", values[i]);
#endif
	  }
      }
#ifdef DEBUG
    fprintf (stdout, " [%i @@@@@@@@@@@@@writing to bayesallfile] \n", myID);
#endif
#ifdef ZNZ
    znzprintf(file,"%s\n",temp);
#else
    fprintf(file,"%s\n",temp);
    fflush(file);
#endif
    myfree(temp);
    //fprintf(stdout,"\n");
    // calculate_parallel_convergence(world, values, size);
}

void handle_message(char *rawmessage,int sender, world_fmt * world)
{
    char *rawptr;
    long  pos=0;
    void *file = (void *) stdout;
    rawptr = rawmessage;
    //fprintf(stderr,"%i> handle_message: %s\n",myID, rawmessage);
    set_filehandle(rawmessage, world, &file, &pos);
    //fprintf(stderr,"%i> handle_message: pos=%li\n",myID,pos);
    //fprintf(stderr,"%i> handle_message: %s\n",myID, rawmessage + pos);
    fprintf((FILE *) file,"%s", rawptr + pos);
    fflush((FILE *) file);
}

void handle_burnin_message(char *rawmessage,int sender, world_fmt * world)
{
  static long       z = 0;
  static boolean done = FALSE;
  static long    *replicates;
  long locus;
  long step;
  MYREAL var;
  MYREAL oldvar;
  MYREAL ess;
  MYREAL acceptance;
  if(!done)
    {
      replicates = mycalloc(world->loci,sizeof(long));
      done = TRUE;
    }
#ifdef USE_MYREAL_FLOAT
  sscanf(rawmessage,"%li%f%f%f%f%li",&locus, &ess, &acceptance, &var, &oldvar, &step);
#else
  sscanf(rawmessage,"%li%lf%lf%lf%lf%li",&locus, &ess, &acceptance, &var, &oldvar, &step);
#endif
  replicates[locus] += 1;
  if(world->options->verbose)
    {
      fprintf(stdout,"%i> Burn-in on node %i (locus=%li, repl=%li) stopped at step %li with variance-ratio=%.2f/%.2f=%.3f\n       min(ess)=%f and avg(acceptance)=%f", myID, sender, locus, replicates[locus], step, var, oldvar, var/oldvar, ess, acceptance);
    }
  if(world->burnin_stops_alloc <= z)
    {
      world->burnin_stops_alloc += z;
      world->burnin_stops = (burnin_record_fmt *) myrealloc(world->burnin_stops, world->burnin_stops_alloc * sizeof(burnin_record_fmt));  
    } 
  world->burnin_stops[z].locus        =  locus; 
  world->burnin_stops[z].replicate =  replicates[locus]; 
  world->burnin_stops[z].stopstep = step; 
  world->burnin_stops[z].ess       = ess;
  world->burnin_stops[z].accept  = acceptance;
  world->burnin_stops[z].variance = var; 
  world->burnin_stops[z].oldvariance = oldvar;
  world->burnin_stops[z].worker = sender;
  z++;
#ifdef DEBUG
  printf("%i> handle burnin_stop: z=%li maxalloc=maxrep*loci\n",myID,z);
#endif
}

///
/// sets up a file database so that the master and worker worker end up writing to the same file
/// the workers send the stuff to the master and the master then figures out (using the db)
/// what file pointer the stuff was intended for.
/// needs globals filedb, and filenum
void setup_filehandle_db(FILE *file, world_fmt *world, option_fmt *options, data_fmt *data)
{
    long filehandle = get_filehandle(file, world, options, data);
    filedb[filenum].file = file;
    filedb[filenum++].handle = filehandle;
#ifdef DEBUG_MPI
    fprintf(stdout,"filedb %li: %p %li\n",filenum, file,filehandle);
#endif
}

long retrieve_filehandle(FILE *file)
{
    long i=0;
    long filehandle = 0;
    while(filedb[i].file != file && i<filenum)
        i++;
    if(i!=filenum)
        filehandle = filedb[i].handle;
    return filehandle;
}

long get_filehandle(void *vfile, world_fmt *world, option_fmt *options, data_fmt *data)
{
#ifdef ZNZ
  if(((znzFile) vfile) == world->bayesmdimfile)
        return BAYESMDIMFILENUM;
    FILE * file = (FILE *) vfile;
#else
    FILE * file = (FILE *) vfile;
    if(file == world->bayesmdimfile)
        return BAYESMDIMFILENUM;
#endif
    if(file == stdout)
        return STDOUTNUM;
    if(file == options->logfile)
        return LOGFILENUM;
    if(file == world->outfile)
        return OUTFILENUM;
    if(file == options->aicfile)
        return AICFILENUM;
    if(file == world->mathfile)
        return MATHFILENUM;
    if(file == world->mighistfile)
        return MIGHISTFILENUM;
    if(file == world->skylinefile)
        return SKYLINEFILENUM;
    if(file == world->bayesfile)
        return BAYESFILENUM;
    if(file == world->pdfoutfile)
        return PDFOUTFILENUM;
    if(file == world->treefile)
        return TREEFILENUM;
    if(file == world->divtimefile)
        return DIVTIMEFILENUM;
    return STDOUTNUM;
}

long get_filehandle2(void *vfile, world_fmt *world)
{
  FILE *file;
#ifdef ZNZ
  if(((znzFile) vfile) == world->bayesmdimfile)
        return BAYESMDIMFILENUM;
#else
    if(((FILE *) vfile) == world->bayesmdimfile)
        return BAYESMDIMFILENUM;
#endif
    file = (FILE*) vfile;
    if(file == stdout)
        return STDOUTNUM;
    if(file == world->options->logfile)
        return LOGFILENUM;
    if(file == world->outfile)
        return OUTFILENUM;
    if(file == world->mighistfile)
        return MIGHISTFILENUM;
    if(file == world->skylinefile)
        return SKYLINEFILENUM;
    if(file == world->options->aicfile)
        return AICFILENUM;
    if(file == world->mathfile)
        return MATHFILENUM;
    if(file == world->bayesfile)
        return BAYESFILENUM;
    if(file == world->treefile)
        return TREEFILENUM;
    if(file == world->divtimefile)
        return DIVTIMEFILENUM;

    return STDOUTNUM;
}

void set_filehandle(char *message, world_fmt *world,
                    void **file, long *msgstart)
{
    char *temp;
    long filenum;
    long i = 1;
    temp = (char *) mycalloc(10,sizeof(char));
    if(message[0] == '\0')
      {
	warning("%i> set_filehandle() the message was NULL",myID);
      }
    temp[0] = message[i];
    while(temp[i-1]!=':' && i < 9 && message[i]!='\0')
      {
#ifdef DEBUG_MPI
	fprintf(stderr,"%i>>>>>> temp     =%s\n",myID,temp);
	fprintf(stderr,"%i>>>>>> temp[%li]=%c\n",myID,i,temp[i]);
	fprintf(stderr,"%i>>>>>> temp     =%s\n",myID,message);
#endif
	i++;
	temp[i-1] = message[i];
      }
    *msgstart = i+1;
    filenum = atol(temp);

    myfree(temp);
    switch(filenum)
      {
      case STDOUTNUM:
	{
	  //		fprintf(stdout,"\n");
	  *file = stdout;
	  return;
	}
      case LOGFILENUM:
	{
	  //	fprintf(stdout," logfile\n");
	  *file = world->options->logfile;
	  return;
	}
      case OUTFILENUM:
	{
	  *file = world->outfile;
	  return;
	}
      case AICFILENUM:
	{
	  *file = world->options->aicfile;
	  return;
	}
      case MATHFILENUM:
	{
	  *file = world->mathfile;
	  return;
	}
      case MIGHISTFILENUM:
	{
	  *file = world->mighistfile;
	  return;
	}
      case SKYLINEFILENUM:
	{
	  *file = world->skylinefile;
	  return;
	}
      case BAYESFILENUM:
	{
	  *file = world->bayesfile;
	  return;
	}
      case BAYESMDIMFILENUM:
	{
	  *file = world->bayesmdimfile;
	  return;
	}
      case TREEFILENUM:
	{
	  *file = world->treefile;
	  return;
	}
      case DIVTIMEFILENUM:
	{
	  *file = world->divtimefile;
	  return;
	}
      case PDFOUTFILENUM:
	{
	  *file = world->pdfoutfile;
	  return;
	}
      }
    *file = stdout;
    return;
}


void
mpi_fprintf(FILE *file, const char *fmt, ...)
{
    char *p1 = NULL;
    char *p  = NULL;
    va_list ap;
    long filehandle = 0;
    long bufsize = 0;
    
    const size_t pallocsize = LONGLINESIZE;
    p  = (char *) calloc(pallocsize,sizeof(char));
    p1 = (char *) calloc(STRSIZE,sizeof(char));
    if(myID != MASTER)
    {
        filehandle = retrieve_filehandle(file);
        bufsize = sprintf(p, "%c%li:",'M',filehandle);
    }
    va_start(ap, fmt);
    bufsize += vsprintf(p+bufsize, fmt, ap);
    if(bufsize >= pallocsize)
      error("failed in mpi_printf(): problem with buffer size!");
    if(myID != MASTER)
    {
        sprintf(p1,"M%li",bufsize);
        MYMPISEND (p1, STRSIZE, MPI_CHAR, (MYINT) MASTER, (MYINT) myID+PRINTTAG, comm_world);
        MYMPISEND (p, bufsize, MPI_CHAR, (MYINT) MASTER, (MYINT) myID+PRINTTAG, comm_world);
    }
    else
      {
        fprintf(file,"%s", p);
      }
    va_end(ap);
    myfree(p);
    myfree(p1);
}

void
mpi_fprintf2(FILE *file, long filesize, const char *fmt, ...)
{
    char *p1;
    char *p;
    va_list ap;
    long filehandle = 0;
    long bufsize = 0;
    
    long pallocsize = filesize+strlen(fmt)+10;//leave room for "M:number"
    
    p  = (char *) mycalloc(pallocsize,sizeof(char));
    p1 = (char *) mycalloc(STRSIZE,sizeof(char));
    if(myID!=MASTER)
    {
        filehandle = retrieve_filehandle(file);
        bufsize = sprintf(p, "%c%li:",'M',filehandle);
    }
    va_start(ap, fmt);
    bufsize += vsprintf(p+bufsize, fmt, ap);
    if(bufsize>=pallocsize)
      {
	warning("Failing because bufsize=%li >= allocsize=%li\n",bufsize,pallocsize);
	error("failed in mpi_printf2()");
      }
    if(myID!=MASTER)
    {
        sprintf(p1,"M%li",bufsize);
        MYMPISEND (p1, STRSIZE, MPI_CHAR, (MYINT) MASTER, (MYINT) myID+PRINTTAG, comm_world);
        MYMPISEND (p, bufsize, MPI_CHAR, (MYINT) MASTER, (MYINT) myID+PRINTTAG, comm_world);
    }
    else
        fprintf(file,"%s", p);
    va_end(ap);
    myfree(p);
    myfree(p1);
}


void
send_divtime(const char *fmt, ...)
{
    char *p1 = NULL;
    char *p  = NULL;
    va_list ap;
    long filehandle = 0;
    long bufsize = 0;
    
    long pallocsize = LINESIZE;
    p  = (char *) mycalloc(pallocsize,sizeof(char));
    p1 = (char *) mycalloc(STRSIZE,sizeof(char));
    if(myID != MASTER)
    {
        filehandle = DIVTIMEFILENUM;
        bufsize = sprintf(p, "%c%li:",'M',filehandle);
    }
    
    va_start(ap, fmt);
    bufsize += vsprintf(p+bufsize, fmt, ap);
    if(bufsize >= pallocsize)
      error("failed in mpi_printf(): problem with buffer size!");
    if(myID != MASTER)
    {
        sprintf(p1,"M%li",bufsize);
        MYMPISEND (p1, STRSIZE, MPI_CHAR, (MYINT) MASTER, (MYINT) myID+PRINTTAG, comm_world);
        MYMPISEND (p, bufsize, MPI_CHAR, (MYINT) MASTER, (MYINT) myID+PRINTTAG, comm_world);
    }
    else
      {
        fprintf(stdout,"%s", p);
      }
    va_end(ap);
    myfree(p);
    myfree(p1);
}


///
/// sends raw bayesian parameters to master, using label 'Z' to match on the master side
#ifdef PARALIO
void mpi_mdim_send(MPI_File *file, char *values, long size)
#else
void mpi_mdim_send(float *values, long size)
#endif
{
#ifdef PARALIO
  MPI_Status status;
  char error_string[LINESIZE];
  int length_of_error_string, error_class;
#else
    char *p1;
    p1 = (char *) mycalloc(SMALLBUFSIZE,sizeof(char));
    if(myID!=MASTER)
    {
#endif /*with paralio the master writes the header using this service*/
#ifdef PARALIO
      my_write_error = MPI_File_write(*file,values,size,MPI_CHAR, &status);
      if (my_write_error != MPI_SUCCESS) 
	{
	  MPI_Error_class(my_write_error, &error_class);
	  MPI_Error_string(error_class, error_string, &length_of_error_string);
	  printf("%i> %s\n", myID, error_string);
	  MPI_Error_string(my_write_error, error_string, &length_of_error_string);
	  printf("%i> %s\n", myID, error_string);
	  my_write_error = TRUE;
	}
#else
      sprintf(p1,"Z%li",size);
      MYMPISEND (p1, SMALLBUFSIZE, MPI_CHAR, (MYINT) MASTER, (MYINT) myID+PRINTTAG, comm_world);
#ifdef DEBUG
      fprintf(stdout,"%i> mdimlast=%f\n",myID,values[size-1]);
#endif
      MYMPISEND (values, size, MPI_FLOAT, (MYINT) MASTER, (MYINT) myID+PRINTTAG, comm_world);
#endif
#ifndef PARALIO /* without paralio the master is not allowed here*/
    }
    else
        error("master sends itself bayesallfile stuff");
    myfree(p1);
#endif
}


///
/// assembles the data from a replicant (receives materials from mpi_send_replicant()
/// 
void 
mpi_receive_replicate(int sender, int tag, long locus, long replicate, world_fmt * world)
{
    MYREAL *buffer;
    long  bufsize=1;    
    MPI_Status status;

    MYMPIRECV (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) tag, comm_world, &status);
#ifdef DEBUG_MPI
    fprintf(stdout,"%i> WORKER: mpi_receive_replicate received bufsize=%li from sender=%i with tag=%i\n",myID, bufsize, status.MPI_SOURCE, status.MPI_TAG);    
    sender = status.MPI_SOURCE;
    tag = status.MPI_TAG;
    fprintf(stdout,"%i> mpi_receive_replicate received bufsize=%li from sender=%i with tag=%i\n",myID, bufsize, sender, tag);    
#endif
    buffer = (MYREAL *) mycalloc (bufsize, sizeof (MYREAL));
    MYMPIRECV (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) tag, comm_world, &status);
    //fprintf(stdout,"%i> received bufsize is really %li and bufsize=%li\n",myID,(long) (long) strlen(buffer),bufsize);
    unpack_single_bayes_buffer(buffer,world->bayes,world,locus);
    if(world->options->mighist)
    {
      MYMPIRECV (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) tag, comm_world, &status);
#ifdef DEBUG_MPI
      fprintf(stdout,"%i> mpi_receive_replicate received mighistogram bufsize=%li from sender=%i with tag=%i\n",
	      myID, bufsize, status.MPI_SOURCE, status.MPI_TAG);    
#endif
      buffer = (MYREAL *) myrealloc (buffer, bufsize * sizeof (MYREAL));
      MYMPIRECV (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) tag, comm_world, &status);
      unpack_mighist_replicate_buffer(buffer, world, locus, world->numpop);
    }
    if(world->options->mighist && world->options->skyline)
    {
      MYMPIRECV (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) tag, comm_world, &status);
      buffer = (MYREAL *) myrealloc (buffer, bufsize * sizeof (MYREAL));
      MYMPIRECV (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) tag, comm_world, &status);
      unpack_skyline_buffer(buffer, world, locus, -1, world->numpop);
    }
  // send best tree if available
    //  if(world->options->treeprint == BEST && world->options->treeinmemory == TRUE)
  if(world->options->treeinmemory == TRUE)
    {
      MYMPIRECV (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) tag, comm_world, &status);
      buffer = (MYREAL *) myrealloc (buffer, bufsize * sizeof (MYREAL));
      //      printf("%i> mpi_receive_buffer() received treespace buffersize: %li %p\n", myID, bufsize,  buffer);
      MYMPIRECV (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) tag, comm_world, &status);
      unpack_treespace_buffer(buffer, world, locus, -1, world->numpop);      
    }

    // BF material
      MYMPIRECV (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) tag, comm_world, &status);
      buffer = (MYREAL *) myrealloc (buffer, bufsize * sizeof (MYREAL));
      MYMPIRECV (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) tag, comm_world, &status);
      if(!world->data->skiploci[locus])
	{
#ifdef DEBUG
	  printf("%i> mpi_receive_buffer() received BF buffersize: %li %p\n", myID, bufsize,  buffer);
#endif
	  bufsize = unpack_BF_buffer(buffer, 0, locus, world);
	}
      // ESS material
      bufsize = unpack_ess_buffer(buffer, bufsize, world);
      // hyper material
      bufsize = unpack_hyper_buffer(buffer,bufsize,world);
  if(world->data->haplotyping_report)
    {
      MYMPIRECV (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) tag, comm_world, &status);
      buffer = (MYREAL *) myrealloc (buffer, bufsize * sizeof (MYREAL));
      MYMPIRECV (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) tag, comm_world, &status);
      unpack_haplotypes_buffer(buffer, world, locus, replicate, world->numpop);
    }
  if(world->has_unassigned)
    {
      MYMPIRECV (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) tag, comm_world, &status);
      buffer = (MYREAL *) myrealloc (buffer, bufsize * sizeof (MYREAL));
      MYMPIRECV (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) tag, comm_world, &status);
      unpack_assign_buffer(buffer,world,locus,replicate, world->numpop);
    }
  if(world->has_estimateseqerror)
    {
      MYMPIRECV (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) tag, comm_world, &status);
      buffer = (MYREAL *) myrealloc (buffer, bufsize * sizeof (MYREAL));
      MYMPIRECV (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) tag, comm_world, &status);
      unpack_seqerror_buffer(buffer,world,locus,replicate, world->numpop);
    }
    myfree(buffer);
}

///
/// replicant sends data to sub-master
/// 
void 
mpi_send_replicate(int sender, long locus,  long replicate, world_fmt * world)
{
  long    allocbufsize   = 2;
  long    bufsize        = 0;
  long    numpop         = world->numpop;
  long    numpop2        = numpop * numpop;
  //long    numpop2plus    = numpop2 + 2 * numpop;
  long    npp            = numpop2 + world->bayes->mu * world->loci;
  MYREAL  *buffer        = NULL;
  //  timearchive_fmt **ta   = world->atl;
  buffer = (MYREAL *) mycalloc (allocbufsize, sizeof (MYREAL));

  bufsize = pack_single_bayes_buffer(&buffer,world->bayes,world,locus);
  allocbufsize = bufsize;

  MYMPISEND (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) (locus+1+ REPTAG), comm_world);
  MYMPISEND (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) (locus+1 + REPTAG), comm_world);
  if(world->options->mighist)
    {
      bufsize = pack_mighist_buffer(&buffer, world, locus, -1, numpop);
      allocbufsize = bufsize;
      MYMPISEND (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) (locus+1+ REPTAG), comm_world);
      MYMPISEND (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) (locus+1 + REPTAG), comm_world);
    }

  if(world->options->mighist && world->options->skyline)
    {
      bufsize = pack_skyline_buffer(&buffer, world, locus, -1, numpop);
      allocbufsize = bufsize;
      MYMPISEND (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) (locus+1+ REPTAG), comm_world);
      MYMPISEND (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) (locus+1 + REPTAG), comm_world);
    }
  // send best tree if available
  if(/*world->options->treeprint == BEST &&*/ world->options->treeinmemory == TRUE)
    {
      //      printf("%i> send treespace buffer to %li", myID,locus+1+REPTAG);
      bufsize = pack_treespace_buffer(&buffer, world, locus, -1, numpop);
      allocbufsize = bufsize;
      MYMPISEND (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) (locus+1+ REPTAG), comm_world);
      MYMPISEND (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) (locus+1 + REPTAG), comm_world);
    }

    // BF material
  bufsize = world->options->heated_chains * world->loci + 5 * world->loci + 20 * (npp+1) + 6 * npp;
  buffer = (MYREAL *) myrealloc (buffer, bufsize *  sizeof (MYREAL));
  if(!world->data->skiploci[locus])
    {
      //fprintf(stderr,"%i> REPLICANT: packed result locus=%li replicate %li\n",myID,locus, replicate);
      bufsize = pack_BF_buffer(&buffer, 0, locus, world);
    }
  bufsize = pack_ess_buffer(&buffer, bufsize, world);
  bufsize = pack_hyper_buffer(&buffer, bufsize, world);
      // send BF and ESS material
  MYMPISEND (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) (locus+1+ REPTAG), comm_world);
  MYMPISEND (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) (locus+1 + REPTAG), comm_world);
  
  if(world->data->haplotyping_report)
    { 
      bufsize = pack_haplotypes_buffer(&buffer, world, locus, replicate /*doesnothing?*/, world->numpop);
      MYMPISEND (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) (locus+1+ REPTAG), comm_world);
      MYMPISEND (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) (locus+1 + REPTAG), comm_world);
      free_haplotypes(world,locus);
    }
  if(world->has_unassigned)
    {
      bufsize = pack_assign_buffer(&buffer, world, locus, replicate,world->numpop);
      MYMPISEND (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) (locus+1+ REPTAG), comm_world);
      MYMPISEND (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) (locus+1 + REPTAG), comm_world);
    }
  if(world->has_estimateseqerror)
    {
      bufsize = pack_seqerror_buffer(&buffer, world, locus, replicate,world->numpop);
      MYMPISEND (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) (locus+1+ REPTAG), comm_world);
      MYMPISEND (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) (locus+1 + REPTAG), comm_world);
    }
  myfree(buffer);
}

#ifdef MPICHECK
void set_memory_limit(rlim_t softsize,rlim_t maxsize)
{
  struct rlimit r;
  r.rlim_cur=softsize;
  r.rlim_max=maxsize;
  setrlimit(RLIMIT_AS, &r);
}

void check_memory_limit()
{
  struct rlimit r;
  getrlimit(RLIMIT_AS, &r);
  fprintf(stdout, "%i> current memory/data usage: %f\n",myID, (double) r.rlim_cur);
}
#endif
#endif /* MPI */
