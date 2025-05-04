/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 R A N D O M   G E N E R A T O R   R O U T I N E S 
 
 creates options structures,
 reads options from parmfile if present
 
 prints options,
 and finally helps to destroy itself.
                                                                                                               
 Peter Beerli 1996, Seattle
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2007 Peter Beerli, Tallahassee FL
 
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
 
 
$Id: random.c 2169 2013-08-24 19:02:04Z beerli $
-------------------------------------------------------*/
/* \file random.c */
#include <stdlib.h>
#include "sighandler.h"
#include "migration.h"
#include "heating.h"
#include "random.h"
#include "tools.h"
#include "priors.h"
#include "migrate_mpi.h"
#ifdef WIN32
#include <sys/timeb.h>
#else /* WIN32 */
#include <sys/time.h>
#endif /* WIN32 */


#ifdef MERSENNE_TWISTER
//#include "../SFMT-src-1.5.1/SFMT.c"
#include "../lib/SFMT-src-1.5.1/SFMT.h"
#endif

#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif


#ifdef GRANDCENTRAL
#include <dispatch/dispatch.h>
#include <dispatch/semaphore.h>
extern dispatch_semaphore_t semaphore;
#endif


#ifdef PTHREADS
extern tpool_t heating_pool;
#endif

extern sfmt_t sfmt;
extern sfmt_t *sfmtp;
extern sfmt_t ** sfmtH;

/* prototypes ----------------------------------------- */
void getseed (option_fmt * options);
void swap_ptr (long **ptr, long **newptr);
unsigned long int random_seed(void);
void set_seed(long autoseed, unsigned long * inseed);
void init_lcg(option_fmt *options);
void getseed_lcg (option_fmt *options);
void getseed_mt (option_fmt *options);
void getseed_quasi (option_fmt *options);
MYREAL failed_random_gamma(MYREAL a);
MYREAL __gamma_rand(MYREAL alpha);

MYREAL randum (void);
#ifdef PTHREADS
MYREAL randum_thread (void);
#endif
long random_integer(long low, long high);

long random_from_freqlist(double *l, long nl);


#ifdef QUASIRANDOM
#define MYINDEX 0

//char[100] generator;

#ifdef CUD
int  nextn=1, kk=1,ii=0,jj=0;
void setup_cud(int thiskk)
{
  strncpy(generator,"Quasi-random number generator: Completely Uniform Distributed numbers", 79);
  kk = thiskk;
}

MYINLINE MYREAL get_cud()
{
  MYREAL rr;
  MYREAL value;
  const MYREAL  magicnumbers[]={2.,2.,4.,15.,65.,315.,1586.,8036.,40448.,200401.,972536.,
	       4609846.,21310545.,96017492., 421606654.,1804551131. };
  //double ps[]={2.,3.,5.,7.,11.,13.,17.,19.,23.,29.,31.,37.,41.,47.,53.,59.,61.,67.,71.,73.,79.,83.,89.,97.,101.};
  const MYREAL lps[] =
  {0.69314718055994530942, 1.0986122886681096914, 
   1.6094379124341003746, 1.9459101490553133051, 
   2.3978952727983705441, 2.5649493574615367361, 
   2.8332133440562160802, 2.9444389791664404600, 
   3.1354942159291496908, 3.3672958299864740272, 
   3.4339872044851462459, 3.6109179126442244444, 
   3.7135720667043078039, 3.8501476017100585868, 
   3.9702919135521218341, 4.0775374439057194506, 
   4.1108738641733112488, 4.2046926193909660597, 
   4.2626798770413154213, 4.2904594411483911291, 
   4.3694478524670214942, 4.4188406077965979235, 
   4.4886363697321398383, 4.5747109785033828221, 
   4.6151205168412594509}; 
  //rr=kk*log(ps[jj]);
#ifdef PTHREADS

  if (heating_pool!=NULL && (pthread_mutex_lock (&(heating_pool->random_lock))) != 0)
        error ("pthread_mutex_lock failed in get_cud()");
#endif

  rr = kk * lps[jj];
  
  if (jj >= ii)
    { 
      kk++; 
      jj=0;
    }
  else 
    jj++;

  if(kk >= magicnumbers[ii])
    {
      kk = 1; 
      ii++;
    }
  value = rr - floor(rr);
#ifdef PTHREADS

    if (heating_pool!=NULL && (pthread_mutex_unlock (&(heating_pool->random_lock))) != 0)
        error ("pthread_mutex_unlock failed in get_cud()");
#endif

    return   value;
  //printf("%f\n",rr);
}

MYINLINE MYREAL get_quasi()
{
    return get_cud();
}

#endif

#ifdef KOROBOV

unsigned int x0=137,ii=0; nextn=1;

void setup_korobov()
{
  strncpy(generator,"Quasi-random number generator: Korobov sequence", 79);
}

MYREAL get_korobov()
{  
unsigned long  a =17364, m=65521;
  unsigned long xn;
  double x;

  xn=(a*x0)%m;
   x=(xn+0.0)/(m+0.0);
   x0=xn; nextn++;
   if (nextn>m){ ii++; nextn=2+ii;}
 return x;
}


MYINLINE MYREAL get_quasi()
{
  return get_korobov();
}

#endif

#ifdef WEYL
static int s, nextn=8;

void setup_weyl()
{
  strncpy(generator,"Quasi-random number generator: Weyl sequence", 79);
}


MYINLINE MYREAL get_weyl()
{
    double r=sqrt(7.0),rr;
    int next;

    next=nextn*nextn;
    rr=r*next;

//    rr= (r*nextn-floor(r*nextn))*nextn;
    rr=rr-floor(rr);

    nextn++;
    return rr;
}

MYINLINE MYREAL get_quasi()
{
  return get_weyl();
}
#endif /*WEYL*/

#ifdef HALTON

#define MAX_D 500 
static int s, nextn=8;

double e, quasi[100];
static int prime[MAX_D];
static double iprime[MAX_D];
static int  primroots[][10]={{1, 2, 3, 3, 8, 11,12,14, 7,18},
    {12,13,17,18,29,14,18,43,41,44},
    {40,30,47,65,71,28,40,60,79,89},
    {56,50,52,61,108,56,66,63,60,66},
    {104,76,111,142,71,154,118,84,127,142},
    {84,105,186,178,188,152,165,159,103,205}, 
    {166,173,188,181,91,233,210,217,153,212},
};

static int warnockOpt[]={1,  2,  2,  5,  3,  7,  3,  10,  18, 11, 
    17, 5, 17,  26, 40, 14, 40, 44, 12, 31,
    45, 70,8,   38, 82, 8,  12, 38, 47, 70,
    29, 57, 97, 110,32, 48, 84, 124,155,26,
    69, 83, 157,171, 8, 22, 112,205, 15, 31,
    61, 105,127,212,12, 57, 109,133,179,210,
    231,34, 161,199,222,255,59, 120,218,237,
    278,341,54, 110,176,218,280,369,17, 97, 
    193,221,331,350,419,21, 85, 173,221,243,
    288,424,45, 78, 173,213,288,426,455,138,
}; 

int gohalt(double *,double *, double *);
double get_halton();
int primes();
int power(int, int, int);
int inhalt(int dimen, int atmost, double tiny, double *quasi);
int power(int a, int b, int m)
{ int i,c=1;
    for(i=0;i<b;i++)
        c=(c*a)%m;
    return c;
} 

void setup_halton()
{
  strncpy(generator,"Quasi-random number generator: Halton sequence", 79);
}

double get_halton()
{ double wq[100],dq[100];
    gohalt(quasi,dq,wq);
    return wq[MYINDEX];
} 

int inhalt(int dimen, int atmost, double tiny, double *quasi)
{
    double delta, f;
    int i,m;
    
    // check dimen
    primes();
    
    s=dimen;
    if (s<1||s>1000)
        return(-1);
    
    // compute and check tolerance
    
    e=0.9*(1.0/(atmost*prime[s-1])-10.0*tiny);
    delta=100*tiny*(double)(atmost+1)*log10((double)atmost);
    if (delta>=0.09*(e-10.0*tiny))
        return(-2);
    
    // now compute first vector
    
    m=1;
    for (i=0;i<s;i++)
    {       
        iprime[i]=1.0/iprime[i];
        quasi[i]=iprime[i]; 
        m=i*prime[i];
    }
    
    printf("largest prime=%d, %f \n",prime[s-1], quasi[1]);
    
    nextn=2;
    
    return 0;
}

int gohalt(double *quasi, double *dq, double *wq)
{
    int i, j, k, ytemp[40],xtemp[40], ztemp, ktemp, ltemp, mtemp;
    double r;
    double t,f,g,h;
    
    // generate quasi one compoment at a time using radix prime[k] for 
    // component k
    
    
    for (i=0;i<s;i++)
    {
        t=iprime[i];
        f=1.0-quasi[i];
        g=1.0;
        h=t;
        while ((f-h)<e)
            // this checks whether q+h>1-e
        {
            g=h;
            h*=t;
        }
        quasi[i]=g+h-f;
    }
    
    for(i=0;i<s;i++)
    {	  
        k=0; mtemp=nextn; 
        ltemp=prime[i]; 
        
        while(mtemp!=0){
            ytemp[k]=mtemp%ltemp;
            mtemp=mtemp/ltemp;
            k++; 
        }
        
        //generating Optimal primitive root 
        for(j=0;j<k;j++)
        {
            // xtemp[j] = (ytemp[j]*power(primroots[i/10][i%10], nextn%ltemp, ltemp))%ltemp;
            if(j>=1) 
                xtemp[j] =(warnockOpt[i]*power(primroots[i/10][i%10], ytemp[j], 
                                               ltemp)+ytemp[j-1])%ltemp;
            else xtemp[j] =(warnockOpt[i]*power(primroots[i/10][i%10], ytemp[j],
                                                ltemp))%ltemp;
            xtemp[j] -= ytemp[j];
        }  
        
        dq[i]=0;t=iprime[i]; 
        for(j=0;j<k; j++)
        { 
            dq[i] += xtemp[j]*t;
            t *= iprime[i];
        }
        
        dq[i] += quasi[i];
        
        
        // generating Warnock Optimal sequences
        for(j=0;j<k;j++)	   
        {      
            if(j>=1)
                xtemp[j]= (ytemp[j]*power(warnockOpt[i],i+1,ltemp)+ytemp[j-1])%ltemp;
            else
                xtemp[j]= (ytemp[j]*power(warnockOpt[i],i+1,ltemp))%ltemp;
            
            xtemp[j] -= ytemp[j];
        }
        
        wq[i]=0;t=iprime[i];
        for(j=0;j<k; j++)
        {
            wq[i] += xtemp[j]*t;
            t *= iprime[i];
        }
        
        wq[i] += quasi[i];
    }
    
    nextn++;
    return(0);
}

int primes()
{
    int i, j, a[MAX_D+1];
    for (a[1] = 0, i = 2; i <= MAX_D; i++)
        a[i] = 1;
    for (i = 2; i <= MAX_D/2; i++)
        for (j = 2; j <= MAX_D/i; j++)
            a[i*j] = 0;
    for (i = 1, j = 0; i <= MAX_D; i++){
        if (a[i]){
            prime[j] =i; 
			iprime[j]=i;
            j++;
        }
    }   
    return j;
}       

MYINLINE MYREAL get_quasi()
{
  return get_halton();
}

#endif /*HALTON */


#endif

/*
#ifdef MERSENNE_TWISTER
void setup_mersennetwister()
{
  strncpy(generator,"Pseudo-random number generator: Mersenne twister", 79);
}
#endif
*/

/// random integer
//=============================================
MYINLINE long random_integer(long low, long high)
{
    //Math.floor(Math.random()*(N-M+1))%(N-M+1)+M
    long r;
    long tt = high-low+1;
    r =  low + (long) (RANDUM() * tt);
    if(r<low || r>high)
      {
        // warning("Random %li integer is out of bounds (%li, %li)\n",r, low, high);
	if(r>high)
	  return high;
	else
	  return low;
	//	error("Stop because this should never happen");
      }
    return r;
}


///
/// better random_seed function taken from a discussion on 
/// http://sourceware.org/ml/gsl-discuss/2004-q1/msg00071.html
/// March 30, 2007, Robert G. Brown (http://www.phy.duke.edu/~rgb/)
/// at Duke University Dept. of Physics, Box 90305  Durham, N.C. 27708-030
/// suggested the code
/// changed December 2012 because on som cluster machines little entropy
/// is generated thus the machine stalls (AMD?) 
unsigned long int random_seed()
{
  
  unsigned int myseed;
  struct timeval tv;
  FILE *devrandom;
  //size_t retval;
  if ((devrandom = fopen("/dev/urandom","r")) == NULL) 
    {
      gettimeofday(&tv,0);
      myseed = (unsigned int) (tv.tv_sec + tv.tv_usec);
    }
  else 
    {
      fread(&myseed,sizeof(myseed),1,devrandom);
      fclose(devrandom);
      if(myseed == 0)
	{
	  gettimeofday(&tv,0);
	  myseed = (unsigned int) (tv.tv_sec + tv.tv_usec);
	}
    }
  return(myseed);
}

void set_seed(long autoseed, unsigned long * inseed)
{
  unsigned long timeseed = 0;
   switch (autoseed)
    {
    case AUTO:
      timeseed = random_seed();
#ifdef LCG
      switch (timeseed % 4)
        {
        case 0:
	  ++timeseed; break;
        case 2:
            ++timeseed;
            break;
        case 1:
        case 3:
            break;
        }
#endif
      *inseed = (unsigned long) labs((long)timeseed);
      break;
    case NOAUTO:
    case NOAUTOSELF:
        break;
    default:
        error ("Error: Seed value not defined");
        //break;
    }
#ifdef GRANDCENTRAL
   semaphore = dispatch_semaphore_create(1);
#endif
}

void init_lcg(option_fmt *options)
{
  long i;
  for (i = 0; i <= 2; i++)
    seed[i] = 0;
  i = 0;
  do
    {
      seed[i] = options->inseed & 2047;
      printf("RANDOM SEEDING: i=%li seed[i]=%li inseed=%li\n",i,seed[i], options->inseed);
      options->inseed /= 2048;
      i++;
    }
  while (options->inseed != 0);
  if(i>3)
    error("lcg random number generator is broken! random.c:507");
}

void getseed_lcg (option_fmt *options)
{
    strncpy(generator,"Pseudo-random number generator: Least Congruental Generator", 79);
    set_seed(options->autoseed, &options->inseed);
    options->saveseed = options->inseed;
    init_lcg(options);
    if(options->randomsubset>0 && options->randomsubsetseed > 0)
      {
#ifdef WIN32
	srand((unsigned int) options->randomsubsetseed);
#else
	srand48((unsigned int) options->randomsubsetseed);
#endif
      }
}

void getseed_mt (option_fmt *options)
{
  long i;
  strncpy(generator,"Pseudo-random number generator: Mersenne-Twister", 79);
  set_seed(options->autoseed, &options->inseed);
  options->saveseed = options->inseed;
  sfmt_init_gen_rand(&sfmt, (unsigned int) options->inseed);
  sfmtp = &sfmt;
  // multiple heated chains will use mutltiple mersenne twister streams
  // remedy for weird crashes in GCD on macs after long runs
  // the splitting of the stream should to my knowledge not lead to
  // problems but will speed up the GCD random number picking because
  // I do not need to block anymore. This change leads to several other
  // changes because the variable sfmt changed to a pointer sfmtp.
  if (options->heated_chains>1)
    {
      sfmtH = (sfmt_t **) mymalloc((size_t) options->heated_chains * sizeof(sfmt_t*));
      for (i=0;i<options->heated_chains;i++)
	{
	  sfmtH[i] = (sfmt_t*) malloc(sizeof(sfmt_t));
	  long mynewseed = RANDINT(0,MAXLONG-1);
	  sfmt_init_gen_rand(sfmtH[i], (unsigned int) mynewseed);
	}
    }
  if(options->randomsubset>0 && options->randomsubsetseed > 0)
    {
#ifdef WIN32
      srand((unsigned int) options->randomsubsetseed);
#else
      srand48((unsigned int) options->randomsubsetseed);
#endif
    }
}
///
/// set up the quasi random number material
/// this function is dependent on running the standard MT or LCG random number generator first.
void getseed_quasi (option_fmt *options)
{
  (void) options;
#ifdef CUD
  setup_cud((int) RANDINT (1, MAXLONG-1));
#endif
#ifdef WEYL
  setup_weyl();
#endif
#ifdef KOROBOV
  setup_korobov();
#endif
#ifdef HALTON
  setup_halton();
#endif

}


#ifdef MPI
/// Random number seed function for MPI usage, the master generates numcpu random numbers
/// and distributes them to the workers.
void getseed(option_fmt *options)
{
  //long test = -1;
  long i;
  long worker;
  long *newseeds;
  long newseed;
  boolean found=FALSE;
  MPI_Status status;
  if (myID==MASTER)
    {
#ifdef MERSENNE_TWISTER
      getseed_mt(options);
#else
      getseed_lcg(options);
#endif
      if(options->progress)
	printf("Master random number seed: %li\n", options->inseed);
      newseeds = (long*) calloc(numcpu,sizeof(long));
      for (worker=1;worker<numcpu;worker++)
	{
	  newseed = RANDINT(0,MAXLONG-1);
	  for (i=1;i<worker;i++)
	    {
	      if (newseed == newseeds[i])
		{
		  found=TRUE;
		  break;
		}
	    }
	  if (found)
	    worker--;
	  else
	    newseeds[worker] = newseed;
	}
      for (worker=1;worker<numcpu;worker++)
	{
	  newseed =  newseeds[worker];
	  MYMPISEND(&newseed,(MYINT) ONE, MPI_LONG, (MYINT) worker, 0, comm_world);
	}
    }
  else /*worker*/
    {
      MYMPIRECV (&options->inseed, (MYINT) ONE, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, comm_world, &status); 
#ifdef MERSENNE_TWISTER
      getseed_mt(options);
#else
      getseed_lcg(options);
#endif
#ifdef QUASIRANDOM
      getseed_quasi(options);
#endif
#ifdef DEBUG
      char nowstr[LINESIZE];
      get_time (nowstr, "%c");
      printf("%i> RANDOM NUMBER SEED:%li at time %s\n",myID,options->inseed,nowstr);
#endif
    }
}
#if 0 /*does not execute this section old getseed function*/
/// Random number seed function for MPI usage, the main seed gets distributed to 
/// all nodes and then the nodes pick a random number from the master seed for
/// their own seed. given that the nodes typically analyze different loci
/// it will be unlikely that even when picked by the off chance the same seed on two
/// nodes we get the same stream of events. the node ID gets the randomnumber[#ID] from the stream.
void getseed(option_fmt *options)
{
  long test = -1;
  long i;
  if (myID==MASTER)
    {
#ifdef MERSENNE_TWISTER
      getseed_mt(options);
#else
      getseed_lcg(options);
#endif
    }
#ifdef DEBUG
  char nowstr[LINESIZE];
  get_time (nowstr, "%c");
  printf("%i> random number seed broadcast start %s\n",myID, nowstr);
#endif
  MYMPIBCAST (&options->inseed, 1, MPI_LONG, MASTER, comm_world);
#ifdef DEBUG
  get_time (nowstr, "%c");
  printf("%i> random number seed broadcast finished %s\n",myID, nowstr);
#endif
  if(myID != MASTER)
    {
      // sets all worker to same seed as master
#ifdef MERSENNE_TWISTER
      getseed_mt(options);
#else
      getseed_lcg(options);
#endif
#ifdef DEBUG_MPI
      printf("WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING\n");
      printf("WARNING all random number seeds on the nodes are identical to the master     WARNING\n");
      printf("WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING\n");
      return;
#endif
      // gets a random value using the master seed
      test = RANDINT (1, MAXLONG-1); 
      options->inseed = test;
      // reset the master seed to the own seed.
#ifdef MERSENNE_TWISTER
      getseed_mt(options);
#else
      getseed_lcd(options);
#endif
#ifdef QUASIRANDOM
      getseed_quasi(options);
#endif
#ifdef DEBUG
      char nowstr[LINESIZE];
      get_time (nowstr, "%c");
      printf("%i> Random number seed: %li at time %s\n",myID,test,nowstr);
#endif
    }
  else
    {
      printf("%i> Master random number seed: %li\n",myID,options->inseed);
    }
}
#endif /*0*/
#else
void getseed(option_fmt *options)
{
#ifdef MERSENNE_TWISTER
  getseed_mt(options);
#else
  getseed_lcg(options);
#endif
#ifdef QUASIRANDOM
  getseed_quasi(options);
#endif
}
#endif


void
swap_ptr (long **ptr, long **newptr)
{
    long *temp;
    temp = *ptr;
    *ptr = *newptr;
    *newptr = temp;
}

#ifdef PTHREADS
MYREAL
randum_thread (void)
/* thread save random number generator */
{
    MYREAL value;
    if (heating_pool!=NULL && (pthread_mutex_lock (&(heating_pool->random_lock))) != 0)
        error ("pthread_mutex_lock failed in random_thread()");
#ifdef MERSENNE_TWISTER
    value = sfmt_genrand_res53(sfmtp);
#else
    newseed[0] = 1549 * seed[0];
    newseed[1] = newseed[0] / 2048;
    newseed[0] &= 2047;
    newseed[2] = newseed[1] / 2048;
    newseed[1] &= 2047;
    newseed[1] += 1549 * seed[1] + 812 * seed[0];
    newseed[2] += newseed[1] / 2048;
    newseed[1] &= 2047;
    newseed[2] += 1549 * seed[2] + 812 * seed[1];
    swap_ptr (&newseed, &seed);
    seed[2] &= 1023;
    value = (((seed[0] / 2048.0 + seed[1]) / 2048.0 + seed[2]) / 1024.0);
#endif
    if (heating_pool!=NULL && (pthread_mutex_unlock (&(heating_pool->random_lock))) != 0)
        error ("pthread_mutex_unlock failed in randum_thread()");
    return value;
}    /* threaded randum */
#else
MYREAL
randum (void)
/* Non thread-safe random number generator (faster) */
{
#ifdef MERSENNE_TWISTER
#ifdef MESS
    MYREAL val = sfmt_genrand_res53(sfmtp);
    fprintf(stdout,"%i> R:%f\n",myID,val);
    return val;
#else
#ifdef GRANDCENTRAL 
    // Wait for a free file descriptor
    dispatch_semaphore_wait(semaphore, DISPATCH_TIME_FOREVER);
    MYREAL r = sfmt_genrand_res53(sfmtp);
    dispatch_semaphore_signal(semaphore);
    return r;
#else
    MYREAL r = sfmt_genrand_res53(sfmtp);
    return r;
#endif
#endif
#else
MYREAL value;
newseed[0] = 1549 * seed[0];
newseed[1] = newseed[0] / 2048;
newseed[0] &= 2047;
newseed[2] = newseed[1] / 2048;
newseed[1] &= 2047;
newseed[1] += 1549 * seed[1] + 812 * seed[0];
newseed[2] += newseed[1] / 2048;
newseed[1] &= 2047;
newseed[2] += 1549 * seed[2] + 812 * seed[1];
swap_ptr (&newseed, &seed);
seed[2] &= 1023;
value = (((seed[0] / 2048.0 + seed[1]) / 2048.0 + seed[2]) / 1024.0);
return value;
#endif /*mersenne_twister*/
}    /* randum */
#endif /*phtreads*/


///////////// this seems to fail, but unclear why ///////////////////////
/// replaced with Marsaglia and Tsang 2000, that seems to work         //
// Hisashi Tanizaki 2008 Economics bulletin vol 3 no 7 pp 1-10
// A simple gamma random number generatior for arbitrary shape parameters
//(i) Given alpha, set n, b1, b2, c1 and c2 as follows:
//  if 0 < alpha <= 0.4   n = alpha^−1 
//  if 0.4 < alpha <= 4.0 n = alpha^-1 + alpha−1(alpha − 0.4)/3.6
//  if 4.0 < alpha        n = alpha^(1/2)
//  b1 = alpha−1/n 
//  b2 = alpha+1/n
//  if 0.0<alpha≤0.4      c1 = 0
//  if 0.4 < alpha        c1 = b1(log b1 − 1)/2
//			  c2 =b2(logb2 −1)/2.
//(ii) Generate v1 and v2 independently from U(0, 1). Set 
//  w1 = c1 + log v1, 
//  w2 = c2 + log v2,
//  y  = n(b1 w2 −b2 w1).
//(iii) Goto(ii)if y<0.
//(iv) Set x = n(w −w ), and take e^x as a gamma random draw 
//     with shape parameter  alpha if 
//     log y >= x and go to (ii) otherwise.
//
// to get gamma(alpha,beta) distributed values use random_gamma(alpha)*beta
//
MYREAL failed_random_gamma(MYREAL a) 
{
  const MYREAL inva = 1./ a;
  const MYREAL n = ((4.0 < a) ? sqrt(a) : ((0.4 < a) ? inva + inva *(a-0.4)/3.6 : inva)); 
  const MYREAL invn = 1. / n;
  const MYREAL b1 = a - inva;
  const MYREAL b2 = a + invn;
  const MYREAL c1 = ((0.4 < a) ? b1 * ( log(b1) - 1.)/2. : 0); 
  const MYREAL c2 = b2 * ( log(b2) - 1.0)/2.;
  
  MYREAL v1 = RANDUM();
  MYREAL v2 = RANDUM();
  MYREAL w1 = c1 + log(v1);
  MYREAL w2 = c2 + log(v2);
  MYREAL y  = n * (b1 * w2 - b2 * w1);
  MYREAL x = (double) -HUGE;
  while(log(y)>=x)
    {  
      while (y<0.0)
	{
	  v1 = RANDUM();
	  v2 = RANDUM();
	  w1 = c1 + log(v1);
	  w2 = c2 + log(v2);
	  y  = n * (b1 * w2 - b2 * w1);
	}
      x = n * (w2 - w1);
      //printf("random gamma: x=%f\n",x);
    }
  return exp(x);
}

///
// Gamma distribution
//
// gamma deviated random number 
// from Marsaglia G., Tsang, W. W. (2000) A simple method for generating Gamma variables. 
// ACM Transactions on Mathematical Software Vol 26. No. 3: 363-372
// requirements are a 
//   - normal random number generator
//   - uniform random number generator
//
// [see Python program rgamma.py to test implementation, a copy is at the end of this file in comments] 
//
// returns a random number for shape parameter alpha>=1 and scale parameter beta=1
MYREAL __gamma_rand(MYREAL alpha)
{
  const MYREAL d = alpha - 1./3.;
  const MYREAL c = 1./sqrt(9. * d);
  MYREAL v,x,xx, u;
  while(1)
    {
      x = normal_rand(0.,1.);
      v = 1. + c * x;
      while(v<=0.0)
	{
	  x = normal_rand(0.,1.);
	  v = 1. + c * x;
	}
      v  = v * v * v;
      u  = UNIF_RANDUM();
      xx = x * x;
      if (u < 1.0-0.0331*xx*xx)
	return d*v;
      if (log(u) < 0.5*xx + d * (1.0 - v + log(v)))
	return d*v;
    }
}

/// returns a random number from a gamma distribution with shape parameter alpha and scale parameter beta
MYREAL gamma_rand(MYREAL alpha, MYREAL beta)
{
  MYREAL aa;
  if (alpha<1.0)
    {
      aa = 1.0 + alpha;
      return __gamma_rand(aa)*pow(UNIF_RANDUM(),(1.0/alpha)) * beta;
    }
  else
    {
      return __gamma_rand(alpha) * beta;
    }
}

/// returns a random number from a truncated gamma distribution
MYREAL trunc_gamma_rand(MYREAL alpha, MYREAL beta, MYREAL lower, MYREAL upper)
{
  MYREAL x = -1;
  long count = 0;
  while(count++ < HUNDRED)
    {
      x = gamma_rand(alpha,beta);
      if (lower < x)
	{
	  if (x < upper)
	    return x;
	}
    }
  if (x<lower)
    return lower;
  else
    {
      if (x>upper)
	return upper;
      else
	return x;
    }
}

MYREAL random_beta(MYREAL a, MYREAL b)
{
  MYREAL x,y; 
  x = gamma_rand(a,1.0); 
  y = gamma_rand(b,1.0); 
  return x/(x+y);
}


void assign_random_startsites(long **random_startsites, long fulllength, long shortsites, long rrepeats)
{
  boolean overlap;
  long r;
  long x;
  long i;
  long l=0;
  long *random_endsites;
  long repeats = shortsites > 0 ? rrepeats : rrepeats-1;
  random_endsites = (long*) mycalloc(rrepeats,sizeof(long));
  if (*random_startsites==NULL)
    *random_startsites = (long*) mycalloc(rrepeats,sizeof(long));
  else
    *random_startsites = (long*) myrealloc(*random_startsites, repeats * sizeof(long));
  for (r=0;r<repeats;r++)
    {
      x=RANDINT(0,fulllength-1-shortsites);
      i=0;
      overlap=FALSE;
      while(i<l)
	{
	  if (x < random_endsites[i] && x >= (*random_startsites)[i] - shortsites)
	    {
	      overlap=TRUE;
	      break;
	    }
	  i++;
	}
      if(overlap==FALSE)
	{
	  random_endsites[l]=x+shortsites;
	  (*random_startsites)[l]=x;
	  l++;
	}
    }
  myfree(random_endsites);
  qsort((void*) (*random_startsites), (size_t) rrepeats, sizeof(long), longcmp);
}

long random_from_freqlist(double *l, long nl)
{
  double *ll = (double*) mymalloc(sizeof(double)*(nl));
  memcpy(ll,l,sizeof(double)*nl);
  double s = 0.0;
  long i;
  for (i=0;i<nl;i++)
    {
      s += ll[i];
    }
  ll[0] = ll[0]/s;
  for (i=1;i<nl;i++)
    {
      ll[i] = ll[i]/s + ll[i-1];
    }
  double r = RANDUM();
  for (i=0;i<nl;i++)
    {
      if (r < ll[i])
	{
	  myfree(ll);
	  return i;
	}
    }
  myfree(ll);
  return nl-1;
}


/* copy of rgamma.py 
#!/usr/bin/env python 
#
# Marsaglia and Tang ? 2000
#
# adapted to python 3 January 2020
#
from __future__ import print_function
import random as ran
import math
import matplotlib.pyplot as plt
import sys

def truncrgama(a,b,l,u):
    while(1):
        x = rgama(a)
        if l < x*b < u:
            return x*b


def rgama(a):
    if a<1.0:
        aa = 1.0+a
        return _rgama(aa)*ran.random()**(1.0/a)
    else:
        return _rgama(a)

def _rgama(aa):
    d = aa - 1./3.
    c = 1./ math.sqrt(9. * d)
    while(1):
        x = ran.gauss(0.,1.)
        v = 1. + c*x
        while(v<=0.0):
            x = ran.gauss(0.,1.)
            v = 1. + c*x
        v = v * v * v
        u = ran.random()
        xx = x*x
        if u<1.-0.0331*xx*xx:
            return d*v
        if math.log(u) < 0.5*xx+ d*(1.-v+math.log(v)):
            return d*v

dist = None        
if len(sys.argv)>4:
    if 'var' in sys.argv:
        variable = True
        if len(sys.argv)==7:
            dist = float(sys.argv[-1])
    else:
        variable = False
else:
    variable = False
    
if len(sys.argv)<5:
    lower = 0.0
    upper = 50.0
else:
    lower = float(sys.argv[3])
    upper = float(sys.argv[4])
if len(sys.argv)<3:
    a = 2.5
    b = 0.5
else:
    a = float(sys.argv[1])
    b = float(sys.argv[2])
data=[]
print(f"Using alpha={a} and beta={b} in range {lower} to {upper}")
if variable:
    oldxx = a * b
    if not dist:
        dist=  (upper - lower)/20.0
    print("variable start with ",oldxx,dist)
    for i in range(100000):
        l = oldxx - dist
        u = oldxx + dist
        if l<lower:
            l = lower
        if u > upper:
            u = upper
        xx=truncrgama(a,b,l,u)
        oldxx = xx
        #print(xx)
        data.append(xx)
else:
    for i in range(100000):
        xx=truncrgama(a,b,lower,upper)
        #print(xx)
        data.append(xx)

plt.figure()
plt.hist(data,50)
plt.show()

*/
