/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 H E L P E R     R O U T I N E S 
 
 some math stuff and 
 string and file manipulation routines
 
 
 Peter Beerli started 1996, Seattle
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2006 Peter Beerli, Tallahassee FL
 
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
 
 
$Id: tools.c 2158 2013-04-29 01:56:20Z beerli $
-------------------------------------------------------*/
/* \file tools.c */
#pragma clang diagnostic ignored "-Wformat-nonliteral"
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#include "migration.h"
#include "sighandler.h"
#include "data.h"
#include "world.h"
#include "tree.h"
#include "random.h"
#include "migrate_mpi.h"
#include "options.h"
#include "laguerre.h"
#if !defined(TOOLSTEST) && !defined(PRIORTEST)
#include "pretty.h"
#endif
#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif

/* external globals*/
extern int myID;
#ifdef CAUTIOUS
extern boolean cautious;
#endif


/* prototypes ------------------------------------------- */
#ifndef HAVE_STRDUP
char *my_strdup(const char *s);
#endif
#ifndef HAVE_STRSEP
char *strsep(char **, const char *);
#endif
MYREAL lengthof (node * p);
MYINLINE node *crawlback (const node * theNode);
/*node *crawl(node * theNode); */
MYINLINE node *showtop (node * theNode);
void adjust_time (node * theNode, MYREAL tyme);
void adjust_time_all (node * theNode, MYREAL tyme);
void set_nodelet_time(node *p, double thetime);
void insert_migr_node (world_fmt * world, node * up, node * down,
                       migr_table_fmt * migr_table, long *migr_table_counter);
long offspring (node * mother, node ** myoffspring);
void children (node * mother, node ** brother, node ** sister);
void siblings (node * mother, node ** brother, node ** sister);
/* math tools */
MYREAL incompletegamma (MYREAL tx, MYREAL talpha);
MYREAL logincompletegamma (MYREAL tx, MYREAL talpha);
MYREAL uppergamma_reg (MYREAL tx, MYREAL talpha);
MYREAL uppergamma (MYREAL tx, MYREAL talpha);

MYREAL polygamma (long n, MYREAL z);
void invert_matrix (MYREAL **a, long nsize);
boolean nrcheck (MYREAL **m, MYREAL **tm, MYREAL *v, long nrows, MYREAL *r1,
                 MYREAL *r2, boolean do_newton);
MYREAL rannor (MYREAL mean, MYREAL sd);
char lowercase (int c);
char uppercase (int c);
MYREAL calc_sum (MYREAL *vector, long n);

void gamma_rates (MYREAL *rate, MYREAL *probcat, long categs, char *input);
void calc_gamma (MYREAL alpha, MYREAL *gama, long categs);

MYREAL mylgamma (MYREAL z);

MYREAL logfac (long n);

void onepass_mean_std_start(MYREAL *mean, MYREAL *std, long *n);
void onepass_mean_std_calc(MYREAL *mean, MYREAL *std, long *n, MYREAL x);
void onepass_mean_std_end(MYREAL *mean, MYREAL *std, long *n);

//MYINLINE double fast_exp(double y);
//MYINLINE float fast_log (float vval);

MYREAL bezier(MYREAL x0, MYREAL y0, MYREAL x1, MYREAL y1);

/* vector initialization */
void doublevec1d (MYREAL **v, long size);
void doublevec2d (MYREAL ***v, long size1, long size2);
void floatvec2d (float ***v, long size1, long size2);
void charvec1d (char **v, long size1);
void charvec2d (char ***v, long size1, long size2);
void intvec2d (long ***v, long size1, long size2);
void free_doublevec2d (MYREAL **v);
void free_floatvec2d (float **v);
void free_charvec2d (char **v);
void free_intvec2d (long **v);
void setdoublevec1d (MYREAL **v, MYREAL *w, long size);
void add_vector (MYREAL *result, MYREAL *v, long size);

/*filemanipulation */
void init_files (world_fmt * world, data_fmt * data, option_fmt * options);
void exit_files (world_fmt * world, data_fmt * data, option_fmt * options);
void openfile (FILE ** fp, char *filename, char *mode, char *perm);
#ifdef ZNZ
void znzopenfile (znzFile * fp, char *filename, char *mode, int use_compressed);
#endif
long read_savesum (world_fmt * world, option_fmt * options, data_fmt * data);
void write_savesum (world_fmt * world);

/* string manipulation */
void translate (char *text, char from, char to);
long count_words (char *text);
long count_char (char *text, char needle);
long locate_char (char *text, char needle);
char * char_position(char *s1, char* s2);
void fprintf2(FILE *file, long filesize, const char *fmt, ...);
void print_line (FILE * outfile, char c, long nn, long flag);
void sprint_line (char *buffer, char c, long nn, long flag);
void add_to_buffer(char *fp, long *bufsize, char **buffer, long *allocbufsize);
void increase_buffer(char **buffer, long *allocbufsize, long *bufsize,long amount);
/* time reporting */
void get_time (char *nowstr, char ts[]);
void get_runtime (char *runtime, time_t start, time_t end);
/*printing aid */
void print_llike (MYREAL llike, char *strllike);
int sprint_tabdigits(double x, char * temp);
/* searching and finding*/
boolean find (long i, long *list, long listlen);

/* conversion between the parameter schemes*/

long mstart (long pop, long numpop);
long mend (long pop, long numpop);
long mmstart (long pop, long numpop);
long mmend (long pop, long numpop);
long mm2m (long frompop, long topop, long numpop);
void m2mm (long i, long numpop, long *frompop, long *topop);
long m2mml (long i, long numpop);
long m2mml2 (long i, long topop, long numpop);

void set_paramstr(char * paramstr, long j, world_fmt * world);
boolean shortcut(long j0, world_fmt *world, long *j);
/* private functions */
MYREAL alnorm (MYREAL x, int up);
void lu_decomp (MYREAL **m, long *indeks, long nrows);
void lu_substitution (MYREAL **m, long *indeks, MYREAL *v, long nrows);
MYREAL d1mach (long i);
long i1mach (long i);
int dpsifn (MYREAL *x, long *n, long kode, long m, MYREAL *ans, long *nz,
            long *ierr);
MYREAL find_chi (long df, MYREAL prob);
MYREAL probchi (long df, double chi);
MYREAL chisquare (long df, MYREAL alpha);
MYREAL chiboundary (long zeros, long nonzeros, MYREAL alpha);
void clearscreen(void);
void stringvec2d (char ****v, long size1, long size2, long size3);
void free_stringvec2d (char ***v,long size1, long size2);
void print_warning2(FILE *file, world_fmt *world);
double myerf1(double x);
double myerfc1(double x);

double wew ( double x, double *en );
double get_time_for_growth(double theta0, double growth, double k, double t0);
double interval_growth(double r, double t0, double theta0, double growth, double k, double tmin, double tmax);
//##

extern long  calculate_newpop_numpop(option_fmt *options, data_fmt *data);

#ifdef MESS
long *check_collection;
long check_collection_count;
extern long unique_id_global;
#endif

/*commandshell utility*/
void clearscreen(void)
{
#ifdef WINDOWS
  system("cls");
#else
  system("clear");
#endif
}

#ifndef HAVE_STRDUP
/* strdup is not ansi-C */
char *my_strdup(const char *s)
{
  char *p = mymalloc(strlen(s) + 1);
  if(p)
    {
      strcpy(p, s);
    }
  else
    p=NULL;
  return p;
}
#endif

#ifndef HAVE_STRSEP
/*-
 * Copyright (c) 1990, 1993
 *      The Regents of the University of California.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *      This product includes software developed by the University of
 *      California, Berkeley and its contributors.
 * 4. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

// #include "config.h"

// #include <sys/types.h>

// #include <string.h>
// #include <stdio.h>

// #ifndef HAVE_STRSEP
char *strsep(char **, const char *);
// #endif

// #if defined(LIBC_SCCS) && !defined(lint)
// static char sccsid[] = "@(#)strsep.c    8.1 (Berkeley) 6/4/93";
// #endif /* LIBC_SCCS and not lint */
// #ifndef lint
// static const char rcsid[] =
//   "$FreeBSD: src/lib/libc/string/strsep.c,v 1.2.12.1 2001/07/09 23:30:07 obrien Exp $";
// #endif

/*
 * Get next token from string *stringp, where tokens are possibly-empty
 * strings separated by characters from delim.
 *
 * Writes NULs into the string at *stringp to end tokens.
 * delim need not remain constant from call to call.
 * On return, *stringp points past the last NUL written (if there might
 * be further tokens), or is NULL (if there are definitely no more tokens).
 *
 * If *stringp is NULL, strsep returns NULL.
 */
char *
strsep(register char **stringp, register const char *delim)
{
    register char *s;
    register const char *spanp;
    register int c, sc;
    char *tok;
    
    if ((s = *stringp) == NULL)
        return (NULL);
    for (tok = s;;) {
        c = *s++;
        spanp = delim;
        do {
            if ((sc = *spanp++) == c) {
                if (c == 0)
                    s = NULL;
                else
                    s[-1] = 0;
                *stringp = s;
                return (tok);
            }
        } while (sc != 0);
    }
    /* NOTREACHED */
}
#endif


/* vector initialization */
void
doublevec1d (MYREAL **v, long size)
{
    *v = (MYREAL *) mycalloc (size, sizeof (MYREAL));
}

/// allocate a 2D array with MYREALS
void
doublevec2d (MYREAL ***v, long size1, long size2)
{
    long i;
    *v = (MYREAL **) mycalloc (size1, sizeof (MYREAL *));
    (*v)[0] = (MYREAL *) mycalloc ((size1 * size2), sizeof (MYREAL));
    for (i = 1; i < size1; i++)
    {
        (*v)[i] = (*v)[0] + size2 * i;
    }
}

/// allocate a 2D array with floats
void
floatvec2d (float ***v, long size1, long size2)
{
    long i;
    *v = (float **) mycalloc (size1, sizeof (float *));
    (*v)[0] = (float *) mycalloc ((size1 * size2), sizeof (float));
    for (i = 1; i < size1; i++)
    {
        (*v)[i] = (*v)[0] + size2 * i;
    }
}

///
/// allocates columns and rows for a matrix of chars
/// use freecharvec2d() to free this.
void
charvec2d (char ***v, long size1, long size2)
{
    long i;
    *v = (char **) mycalloc (size1, sizeof (char *));
    (*v)[0] = (char *) mycalloc ((size1 * size2), sizeof (char));
    for (i = 1; i < size1; i++)
    {
        (*v)[i] = (*v)[0] + size2 * i;
    }
}

///
/// allocates columns and rows for a matrix of chars
/// use freecharvec2d() to free this.
void stringvec2d (char ****v, long size1, long size2, long size3)
{
  long i,j;
    *v = (char ***) mycalloc (size1, sizeof (char **));
    (*v)[0] = (char **) mycalloc ((size1 * size2), sizeof (char *));
    for (i = 1; i < size1; i++)
    {
        (*v)[i] = (*v)[0] + size2 * i;
	for(j=0; j < size2; j++)
	  {
	    (*v)[i][j] = (char *) mycalloc(size3,sizeof(char));
	  }
    }
}

///
/// allocates columns and rows for a matrix of strings
/// use free_intvec2d() to free this.
void
intvec2d (long ***v, long size1, long size2)
{
    long i;
    *v = (long **) mycalloc (size1, sizeof (long *));
    (*v)[0] = (long *) mycalloc ((size1 * size2), sizeof (long));
    for (i = 1; i < size1; i++)
    {
        (*v)[i] = (*v)[0] + size2 * i;
    }
}


/// free 2D vector of MYREAL
void
free_doublevec2d (MYREAL **v)
{
    myfree(v[0]);
    myfree(v);
}

/// free 2D vector of floats
void
free_floatvec2d (float **v)
{
    myfree(v[0]);
    myfree(v);
}

void
free_charvec2d (char **v)
{
    myfree(v[0]);
    myfree(v);
}
void
free_intvec2d (long **v)
{
    myfree(v[0]);
    myfree(v);
}

void free_stringvec2d (char ***v,long size1, long size2)
{
  long j,z;
  for(j=0;j<size1;j++)
    for(z=0;z<size2;z++)
      myfree(v[j][z]);
  myfree(v[0]);
  myfree(v);
}

void
setdoublevec1d (MYREAL **v, MYREAL *w, long size)
{
    doublevec1d (v, size);
    memcpy (*v, w, sizeof (MYREAL) * (size_t) size);
}

void
add_vector (MYREAL *result, MYREAL *v, long size)
{
    long i;
    for (i = 0; i < size; i++)
        result[i] += v[i];
}


MYREAL
logfac (long n)
{
    /* log(n!) values were calculated with Mathematica
       with a precision of 30 digits */
    switch (n)
    {
    case 0:
        return 0.;
    case 1:
      return 0.;
    case 2:
      return 0.693147180559945309417232121458;
    case 3:
      return 1.791759469228055000812477358381;
    case 4:
      return 3.1780538303479456196469416013;
    case 5:
      return 4.78749174278204599424770093452;
    case 6:
      return 6.5792512120101009950601782929;
    case 7:
      return 8.52516136106541430016553103635;
    case 8:
      return 10.60460290274525022841722740072;
    case 9:
        return 12.80182748008146961120771787457;
    case 10:
      return 15.10441257307551529522570932925;
    case 11:
      return 17.50230784587388583928765290722;
    case 12:
      return 19.98721449566188614951736238706;
    case 13:
      return 22.5521638531234228855708498286;
    case 14:
      return 25.1912211827386815000934346935;
    case 15:
      return 27.8992713838408915660894392637;
    case 16:
      return 30.6718601060806728037583677495;
    case 17:
      return 33.5050734501368888840079023674;
    case 18:
      return 36.3954452080330535762156249627;
    case 19:
      return 39.3398841871994940362246523946;
    case 20:
      return 42.3356164607534850296598759707;
    case 21:
      return 45.3801388984769080261604739511;
    case 22:
      return 48.4711813518352238796396496505;
    case 23:
      return 51.6066755677643735704464024823;
    case 24:
      return 54.7847293981123191900933440836;
    case 25:
      return 58.0036052229805199392948627501;
    case 26:
      return 61.2617017610020019847655823131;
    case 27:
      return 64.5575386270063310589513180238;
    case 28:
      return 67.8897431371815349828911350102;
    case 29:
      return 71.2570389671680090100744070426;
    case 30:
      return 74.6582363488301643854876437342;
    default:
        return LGAMMA (n + 1.);
    }
}

///
/// initialize the mean and standard deviation calculations 
/// for the one-pass calculations that allows to throw away 
/// the intermediate results
void onepass_mean_std_start(MYREAL *mean, MYREAL *std, long *n)
{
  (*mean) = 0.0;
  (*std)  = 0.0;
  (*n)    = 0;
}

///
/// calculates the mean and standard deviation using a one-pass filter 
/// that allows to throw away the intermediate results
void onepass_mean_std_calc(MYREAL *mean, MYREAL *std, long *n, MYREAL x)
{
  
  long   nn    = (*n)+1;
  MYREAL mm    = (*mean);
  MYREAL delta = x - mm;
  mm           = mm + delta / nn;
  (*std)       = (*std) + delta * (x - mm);
  (*mean)      = mm;
  (*n)         = nn;
}

///
/// finishes one-pass mean and std calculations
void onepass_mean_std_end(MYREAL *mean, MYREAL *std, long *n)
{
  (void) mean;
  if(*n > 1)
    (*std) = sqrt((*std) / ((*n) - 1));
  else
    (*std) = -999;
}
 
void combine_meanstd(MYREAL *mean1, MYREAL *sd1, long *n1, MYREAL mean2, MYREAL sd2, long n2)
{
  MYREAL _n1 = (MYREAL) *n1;
  MYREAL _n2 = (MYREAL) n2;
  MYREAL ssd1 = (*sd1) * (*sd1);
  MYREAL ssd2 = sd2 * sd2;
  MYREAL mm1 = (*mean1) *(*mean1);
  MYREAL mm2 = mean2 * mean2;
  MYREAL denom = (_n1-1.) * ssd1 + (_n2-1.)*ssd2 + _n1*_n2/(_n1+_n2) * (mm1 + mm2 - 2.0*(*mean1)*mean2);
  *sd1 = sqrt(denom / (_n1 + _n2 - 1.0)); 
  *mean1 = (_n1 * (*mean1) + _n2*mean2) / (_n1+_n2);
  *n1 = *n1 + n2;
}

    


#ifndef PRIORTEST
/*garbage collection and recycle*/
#ifdef MESS
boolean find_in_list(long id)
{
  long i;
  for(i=0; i < check_collection_count; i++)
    {
      if(id == check_collection[i])
	return TRUE;
    }
  return FALSE;
}

/// tags the lineage that was added the last time the tree was changed.
///
boolean traverse_check_id(node *theNode, long id)
{
  boolean found = FALSE;
  if(theNode->type=='r')
    { 
      if(theNode->id == id || theNode->next->id == id || theNode->next->next->id == id)
	return TRUE;	
      else
	{
	  found = traverse_check_id(theNode->next->back, id);
	  return found;
	}
    }
  if(theNode->type == 'i')
    {
      if(theNode->id == id || theNode->next->id == id || theNode->next->next->id == id)
	return TRUE;	
      else
	{
	  found = traverse_check_id(theNode->next->back, id);
	  if(found == TRUE)
	    return found;
	  found = traverse_check_id(theNode->next->next->back, id);
	  return found;
	}
    }
  if(theNode->type == 'm' || theNode->type == 'd')
    {
      if(theNode->id == id || theNode->next->id == id)
	return TRUE;	
      else
	{
	  found = traverse_check_id(theNode->next->back, id);
	  return found;
	}
    }
  if(theNode->type == 't')
    {
      if(theNode->id == id)
	return TRUE;	
      else
	return FALSE;
    }
  error("node with no type");
  return FALSE;
}

void traverse_checker(node *root)
{
  long i;
  long x;
  for(i=0;i<check_collection_count; i++)
    {
      x = check_collection[i];
      if(!traverse_check_id(root,x))
	printf ("Node id %li not found\n",x);
    }
}
#endif /*MESS*/

void start_node_collection(world_fmt *world)
{
  long i;
#ifdef MESS
  unique_id_global=0;
#endif
  world->node_collection_allocated = 2 * HUNDRED  + 2 * world->numpop * MIGRATION_LIMIT;
  world->node_collection = (node **) mycalloc(world->node_collection_allocated,sizeof(node*));
#ifdef MESS
  check_collection = (long *) calloc(world->node_collection_allocated, sizeof(long));
#endif
  for(i = 0; i < world->node_collection_allocated; i++)
    {
      world->node_collection[i] = (node *) calloc(1, sizeof(node));
#ifdef MESS
      (*world->node_collection[i]).id = unique_id_global++;
#endif
    }
  world->node_collection_count = world->node_collection_allocated;
}

void stop_node_collection(world_fmt *world)
{
  long i;
  for(i=0;i<world->node_collection_count;i++)
    {
      if(world->node_collection[i] != NULL)
	{
	  free(world->node_collection[i]);
	}
    }
  myfree(world->node_collection);
}

void set_nodelet_time(node *p, double thetime)
{
  //multi-branch safe
  node * q = p;
  if (q->top)
    {
      q->tyme= thetime;
      q = q->next;
      while(q!=NULL && !q->top)
	{
	  q->tyme = thetime;
	  q = q->next;
	}
    }
  else
    error("set_nodelet_time() failed");
}


void collect_nodelet(world_fmt *world, node *p)
{
#ifdef MESS
  long i;
#endif
  long oldalloc;

  if(world->node_collection_count < world->node_collection_allocated)
    {
      world->node_collection[world->node_collection_count] = p;
#ifdef MESS
      for(i=0;i<check_collection_count;i++)
	{
	  if(p->id == check_collection[i])
	    break;
	}
      if(i==check_collection_count)
	{
	  warning("found node not accounted for: %li",p->id);
	  error("");
	}
      else
	{
	  check_collection_count--;
	  check_collection[i] = check_collection[check_collection_count];
	}
#endif      
      world->node_collection_count++;
      // printf("d%li %c\n",p->id, p->type);
    }
  else
    {
#ifdef DEBUG
      warning("in node collection phase we should not need to allocate new space\n");
#endif
      oldalloc = world->node_collection_allocated;
      world->node_collection_allocated += 300;
      world->node_collection = (node **) realloc(world->node_collection, (size_t) world->node_collection_allocated * sizeof(node*));
      memset(world->node_collection+oldalloc,0, 300 * sizeof(node*));
      world->node_collection[world->node_collection_count] = p;
#ifdef MESS
      for(i=0;i<check_collection_count;i++)
	{
	  if(p->id == check_collection[i])
	    break;
	}
      if(i==check_collection_count)
	{
	  warning("found node not accounted for: %li",p->id);
	  error("");
	}
      else
	{
	  check_collection_count--;
	  check_collection[i] = check_collection[check_collection_count];
	}
#endif
      world->node_collection_count++;
#ifdef DEBUG
      FPRINTF(stdout, "%i> Heat=%f: Collect_nodelet: recycler increased size to: %li/%li\n", myID, 1./world->heat, world->node_collection_count, world->node_collection_allocated);
#endif
    }
  //printf("\n");
}

///
/// dispense memory for nodes from a list of avialable nodes, allocates more memory 
/// no nodes are available in the list if nodes
node *  dispense_nodelet(world_fmt *world)
{
  long i;
  node *p = NULL;
  if(world->node_collection_count > 0)
    {
      world->node_collection_count--;
      p = world->node_collection[world->node_collection_count];
#ifdef MESS
      check_collection[check_collection_count++] = p->id;
#endif
      world->node_collection[world->node_collection_count] = NULL;
    }
  else
    {
      //all node memory is checked out, we
      //increase available nodes in chunks of 300
      world->node_collection_allocated += 300;
      world->node_collection = (node **) myrealloc(world->node_collection,sizeof(node *)* (size_t) world->node_collection_allocated);
      // set new slots to NULL
      memset(world->node_collection + world->node_collection_allocated - 300, 0,sizeof(node *)*300);
      // we allocate memory of 300 nodes and start at the 0 element because all
      // elements allocated earlier are checked out and the node_collection array
      // should contain only NULLs
      for(i = 0; i < 300; i++)
	{
	  world->node_collection[i] = (node *) mycalloc(1, sizeof(node));
#ifdef MESS
	  world->node_collection[i]->id = unique_id_global++;
#endif
	}
#ifdef MESS
      check_collection = (long *) realloc(check_collection, world->node_collection_allocated * sizeof(long));
#endif
      world->node_collection_count += 300;
      world->node_collection_count--;
      p = world->node_collection[world->node_collection_count]; // dispense
      world->node_collection[world->node_collection_count] = NULL; // remove from list
#ifdef MESS
      check_collection[check_collection_count++] = p->id;
#endif
    }
  return p;
}

///
/// swaps the nodelet recycler between different temperatured chains
/// 
void swap_node_collection(world_fmt * tthis, world_fmt * tthat)
{
#ifdef DEBUG
  fprintf(stdout,"%i> swapping node_collection parts between temps: %f %f\n",myID,tthis->heat, tthat->heat);
#endif
  node **temp;
  //long tempalloc;
  long tempcount;
  long tempallocated;
  temp = tthis->node_collection;
  tempcount = tthis->node_collection_count;
  tempallocated = tthis->node_collection_allocated;

  tthis->node_collection = tthat->node_collection;
  tthis->node_collection_count = tthat->node_collection_count;
  tthis->node_collection_allocated = tthat->node_collection_allocated;

  tthat->node_collection = temp;
  tthat->node_collection_count = tempcount;
  tthat->node_collection_allocated = tempallocated;
}

/*FILEMANIPULATION======================================================= */
void
init_files (world_fmt * world, data_fmt * data, option_fmt * options)
{
  long *tempdb;
  long maxfilenum = MYMAXFILENUM; //see definitions.h all files have a number and a filename;
#ifdef MPI
  long i,j;
#ifdef PARALIO
  long l;
  boolean my_file_open_error;
  char error_string[LINESIZE];
  int length_of_error_string, error_class;
  if(options->has_bayesmdimfile)
    {
      l = strlen(options->bayesmdimfilename);
      sprintf(options->bayesmdimfilename+l,".%d",myID);
    }
#endif
#endif
  
  tempdb = (long *) mycalloc(maxfilenum,sizeof(long));
  
  if (myID == MASTER)
    {
#ifdef MPI
      setup_filehandle_db((void *) stdout, world, options,data);
#endif
      if (!options->readsum || options->checkpointing)
	{
	  openfile (&data->infile, options->infilename, "r+",  NULL);
	  if (options->usertree)
	    openfile (&data->utreefile, options->utreefilename, "r+",  NULL);
	  if (options->weights)
	    openfile (&data->weightfile, options->weightfilename, "r+",  NULL);
	  if (options->categs > 1)
	    openfile (&data->catfile, options->catfilename, "r+",  NULL);
	  if (options->dist)
	    openfile (&data->distfile, options->distfilename, "r+",  NULL);
	  if (options->writesum)
	    openfile (&data->sumfile, options->sumfilename, "w+",  NULL);
	}
      else
	{
	  if(!options->bayes_infer)
	    openfile (&data->sumfile, options->sumfilename, "r+",  NULL);
	}
      
      openfile (&world->outfile, options->outfilename, "w+",  NULL);
#ifdef MPI
      setup_filehandle_db((void *) world->outfile, world, options,data);
#endif
      
      if (options->treeprint > 0 && (!options->readsum))
	{
	  openfile (&world->treefile, options->treefilename, "w+",  NULL);
#ifdef MPI
	  setup_filehandle_db((void *) world->treefile, world, options,data);
#endif
	}
      
      if (options->mighist)
	{
	  if(options->datatype != 'g')
	    {
	      if(options->checkpointing)
		openfile (&world->mighistfile, options->mighistfilename, "a+", NULL);
	      else
		openfile (&world->mighistfile, options->mighistfilename, "w+", NULL);
	    }
	  else
	    {
	      openfile (&world->mighistfile, options->mighistfilename, "r+", NULL);
	    }
#ifdef MPI
	  setup_filehandle_db((void *) world->mighistfile, world, options, data);
#endif			
	  if(options->skyline)
	    {
	      if(options->datatype != 'g')
		{
		  if(options->checkpointing)
		    openfile (&world->skylinefile, options->skylinefilename, "a+", NULL);
		  else
		    openfile (&world->skylinefile, options->skylinefilename, "w+", NULL);
		}
	      else
		{
		  openfile (&world->skylinefile, options->skylinefilename, "r+", NULL);
		}
	    }
#ifdef MPI
	  setup_filehandle_db((void *) world->skylinefile, world, options, data);
#endif			
	}
      
      if (options->writelog)
	{
	  openfile (&options->logfile, options->logfilename, "w+",  NULL);
#ifdef MPI
	  setup_filehandle_db((void *) options->logfile, world, options, data);
#endif			
	}
      
      if (options->mixplot)
	{
	  openfile (&options->mixfile, options->mixfilename, "w+",  NULL);
#ifdef MPI
	  setup_filehandle_db((void *) options->mixfile, world, options, data);
#endif			
	}

      if (options->recorddivtime)
	{
	  openfile (&world->divtimefile, options->divtimefilename, "w+",  NULL);
#ifdef MPI
	  setup_filehandle_db((void *) world->divtimefile, world, options, data);
#endif			
	}


      if (options->aic)
	{
	  openfile (&options->aicfile, options->aicfilename, "w+",  NULL);
#ifdef MPI
	  setup_filehandle_db((void *) options->aicfile, world, options, data);
#endif
	}			
      if (options->plot)
	{
	  switch (options->plotmethod)
	    {
	    case PLOTALL:
	      openfile (&world->mathfile, options->mathfilename, "w+", 
			NULL);
#ifdef MPI
	      setup_filehandle_db((void *) world->mathfile, world, options, data);
#endif
	      break;
	    default:  /*e.g. 0 this create just the plots in outfile */
	      break;
	    }
	}
      if (options->geo)
	{
	  openfile (&data->geofile, options->geofilename, "r+",  NULL);
	}
      if (options->has_datefile)
	{
	  openfile (&data->datefile, options->datefilename, "r+",  NULL);
	}
#ifdef UEP
      
      if (options->uep)
	openfile (&data->uepfile, options->uepfilename, "r+",  NULL);
#endif
      if (options->bayes_infer)
	{
	  if(options->has_bayesfile)
	    openfile (&world->bayesfile, options->bayesfilename, "w+",  NULL);
	  
	  if(options->has_bayesmdimfile)
	    {
	      if(options->datatype != 'g' && !options->checkpointing)
		{
#ifdef MPI
#ifndef PARALIO
#ifdef ZNZ
		  znzopenfile (&world->bayesmdimfile, options->bayesmdimfilename, "w", options->use_compressed);
#else
		  openfile (&world->bayesmdimfile, options->bayesmdimfilename, "w+",  NULL);
#endif /*znz*/
#endif /*not paralio*/			  
#else
		  // not MPI
#ifdef ZNZ
		  znzopenfile (&world->bayesmdimfile, options->bayesmdimfilename, "w",  options->use_compressed);
#else
		  openfile (&world->bayesmdimfile, options->bayesmdimfilename, "w+",  NULL);
#endif
#endif /*mpi*/
		}
	      else
		{
#ifdef ZNZ
		  znzopenfile (&world->bayesmdimfile, options->bayesmdimfilename, "r",  options->use_compressed);
#else
		  openfile(&world->bayesmdimfile, options->bayesmdimfilename,"r+",NULL);
#endif
		}
	    }
#ifdef MPI
	  if(options->has_bayesfile)
	    setup_filehandle_db((void *) world->bayesfile, world, options,data);
	  if(options->has_bayesmdimfile)
	    {
#ifdef PARALIO
	      fprintf(stdout,"%i> before open %s\n",myID,options->bayesmdimfilename);
	      my_file_open_error = MPI_File_open(MPI_COMM_SELF, options->bayesmdimfilename, 
						 MPI_MODE_CREATE | MPI_MODE_WRONLY ,
						 MPI_INFO_NULL, &world->mpi_bayesmdimfile);
	      if (my_file_open_error != MPI_SUCCESS) 
		{
		  MPI_Error_class(my_file_open_error, &error_class);
		  MPI_Error_string(error_class, error_string, &length_of_error_string);
		  printf("%i> %s %s\n", myID, error_string, options->bayesmdimfilename);
		  MPI_Error_string(my_file_open_error, error_string,
				   &length_of_error_string);
		  printf("%i> %s\n", myID, error_string);
		  my_file_open_error = TRUE;
		}
	      MPI_File_set_size(world->mpi_bayesmdimfile, 0);		   
#else
	      setup_filehandle_db((void *) world->bayesmdimfile, world, options,data);
#endif
	    } //has bayesmdimfile
#endif
	}			
#ifdef MPI
      tempdb[0]=filenum;
      for(i=0, j=1;i<filenum;i++)
	{
	  tempdb[j++] = (long) filedb[i].file;
	  tempdb[j++] = (long) filedb[i].handle;
	}
      MYMPIBCAST (tempdb, maxfilenum, MPI_LONG, MASTER, comm_world);
#endif
    }
  else
    {
      // all non myID==MASTER
#ifdef MPI
#ifdef PARALIO
      if(options->has_bayesmdimfile)
	{
	  MPI_Info info;
	  fprintf(stdout,"%i> before open %s\n",myID,options->bayesmdimfilename);
	  my_file_open_error = MPI_File_open(MPI_COMM_SELF, options->bayesmdimfilename, 
					     MPI_MODE_CREATE | MPI_MODE_WRONLY ,
					     MPI_INFO_NULL, &world->mpi_bayesmdimfile);
	  if (my_file_open_error != MPI_SUCCESS) 
	    {
	      MPI_Error_class(my_file_open_error, &error_class);
	      MPI_Error_string(error_class, error_string, &length_of_error_string);
	      printf("%i> @%s %s\n", myID, error_string, options->bayesmdimfilename);
	      MPI_Error_string(my_file_open_error, error_string,
			       &length_of_error_string);
	      printf("%i> @%s\n", myID, error_string);
	      my_file_open_error = TRUE;
	    }
	  MPI_File_set_size(world->mpi_bayesmdimfile, 0);
	}
#endif
      MYMPIBCAST (tempdb, maxfilenum, MPI_LONG, MASTER, comm_world);
      filenum = tempdb[0]; 
      for(i=0, j=1;i<filenum;i++)
	{
	  filedb[i].file = (FILE *) tempdb[j++];
	  filedb[i].handle =  tempdb[j++];
	}
      //assumess same order as in the master
      j=1; //zero is stdout
      world->outfile = filedb[j++].file;
      if (options->treeprint > 0 && (!options->readsum))
	world->treefile = filedb[j++].file;
      if (options->mighist)
	world->mighistfile = filedb[j++].file;
      if (options->skyline)
	world->skylinefile = filedb[j++].file;
      if (options->writelog)
	options->logfile = filedb[j++].file;
      if (options->recorddivtime)
	options->divtimefile = filedb[j++].file;
      if (options->aic)
	options->aicfile = filedb[j++].file;
      if (options->plot && options->plotmethod == PLOTALL)
	world->mathfile = filedb[j++].file;
      if (options->bayes_infer)
	{
	  if(options->has_bayesfile)
	    world->bayesfile = filedb[j++].file;
	  // TEST: experiment with parallele I/O
	  if(options->has_bayesmdimfile)
	    {
#ifdef PARALIO
	      //has to be  dealt with before the broadcast	      MPI_File_open(comm_world,options->bayesmdimfilename,MPI_MODE_CREATE,MPI_INFO_NULL, &world->mpi_bayesmdimfile);
#else
#ifdef ZNZ
	      world->bayesmdimfile = (znzFile) filedb[j++].file;
#else
	      world->bayesmdimfile = filedb[j++].file;
#endif
#endif
	    }
	}
#endif /*MPI*/
    }
  myfree(tempdb);
}

void
exit_files (world_fmt * world, data_fmt * data, option_fmt * options)
{
    if(myID==MASTER)
    {
      if (!options->readsum)
	{
	  FClose (data->infile);
	  //  if (options->treeprint > 0)
	  //  FClose (world->treefile);
	  if(options->usertree)
	    FClose(data->utreefile);
	  if (options->weights)
	    FClose (data->weightfile);
	  if (options->categs > 1)
	    FClose (data->catfile);
	  if (options->dist)
	    FClose (data->distfile);
	}
      
      if (options->writesum || options->readsum)
	FClose (data->sumfile);
      
      FClose (world->outfile);
      
      if (options->writelog)
	FClose (options->logfile);
      if (options->mighist)
	FClose (world->mighistfile);
      if (options->geo)
	FClose (data->geofile);
      if (options->has_datefile)
	FClose (data->datefile);
      if(options->recorddivtime)
	 FClose(world->divtimefile);
      if (options->aic)
	FClose (options->aicfile);
      if (options->plot && options->plotmethod == PLOTALL)
	FClose (world->mathfile);
#ifdef UEP
      
      if (options->uep)
	FClose (data->uepfile);
#endif
      if (options->bayes_infer)
	{
	  if(options->has_bayesfile)
	    FClose (world->bayesfile);
	  if(options->has_bayesmdimfile)
#ifdef ZNZ
	    znzClose (world->bayesmdimfile);
#else
	    FClose (world->bayesmdimfile);
#endif
	}
    }
}

/* string manipulation ================================== */
/* Converts any character from to character to in string text */
void
translate (char *text, char from, char to)
{
    int i, j, gap = 0;
    while (text[gap] == from)
        gap++;
    for (i = gap, j = 0; text[i] != '\0'; i++)
    {
        if (text[i] != from)
        {
            text[j++] = text[i];
        }
        else
        {
            if (text[i - 1] != from)
            {
                text[j++] = to;
            }
        }
    }
    text[j] = '\0';
}
/// unpad from the end of the word until there is none of the
/// specified character left
void
unpad (char *text, char removechars[])
{
  long gap = (long) strlen(text)-1;
  while (gap >= 0 && strchr(removechars,text[gap]))
        gap--;
    text[gap+1] = '\0';
}

///
/// get next word from list with list of delimiters
/// uses strsep
void
get_next_word(char **instring, char *delimiters, char **nextword)
{
  *nextword = strsep(instring, delimiters);
   while(*nextword != NULL && strchr(delimiters,(*nextword)[0]))
  {
      *nextword = strsep(instring,delimiters);
  }
}


/*===============================================
  count words in a string delimited by delimiter
*/
long
count_words (char *text)
{
    long counts = 0;
    char *pt = text;
    while (isspace (*pt) && *pt != '\0')
        pt++;
    while (*pt != '\0')
    {
        while (!isspace (*pt) && *pt != '\0')
            pt++;
        while (isspace (*pt) && *pt != '\0')
            pt++;
        counts++;
    }
    return counts;
}

long count_char (char *text, char needle)
{
  char nn = needle;
  long count = 0;
  char *i = text;
  while(*i !='\0')
    {
      if (nn == *i)
	count++;
      i++;
    }
  return count;
}

long locate_char (char *text, char needle)
{
  char nn = needle;
  long count = 0;
  char *i = text;
  while(*i !='\0')
    {
      if (nn != *i)
	{
	  count++;
	}
      else
	return count;
      i++;
    }
  return count;
}

char * char_position(char *s1, char* s2)
{
  char *val;
  boolean found=FALSE;
  while (*s1!='\0' && !found)
    {
      char c=*s1++;
      val = strchr(s2,c);
      if (val!=NULL)
	return val;
    }
  return NULL;
}


/*===============================================
 timer utility
 
 ts = "%c" -> time + full date (see man strftime)
      = "%H:%M:%S" -> time hours:minutes:seconds */

void
get_time (char *nowstr, char ts[])
{
#ifdef NOTIME_FUNC
    switch (strlen (ts))
    {
    case 2:
        strcpy (nowstr, " ");
        break;
    case 3:
        strcpy (nowstr, "  ");
        break;
    case 8:
        strcpy (nowstr, "        ");
        break;
    default:
        strcpy (nowstr, " ");
        break;
    }
#else
    time_t nowbin;
    struct tm *nowstruct;
    if (time (&nowbin) != (time_t) - 1)
    {
        nowstruct = localtime (&nowbin);
#ifdef WIN32
	if (!strcmp(ts,"%s"))
	  {
	    //windows 7 does not have the "%s" format for strftime()
	    sprintf(nowstr,"%i",nowbin);
	  }
	else
	  strftime (nowstr, LINESIZE, ts, nowstruct);
#else
	strftime (nowstr, LINESIZE, ts, nowstruct);
#endif
    }
#endif
}

void get_runtime (char *runtime, time_t start, time_t end)
{
  time_t delta = end - start;
  time_t seconds = delta % 60;
  time_t minutes = ((time_t) delta/60) % 60;
  time_t hours   = ((time_t) delta/3600) % 24;
  time_t days    = ((time_t) delta/86400) % 365;
  sprintf(runtime, "Runtime:%04li:%02li:%02li:%02li",days,hours,minutes,seconds);
}

/*===============================================
 printer utility
 */
void
print_llike (MYREAL llike, char *strllike)
{
    if (fabs (llike) > 10e20)
    {
        sprintf (strllike, "%cInfinity ", llike < 0 ? '-' : ' ');
    }
    else
        sprintf (strllike, "%-10.5f", llike);
}

int sprint_tabdigits(double x, char * temp)
{
  int fmt=5;
  long digits = (long) (floor(log10(x)));
  switch(digits)
    {
    case -8:
    case -6:
    case -5:
      fmt = 10;
      break;
    case -4:
    case -3:
      fmt = 8;
      break;
    case -2:
    case -1:
      fmt= 5;
      break;
    case 0:
    case 1:
      fmt = 4;
      break;
    case 2:
      fmt = 2;
      break;
    case 3:
      fmt = 1;
      break;
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
      fmt = 0;
      break;
    default:
      if(digits<-8)
	fmt=20;
      else
	fmt = 5;
      break;
    }
  int c = sprintf(temp,"\t%.*f",fmt, x);
  return c;
}


///
/// remove traling whitespace from a *string
void remove_trailing_blanks(char **filename)
{
  long pos;
  pos = (long) strlen(*filename);
  while (pos>0 && isspace( (*filename)[pos-1]))
    pos--;
  (*filename)[pos] = '\0';   
}

///
/// trim white space on both sides, keep pointer
/// example "   this is a test      "
/// find total number of characters:23, start of non-whitespace: 3 (start at zero)
/// end is at 13
void trim(char **line)
{
  char *tmp = *line;
  unsigned long start=0;
  while (isspace(*tmp)) 
    {
      start++;
      tmp++;
    }
  unsigned long total = strlen(*line);
  unsigned long end = total;
  while (isspace(*(*line + end - 1))) 
    {
      end--;
    }
  memmove((*line),(*line)+start,end-start);
  (*line)[end-start]='\0';
}
///
/// trim white space on both sides, keep pointer
/// example "   this is a test      "
/// find total number of characters:23, start of non-whitespace: 3 (start at zero)
/// end is at 1
/// the start of the line stays where it is, the end is fixed with a '\0'
/// the start index is returned.
long nondestruct_trim(char **line)
{
  char *tmp = *line;
  long start=0;
  while (isspace(*tmp)) 
    {
      start++;
      tmp++;
    }
  long total = (long) strlen(*line);
  long end = total;
  while (isspace(*(*line + end - 1))) 
    {
      end--;
    }
  (*line)[end-start]='\0';
  return start;
}

///
/// read a single word from a string up to the position of 
/// the breakchar and return that position 
// the ignore_repeats flag jumps over repeated delimiters and whitespace
// returns last character read
long read_word_delim(char *input, char *word, char *delim, boolean ignore_repeats)
{
  long i = 1;
  int ch;
  long count = 0;
  ch = input[0];
  if(ch=='\0')
    return 0;
  if(ignore_repeats)
    {
      while (isspace (ch) && strchr(delim,ch)!=NULL)
	ch = input[i++];
    }
  else
    {
      while (isspace (ch) && strchr(delim,ch)==NULL)
	ch = input[i++];
    }
  while(strchr(delim, ch)==NULL && ch!=EOF)
    {
      word[count++] = (char) ch;
      ch = input[i++];
    }
  word[count]='\0';
  return i;
}
///
/// read a single word from a file
/// words are delimited by delim-string chars
// or if set to NULL by whitespace
char read_word(FILE *infile, char *word, char *delim)
{
  char standard[5] = " \t\r\n";
  char *localdelim;
  int ch;
  long count = 0;
  if(delim==NULL)
    localdelim = standard;
  else
    localdelim = delim;
  ch = getc(infile);
  while (isspace (ch))
    ch = getc(infile);
  while(strchr(localdelim, ch)==NULL && ch!=EOF)
    {
      word[count++] = (char) ch;
      ch = getc(infile);
    }
    word[count]='\0';
    if(ch==EOF)
      return EOF;
    else
      {
	//	ungetc(ch,infile);
	//if(count!=0)
	//  return word[count-1];
	return (char) ch;
      }
    //return (char) '0';
}
///
/// read a single word from a file
/// words are delimited by withespace 
/*char read_word(FILE *infile, char *word)
{
    int ch;
    long count = 0;
    ch = getc(infile);
    while (isspace (ch))
        ch = getc(infile);
    while(strchr(" \t\n\r", ch)==NULL && ch!=EOF)
    {
        word[count++] = ch;
        ch = getc(infile);
    }
    word[count]='\0';
    if(ch==EOF)
      return EOF;
    else
      {
	ungetc(ch,infile);
	return 0;
      }
    return word[count-1];
    }*/

//========================================
// prepends a string to a file
// a blank will be inserted between the string and the
// text in the file
void unread_word(FILE *infile, char *word)
{
  char c = ' ';
  long count = (long) strlen(word);
  if(word[0] != '\n' && word[0] != ' ')
    ungetc(c, infile);
  while(count > 0)
    {
      count--;
      ungetc(word[count], infile);
    }
}

#ifdef STRANGEDEBUG
extern int errno;
#endif

#ifdef ZNZ
/// opens a zip file for reading or writing
void
znzopenfile (znzFile * fp, char *filename, char *mode, int use_compressed)
{
  //int trials = 0;
    znzFile of;
    char *file;
    char *p;
    file = (char *) mycalloc(LINESIZE,sizeof(char));
    if ((p = strpbrk (filename, CRLF)) != NULL)
        *p = '\0';
    remove_trailing_blanks(&filename);
    sprintf(file, "%-.*s", (int) LINESIZE, filename);
    of = znzopen (file, mode,use_compressed);
    if (znz_isnull(of))
      error("Could not open zipped bayesallfile");
    *fp = of;
    myfree(file);
}
#endif


/// opens a file for reading or writing
void
openfile (FILE ** fp, char *filename, char *mode, char *perm)
{
    int trials = 0;
    FILE *of = NULL;
    char *file;
    char *p;
#ifdef CAUTIOUS
    struct stat sb;
#endif
    file = (char *) mycalloc(LINESIZE,sizeof(char));
    if ((p = strpbrk (filename, CRLF)) != NULL)
        *p = '\0';
    remove_trailing_blanks(&filename);
    sprintf(file, "%-.*s",(int) LINESIZE, filename);
#ifdef CAUTIOUS
    if(cautious)
      {
	//	sb = (stat *) mycalloc(1,sizeof(stat));
	if(strchr(mode,'w') && (0 == stat(file, &sb)))
	  {
	    file[0] = '\0';
	    while (file[0] == '\0' && trials++ < 10)
	      {
		printf ("File %s exists, do you want to overwrite it? [YES/No]\n===>",filename);
		FGETS (file, LINESIZE, stdin);
		if(file[0]=='Y' || file[0] == 'y')
		  {
		    file[0] = '\0';
		    while (file[0] == '\0' && trials++ < 10)
		      {
			printf ("Please enter a new filename>");
			FGETS (file, LINESIZE, stdin);
		      }
		    
		  }
	      }
	  }
	//	myfree(sb);
      }
#endif
    while (trials++ < 10)
    {
      of = fopen (file, mode);
#ifdef STRANGEDEBUG
      fprintf(stdout,"original filename -->%s<--\n",filename);
      fprintf(stdout,"used filename     -->%s<--\n",file);
#endif
        if (of!=NULL)
            break;
        else
        {
	  //	  extension = strrchr(file,'.');
	  //	  if(extension!=NULL && !strncmp(extension,".txt",4))
#ifdef STRANGEDEBUG
	  fprintf(stdout,"mode              -->%s<--\n", mode);
	  fprintf(stdout,"errorcode         -->%i<--\n",errno);
#endif
            switch (*mode)
	      {
	      case 'r':
#ifdef MPI

                printf("%i> ",myID);
#endif

                printf ("Cannot read from file \"%s\"\n", file);
                file[0] = '\0';
                while (file[0] == '\0' && trials++ < 10)
                {
#ifdef MPI
                    printf("%i> ",myID);
#endif

                    printf ("Please enter a new filename for reading>");
                    FGETS (file, LINESIZE, stdin);
                }
                break;
            case 'w':
#ifdef MPI

                printf("%i> ",myID);
#endif

                printf ("Cannot write to file %s\n", file);
                file[0] = '\0';
                while (file[0] == '\0' && trials++ < 10)
                {
#ifdef MPI
                    printf("%i> ",myID);
#endif

                    printf ("Please enter a new filename for writing>");
                    FGETS (file, LINESIZE, stdin);
                }
                break;
            }
        }
    }
    if (trials >= 10)
    {
#ifdef MPI
        printf("%i> ",myID);
#endif

        printf ("You cannot find your file either, so I stop\n\n");
        exit (0);
    }
    *fp = of;
    if (perm != NULL)
        strcpy (perm, file);
    strcpy (filename, file);
    myfree(file);
}

void get_filename(char **store, char * value)
{
  long ls;
  long lv;
  remove_trailing_blanks(&value);
  ls = (long) strlen(*store)+1;
  lv = (long) strlen(value)+1;
      if (ls < lv)
	{
	  (*store) = (char *) myrealloc(*store,sizeof(char)*(size_t) lv);
	  memset(*store,0,sizeof(char)*(size_t) lv);
	}
      strcpy (*store, value);
}




/*=======================================================*/




/*--------------------------------
creates the length value in a node
*/
MYREAL
lengthof (node * p)
{
    if (p->type == 'm' || p->type == 'd')
        error ("A migration node was feed into lengthof");
    return fabs (p->tyme - showtop(crawlback (p))->tyme);
}    /* length */


/*------------------------------------------------
Find the next non-migration node starting
with the theNode, returns to backnode which is not 
a migration, does NOT return always a top-node!
multibranch node safe
*/
MYINLINE node *
crawlback (const node * theNode)
{
    node *tmp = theNode->back;

    while (tmp->type == 'm' || tmp->type == 'd')
    {
        tmp = tmp->next->back;
    }
    return tmp;
}

 
MYINLINE node *
showtop (node * theNode)
{
  node *q;
  q = theNode;
  while(!q->top)
    {
      if (q->next !=NULL)
	q=q->next;
      else
	error("showtop() crashed with q->next==NULL");
    }
  return q;
}

/* adjust the time in a node to time */
void
adjust_time (node * theNode, MYREAL tyme)
{
    switch (theNode->type)
    {
    case 'm':
    case 'd':
        theNode->tyme = theNode->next->tyme = tyme;
        break;
    case 'i':
      set_nodelet_time(theNode,tyme);
      //        theNode->tyme = theNode->next->tyme = theNode->next->next->tyme = tyme;
        break;
    case 'r':
      break;
    case 't':
      break;
    default:
      error ("Wrong node type");
      //break;
    }
}
/* adjust the time in any node to time */
void
adjust_time_all (node * theNode, MYREAL tyme)
{
  switch (theNode->type)
    {
    case 'm':
    case 'd':
      theNode->tyme = theNode->next->tyme = tyme;
      //set_dirty(theNode);
      break;
    case 'i':
      set_nodelet_time(theNode,tyme);
      //theNode->tyme = theNode->next->tyme = theNode->next->next->tyme = tyme;
      set_dirty(theNode);
      break;
    case 'r':
      warning("root adjusted time");
      set_nodelet_time(theNode,tyme);
      //theNode->tyme = theNode->next->tyme = theNode->next->next->tyme = tyme;
      set_dirty(theNode);
      break;
    case 't':
      theNode->tyme = tyme;
      break;
    default:
      error ("Wrong node type");
      //break;
    }
}

// inserts a migration node or 
// 2013: inserts a speciation node
void
insert_migr_node (world_fmt * world, node * up, node * down,
                  migr_table_fmt * migr_table, long *migr_table_counter)
{
  long i;
  node *theNode;
  if (!up->top)
    error ("up has to be a top-node");
  //xcode theNode = showtop (up)->back;
  if (*migr_table_counter > 0 && up->tyme > migr_table[0].time)
    {
      printf("%i> up->type=%c tyme=%f showtop(up)time=%f mig[0]tabletyme=%f\n", 
	     myID, up->type, up->tyme, showtop(up)->tyme, migr_table[0].time);
      error ("insert_migr_node: the first migration/speciation node has a wrong time for up");
    }
  if (migr_table[(*migr_table_counter) - 1].from != down->actualpop)
    {
      error ("this should never happen -> wrong choice of nodes\n");
    }
  if (((*migr_table_counter) > 0) && (migr_table[(*migr_table_counter) - 1].from != down->actualpop))
    {
      error ("problem caught in inser_migr_table");
    }
  theNode=up;
  long mm = (*migr_table_counter) - 1;
  if((down->tyme - migr_table[mm].time) < DBL_EPSILON)
    {
      warning("%i> adjusted time of migration because it was identical to internal node\n",myID);
      migr_table[mm].time -= DBL_EPSILON;
      while(mm>0 && migr_table[mm-1].time >= migr_table[mm].time)
	{
	  warning("%i> adjusted time of migration because it was identical to internal node\n",myID);
	  mm--;
	  migr_table[mm].time -= DBL_EPSILON;
	}
    }
  
  for (i = 0; i < (*migr_table_counter); i++)
    {
      theNode = add_migration(world, theNode, migr_table[i].event,
			      migr_table[i].from, 
			      migr_table[i].to,
			      migr_table[i].time - theNode->tyme);
    }
    if(down->tyme < migr_table[i-1].time)
      {
	error("Problem with migration addition in coalesce1p()");
      }
    down->back = theNode;
    theNode->back = down;
}

long offspring (node * mother, node ** myoffspring)
{
    node *m;
    long z=0;
    node * q;
    m = showtop (mother);

    if (m->type == 't')
    {
        error ("this is a tip, so there are no more child nodes\n");
    }
    else
    {
      z=0;
      q = m->next;
      while(!q->top)
	{
	  myoffspring[z] = crawlback (q);
	  q = q->next;
	  z++;
	}
    }
    return z;
}

void
children (node * mother, node ** brother, node ** sister)
{
    node *m;

    m = showtop (mother);

    if (m->type == 't')
    {
        error ("this is a tip, so there are no more child nodes\n");
    }
    else
    {
        (*brother) = crawlback (m->next);
        (*sister) = crawlback (m->next->next);
    }
}

// try to be same syntax as children(), but sister is an array with 1 or more entries
void siblings (node * mother, node ** brother, node ** sister)
{
  node * sisters[100];
  size_t z=0;
  node *m;
  //node *circle=NULL;
  node * sis=NULL;
  m = showtop (mother);
  
  if (m->type == 't')
    {
      error ("this is a tip, so there are no more child nodes\n");
    }
  else
    {
      m = m->next;
      z=0;
      while (!m->top)
	{
	  sis = crawlback(m->back);
	  m = m->next;
	  sisters[z++] = sis;
	}
      (*brother) = sisters[0];
      memcpy(sister,&sisters[1],sizeof(node*)*(z-1));
      //(*sister) = crawlback (m->next->next);
    }
}

#endif /*not PRIORTEST*/
/*       Uses Lanczos-type approximation to ln(gamma) for z > 0. */
/*       Reference: */
/*            Lanczos, C. 'A precision approximation of the gamma */
/*                    function', J. SIAM Numer. Anal., B, 1, 86-96, 1964. */
/*       Accuracy: About 14 significant digits except for small regions */
/*                 in the vicinity of 1 and 2. */
/*       Programmer: Alan Miller */
/*                   CSIRO Division of Mathematics & Statistics */
/*       Latest revision - 17 April 1988 */
/* translated and modified into C by Peter Beerli 1997 */
MYREAL
mylgamma (MYREAL z)
{
    MYREAL a[9] = { 0.9999999999995183, 676.5203681218835,
                    -1259.139216722289, 771.3234287757674, -176.6150291498386,
                    12.50734324009056, -0.1385710331296526, 9.934937113930748e-6,
                    1.659470187408462e-7
                  };
    MYREAL lnsqrt2pi = 0.9189385332046727;
    MYREAL result;
    long j;
    MYREAL tmp;
    if (z <= 0.)
    {
        return MYREAL_MAX;  /*this will kill the receiving calculation */
    }
    result = 0.;
    tmp = z + 7.;
    for (j = 9; j >= 2; --j)
    {
        result += a[j - 1] / tmp;
        tmp -= 1.;
    }
    result += a[0];
    result = log (result) + lnsqrt2pi - (z + 6.5) + (z - 0.5) * log (z + 6.5);
    return result;
}    /* lgamma */

/* ALGORITHM AS239  APPL. STATIST. (1988) VOL. 37, NO. 3
   Computation of the Incomplete Gamma Integral 
   Auxiliary functions required: lgamma() = logarithm of the gamma 
   function, and alnorm() = algorithm AS66 
   in Mathematica this is GammaRegularized[a,0,x] === Gamma[a,0,x]/Gamma[a]
   
This is the integral from 0 to x (the lower gamma function) divided by gamma(a)
the integral from x to infinity is in uppergamma_reg() and uppergamma()
 */
MYREAL
incompletegamma (MYREAL tx, MYREAL talpha)
{
  const  MYREAL eps = DBL_EPSILON;
  double gama, d_1, d_2, d_3;
  /*static */
  double  a, b, c, an, rn;
  /*static */
  double  pn1, pn2, pn3, pn4, pn5, pn6, arg;
  double x = (double) tx;
  double alpha = (double) talpha;

  gama = 0.;
  /*  Check that we have valid values for X and P */
  if (alpha <= 0. || x < 0.)
    error ("failed in imcompletegamma(): wrong alpha or x\n");
  if (fabs (x) < eps)
        return (MYREAL) gama;

    /*  Use a normal approximation if P > PLIMIT */
    if (alpha > 1e3)
    {
        pn1 =
            sqrt (alpha) * 3. * (pow (x / alpha, (1. / 3.)) + 1. / (alpha * 9.) -
                                 1.);
        gama = alnorm (pn1, FALSE);
        return (MYREAL) gama;
    }

   /*  If X is extremely large compared to P then set GAMMAD = 1 */
    if (x > 1e8)
    {
        gama = 1.;
        return (MYREAL) gama;
    }

    if (x <= 1. || x < alpha)
    {
        /*  Use Pearson's series expansion. */
        /*  (Note that P is not large enough to force overflow in lgamma()). */
        arg = alpha * LOG (x) - x - LGAMMA (alpha + 1.);
        c = 1.;
        gama = 1.;
        a = alpha;
        while (c > 1e-14)
        {
            a += 1.;
            c = c * x / a;
            gama += c;
        }
        arg += LOG (gama);
        gama = 0.;
        if (arg >= -88.)
        {
            gama = EXP (arg);
        }

    }
    else
    {
        /*  Use a continued fraction expansion */
        arg = alpha * LOG (x) - x - LGAMMA (alpha);
        a = 1. - alpha;
        b = a + x + 1.;
        c = 0.;
        pn1 = 1.;
        pn2 = x;
        pn3 = x + 1.;
        pn4 = x * b;
        gama = pn3 / pn4;
        for (;;)
        {
            a += 1.;
            b += 2.;
            c += 1.;
            an = a * c;
            pn5 = b * pn3 - an * pn1;
            pn6 = b * pn4 - an * pn2;
            if (fabs (pn6) > 0.)
            {
                rn = pn5 / pn6;
                /* Computing MIN */
                d_2 = 1e-14;
		d_3 = rn * 1e-14;
		d_1 = gama - rn;
                if (fabs (d_1) <= MIN (d_2, d_3))
                {
                    arg += LOG (gama);
                    gama = 1.;
                    if (arg >= -88.)
                    {
                        gama = 1. - EXP (arg);
                    }
                    return (MYREAL) gama;
                }
                gama = rn;
            }
            pn1 = pn3;
            pn2 = pn4;
            pn3 = pn5;
            pn4 = pn6;
            if (fabs (pn5) >= 1e37)
            {
                /*  Re-scale terms in continued fraction if terms are large */
                pn1 /= 1e37;
                pn2 /= 1e37;
                pn3 /= 1e37;
                pn4 /= 1e37;
            }
        }
    }
    return (MYREAL) gama;
}    /* incompletegamma() */

/* returns the log of incomplete gamma, copied from above and edited so that we do only use this when needed*/
MYREAL
logincompletegamma (MYREAL tx, MYREAL talpha)
{
  const  MYREAL eps = DBL_EPSILON;
  double gama, d_1, d_2, d_3;
  /*static */
  double  a, b, c, an, rn;
  /*static */
  double  pn1, pn2, pn3, pn4, pn5, pn6, arg;
  double x = (double) tx;
  double alpha = (double) talpha;

  //gama = 0.;
  /*  Check that we have valid values for X and P */
  if (alpha <= 0. || x < 0.)
    error ("failed in imcompletegamma(): wrong alpha or x\n");
  if (fabs (x) < eps)
    return (MYREAL) -HUGE;

    /*  Use a normal approximation if P > PLIMIT */
    if (alpha > 1e3)
    {
        pn1 =
            sqrt (alpha) * 3. * (pow (x / alpha, (1. / 3.)) + 1. / (alpha * 9.) -
                                 1.);
        gama = alnorm (pn1, FALSE);
	if (gama>0.)
	  return (MYREAL) log(gama);
	else
	  return (MYREAL) -HUGE;
    }

   /*  If X is extremely large compared to P then set GAMMAD = 1 */
    if (x > 1e8)
    {
      //gama = 1.;
        return (MYREAL) 0.0;
    }

    if (x <= 1. || x < alpha)
    {
        /*  Use Pearson's series expansion. */
        /*  (Note that P is not large enough to force overflow in lgamma()). */
        arg = alpha * LOG (x) - x - LGAMMA (alpha + 1.);
        c = 1.;
        gama = 1.;
        a = alpha;
        while (c > 1e-14)
        {
            a += 1.;
            c = c * x / a;
            gama += c;
        }
        arg += LOG (gama);
	gama = arg;
    }
    else
    {
        /*  Use a continued fraction expansion */
        arg = alpha * LOG (x) - x - LGAMMA (alpha);
        a = 1. - alpha;
        b = a + x + 1.;
        c = 0.;
        pn1 = 1.;
        pn2 = x;
        pn3 = x + 1.;
        pn4 = x * b;
        gama = pn3 / pn4;
        for (;;)
        {
            a += 1.;
            b += 2.;
            c += 1.;
            an = a * c;
            pn5 = b * pn3 - an * pn1;
            pn6 = b * pn4 - an * pn2;
            if (fabs (pn6) > 0.)
            {
                rn = pn5 / pn6;
                /* Computing MIN */
                d_2 = 1e-14;
		d_3 = rn * 1e-14;
		d_1 = gama - rn;
                if (fabs (d_1) <= MIN (d_2, d_3))
                {
                    arg += LOG (gama);
		    gama = LOG(1. - EXP (arg));
                    return (MYREAL) gama;
                }
                gama = rn;
            }
            pn1 = pn3;
            pn2 = pn4;
            pn3 = pn5;
            pn4 = pn6;
            if (fabs (pn5) >= 1e37)
            {
                /*  Re-scale terms in continued fraction if terms are large */
                pn1 /= 1e37;
                pn2 /= 1e37;
                pn3 /= 1e37;
                pn4 /= 1e37;
            }
        }
    }
    return (MYREAL) gama;
}    /* incompletegamma() */


MYREAL
uppergamma_reg (MYREAL tx, MYREAL talpha)
{
  double val;
  val = 1.0 - incompletegamma(tx, talpha);
  return val;
}

MYREAL
uppergamma (MYREAL tx, MYREAL talpha)
{
  double val;
  val = tgamma(talpha) * (1.0 - incompletegamma(tx, talpha));
  return val;
}

/* calculation is replaced by the correct function in
   polygamma.c (which is a translation of a fortran program by amos
 
   driver for the polygamma calculation */
MYREAL
polygamma (long n, MYREAL z)
{
    MYREAL ans;
    long nz, ierr;
    dpsifn (&z, &n, 1, 1, &ans, &nz, &ierr);
    if (n == 0)
        return -ans;
    else
        return ans;
}

/*-------------------------------------------------------*/
/* nrcheck subroutine (used in damped newton raphson proc */
/* syntax: nrcheck(matrix,inversematrix,ncols=nrows,returnval1,returnval2) */
/* mai 95 PB                                             */
boolean
nrcheck (MYREAL **m, MYREAL **tm, MYREAL *v, long nrows, MYREAL *r1,
         MYREAL *r2, boolean do_newton)
{
    long i, j, k;
    MYREAL *tmp, *tmp2, tmp3 = 0.0, tmp4 = 0.0;
    tmp = (MYREAL *) mycalloc (nrows, sizeof (MYREAL));
    tmp2 = (MYREAL *) mycalloc (nrows, sizeof (MYREAL));
    /*first evaluate r1 */
    (*r1) = (*r2) = 0.0;
    for (i = 0; i < nrows; i++)
    {
        (*r1) += v[i] * v[i];
    }
    /*                                       T    */
    for (j = 0; j < nrows; j++)
    {    /* g . G */
        for (k = 0; k < nrows; k++)
        {
            tmp[j] += v[k] * m[j][k];
            tmp2[j] += v[k] * tm[j][k];
        }
    }
    /*                                       T        */
    for (i = 0; i < nrows; i++)
    {    /* g . G . g */
        (*r2) += tmp[i] * v[i];
        tmp3 += tmp2[i] * v[i];
    }
    tmp4 = LOG (fabs ((*r1)));
    tmp4 = tmp4 + tmp4 - LOG (fabs ((*r2)));
    tmp4 = ((*r2) < 0 ? -1 : 1) * EXP (tmp4);
    myfree(tmp);
    if (do_newton && (tmp3 > (tmp4 > 0 ? tmp4 : 0)))
    {
      memcpy (v, tmp2, sizeof (MYREAL) * (size_t) nrows);
        myfree(tmp2);
        return TRUE;
    }
    myfree(tmp2);
    return FALSE;
}


/*-------------------------------------------------------*/
/* Matrix inversion subroutine                           */
/* The passed matrix will be replaced by its inverse!!!!! */
/* Gauss-Jordan reduction -- invert matrix a in place,   */
/* overwriting previous contents of a.  On exit, matrix a */
/* contains the inverse.                                 */
void
invert_matrix (MYREAL **a, long nsize)
{
    long i, j;
    long *indeks;
    MYREAL *column, **result;
    indeks = (long *) mymalloc (sizeof (long) * (size_t) nsize);
    column = (MYREAL *) mymalloc (sizeof (MYREAL) * (size_t) nsize);
    result = (MYREAL **) mymalloc (sizeof (MYREAL *) * (size_t) nsize);
    for (i = 0; i < nsize; i++)
    {
        result[i] = (MYREAL *) mymalloc (sizeof (MYREAL) * (size_t) nsize);
    }
    lu_decomp (a, indeks, nsize);
    for (j = 0; j < nsize; j++)
    {
        memset (column, 0, sizeof (MYREAL) * (size_t) nsize);
        column[j] = 1.0;
        lu_substitution (a, indeks, column, nsize);
        for (i = 0; i < nsize; i++)
            result[i][j] = column[i];
    }
    for (i = 0; i < nsize; i++)
    {
        memcpy (a[i], result[i], sizeof (MYREAL) * (size_t) nsize);
        myfree(result[i]);
    }
    myfree(result);
    myfree(column);
    myfree(indeks);
}

/*=======================================================*/

/*-------------------------------------------------------*/
/* LU decomposition                                      */
/* after Dahlquist et al. 1974 and Press et al. 1988     */
/* the method's uses Crout's procedure and the pivoting  */
/* described in Press et al.                             */
/* Syntax: lu_decomp(matrix, indeks, nrows)               */
/* matrix will be destroyed and filled with the two      */
/* triangular matrices, indeks is the index vector for the */
/* pivoting and the row change in case of 0 pivot values */
/* nrows is the number of rows and columns in matrix     */
/* april 95 PB                                           */
void
lu_decomp (MYREAL **m, long *indeks, long nrows)
{
    long i, j, k, p, kmax = -1;
    MYREAL *max_row_vals, big, summ, pivot, bigt;
    max_row_vals = (MYREAL *) mycalloc (1, sizeof (MYREAL) * (size_t) nrows);
    for (i = 0; i < nrows; i++)
    {
        big = 0.0;
        for (j = 0; j < nrows; j++)
        {
            if ((bigt = fabs (m[i][j])) > big)
                big = bigt;
        }
        max_row_vals[i] = 1.0 / big;
        if (big == 0.0)
        {
            error ("Singular matrix detected in lu_decomp\n");
        }
    }
    for (i = 0; i < nrows; i++)
    {
        for (k = 0; k < i; k++)
        {   /* upper half of matrix */
            summ = m[k][i];
            for (p = 0; p < k; p++)
                summ -= m[k][p] * m[p][i];
            m[k][i] = summ;
        }
        big = 0.0;
        for (k = i; k < nrows; k++)
        {   /* lower half of matrix */
            summ = m[k][i];
            for (p = 0; p < i; p++)
                summ -= m[k][p] * m[p][i];
            m[k][i] = summ;
            pivot = fabs (summ) /**max_row_vals[k]*/;
            /*  printf(stdout,"i=%li,pivot=%f,big=%f\n",i,pivot,big); */
            if (pivot >= big)
            {
                big = pivot;
                kmax = k;
            }
        }
        if (i != kmax)
        {
            for (p = 0; p < nrows; p++)
            {
                pivot = m[kmax][p];
                m[kmax][p] = m[i][p];
                m[i][p] = pivot;
            }
            max_row_vals[kmax] = max_row_vals[i];
        }
        indeks[i] = kmax;
        if (m[i][i] == 0.0)
            m[i][i] = SMALL_VALUE;
        if (i != nrows - 1)
        {
            pivot = 1. / m[i][i];
            for (k = i + 1; k < nrows; k++)
                m[k][i] *= pivot;
        }
    }
    myfree(max_row_vals);
}    /* end of lu_decomp */

/*-------------------------------------------------------*/
/* LU substitution                                       */
/* after Dahlquist et al. 1974 and Press et al. 1988     */
/* needs first the evaluation LU decomposition           */
/* Syntax: lu_substition(matrix, indeks, vector, nrows)   */
/* matrix = LU decomposed matrix, indeks = order of matrix */
/* vector = value vector, nrows = number of rows/columns */
/* april 95 PB                                           */
void
lu_substitution (MYREAL **m, long *indeks, MYREAL *v, long nrows)
{
    long i, j;
    MYREAL summ;
    for (i = 0; i < nrows; i++)
    {
        summ = v[indeks[i]];
        v[indeks[i]] = v[i];
        for (j = 0; j < i; j++)
            summ -= m[i][j] * v[j];
        v[i] = summ;
    }
    for (i = nrows - 1; i >= 0; i--)
    {
        summ = v[i];
        for (j = i + 1; j < nrows; j++)
            summ -= m[i][j] * v[j];
        v[i] = summ / m[i][i];
    }
}


/* Algorithm AS66 Applied Statistics (1973) vol22 no.3
   Evaluates the tail area of the standardised normal curve
   from x to infinity if upper is .true. or
   from minus infinity to x if upper is .false. */
MYREAL
alnorm (MYREAL x, int up)
{
    /* Initialized data */
    /* *** machine dependent constants ????????????? */
    /*static */ MYREAL zero = 0.;
    /*static */
    MYREAL a1 = 5.75885480458;
    /*static */
    MYREAL a2 = 2.62433121679;
    /*static */
    MYREAL a3 = 5.92885724438;
    /*static */
    MYREAL b1 = -29.8213557807;
    /*static */
    MYREAL b2 = 48.6959930692;
    /*static */
    MYREAL c1 = -3.8052e-8;
    /*static */
    MYREAL c2 = 3.98064794e-4;
    /*static */
    MYREAL c3 = -.151679116635;
    /*static */
    MYREAL c4 = 4.8385912808;
    /*static */
    MYREAL c5 = .742380924027;
    /*static */
    MYREAL one = 1.;
    /*static */
    MYREAL c6 = 3.99019417011;
    /*static */
    MYREAL d1 = 1.00000615302;
    /*static */
    MYREAL d2 = 1.98615381364;
    /*static */
    MYREAL d3 = 5.29330324926;
    /*static */
    MYREAL d4 = -15.1508972451;
    /*static */
    MYREAL d5 = 30.789933034;
    /*static */
    MYREAL half = .5;
    /*static */
    MYREAL ltone = 7.;
    /*static */
    MYREAL utzero = 18.66;
    /*static */
    MYREAL con = 1.28;
    /*static */
    MYREAL p = .398942280444;
    /*static */
    MYREAL q = .39990348504;
    /*static */
    MYREAL r = .398942280385;

    /*static */
    MYREAL y, result;

    if (x < zero)
    {
        up = !up;
        x = -x;
    }
    if (x <= ltone || (up && x <= utzero))
    {
        y = half * x * x;
        if (x > con)
        {
            result =
                r * EXP (-y) / (x + c1 +
                                d1 / (x + c2 +
                                      d2 / (x + c3 +
                                            d3 / (x + c4 +
                                                  d4 / (x + c5 +
                                                        d5 / (x + c6))))));
            return ((!up) ? one - result : result);
        }
        result =
            half - x * (p - q * y / (y + a1 + b1 / (y + a2 + b2 / (y + a3))));
        return ((!up) ? one - result : result);
    }
    else
    {
        return ((!up) ? 1.0 : 0.);
    }
    /*fake */ //return -99;
}    /* alnorm */

/* dpsifn.c -- translated by f2c (version 19950808).
   and hand-patched by Peter Beerli Seattle, 1996
   SUBROUTINE DPSIFN (X, N, KODE, M, ANS, NZ, IERR)
 
   C***BEGIN PROLOGUE  DPSIFN
   C***PURPOSE  Compute derivatives of the Psi function.
   C***LIBRARY   SLATEC
   C***CATEGORY  C7C
   C***TYPE      MYREAL PRECISION (PSIFN-S, DPSIFN-D)
   C***KEYWORDS  DERIVATIVES OF THE GAMMA FUNCTION, POLYGAMMA FUNCTION,
   C             PSI FUNCTION
   C***AUTHOR  Amos, D. E., (SNLA)
   C***DESCRIPTION
   C
   C         The following definitions are used in DPSIFN:
   C
   C      Definition 1
   C         PSI(X) = d/dx (ln(GAMMA(X)), the first derivative of
   C                  the log GAMMA function.
   C      Definition 2
   C                     K   K
   C         PSI(K,X) = d /dx (PSI(X)), the K-th derivative of PSI(X).
   C   ___________________________________________________________________
   C      DPSIFN computes a sequence of SCALED derivatives of
   C      the PSI function; i.e. for fixed X and M it computes
   C      the M-member sequence
   C
   C                    ((-1)**(K+1)/GAMMA(K+1))*PSI(K,X)
   C                       for K = N,...,N+M-1
   C
   C      where PSI(K,X) is as defined above.   For KODE=1, DPSIFN returns
   C      the scaled derivatives as described.  KODE=2 is operative only
   C      when K=0 and in that case DPSIFN returns -PSI(X) + LN(X).  That
   C      is, the logarithmic behavior for large X is removed when KODE=2
   C      and K=0.  When sums or differences of PSI functions are computed
   C      the logarithmic terms can be combined analytically and computed
   C      separately to help retain significant digits.
   C
   C         Note that CALL DPSIFN(X,0,1,1,ANS) results in
   C                   ANS = -PSI(X)
   C
   C     Input      X is MYREAL PRECISION
   C           X      - Argument, X .gt. 0.0D0
   C           N      - First member of the sequence, 0 .le. N .le. 100
   C                    N=0 gives ANS(1) = -PSI(X)       for KODE=1
   C                                       -PSI(X)+LN(X) for KODE=2
   C           KODE   - Selection parameter
   C                    KODE=1 returns scaled derivatives of the PSI
   C                    function.
   C                    KODE=2 returns scaled derivatives of the PSI
   C                    function EXCEPT when N=0. In this case,
   C                    ANS(1) = -PSI(X) + LN(X) is returned.
   C           M      - Number of members of the sequence, M.ge.1
   C
   C    Output     ANS is MYREAL PRECISION
   C           ANS    - A vector of length at least M whose first M
   C                    components contain the sequence of derivatives
   C                    scaled according to KODE.
   C           NZ     - Underflow flag
   C                    NZ.eq.0, A normal return
   C                    NZ.ne.0, Underflow, last NZ components of ANS are
   C                             set to zero, ANS(M-K+1)=0.0, K=1,...,NZ
   C           IERR   - Error flag
   C                    IERR=0, A normal return, computation completed
   C                    IERR=1, Input error,     no computation
   C                    IERR=2, Overflow,        X too small or N+M-1 too
   C                            large or both
   C                    IERR=3, Error,           N too large. Dimensioned
   C                            array TRMR(NMAX) is not large enough for N
   C
   C         The nominal computational accuracy is the maximum of unit
   C         roundoff (=D1MACH(4)) and 1.0D-18 since critical constants
   C         are given to only 18 digits.
   C
   C         PSIFN is the single precision version of DPSIFN.
   C
   C *Long Description:
   C
   C         The basic method of evaluation is the asymptotic expansion
   C         for large X.ge.XMIN followed by backward recursion on a two
   C         term recursion relation
   C
   C                  W(X+1) + X**(-N-1) = W(X).
   C
   C         This is supplemented by a series
   C
   C                  SUM( (X+K)**(-N-1) , K=0,1,2,... )
   C
   C         which converges rapidly for large N. Both XMIN and the
   C         number of terms of the series are calculated from the unit
   C         roundoff of the machine environment.
   C
   C***REFERENCES  Handbook of Mathematical Functions, National Bureau
   C                 of Standards Applied Mathematics Series 55, edited
   C                 by M. Abramowitz and I. A. Stegun, equations 6.3.5,
   C                 6.3.18, 6.4.6, 6.4.9 and 6.4.10, pp.258-260, 1964.
   C               D. E. Amos, A portable Fortran subroutine for
   C                 derivatives of the Psi function, Algorithm 610, ACM
   C                 Transactions on Mathematical Software 9, 4 (1983),
   C                 pp. 494-502.
   C***ROUTINES CALLED  D1MACH, I1MACH
   C***REVISION HISTORY  (YYMMDD)
   C   820601  DATE WRITTEN
   C   890531  Changed all specific intrinsics to generic.  (WRB)
   C   890911  Removed unnecessary intrinsics.  (WRB)
   C   891006  Cosmetic changes to prologue.  (WRB)
   C   891006  REVISION DATE from Version 3.2
   C   891214  Prologue converted to Version 4.0 format.  (BAB)
   C   920501  Reformatted the REFERENCES section.  (WRB)
   C***END PROLOGUE  DPSIFN
 
 
 */

static long fifteen = 15;
static long sixteen = 16;
static long five = 5;
static long four = 4;
static long fourteen = 14;

MYREAL
d1mach (long i)
{
  //#ifdef MYREAL == float
  //const  MYREAL eps = FLT_EPSILON;
  //const  MYREAL numbermax = FLT_MAX;
  //const  MYREAL numbermin = FLT_MIN;
  //#else
  const  MYREAL eps = DBL_EPSILON;
  const  MYREAL numbermax = MYREAL_MAX;
  const  MYREAL numbermin = DBL_MIN;
  //#endif

    switch (i)
    {
    case 1:
        return numbermin;
    case 2:
        return numbermax;
    case 3:
      return eps / (double) FLT_RADIX;
    case 4:
        return eps;
    case 5:
        return log10 ((MYREAL) FLT_RADIX);
    }
    usererror ("invalid argument: d1mach(%ld)\n", i);
    //return 0;   /* for compilers that complain of missing return values */
}

long
i1mach (long i)
{
    switch (i)
    {
    case 1:
        return 5;   /* standard input */
    case 2:
        return 6;   /* standard output */
    case 3:
        return 7;   /* standard punch */
    case 4:
        return 0;   /* standard error */
    case 5:
        return 32;  /* bits per integer */
    case 6:
        return 1;   /* Fortran 77 value */
    case 7:
        return 2;   /* base for integers */
    case 8:
        return 31;  /* digits of integer base */
    case 9:
        return LONG_MAX;
    case 10:
        return FLT_RADIX;
    case 11:
        return FLT_MANT_DIG;
    case 12:
        return FLT_MIN_EXP;
    case 13:
        return FLT_MAX_EXP;
    case 14:
        return DBL_MANT_DIG;
    case 15:
        return DBL_MIN_EXP;
    case 16:
        return DBL_MAX_EXP;
    }
    usererror ("invalid argument: i1mach(%ld)\n", i);
    //return 0;   /* for compilers that complain of missing return values */
}

int
dpsifn (MYREAL *x, long *n, long kode, long m, MYREAL *ans, long *nz,
        long *ierr)
{
    /* Initialized data */

    /*static */ long nmax = 100;
    /*static */
    MYREAL b[22] = { 1., -.5, .166666666666666667,
                     -.0333333333333333333, .0238095238095238095, -.0333333333333333333,
                     .0757575757575757576, -.253113553113553114, 1.16666666666666667,
                     -7.09215686274509804, 54.9711779448621554, -529.124242424242424,
                     6192.1231884057971, -86580.2531135531136, 1425517.16666666667,
                     -27298231.067816092, 601580873.900642368, -15116315767.0921569,
                     429614643061.166667, -13711655205088.3328, 488332318973593.167,
                     -19296579341940068.1
                   };

    /* System generated locals */
    long i1, i2;
    MYREAL d1, d2;


    /* Local variables */
    /*static */
    MYREAL elim, xinc, xmin, tols, xdmy, yint, trmr[100], rxsq;
    /*static */
    long i__, j, k;
    /*static */
    MYREAL s, t, slope, xdmln, wdtol;
    /*static */
    MYREAL t1, t2;
    /*static */
    long fn;
    /*static */
    MYREAL ta;
    /*static */
    long mm, nn, np;
    /*static */
    MYREAL fx, tk;
    /*static */
    long mx, nx;
    /*static */
    MYREAL xm, tt, xq, den, arg, fln, r1m4, r1m5, eps, rln, tol,
    xln, trm[22], tss, tst;
    int i;
	for(i=0;i<22;i++)
        trm[i]=0.0;
    for(i=0;i<100;i++)
        trmr[i]=0.0;
    
    /* Parameter adjustments */
    --ans;

    /* Function Body */
    /* ----------------------------------------------------------------------- */
    /*             BERNOULLI NUMBERS */
    /* ----------------------------------------------------------------------- */

    /* ***FIRST EXECUTABLE STATEMENT  DPSIFN */
    *ierr = 0;
    *nz = 0;
    if (*x <= 0.)
    {
        *ierr = 1;
    }
    if (*n < 0)
    {
        *ierr = 1;
    }
    if (kode < 1 || kode > 2)
    {
        *ierr = 1;
    }
    if (m < 1)
    {
        *ierr = 1;
    }
    if (*ierr != 0)
    {
        return 0;
    }
    mm = m;
    /* Computing MIN */
    //xcode i1 = -fifteen;
    //xcode i2 = sixteen;
    nx = MIN (-i1mach (fifteen), i1mach (sixteen));
    r1m5 = d1mach (five);
    r1m4 = d1mach (four) * .5;
    wdtol = MAX (r1m4, 5e-19);
    /* ----------------------------------------------------------------------- */
    /*     ELIM = APPROXIMATE EXPONENTIAL OVER AND UNDERFLOW LIMIT */
    /* ----------------------------------------------------------------------- */
    elim = (nx * r1m5 - 3.) * 2.302;
    xln = LOG (*x);
L41:
    nn = *n + mm - 1;
    fn = nn;
    t = (fn + 1) * xln;
    /* ----------------------------------------------------------------------- */
    /*     OVERFLOW AND UNDERFLOW TEST FOR SMALL AND LARGE X */
    /* ----------------------------------------------------------------------- */
    if (fabs (t) > elim)
    {
        goto L290;
    }
    if (*x < wdtol)
    {
        goto L260;
    }
    /* ----------------------------------------------------------------------- */
    /*     COMPUTE XMIN AND THE NUMBER OF TERMS OF THE SERIES, FLN+1 */
    /* ----------------------------------------------------------------------- */
    rln = r1m5 * i1mach (fourteen);
    rln = MIN (rln, 18.06);
    fln = MAX (rln, 3.) - 3.;
    yint = fln * .4 + 3.5;
    slope = fln * (fln * 6.038e-4 + .008677) + .21;
    xm = yint + slope * fn;
    mx = (long) xm + 1;
    xmin = (MYREAL) mx;
    if (*n == 0)
    {
        goto L50;
    }
    xm = rln * -2.302 - MIN (0., xln);
    arg = xm / *n;
    arg = MIN (0., arg);
    eps = EXP (arg);
    xm = 1. - eps;
    if (fabs (arg) < .001)
    {
        xm = -arg;
    }
    fln = *x * xm / eps;
    xm = xmin - *x;
    if (xm > 7. && fln < 15.)
    {
        goto L200;
    }
L50:
    xdmy = *x;
    xdmln = xln;
    xinc = 0.;
    if (*x >= xmin)
    {
        goto L60;
    }
    nx = (long) (*x);
    xinc = xmin - nx;
    xdmy = *x + xinc;
    xdmln = LOG (xdmy);
L60:
    /* ----------------------------------------------------------------------- */
    /*     GENERATE W(N+MM-1,X) BY THE ASYMPTOTIC EXPANSION */
    /* ----------------------------------------------------------------------- */
    t = fn * xdmln;
    t1 = xdmln + xdmln;
    t2 = t + xdmln;
    /* Computing MAX */
    d1 = fabs (t);
    d2 = fabs (t1);
    d1 = MAX (d1, d2);
    d2 = fabs (t2);
    tk = MAX (d1, d2);
    if (tk > elim)
    {
        goto L380;
    }
    tss = EXP (-t);
    tt = .5 / xdmy;
    t1 = tt;
    tst = wdtol * tt;
    if (nn != 0)
    {
        t1 = tt + 1. / fn;
    }
    rxsq = 1. / (xdmy * xdmy);
    ta = rxsq * .5;
    t = (fn + 1) * ta;
    s = t * b[2];
    if (fabs (s) < tst)
    {
        goto L80;
    }
    tk = 2.;
    for (k = 4; k <= 22; ++k)
    {
        t = t * ((tk + fn + 1) / (tk + 1.)) * ((tk + fn) / (tk + 2.)) * rxsq;
        trm[k - 1] = t * b[k - 1];
	d1 = trm[k - 1];
        if ( fabs (d1) < tst)
        {
            goto L80;
        }
        s += trm[k - 1];
        tk += 2.;
        /* L70: */
    }
L80:
    s = (s + t1) * tss;
    if (xinc == 0.)
    {
        goto L100;
    }
    /* ----------------------------------------------------------------------- */
    /*     BACKWARD RECUR FROM XDMY TO X */
    /* ----------------------------------------------------------------------- */
    nx = (long) xinc;
    np = nn + 1;
    if (nx > nmax)
    {
        goto L390;
    }
    if (nn == 0)
    {
        goto L160;
    }
    xm = xinc - 1.;
    fx = *x + xm;
    /* ----------------------------------------------------------------------- */
    /*     THIS LOOP SHOULD NOT BE CHANGED. FX IS ACCURATE WHEN X IS SMALL */
    /* ----------------------------------------------------------------------- */
    i1 = nx;
    for (i__ = 1; i__ <= i1; ++i__)
    {
        i2 = -np;
        trmr[i__ - 1] = pow (fx, (MYREAL) i2);
        s += trmr[i__ - 1];
        xm += -1.;
        fx = *x + xm;
        /* L90: */
    }
L100:
    ans[mm] = s;
    if (fn == 0)
    {
        goto L180;
    }
    /* ----------------------------------------------------------------------- */
    /*     GENERATE LOWER DERIVATIVES, J.LT.N+MM-1 */
    /* ----------------------------------------------------------------------- */
    if (mm == 1)
    {
        return 0;
    }
    i1 = mm;
    for (j = 2; j <= i1; ++j)
    {
        --fn;
        tss *= xdmy;
        t1 = tt;
        if (fn != 0)
        {
            t1 = tt + 1. / fn;
        }
        t = (fn + 1) * ta;
        s = t * b[2];
        if (fabs (s) < tst)
        {
            goto L120;
        }
        tk = (MYREAL) (fn + 4);
        for (k = 4; k <= 22; ++k)
        {
            trm[k - 1] = trm[k - 1] * (fn + 1) / tk;
	    d1 = trm[k - 1];
            if (fabs (d1) < tst)
            {
                goto L120;
            }
            s += trm[k - 1];
            tk += 2.;
            /* L110: */
        }
L120:
        s = (s + t1) * tss;
        if (xinc == 0.)
        {
            goto L140;
        }
        if (fn == 0)
        {
            goto L160;
        }
        xm = xinc - 1.;
        fx = *x + xm;
        i2 = nx;
        for (i__ = 1; i__ <= i2; ++i__)
        {
            trmr[i__ - 1] *= fx;
            s += trmr[i__ - 1];
            xm += -1.;
            fx = *x + xm;
            /* L130: */
        }
L140:
        mx = mm - j + 1;
        ans[mx] = s;
        if (fn == 0)
        {
            goto L180;
        }
        /* L150: */
    }
    return 0;
    /* ----------------------------------------------------------------------- */
    /*     RECURSION FOR N = 0 */
    /* ----------------------------------------------------------------------- */
L160:
    i1 = nx;
    for (i__ = 1; i__ <= i1; ++i__)
    {
        s += 1. / (*x + nx - i__);
        /* L170: */
    }
L180:
    if (kode == 2)
    {
        goto L190;
    }
    ans[1] = s - xdmln;
    return 0;
L190:
    if (fabs(xdmy - *x) <= (double) FLT_EPSILON)
    {
        return 0;
    }
    xq = xdmy / *x;
    ans[1] = s - LOG (xq);
    return 0;
    /* ----------------------------------------------------------------------- */
    /*     COMPUTE BY SERIES (X+K)**(-(N+1)) , K=0,1,2,... */
    /* ----------------------------------------------------------------------- */
L200:
    nn = (long) fln + 1;
    np = *n + 1;
    t1 = (*n + 1) * xln;
    t = EXP (-t1);
    s = t;
    den = *x;
    i1 = nn;
    for (i__ = 1; i__ <= i1; ++i__)
    {
        den += 1.;
        i2 = -np;
        trm[i__ - 1] = pow (den, (MYREAL) i2);
        s += trm[i__ - 1];
        /* L210: */
    }
    ans[1] = s;
    if (*n != 0)
    {
        goto L220;
    }
    if (kode == 2)
    {
        ans[1] = s + xln;
    }
L220:
    if (mm == 1)
    {
        return 0;
    }
    /* ----------------------------------------------------------------------- */
    /*     GENERATE HIGHER DERIVATIVES, J.GT.N */
    /* ----------------------------------------------------------------------- */
    tol = wdtol / 5.;
    i1 = mm;
    for (j = 2; j <= i1; ++j)
    {
        t /= *x;
        s = t;
        tols = t * tol;
        den = *x;
        i2 = nn;
        for (i__ = 1; i__ <= i2; ++i__)
        {
            den += 1.;
            trm[i__ - 1] /= den;
            s += trm[i__ - 1];
            if (trm[i__ - 1] < tols)
            {
                goto L240;
            }
            /* L230: */
        }
L240:
        ans[j] = s;
        /* L250: */
    }
    return 0;
    /* ----------------------------------------------------------------------- */
    /*     SMALL X.LT.UNIT ROUND OFF */
    /* ----------------------------------------------------------------------- */
L260:
    i1 = -(*n) - 1;
    ans[1] = pow (*x, (MYREAL) i1);
    if (mm == 1)
    {
        goto L280;
    }
    k = 1;
    i1 = mm;
    for (i__ = 2; i__ <= i1; ++i__)
    {
        ans[k + 1] = ans[k] / *x;
        ++k;
        /* L270: */
    }
L280:
    if (*n != 0)
    {
        return 0;
    }
    if (kode == 2)
    {
        ans[1] += xln;
    }
    return 0;
L290:
    if (t > 0.)
    {
        goto L380;
    }
    *nz = 0;
    *ierr = 2;
    return 0;
L380:
    ++(*nz);
    ans[mm] = 0.;
    --mm;
    if (mm == 0)
    {
        return 0;
    }
    goto L41;
L390:
    *nz = 0;
    *ierr = 3;
    return 0;
}    /* dpsifn_ */




MYREAL
rannor (MYREAL mean, MYREAL sd)
{
    MYREAL r1, r2;
    r1 = RANDUM ();
    r2 = RANDUM ();
    return sd * sqrt (-2. * LOG (r1)) * cos (TWOPI * r2) + mean;
}


char
lowercase (int c)
{
    return (char) tolower (c);
}

char
uppercase (int c)
{
    return (char) toupper (c);
}

void upper(char *from, char **to)
{
    long i=0;
    while(from[i] != '\0')
    {
       (*to)[i] = from[i];
        ++i;
    }
    (*to)[++i] = '\0';
}

MYREAL
find_chi (long df, MYREAL prob)
{
    double a, b, m;
    double xb = 200.0;
    double xa = 0.0;
    double xm = 5.;
    double dprob = (double) prob;
    a = probchi (df, xa);
    m = probchi (df, xm);
    b = probchi (df, xb);
    while (fabs (m - dprob) > EPSILON)
    {
        if (m < dprob)
        {
            b = m;
            xb = xm;
        }
        else
        {
            a = m;
            xa = xm;
        }
        xm = (-(b * xa) + prob * xa + a * xb - dprob * xb) / (a - b); //(xa + xb)/2.;

        m = probchi (df, xm);
    }
    return (MYREAL) xm;
}


MYREAL
probchi (long df, double chi)
{
  const  MYREAL eps = DBL_EPSILON;
    double prob;
    double v = ((MYREAL) df) / 2.;
    
    if (chi > eps && v > eps)
    {
        //lg = EXP (LGAMMA (v));
        prob = 1. - incompletegamma (chi / 2., v);
    }
    else
        prob = 1.0;
    //  printf("prob=%f v=%f chi=%f lg(v/2)=%f  ig(chi/2,v/2)=%f\n",
    //  prob,v,chi,lg, incompletegamma(chi/2.,v/2.));

    return prob;
}

MYREAL
probchiboundary (MYREAL chi, long zeros, long all)
{
    long nonzeros = all - zeros;

//xcode    MYREAL a, b, m;
    MYREAL m;
    MYREAL xb = 1.0;  //1.0-EPSILON/1000.;
    MYREAL xa = 0.0;  //EPSILON/1000.;
    MYREAL xm = 0.51;
    if(all==0)
        return 1.;
    if (zeros == 0)
    {
        return probchi (all, chi);
    }
    //xcode a = chiboundary (zeros, nonzeros, xa);
    m = chiboundary (zeros, nonzeros, xm);
    //xcode b = chiboundary (zeros, nonzeros, xb);
    while (fabs (m - chi) > EPSILON
            && (fabs (xa - xm) > EPSILON && fabs (xb - xm) > EPSILON))
    {
        if (m < chi)
        {
            //xcode b = m;
            xb = xm;
        }
        else
        {
            //xcode a = m;
            xa = xm;
        }
        xm =   /*(-(b * xa) + chi * xa + a * xb - chi * xb) / (a - b);      */
            (xa + xb) / 2.;

        m = chiboundary (zeros, nonzeros, xm);
    }
    return xm;
}

MYREAL
chiboundary (long zeros, long nonzeros, MYREAL alpha)
{
    //  MYREAL prob;
    MYREAL sum = 0.;
    long i;
    long k = zeros;
    MYREAL freq;
    MYREAL summ;
    //printf("z=%li a=%4.2f ", k,alpha);
    for (i = 0; i <= zeros; i++)
    {
        //      sum = zerovec[i]  * (i+nonzeros) == 0 ? 0. : find_chi(i+nonzeros,alpha);
        freq = EXP (logfac (k) - logfac (k - i) - logfac (i) - LOG2 * k);
        summ = (i + nonzeros) == 0 ? 0. : find_chi (i + nonzeros, alpha);
        //      printf(" %.2f(%.4f)",freq*pow(2.,k),summ);
        sum += freq * summ;
    }
    //printf("\n");
    return sum;
}


MYREAL
chisquare (long df, MYREAL alpha)
{
    const MYREAL table05[] =
        {
            3.84146, 5.99147, 7.81473, 9.48773, 11.0705, 12.5916
        };
    const MYREAL table01[] =
        {
            6.63490, 9.21034, 11.3449, 13.2767, 15.0863, 16.8119
        };
    if( df<=6)
      {
	if (fabs(alpha - 0.05) <= (double) FLT_EPSILON)
	  return table05[df - 1];
	if (fabs(alpha - 0.01) <= (double) FLT_EPSILON)
	  return table01[df - 1];
      }
    return -100000000;
    //    return generalized_gamma(0.5,df/2.0,1.0,alpha);
}

MYREAL
calc_sum (MYREAL *vector, long n)
{
    long i;
    MYREAL summ = 0.0;
    for (i = 0; i < n; i++)
        summ += vector[i];
    return summ;
}

//==========================================
// searching and finding

boolean
find (long i, long *list, long listlen)
{
    long j;
    for (j = 0; j < listlen; j++)
    {
        if (i == list[j])
            return TRUE;
    }
    return FALSE;
}

//====================================================
// conversion between the different parameter schemes
// returns the begining of  mig_.i
long
mstart (long pop, long numpop)
{
    return numpop + pop * numpop - pop;
}

// returns the end of  mig_.i
long
mend (long pop, long numpop)
{
    return numpop + pop * numpop - pop + numpop - 1;
}

///
/// return the first element of population pop
long
mmstart (long pop, long numpop)
{
    return pop * (numpop);
}


///
/// return the element past the last element in of population pop
long
mmend (long pop, long numpop)
{
    return pop * numpop + numpop;
}

/// 
/// Returns the location in a full matrix given the abbreviated matrix
/// given i,j coordinates in an abbreviated matrix {d,d,...,d,a,a,a,...,a,b,b,...,b, ...}
/// it returns the position j* in a linearized full matrix with diagonal element "-"
/// the linearized matrix would look like this {-, a,a,a,...,a, b, - , b,  .... ,b, c, c, -, c, ....},
/// the position of i is the same for the abbreviated and the full matrix, and is not returned because
/// the calling function knows this already.
long
mm2m (long frompop, long topop, long numpop)
{
    if (frompop == topop)
        return (frompop);
    if (frompop < topop)
        return numpop + topop * (numpop - 1) + frompop;
    else
        return numpop + topop * (numpop - 1) + (frompop - 1);
}

///
/// Calulates the j and i from a linear abbreviated matrix
/// position z in {d,d,d,...,d,a,a,...,a,b,b,...,b, c, ...} we know how many populations c (columns or rows) 
/// exist and wnat to calculate the i,j coordinates in the c x c matrix {{d,a,a,...,a},{b,d,b,...},...}
/// in the function z above is im and i = frompop, j is topop
/// frompop and topop contain the result
void
m2mm (long i, long numpop, long *frompop, long *topop)
{
    if (i < numpop || numpop==1)
    {
        *frompop = i;
        *topop = i;
        return;
    }
    else
    {
        (*topop) = (long) (i - numpop) / (numpop - 1);
        (*frompop) = i - numpop - (*topop) * (numpop - 1);
        if (*frompop >= *topop)
            *frompop += 1;
    }
}

void d2mm(long pop,world_fmt *world, long *frompop, long *topop)
{
  species_fmt *s;
  long i = (pop-world->numpop2-world->bayes->mu)/2;
  if (pop<world->numpop2)
    error("d2mm called in the wrong context");
  s = &world->species_model[i]; 
  *frompop = s->from;
  *topop = s->to; 
}

long m2mmm(long frompop, long topop, long numpop)
{
  // from=1 to=0
  long v = numpop + topop * (numpop - 1);
  // v = 3 + 0*(3-1)=3
  if (frompop < topop)
    return v + frompop;
  else
    return v + frompop - 1; //3+1-1
}

// return linearize entry for speciation given the from and to
long mm2d(long frompop, long topop,world_fmt *world)
{
  long i;
  species_fmt *s;
  for (i=0;i<world->species_model_size;i++)
    {
      s = &world->species_model[i]; 
      if (frompop == s->from && topop == s->to)
	return world->numpop2 + world->bayes->mu + 2 * s->id;
    }
  error("mm2d called in the wrong context ");
  //return -1; 
}


///
/// position i in linear abbreviated array
long
m2mml (long i, long numpop)
{
    long topop, frompop;

    if (i < numpop || numpop==1)
    {
        return i * numpop + i;
    }
    else
    {
        topop = (long) (i - numpop) / (numpop - 1);
        frompop = i - numpop - (topop) * (numpop - 1);
        if (frompop >= topop)
            frompop += 1;
        return numpop * topop + frompop;
    }
}


long
mml2m (long pos, long numpop)
{
    long topop = 0, frompop = 0, i = 1;
    if (pos == 0)
        return 0;
    while (pos > numpop * (i++))
        topop++;
    frompop = pos - topop * numpop;
    return mm2m (frompop, topop, numpop);
}


long
m2mml2 (long i, long topop, long numpop)
{
    long frompop;

    if (i < numpop || numpop == 1)
    {
        return i * numpop + i;
    }
    else
    {
        frompop = i - numpop - (topop) * (numpop - 1);
        if (frompop >= topop)
            frompop += 1;
        return numpop * topop + frompop;
    }
}

#ifndef PRIORTEST
void  set_paramstr(char *paramstr, long j, world_fmt *world)
{
  long frompop;
  long topop;
  const long numpop2 = world->numpop2;
  const long numpop = world->numpop;
  const boolean usem = world->options->usem;
  const boolean rate = world->bayes->mu;
  //const long npx = numpop2 + rate + world->species_model_size * 2;
  species_fmt *s;
  long x;
  long remainder;
  if (j < numpop)
    sprintf(paramstr,"Theta_%-3li",j+1);
  else
    {
      if (j < numpop2)
	{
	  m2mm (j, numpop, &frompop, &topop);
	  if(usem)
	    {
		sprintf(paramstr, "M_%li->%li", frompop+1, topop+1);
	    }
	  else
	    {
	      sprintf(paramstr, "xN_%lim_%li->%li", topop+1, frompop+1, topop+1);
	    }
	}
      else
	{
	  if (j==numpop2 && rate)
	    sprintf(paramstr, "Rate");
	  else
	    {
	      x = j-numpop2-rate;
	      remainder = x % 2;
	      x = (long) x/2;
	      s = &world->species_model[x];
	      if (remainder == 0)
		sprintf(paramstr, "D_%li->%li",s->from+1,s->to+1);
	      else
		sprintf(paramstr, "S_%li->%li",s->from+1,s->to+1);
	    }
	}
    }
}

void
gamma_rates (MYREAL *rate, MYREAL *probcat, long categs, char *input)
{
    long i;
    MYREAL alpha = MYREAL_MAX;
    MYREAL value;
    while (!isdigit (*input) && *input != '\0')
        input++;
    if ((value = strtod (input, (char **) NULL)) > 0)
        alpha = value;
    initgammacat (categs, alpha, 1., rate, probcat);
    for (i = 0; i < categs; i++)
    {
        probcat[i] = EXP (probcat[i]);
    }
    //  calc_gamma (alpha, rate, categs);
}

/* calculation of rate values following a gamma distribution for
   given probability values */
void
calc_gamma (MYREAL alpha, MYREAL *gama, long categs)
{
    long i, panic;
    MYREAL low, mid, high, xlow, xhigh, tmp, freq = 0, x = 10, elements =
                (MYREAL) categs;
    freq = -(0.5 / elements); /*so we have midpoints instead of endpoints */
    for (i = 0; i < categs; i++)
    {
        low = 0;
        mid = incompletegamma (10., alpha);
        high = 1.;
        freq += 1. / (elements);
        if (freq < mid)
        {
            high = mid;
            xlow = 0;
            xhigh = 10.;
            x = 5.;
        }
        else
        {
            low = mid;
            xhigh = 1e10;
            xlow = 10.;
            x = 1e5;
        }
        panic = 0;
        while (panic++ < 1000 && fabs (low - high) > 0.0001 && x > 0.000000001)
        {
            mid = incompletegamma (x, alpha);
            if (freq < mid)
            {
                high = mid;
                tmp = x;
                x = (x + xlow) / 2.;
                xhigh = tmp;
            }
            else
            {
                low = mid;
                tmp = x;
                x = (x + xhigh) / 2.;
                xlow = tmp;
            }
        }
        gama[i] = x / alpha;
        //Debug
        //      printf (stderr, "  %li> %f\n", i, gama[i]);

        if (x >= 10e10)
        {
            error ("calc_gamma(): x is too big");
        }
    }
}


void
fprintf2(FILE *file, long filesize, const char *fmt, ...)
{
    char *p;
    va_list ap;
    long bufsize;
    long pallocsize = ((long) (strlen(fmt))+1+filesize);
    p  = (char *) mycalloc(pallocsize,sizeof(char));
    va_start(ap, fmt);
    bufsize = vsprintf(p, fmt, ap);
    if(bufsize>=pallocsize)
      error("failed in printf2()");
    fprintf(file,"%s", p);
    va_end(ap);
    myfree(p);
}


void
print_line (FILE * outfile, char c, long nn, long flag)
{
    long i, start = 0;
    switch (flag)
    {
    case START:
        start = 2;
        FPRINTF (outfile, "=--");
        break;
    case STOP:
        start = 2;
        FPRINTF (outfile, "==-");
        break;
    default:
        start = 0;
    }
    for (i = start; i < nn; i++)
    {
        FPRINTF (outfile, "%c",c);
    }
    FPRINTF(outfile, "\n");
}


void
sprint_line
(char *buffer, char c, long nn, long flag)
{
    char ch[2];
    long i, start = 0;
    char fp[LINESIZE];
    buffer[0] = '\0';
    ch[0] = c;
    ch[1] = '\0';
    switch (flag)
    {
    case START:
        start = 2;
        sprintf (fp, "=%c%c", c, c);
        strcat (buffer, fp);
        break;
    case STOP:
        start = 2;
        sprintf (fp, "==%c", c);
        strcat (buffer, fp);
        break;
    default:
        start = 0;
    }
    for (i = start; i < nn; i++)
    {
        strcat (buffer, ch);
    }
    strcat (buffer, "\n");
}

char
sgetc (char **buffer)
{
    char ch;
    ch = **buffer;
    (*buffer)++;
    return ch;
}

/// line-end transparent string gets command
char *
sgets (char *s, int size, char **stream)
{
    long ch = '\0';
    long counter = 0;
    while (counter < size - 1)
    {
        ch = **stream;
        (*stream)++;
        switch (ch)
        {
        case '\0':
        case '\r':
        case '\n':
            s[counter] = '\0';
            return s;
        default:
            s[counter] = (char) ch;
            break;
        }
        counter++;
    }
    return s;
}
/// line-end transparent string gets command,
/// this needs and end on the stream but allocates memory for s
char *
sgets_safe (char **s, long *size, char **stream)
{
  boolean notdone=TRUE;
  long ch = '\0';
  long counter = 0;
  memset(*s,0,sizeof(char)* (size_t) (*size));
  while (notdone)
    {
	if(counter >= *size)
	  {
	    *size *= 2;
	    *s = (char*) myrealloc(*s, (size_t)*size * sizeof(char));
	  }
        ch = **stream;
	(*stream)++;        
        switch (ch)
        {
        case '\0':
        case '\r':
        case '\n':
	  (*s)[counter] = '\0';
            return *s;
        default:
	  (*s)[counter] = (char) ch;
            break;
        }
        counter++;
    }
    return *s;
}

/// adds memory to buffer if needed to satisfy the expected increase of characters
extern void increase_buffer(char **buffer, long *allocbufsize, long *bufsize,long amount)
{
  if(*allocbufsize <= (*bufsize + amount))
    {
      *allocbufsize += amount;
      (*buffer) = (char *) myrealloc (*buffer, (*allocbufsize) * sizeof (char));
    }

}

/// sticks a text into the printing buffer
void add_to_buffer(char *fp, long *bufsize, char **buffer, long *allocbufsize)
{
  long fpsize = (long) strlen (fp) + 1;
  if(*allocbufsize <= (*bufsize + fpsize))
    {
      *allocbufsize += 100 * fpsize;
      (*buffer) = (char *) myrealloc (*buffer, (*allocbufsize) * sizeof (char));
    }
  (*bufsize) += sprintf((*buffer) + (*bufsize),"%s",fp);
}


/// sticks a text into the printing buffer
long print_to_buffer(char **buffer, long *maxbufsize, char *tempbuffer, long *pos, const char *fmt, ...)
{
  long mypos=0;
  char *p = tempbuffer;
  va_list ap;
  //p = (char *) mycalloc(1024,sizeof(char));
  va_start(ap, fmt);
  mypos = vsprintf(p, fmt, ap);
  va_end(ap);
  if((*pos + mypos) < (*maxbufsize))
      {
	(*pos) += sprintf((*buffer) + (*pos), "%s",p);
	//	if(*pos + mypos >= *maxbufsize )
	//  {
	//    printf("%i> pos=%li + mypos=%li >  maxbufsize=%li \n",myID, *pos, mypos, *maxbufsize); 
	//  }
      }
    else
      {
	*maxbufsize = *pos + 4 * mypos; // add some extra space
	(*buffer) = (char *) myrealloc ((*buffer), (*maxbufsize) * sizeof (char));
	(*pos) += sprintf((*buffer) + (*pos), "%s",p);
      }
  //  myfree(p);
  return (*pos);
}

/// sticks a text into the warning buffer that is printed at the end of the PDF and the end of the TEXT file
void record_warnings(world_fmt * world, const char *fmt, ...)
{
  long mypos=0;
  char *p = (char *) calloc(LINESIZE,sizeof(char));
  va_list ap;
  if(myID==MASTER && world->cold)
    {
      va_start(ap, fmt);
      mypos = vsprintf(p, fmt, ap);
      va_end(ap);
      if((world->warningsize + mypos) < world->warningallocsize)
	{
	  world->warningsize += sprintf(world->warning + world->warningsize, "%s\n",p);
	}
      else
	{
	  world->warningallocsize = world->warningsize + 4 * mypos; // add some extra space
	  if(world->warning != NULL)
	    world->warning = (char *) myrealloc (world->warning, world->warningallocsize * sizeof (char));
	  else
	    world->warning = (char *) mycalloc (world->warningallocsize, sizeof (char));
	  world->warningsize += sprintf(world->warning + world->warningsize, "%s\n",p);
	}
    }
  myfree(p);
}




void print_warning2(FILE *file, world_fmt *world)
{
  char paragraph[] = "This section reports potential problems with your run, but such reporting is often not very accurate. Whith many parameters in a multilocus analysi\
s, it is very common that some parameters for some loci will not be very informative, triggering suggestions (for example to increase the prior ran\
ge) that are not sensible. This suggestion tool will improve with time, therefore do not blindly follow its suggestions. If some parameters are fla\
gged, inspect the tables carefully and judge wether an action is required. For example, if you run a Bayesian inference with sequence data, for mac\
roscopic species there is rarely the need to increase the prior for Theta beyond 0.1; but if you use microsatellites it is rather common that your \
prior distribution for Theta should have a range from 0.0 to 100 or more. With many populations (>3) it is also very common that some migration rou\
tes are estimated poorly because the data contains little or no information for that route. Increasing the range will not help in such situations, \
reducing number of parameters may help in such situations.\0";
  long i=0;
  long z=0;
  char *section;
  section = (char*) mycalloc(LINESIZE,sizeof(char));
  fprintf(file,"\nPOTENTIAL PROBLEMS\n");
  print_line(file,'-',90,10);
  while (paragraph[i] != '\0')
    {
      section[z]=paragraph[i];
      section[z+1] = '\0';
      z++;
      i++;
      if(z >= 90)
	{
	  while(section[z-1]!=' ')
	    {
	      z--;
	      i--;
	    }
	  section[z]='\0';
	  fprintf(file,"%s\n",section);
	  section[0]='\0';
	  z=0;
	}
    }
  fprintf(file,"%s\n",section);
  print_line(file,'-',90,10);
  if (world->warning!=NULL && world->warning[0]=='\0')
    fprintf(file,"No warning was recorded during the run\n\n");
  else
    fprintf(file,"%s", world->warning);
  print_line(file,'-',90,10);
  myfree(section);
}


void print_stored_warnings(world_fmt *world)
{
  if(world->warningsize>0)
    {  
      print_warning2(world->outfile,world);
      if(world->options->progress)
	{
	  print_warning2(stdout,world);
	}
    }
}


// true should continue/shortcut loop
// false should do the rest of the loop, and reset j0 to j
// selects variable to work on, used on Bayesian context
// if j0 is bigger than the map return j0 and False
// except when there is growth then one gets j for growth, and after that
// return j0 and FALSE
boolean shortcut(long j0, world_fmt *world, long *j)
{
  bayes_fmt *bayes = world->bayes;
  *j = -1;
  if(j0<bayes->mapsize)
    {
      if(bayes->map[j0][1] == INVALID)
	return TRUE;
      else
	{
	  *j = bayes->map[j0][1];
	}
      if(*j < j0)
	{
	  return TRUE;
	}
      else
	{
	  return FALSE;
	}
    }
  else
    {
      if (world->has_growth)
	{
	  *j = j0;
	  long np = world->numpop2 + bayes->mu + 2 * world->species_model_size;
	  if ((j0-np<0) || (j0-np >= world->options->growpops_numalloc  ))
	    return TRUE;
	  long pick = world->options->growpops[j0-np];
	  if (pick == 0)
	    return TRUE;
      else
	{
	  *j = pick + np - 1;
	  return FALSE;
	}
    }
  else
    {
      *j = j0;
      return FALSE;
    }
}
}
#endif

// find the number of digits of a number, needs to be transformed to long before use
// only works for numbers smaller than 10^6 returns 8 otherwise
int finddigits(long number)
{
  if (number < 10000000)
    { //includes 1..7
      if (number < 10000)
	{  //includes 1..4
	  if (number < 100)
	    {   //includes 1..2
	      if (number < 10)
		{
		  return 1;
		}
	      else
		{
		  return 2;
		}
	    }
	  else
	    { //includes 3..4
	      if (number < 1000)
		{
		  return 3;
		}
	      else
		{
		  return 4;
		}
	    }
	}
      else
	{  //includes 5..7
	  if (number < 1000000)
	    { //includes 5,6
	      if (number < 100000)
		{
		  return 5;
		}
	      else
		{
		  return 6;
		}
	    }
	  else
	    { //includes 7
	      return 7;
	    }
	}
    }
  else
    {  //includes 8..inf
      return 8;
    }
  //return 8;
}

/// erf function and related
/// moved from speciate.c and added inverse_erf on July 10 2015
/// 
/// Erf approximations and hardware assistance (AVX erf does not work on my laptop macbook 2013late)
#ifdef AVX_not_working
#include <immintrin.h>
void myerf(double *x, double *result)
{
  __m256d v1 = _mm256_loadu_pd(x);
  __m256d yy = _mm256_erf_pd(v1);
  _mm256_storeu_pd(result,yy);
}
#endif

// Abramovitz and Stegun 7.1.26
// gain on Mac seems about 1/3
double myerf(double x)
{
  const double p = 0.3275911;
  const double a1= 0.254829592;
  const double a2= -0.284496736;
  const double a3=1.421413741;
  const double a4 = -1.453152027;
  const double a5 = 1.061405429;
  const double t = 1.0 / (1.0 + p*x);
  const double t2 = t * t;
  const double t3 = t * t2;
  const double t4 = t2 * t2;
  const double t5 = t2 * t3;
  return 1.0 - (a1*t + a2*t2 + a3*t3 + a4*t4 + a5*t5)*exp(-(x*x));
}
double myerfc(double x)
{
  const double p = 0.3275911;
  const double a1= 0.254829592;
  const double a2= -0.284496736;
  const double a3=1.421413741;
  const double a4 = -1.453152027;
  const double a5 = 1.061405429;
  const double t = 1.0 / (1.0 + p*x);
  const double t2 = t * t;
  const double t3 = t * t2;
  const double t4 = t2 * t2;
  const double t5 = t2 * t3;
  return (a1*t + a2*t2 + a3*t3 + a4*t4 + a5*t5)*exp(-(x*x));
}


double myerf1(double x)
{
  const double p = 0.47047;
  const double a1= 0.3480242;
  const double a2= -0.0958798;
  const double a3=0.7478556;
  double t = 1.0 / (1.0 + p*x);
  double t2 = t * t;
  double t3 = t * t2;
  return 1.0 - (a1*t + a2*t2 + a3*t3)*exp(-(x*x));
}

//alternative implementation (not exposed to other files)
double myerfc1(double x)
{
  const double p = 0.47047;
  const double a1= 0.3480242;
  const double a2= -0.0958798;
  const double a3=0.7478556;
  double t = 1.0 / (1.0 + p*x);
  double t2 = t * t;
  double t3 = t * t2;
  return (a1*t + a2*t2 + a3*t3)*exp(-(x*x));
}

//
// Lower tail quantile for standard normal distribution function.
//
// This function returns an approximation of the inverse cumulative
// standard normal distribution function.  I.e., given P, it returns
// an approximation to the X satisfying P = Pr{Z <= X} where Z is a
// random variable from the standard normal distribution.
//
// The algorithm uses a minimax approximation by rational functions
// and the result has a relative error whose absolute value is less
// than 1.15e-9.
//
// Author:      Peter John Acklam
// (Javascript version by Alankar Misra @ Digital Sutras (alankar@digitalsutras.com))
// Time-stamp:  2003-05-05 05:15:14
// E-mail:      pjacklam@online.no
// WWW URL:     http://home.online.no/~pjacklam

// An algorithm with a relative error less than 1.15*10-9 in the entire region.
// converted from a javascript function to C [pbeerli]
double inverse_cumstd_normal(double p)
{
  // Coefficients in rational approximations
  const double a[] = {-3.969683028665376e+01,  2.209460984245205e+02,
		      -2.759285104469687e+02,  1.383577518672690e+02,
		      -3.066479806614716e+01,  2.506628277459239e+00};

  const double b[] = {-5.447609879822406e+01,  1.615858368580409e+02,
		      -1.556989798598866e+02,  6.680131188771972e+01,
		      -1.328068155288572e+01 };

  const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01,
		      -2.400758277161838e+00, -2.549732539343734e+00,
		      4.374664141464968e+00,  2.938163982698783e+00};

  const double d[] = {7.784695709041462e-03, 3.224671290700398e-01,
		     2.445134137142996e+00,  3.754408661907416e+00};

  // Define break-points.
  const double plow  = 0.02425;
  const double phigh = 1 - plow;

  // Rational approximation for lower region:
  if ( p < plow ) {
    double q  = sqrt(-2.*log(p));
    return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
      ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }

  // Rational approximation for upper region:
  if ( phigh < p ) {
    double q  = sqrt(-2.*log(1-p));
    return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
      ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }

  // Rational approximation for central region:
  double q = p - 0.5;
  double r = q*q;
  return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
    (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
}

// z = erfinv(x) where x is between -1 and 1 returns a value z so that x=erf(z)
double erfinv(double x)
{
  if (x < -1.0 || x > 1.0)
    {
      warning("erfinv received illegal value");
      return (double) -HUGE;
    }
  //  printf("%f %f\n", x, inverse_cumstd_normal((x+1.0)/2.0)/SQRT2);
  return inverse_cumstd_normal((x+1.0)/2.0)/SQRT2;
}


/******************************************************************************/
/*
  double wew(double, *double)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 February 2005
  Author:
    John Burkardt
  Purpose:
    WEW estimates Lambert's W function.

  Discussion:
    For a given X, this routine estimates the solution W of Lambert's 
    equation:
      X = W * EXP ( W )
    This routine has higher accuracy than WEW_B.
  Modified:
    11 June 2014
  Reference:
    Fred Fritsch, R Shafer, W Crowley,
    Algorithm 443: Solution of the transcendental equation w e^w = x,
    Communications of the ACM,
    October 1973, Volume 16, Number 2, pages 123-124.
  Parameters:
    Input, double X, the argument of W(X)
    Output, double *EN, the last relative correction to W(X).
    Output, double WEW_A, the estimated value of W(X).
*/
double wew( double x, double *en )
{
  const double c1 = 4.0 / 3.0;
  const double c2 = 7.0 / 3.0;
  const double c3 = 5.0 / 6.0;
  const double c4 = 2.0 / 3.0;
  double f;
  double temp;
  double temp2;
  double wn;
  double y;
  double zn;
/*
  Initial guess.
*/
  f = log ( x );

  if ( x <= 6.46 )
  {
    wn = x * ( 1.0 + c1 * x ) / ( 1.0 + x * ( c2 + c3 * x ) );
    zn = f - wn - log ( wn );
  }
  else
  {
    wn = f;
    zn = - log ( wn );
  }
/*
  Iteration 1.
*/
  temp = 1.0 + wn;
  y = 2.0 * temp * ( temp + c4 * zn ) - zn;
  wn = wn * ( 1.0 + zn * y / ( temp * ( y - zn ) ) );
/*
  Iteration 2.
*/
  zn = f - wn - log ( wn );
  temp = 1.0 + wn;
  temp2 = temp + c4 * zn;
  *en = zn * temp2 / ( temp * temp2 - 0.5 * zn );
  wn = wn * ( 1.0 + *en );

  return wn;
}

// returns t for growth
//ProductLog[-((g theta0 Log[r])/(2 k))]/g
// assumes that growth != 0.0
double interval_growth(double r, double t0, double theta0, double growth, double k, double tmin, double tmax)
{
  double minu=tmin;
  double maxu=tmax;
  double u=(maxu+minu)/2.0;
  double prob;
  long maxcount=200;
  double logr = log(r);
  while((maxu-minu)>EPSILON && maxcount-- > 0)
    {
      prob =  - k * (exp(growth * (t0+u)) - exp(growth * t0))/(theta0 * growth);
      //printf("@ %f %f %f prob=%f r=%f \n",minu,u,maxu, prob, logr);
      if (logr>prob)
	{
	  maxu = u;
	}
      else
	{
	  minu = u;
	}
      u = (maxu+minu)/2.0;
    }
  return u;
}



double get_time_for_growth(double theta0, double growth, double k, double t0)
{
  // revision using my own
  // logr <= 0
  double logr = LOG(UNIF_RANDUM());
  //double interval =  (1./(t0 * growth) * log(-1/k *  logr * (theta0 * growth) - exp(growth * t0)))/growth;
#ifdef WINDOWS
  MYCOMPLEX tmp = {-1.0 + (growth * theta0 * logr)/(k * exp(growth*t0)), 0.0};
  MYCOMPLEX tmp2 = clog(tmp);
  MYCOMPLEX tmp3 = _Cmulcr(tmp2, 1.0/growth);
  double interval = creal(tmp3);
#else
  double interval = creal(clog(-1.0 + (growth * theta0 * logr)/(k * exp(growth*t0)))/growth);
#endif
  /////logr =  - log(k * (exp(growth * (t0+u)) - exp(growth * t0))/(theta0 * growth));
  //logr =  - [log(k) + log(exp(growth * (t0+u)) - exp(growth * t0))-log(theta0 * growth))];
  
  
  //double tk = - LOG(UNIF_RANDUM()) * theta0 / k;
  //double interval = (-growth * t0 + LOG(exp(growth * t0) + growth * tk))/growth;
  if (interval<0.0)
    return 0.0;
  else
    return interval;
}





// calculates a rough approximation to log(val)
// max error = 0.0007
// Submitted by Laurent de Soras, posted on 30 March 2001
// http://www.flipcode.com/cgi-bin/msg.cgi?showThread=Tip-Fastlogfunction&forum=totd&id=-1
// Fast log() Function, by Laurent de Soras:
// Here is a code snippet to replace the slow log() function...
// It just performs an approximation, but the maximum error is below 0.007.
// Speed gain is about x5, and probably could be increased by tweaking the assembly code.
// The function is based on floating point coding.
// It's easy to get floor (log2(N)) by isolating exponent part.
// We can refine the approximation by using the mantissa. This function returns log2(N)
/*
MYINLINE float fast_log2 (float vval)
{
    float vv;
    int * const exp_ptr =  (int *) (&vval);
    int            x = *exp_ptr;
    const int      log_2 = ((x >> 23) & 255) - 128;
    x &= ~(255 << 23);
    x += 127 << 23;
    *exp_ptr = x;

    vval = ((-1.0/3) * vval + 2) * vval - 2.0/3;   // (1)
    vv = vval + log_2;
    return (vv);
}

// The line (1) computes 1+log2(m), m ranging from 1 to 2. The proposed
// formula is a 3rd degree polynomial keeping first derivate
// continuity. Higher degree could be used for more accuracy. For faster
// results, one can remove this line, if accuracy is not the matter (it
// gives some linear interpolation between powers of 2).
//Now we got log2(N), we have to multiply it by ln(2) to get the natural log :


MYINLINE float fast_log (float vval)
{
    float v = fast_log2 (vval) * 0.69314718f;
    //fprintf(stdout,"val= %f fast_log=%f log=%f\n", (float) vval, (float) v, log(vval));
    return v;
}

//#define myEXPA (1048576 / 0.693147180559945309417232121458)
#define myEXPA 1512775.39519518569383584038231 
#define myEXPC 60801 

MYINLINE double fast_exp(double y) 
{ 
    union 
    { 
        double d; 
      // weird the implementation says ifdef but that does no work at all
#ifndef LITTLE_ENDIAN 
        struct { int j, i; } n; 
#else 
        struct { int i, j; } n; 
#endif 
    } eco;
    eco.n.i = (int)(myEXPA*(y)) + (1072693248 - myEXPC);
    eco.n.j = 0; 
    return eco.d; 
} 

*/

#ifdef TOOLSTEST
int main()
{
  // testing erfinv
  double a1 = erfinv_approx(0.5);
  double a2 = erfinv_maclaurin(0.5);
  double b1 = erf(a1);
  double b2 = erf(a2);
  printf("a1(%f)=erfinv(%f),a2(%f)=erfinv(%f),\n",a1,0.5,a2,0.5);
  printf("b1(%f)=erf(%f),b2(%f)=erf(%f))\n",b1,a1,b2,a2);
  printf("testing trim()\n");
  char * line = (char *) calloc(100, sizeof(char));
  strcpy(line,"   this has trailing and leading blanks   ");
  printf("priginal   :%s\n",line);
  trim(&line);
  printf("manipulated:%s\n",line);
  free(line);
}
#endif
