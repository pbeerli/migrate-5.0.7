/* -----------------------------------------------------
   sighandler.c                         
   handels to following signals:                
   SIGIOT         Input  Output problems             
   SIGIRAP        Over  Underflow, 0 divide          
   SIGFPE         Floating point exceptions         
   SIGBUS         Bus error                 
   SIGSEGV        Segmentation fault                
   SIGXCPU        CPU time limit exceeded           
   SIGXFSZ        File size limit exceeded          
   SIGILL         Illegal instruction                       
   SIGUSR1        User signal 1                             
   if most of those exception are encountered the system     
   tries to exit gracefully, but with some it dies          
   anyway, but tries to say why in a way which is      
   for humans better understandable  (I hope)                       
   -----------------------------------------------------     
   part of the migrate                   
 
   Copyright 1996-2005 Peter Beerli
 
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

   $Id: sighandler.c 2158 2013-04-29 01:56:20Z beerli $                  
   ----------------------------------------------------- */
/*! \file sighandler.c */
#ifndef WINDOWS
#pragma clang diagnostic ignored "-Wformat-nonliteral"
#endif

//#ifndef CLANG_ANALYZER_NORETURN
//#if __has_feature(attribute_analyzer_noreturn)
//#define CLANG_ANALYZER_NORETURN __attribute__((analyzer_noreturn))
//#else
//#define CLANG_ANALYZER_NORETURN
//#endif
//#endif


#include "sighandler.h"
#include <wchar.h>
#include "migrate_mpi.h"
#ifdef ZNZ
#include "znzlib.h"
#endif
#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif
#ifdef WATCOM
#define __MWERKS__
#endif
#undef debug

extern int myID;
#ifdef MEMDEBUG
#include <sys/time.h>
extern FILE *memfile;
extern struct timeval memt_start, memt_finish;
extern double memelapsed;
extern long totalsize;
#endif

void usererror (char string[], ...) __attribute__((noreturn));
void sig_error (char string[], char filename[], long line) __attribute__((noreturn));


void
signalhandling (long switcher)
{
    if (switcher == ON)
    {
#ifndef WIN32
#ifndef __MWERKS__
#ifdef SIGIOT
        signal (SIGIOT, signalhandler);
#endif
        signal (SIGTRAP, signalhandler);
        signal (SIGBUS, signalhandler);
	//              signal (SIGUSR1, signalhandler);
#ifndef SYSTEM_V

        signal (SIGXCPU, signalhandler);
        signal (SIGXFSZ, signalhandler);
#endif
#endif
#endif

        signal (SIGFPE, signalhandler);
        signal (SIGSEGV, signalhandler);
        signal (SIGILL, signalhandler);
    }
    else
    {
#ifndef WIN32
#ifndef __MWERKS__
#ifdef SIGIOT
        signal (SIGIOT, SIG_DFL);
#endif
        signal (SIGTRAP, SIG_DFL);
        signal (SIGBUS, SIG_DFL);
	//	      signal (SIGUSR1, SIG_DFL);
#ifndef SYSTEM_V

        signal (SIGXCPU, SIG_DFL);
        signal (SIGXFSZ, SIG_DFL);
#endif
#endif
#endif
        signal (SIGFPE, SIG_DFL);
        signal (SIGSEGV, SIG_DFL);
        signal (SIGILL, SIG_DFL);
    }
}

void
signalhandler (int sig)
{
  //fputc ('\040', stderr);
    fflush (NULL);

    switch (sig)
    {
#ifndef WIN32
#ifndef __MWERKS__
#ifdef SIGIOT
    case SIGIOT:
#endif
#ifdef MPI

        fprintf(stderr, "On node %i happened a", myID);
#endif

        fprintf (stderr, "\nSEVERE ERROR: Input/Output problem!\n");
        fprintf (stderr,
                 "              Most likely a fatal error in the infile or parmfile\n\n");
        fprintf (stderr,
                 "              Please report error with as much detail as possible to\n");
        fprintf (stderr, "              %s\n\n", MAINTAINER);
        exit (sig);
        //break;
    case SIGTRAP:
#ifdef MPI

        fprintf(stderr, "On node %i happened a", myID);
#endif

        fprintf (stderr, "SEVERE ERROR: Trace trap\n");
        fprintf (stderr,
                 "              there was an an overflow/underflow/0 divide problem\n");
        fprintf (stderr, "              or some other problems.\n");
        fprintf (stderr,
                 "              But check your infile for errors, too.\n");
        fprintf (stderr,
                 "              Please report error with as much detail as possible to\n");
        fprintf (stderr, "              %s\n\n", MAINTAINER);
        exit (sig);
        //break;
    case SIGBUS:
#ifdef MPI

        fprintf(stderr, "On node %i happened a", myID);
#endif

        fprintf (stderr, "SEVERE ERROR: Bus problems\n");
        fprintf (stderr,
                 "              this results in an non recoverable crash\n");
        fprintf (stderr,
                 "              But check your infile for errors, too.\n");
        fprintf (stderr,
                 "              Please report error with as much detail as possible to\n");
        fprintf (stderr, "              %s\n\n", MAINTAINER);
        exit (sig);
        //break;
    case SIGXCPU:
#ifdef MPI

        fprintf(stderr, "On node %i happened a", myID);
#endif

        fprintf (stderr, "SEVERE ERROR: This program has a time limit?\n");
        fprintf (stderr, "              We didn' program in a time limit,\n");
        fprintf (stderr,
                 "              report this error to your system administrator, please.\n");
        exit (sig);
        //break;
    case SIGXFSZ:
#ifdef MPI

        fprintf(stderr, "On node %i happened a", myID);
#endif

        fprintf (stderr, "SEVERE ERROR: This program has a file size limit?\n");
        fprintf (stderr,
                 "              We didn't program in a file size limit,\n");
        fprintf (stderr,
                 "              report this error to your system administrator, please.\n");
        exit (sig);
        //break;
	////    case SIGUSR1:
#ifdef MPI
	  ////fprintf(stderr, "On node %i happened a USER (%i) SIGNAL\n", myID, sig);
	  ////MPI_Finalize();
	//	exit(sig);
#endif
        //      warning ("\nUser signal received and ignored\n\n");
        ////break;

#endif /**__MWERKS__*/
#endif /**WIN32*/
    case SIGFPE:
#ifdef MPI

        fprintf(stderr, "On node %i happened a", myID);
#endif

        fprintf (stderr, "SEVERE ERROR: Floating point exception\n");
        fprintf (stderr,
                 "              There was an integer/floating point problem\n");
        fprintf (stderr,
                 "              Often this is a division by zero. If you dataset is moderately sized\n");
        fprintf (stderr,
                 "              this is most likely an error in you data (infile).\n");
        fprintf (stderr,
                 "              Please report error with as much detail as possible to\n");
        fprintf (stderr, "              %s\n\n", MAINTAINER);
        exit (sig);
        //break;
    case SIGSEGV:
#ifdef MPI

        fprintf (stderr, "SEVERE ERROR: Segmentation fault on Node **%i**\n",
                 myID);
#else

        fprintf (stderr, "SEVERE ERROR: Segmentation fault\n");
#endif

        fprintf (stderr,
                 "              this results in an non recoverable crash.\n");
        fprintf (stderr,
                 "              But check the datatype and your infile for errors, too.\n");
        fprintf (stderr,
                 "              Please report error with as much detail as possible to\n");
        fprintf (stderr, "              %s\n\n", MAINTAINER);
        exit (sig);
        //break;
    case SIGILL:
#ifdef MPI

        fprintf(stderr, "On node %i happened a", myID);
#endif

        fprintf (stderr, "SEVERE ERROR: Illegal instruction!\n");
        fprintf (stderr, "              this is maybe a programming error.\n");
        fprintf (stderr,
                 "              But check your infile for errors, too.\n");
        fprintf (stderr,
                 "              Please report error with as much detail as possible to\n");
        fprintf (stderr, "              %s\n\n", MAINTAINER);
        exit (sig);
        //break;
    default:
        break;
    }
#ifndef WIN32
#ifndef __MWERKS__
#ifdef SIGIOT
    signal (SIGIOT, signalhandler);
#endif
    signal (SIGTRAP, signalhandler);
    //    signal (SIGUSR1, signalhandler);
    signal (SIGBUS, signalhandler);
#ifndef SYSTEM_V

    signal (SIGXCPU, signalhandler);
    signal (SIGXFSZ, signalhandler);
#endif
#endif
#endif
    signal (SIGFPE, signalhandler);
    signal (SIGSEGV, signalhandler);
    signal (SIGILL, signalhandler);

}

void
warning (char string[], ...)
{
    va_list args;
    
    fprintf (stdout, "WARNING: ");
    va_start (args, string);
    vfprintf (stdout, string, args);
    va_end (args);
}


void usererror (char string[], ...) 
{
    va_list args;
    fprintf (stdout, "\nERROR: ");
    va_start (args, string);
    vfprintf (stdout, string, args);
    va_end (args);
    exit (EXIT_FAILURE);
}

#ifndef WINDOWS
__attribute((analyzer_noreturn))
#endif
void sig_error (char string[], char filename[], long line)
{
  fprintf (stdout, "\n%i> ERROR in file %s on line %li\nERROR %s\n", myID, filename, line, string);
    exit (EXIT_FAILURE);
}

///
/// replacement for fgets so that migrate can read 
/// files from different computers with different line-endings
/// should be able to read a file with
/// - windows
/// - mac
/// - unix
/// line-endings.
#ifdef ZNZ
char 
myznzfgets (char *buffer, long bufsize, znzFile infile)
{
  long count = 0;
  char ch = '\0';
  char oldch;
  bufsize--;
  ch = (char) znzgetc (infile);
  while (!(ch == '\r' || ch == '\n' || ch == EOF) && count < bufsize)
    {
      buffer[count] = ch;
      count++;
      ch = (char) znzgetc (infile);
    }
  if (ch == '\r')
    {
      oldch = ch;
      ch = (char)  znzgetc (infile);
      if(!((ch == EOF) || (ch == '\n')))
	{
	  znzungetc (ch, infile);
	  buffer[count] = '\0';
	   return oldch;
	}
    }
  buffer[count] = '\0';
  return ch;
}
#endif
///
/// replacement for fgets so that migrate can read 
/// files from different computers with different line-endings
/// should be able to read a file with
/// - windows
/// - mac
/// - unix
/// line-endings.

char 
myfgets (char *buffer, long bufsize, FILE * infile)
{
  long count = 0;
  char ch = '\0';
  char oldch;
  bufsize--;
  ch = (char)  fgetc (infile);
  while (!(ch == '\r' || ch == '\n' || ch == EOF) && count < bufsize)
    {
      //ch = fgetc (infile);
      buffer[count] = ch;
      count++;
      ch = (char) fgetc (infile);
    }
  if (ch == '\r')
    {
      oldch = ch;
      ch = (char) fgetc (infile);
      if(!((ch == EOF) || (ch == '\n')))
	{
	  ungetc (ch, infile);
	  buffer[count] = '\0';
	   return oldch;
	}
    }
  buffer[count] = '\0';
  return ch;
}

///
/// replacement for fgets that reallocates if things are longer than buffer
//  so that migrate can read 
/// files from different computers with different line-endings
/// should be able to read a file with
/// - windows
/// - mac
/// - unix
/// line-endings.
char 
myfgetssafe (char **buffer, long *bufsize, FILE * infile)
{
  long count = 0;
  char ch = '\0';
  char oldch;
  long bs = *bufsize;
  bs--;
  ch = (char) fgetc (infile);
  while (!(ch == '\r' || ch == '\n' || ch == EOF) && count < bs)
    {
      (*buffer)[count] = ch;
      count++;
      if(count>=bs)
	{
	  bs += STRSIZE;
	  *bufsize += STRSIZE;
	  (*buffer) = (char *) realloc(*buffer, sizeof(char) * (size_t) (*bufsize));
	}
      ch = (char) fgetc (infile);
    }
  if (ch == '\r')
    {
      oldch = ch;
      ch = (char) fgetc (infile);
      if(!((ch == EOF) || (ch == '\n')))
	{
	  ungetc (ch, infile);
	  (*buffer)[count] = '\0';
	   return oldch;
	}
    }
  (*buffer)[count] = '\0';
  return ch;
}

#ifdef NEXTAPP
void
malloc_error_found (int value)
{
    switch (value)
    {
    case 0:
        fprintf (stderr, "vm_allocate failed\n");
        break;
    case 1:
        fprintf (stderr, "vm_deallocate failed\n");
        break;
    case 2:
        fprintf (stderr, "vm_copy failed\n");
        break;
    case 3:
        fprintf (stderr,
                 "I tried to reallocate or free space which was already freed\n");
        break;
    case 4:
        fprintf (stderr, "Internal memory verification in heap failed\n");
        break;
    case 5:
        fprintf (stderr,
                 "I tried to reallocate or free space which was never allocated\n");
        break;
    default:
        fprintf (stderr, "huh, what mymalloc error was this?\n");
        break;
    }
    errno = ENOMEM;
}
#endif
#ifdef MEMDEBUG
///*memdebug routine for free()*/
void memdebug_free(void *ptr, const char file[], const long line)
{
  long oldsize;
  gettimeofday(&memt_finish, NULL);  
  memelapsed = memt_finish.tv_sec - memt_start.tv_sec +
    (memt_finish.tv_usec - memt_start.tv_usec) / 1.e6;
  oldsize = malloc_size(ptr);
  totalsize -= oldsize;
  fprintf(memfile,"%li %li %li %s\n",(long) oldsize, totalsize, line, file);
  //  fprintf(memfile,"%i> %p %.9f . free (%li_%li) in %s at %li\n",myID, ptr, memelapsed, oldsize,totalsize, file, line);
  free(ptr);
}
#endif

/* save memory routines, this routines are used
   with defines in migration.h, DO NOT move this part somewhere else */
#ifdef MIGRATE_MALLOC
//#undef mymalloc
//#undef realloc
//#undef calloc
//#undef free

void *
MIGRATE_malloc (size_t size, const char file[], const long line)
{
    void *x;
    size_t newsize;
    if(size == 0)
      {
        fprintf (stderr,
                 "%i> WARNING in memory allocation (malloc(%li)) in file %s on line %li\n",
                 myID, (long) size, file, line);
        newsize = 1;
      }
    else
      newsize = size;

    x = malloc (newsize);
#ifdef MEMDEBUG
  gettimeofday(&memt_finish, NULL);  
  memelapsed = memt_finish.tv_sec - memt_start.tv_sec +
    (memt_finish.tv_usec - memt_start.tv_usec) / 1.e6;
  totalsize += malloc_size(x);
  fprintf(memfile,"%li %li %li %s\n",(long) newsize, totalsize, line, file);
  //  fprintf(memfile,"%i> %p %.9f . malloc size=%li_%li in %s at %li\n",myID, x, memelapsed, (long) size, totalsize, file, line);
#endif    
    if (x == NULL)
    {
#ifdef MPI
        fprintf (stderr,
                 "%i> ERROR in memory allocation (malloc(%li)) in file %s on line %li\n",
                 myID, (long) size, file, line);
#else
        fprintf (stderr,
                 "ERROR in memory allocation (malloc(%li)) in file %s on line %li\n",
                 (long) size, file, line);
#endif
        exit (-1);
    }
    return x;
}


void *
MIGRATE_calloc (size_t repeats, size_t size,
               const char file[], const long line)
{
    void *x;
    size_t newsize;
    if(size == 0)
            {
        fprintf (stderr,
                 "%i> WARNING in memory allocation (malloc(%li)) in file %s on line %li\n",
                 myID, (long) size, file, line);
	newsize = 1;
      }
    else
      newsize = size;

    x = calloc (repeats, newsize);
#ifdef MEMDEBUG
  gettimeofday(&memt_finish, NULL);  
  memelapsed = memt_finish.tv_sec - memt_start.tv_sec +
    (memt_finish.tv_usec - memt_start.tv_usec) / 1.e6;
  totalsize += malloc_size(x);
  fprintf(memfile,"%li %li %li %s\n",(long) repeats*newsize, totalsize, line, file);
  //    fprintf(memfile,"%i> %p %.9f . calloc:%lix%li  size=%li_%li in %s at %li\n", myID, x, memelapsed, (long) repeats, (long) size, (long)(repeats * size), totalsize, file, line);
#endif    
    if (x == NULL)
    {
#ifdef MPI
        fprintf (stderr,
                 "%i> ERROR in memory allocation (calloc(%li,%li)) in file %s on line %li\n",
                 myID, (long) repeats, (long) size, file, line);
#else
        fprintf (stderr,
                 "ERROR in memory allocation (calloc(%li,%li)) in file %s on line %li\n",
                 (long) repeats, (long) size, file, line);
#endif
        exit (-1);
    }
    return x;
}

void *
MIGRATE_realloc (void *ptr, size_t size, const char file[],
                const long line)
{
    void *x;
#ifdef MEMDEBUG
    long oldsize;
#endif
    unsigned long newsize;
    //if(size == 0)
    //        {
    //    fprintf (stderr,
    //             "%i> WARNING in memory allocation (malloc(%li)) in file %s on line %li\n",
    //             myID, (long) size, file, line);
    //	newsize = 1;
    //  }
    //else
    newsize = (unsigned long) size;
#ifdef MEMDEBUG
    oldsize = malloc_size(ptr);
#endif
    x = realloc (ptr, newsize);
#ifdef MEMDEBUG
    gettimeofday(&memt_finish, NULL);  
    memelapsed = memt_finish.tv_sec - memt_start.tv_sec +
    (memt_finish.tv_usec - memt_start.tv_usec) / 1.e6;
    newsize = malloc_size(x);
    totalsize += newsize - oldsize;
    fprintf(memfile,"%li %li %li %s\n",(long) newsize-oldsize, totalsize, line, file);
    //fprintf(memfile,"%i> %p %.9f %p realloc size=%li_%li in %s at %li\n",myID, x, memelapsed, ptr, (long) newsize, totalsize, file, line);
#endif    
    
    if (x == NULL)
    {
#ifdef MPI
        fprintf (stderr,
                 "%i> ERROR in memory allocation (realloc(ptr,%li)) in file %s on line %li\n",
                 myID, (long) size, file, line);
#else
        fprintf (stderr,
                 "ERROR in memory allocation (realloc(ptr,%li)) in file %s on line %li\n",
                  (long) size, file, line);
#endif
        exit (-1);
    }
    return x;
}


#endif /*MIGRATE_MALLOC */
