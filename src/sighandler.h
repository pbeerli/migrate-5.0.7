#ifndef SIGNAL_HANDLER
#define SIGNAL_HANDLER
/* -----------------------------------------------------
   sighandler.h                                           
   handels to following signals:                          
   SIGIOT         Input/Output problems                   
   SIGFPE         Floating point exceptions               
   SIGBUS         Bus error                               
   SIGSEGV        Segmentation fault                      
   SIGXCPU        CPU time limit exceeded                 
   SIGXFSZ        File size limit exceeded                
   if most of those exception are encountered the system  
   tries to exit gracefully, but with some it dies        
   anyway, but tries to say it why in a way which is     
   for humans better understandable                       
   -----------------------------------------------------  
 (c) Peter Beerli 2013
 
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
 
 
   ----------------------------------------------------- */
#ifndef __MWERKS__
#include <sys/types.h>
#endif
#include <stdarg.h>
#include <signal.h>
#include "definitions.h"
#define ON 1L
#define OFF 0L

extern void signalhandling (long switcher);
extern void signalhandler (int sig);
extern void warning (char string[], ...);
extern void usererror (char string[], ...) __attribute__((noreturn)) ;
extern void sig_error (char string[], char filename[], long line) __attribute__((noreturn)) ;
#define error(x) sig_error(x, __FILE__,(long) __LINE__) 

extern char myfgets (char *buffer, long bufsize, FILE * infile);
extern char myfgetssafe (char **buffer, long *bufsize, FILE * infile);
#ifdef ZNZ
#include "znzlib.h"
extern char myznzfgets (char *buffer, long bufsize, znzFile infile);
#endif
#ifdef NEXTAPP
extern void malloc_error_found (int value);
#endif

#ifdef MEMDEBUG
extern void memdebug_free(void *ptr, const char file[], const long line);
#define myfree(a) memdebug_free(a,__FILE__,__LINE__);
#else
#define myfree(a)  free(a)
#endif

#ifdef MIGRATE_MALLOC
extern void *MIGRATE_malloc (size_t size, const char file[],
                                const long line);
extern void *MIGRATE_calloc (size_t repeats,
                                size_t size, const char file[],
                                const long line);
extern void *MIGRATE_realloc (void *ptr, size_t,
                                 const char file[], const long line);
#define mymalloc(a) MIGRATE_malloc((size_t) a,__FILE__,(const long)__LINE__)
#define mycalloc(b,a) MIGRATE_calloc((size_t) b, (size_t) a,__FILE__,(const long)__LINE__)
#define myrealloc(c,a) MIGRATE_realloc(c, (size_t) a,__FILE__,(const long) __LINE__)
#else
#define mymalloc(a) malloc((size_t) a)
#define mycalloc(b,a) calloc((size_t) b, (size_t) a)
#define myrealloc(c,a) realloc(c, (size_t) a)
#endif
#ifdef MPI
#ifdef DEBUG_MPI
#define MYMPISEND(a,b,c,d,e,f) printf("%i>before ssend of %li bytes to %i with tag %i (%s,%i)\n",myID, b* sizeof(c), d, e,__FILE__,__LINE__);fflush(stdout);MPI_Ssend(a,b,c, d,  e,f);printf("%i>after send (%s,%i)\n",myID, __FILE__,__LINE__);fflush(stdout)
//#define MYMPIISEND(a,b,c,d,e,f,g) printf("%i>before isend of %li bytes to %i with tag %i in %i (%s,%i)\n",myID, b* sizeof(c), d, e, (int) f, __FILE__,__LINE__);fflush(stdout);MPI_Isend(a,b,c,d,  e,f,g);printf("%i>after send (%s,%i)\n",myID, __FILE__,__LINE__);fflush(stdout)
#define MYMPIISEND(a,b,c,d,e,f,g) printf("%i>before isend of %li bytes to %i with tag %i in %i (%s,%i)\n",myID, b* sizeof(c), d, e, (int) f, __FILE__,__LINE__);fflush(stdout);MPI_Send(a,b,c,d,  e,f);printf("%i>after send (%s,%i)\n",myID, __FILE__,__LINE__);fflush(stdout)
#define MYMPIRECV(a,b,c,d,e,f,g) printf("%i>before receive of %li bytes from %i with tag %i in %i (%s,%i)\n",myID, b * sizeof(c),d,e, (int) f, __FILE__,__LINE__);fflush(stdout);MPI_Recv(a,b,c,d,e,f,g);printf("%i>after receive (%s,%i)\n",myID, __FILE__,__LINE__);fflush(stdout)
#define MYMPIBCAST(a,b,c,d,e) printf("%i>before broadcast of %li bytes (%s,%i)\n",myID, ((long) b) * ((long) sizeof(c)), __FILE__,__LINE__);fflush(stdout);MPI_Bcast(a,b,c,d,e);printf("%i>after broadcast (%s,%i)\n",myID, __FILE__,__LINE__);fflush(stdout)
#define MYMPIBARRIER(a) printf("%i>before barrier (%s,%i)\n",myID, __FILE__,__LINE__);fflush(stdout);MPI_Barrier(a);printf("%i>after barrier (%s,%i)\n",myID, __FILE__,__LINE__);fflush(stdout)
//#define MYMPIWAITALL(a,b,c) printf("%i>before waitall (%s,%i)\n",myID, __FILE__,__LINE__);fflush(stdout);MPI_Waitall(a,b,c);printf("%i>after waitall (%s,%i)\n",myID, __FILE__,__LINE__);fflush(stdout)
#define MYMPIWAITALL(a,b,c) /*do nothing*/
#else
#define MYMPISEND(a,b,c,d,e,f) MPI_Send(a,(int) b,c,d,e,f)
//#define MYMPIISEND(a,b,c,d,e,f,g) MPI_Isend(a,(int) b,c,d,e,f,g)
#define MYMPIISEND(a,b,c,d,e,f,g) MPI_Send(a,(int) b,c,d,e,f)
#define MYMPIRECV(a,b,c,d,e,f,g) MPI_Recv(a,(int) b,c,d,e,f,g)
#define MYMPIBCAST(a,b,c,d,e) MPI_Bcast(a,(int) b,c,d,e)
#define MYMPIBARRIER(a) MPI_Barrier(a)
//#define MYMPIWAITALL(a,b,c) MPI_Waitall(a,b,c)
#define MYMPIWAITALL(a,b,c) /*do nothing*/
#endif
#endif
#endif

