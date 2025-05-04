#ifndef TOOLS_INCLUDE
#define TOOLS_INCLUDE
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                       
 
 H E L P E R     R O U T I N E S 
 
 Peter Beerli 1996, Seattle
 beerli@fsu.edu
 
(c) 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
(c) 2003-2006 Peter Beerli, Tallahassee FL
 
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

$Id: tools.h 2157 2013-04-16 22:27:40Z beerli $
-------------------------------------------------------*/

#include "migration.h"
// selects variable to work on, used on Bayesian context
// true should continue/shortcut loop
// false should do the rest of the loop, and reset j0 to j
							  
extern void clearscreen(void);
extern char *my_strdup(const char *s);
#ifndef HAVE_STRDUP
#define strdup(x) my_strdup(x)
#endif
#ifndef HAVE_STRSEP
extern char *strsep(char **, const char *);
#endif
extern boolean shortcut(long j0, world_fmt *world, long *j);
extern int finddigits(long number);
extern MYREAL lengthof (node * p);
extern MYINLINE node *crawlback (const node * theNode);
extern node *crawl (node * theNode);
extern MYINLINE node *showtop (node * theNode);
extern void adjust_time (node * theNode, MYREAL tyme);
extern void adjust_time_all (node * theNode, MYREAL tyme);
extern void insert_migr_node (world_fmt * world, node * up, node * down,
                                  migr_table_fmt * migr_table,
                                  long *migr_table_counter);
extern long offspring (node * mother, node ** myoffspring);
extern void children (node * mother, node ** brother, node ** sister);
extern void siblings (node * mother, node ** brother, node ** sister);
/* math tools */
extern MYREAL incompletegamma (MYREAL tx, MYREAL talpha);
extern MYREAL logincompletegamma (MYREAL tx, MYREAL talpha);
extern MYREAL uppergamma_reg (MYREAL tx, MYREAL talpha);
extern MYREAL uppergamma (MYREAL tx, MYREAL talpha);
extern MYREAL polygamma (long n, MYREAL z);
extern void invert_matrix (MYREAL **a, long nsize);
extern boolean nrcheck (MYREAL **m, MYREAL **tm, MYREAL *v, long nrows,
                            MYREAL *r1, MYREAL *r2, boolean do_newton);
extern MYREAL rannor (MYREAL mean, MYREAL sd);
extern MYREAL find_chi (long df, MYREAL prob);
extern MYREAL probchiboundary (MYREAL chi, long zeros, long all);
extern MYREAL chiboundary (long zeros, long nonzeros, MYREAL alpha);
extern MYREAL probchi (long df, MYREAL chi);
extern MYREAL chisquare (long df, MYREAL alpha);
extern void gamma_rates (MYREAL *rate, MYREAL *probcat, long categs,
                             char *input);
#ifndef HAVE_LGAMMA
extern MYREAL mylgamma (MYREAL z);
#endif
extern MYREAL calc_sum (MYREAL *vector, long n);
extern MYREAL logfac (long n);

extern void onepass_mean_std_start(MYREAL *mean, MYREAL *std, long *n);
extern void onepass_mean_std_calc(MYREAL *mean, MYREAL *std, long *n, MYREAL x);
extern void onepass_mean_std_end(MYREAL *mean, MYREAL *std, long *n);
extern void combine_meanstd(MYREAL *mean1, MYREAL *sd1, long *n1, MYREAL mean2, MYREAL sd2, long n2);
extern boolean traverse_check_id(node *theNode, long id);
extern  void start_node_collection(world_fmt *world);
extern void stop_node_collection(world_fmt *world);
extern void collect_nodelet(world_fmt *world, node *p);
extern node *  dispense_nodelet(world_fmt *world);
extern void swap_node_collection(world_fmt * tthis, world_fmt * that);


extern float fast_log(float vval);
extern MYINLINE double fast_exp(double y); 


extern char lowercase (int c);
extern char uppercase (int c);
extern void upper(char *from, char **to);
extern long read_word_delim(char *input, char *word, char * delim, boolean ignore_repeats);

/* vector initialization and calc*/
extern void doublevec1d (MYREAL **v, long size);
extern void doublevec2d (MYREAL ***v, long size1, long size2);
extern void floatvec2d (float ***v, long size1, long size2);
extern void charvec2d (char ***v, long size1, long size2);
extern void intvec2d (long ***v, long size1, long size2);
extern void free_doublevec2d (MYREAL **v);
extern void free_floatvec2d (float **v);
extern void free_charvec2d (char **v);
extern void free_intvec2d (long **v);
extern void setdoublevec1d (MYREAL **v, MYREAL *w, long size);
extern void add_vector (MYREAL *result, MYREAL *v, long size);

/*filemanipulation */
extern void init_files (world_fmt * world, data_fmt * data,
                            option_fmt * options);
extern void exit_files (world_fmt * world, data_fmt * data,
                            option_fmt * options);
extern void openfile (FILE ** fp, char *filename, char *mode, char *perm);
#ifdef ZNZ
extern void znzopenfile (znzFile * fp, char *filename, char *mode, int use_compressed);
#endif
extern void get_filename(char **store, char * value);
extern long read_savesum (world_fmt * world, option_fmt * options,
                              data_fmt * data);
extern void write_savesum (world_fmt * world);

/* string manipulation */
extern void translate (char *text, char from, char to);
extern void unpad (char *text, char removechars[]);
extern void get_next_word(char **instring, char *delimiters, char **nextword);
extern long count_words (char *text);
extern long count_char (char *text, char needle);
extern long locate_char (char *text, char needle);
extern char * char_position(char *s1, char* s2);
extern char read_word(FILE *infile, char *word, char *delim);
extern void unread_word(FILE *infile, char *word);
extern void increase_buffer(char **buffer, long *allocbufsize, long *bufsize,long amount);
extern void add_to_buffer(char *fp, long *bufsize, char **buffer, long *allocbufsize);
extern long print_to_buffer(char **buffer, long *maxbufsize, char *tempbuffer, long *pos, const char *fmt, ...);
extern void record_warnings(world_fmt * world, const char *fmt, ...);
extern void print_stored_warnings(world_fmt *world);
extern void remove_trailing_blanks(char **filename);
extern void trim(char **line);
extern long nondestruct_trim(char **line);
extern void set_paramstr(char *paramstr, long j, world_fmt *world);
/* time reporting */
extern void get_time (char *nowstr, char ts[]);
extern void get_runtime (char *runtime, time_t start, time_t end);
/* printing aid */
extern void print_llike (MYREAL llike, char *strllike);
extern int sprint_tabdigits(double x, char * temp);
/* searching and finding*/
extern boolean find (long i, long *list, long listlen);
/* conversion between the different parameter schemes*/
extern long mstart (long pop, long numpop);
extern long mend (long pop, long numpop);
extern long mm2m (long frompop, long topop, long numpop);
extern void m2mm (long i, long numpop, long *frompop, long *topop);
extern void d2mm(long pop,world_fmt *world, long *frompop, long *topop);
extern long mm2d(long frompop, long topop,world_fmt *world);
extern long m2mmm(long frompop, long topop, long numpop);
extern long m2mml (long i, long numpop);
extern long m2mml2 (long i, long topop, long numpop);
extern long mmstart (long pop, long numpop);
extern long mmend (long pop, long numpop);
extern long mml2m (long pos, long numpop);
extern void print_line (FILE * outfile, char c, long nn, long flag);
extern void sprint_line (char *buffer, char c, long nn, long flag);
extern void fprintf2(FILE *file, long filesize, const char *fmt, ...);

/* reading from char * buffer */
extern char sgetc (char **buffer);
extern char *sgets (char *s, int size, char **stream);
extern char *sgets_safe (char **s, long *size, char **stream);

/* non-system erf functions */
extern double myerf(double x);
extern double myerfc(double x);
extern double erfinv(double x);
extern double inverse_cumstd_normal(double p);

extern double wew ( double x, double *en );

extern double get_time_for_growth(double theta0, double growth, double k, double t0);
extern double interval_growth(double r, double t0, double theta0, double growth, double k, double tmin, double tmax);
#ifdef CAUTIOUS
extern boolean cautious;
#endif

#endif /*TOOLS_INCLUDE */
