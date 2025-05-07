/*
 *  pretty.h
 *  part of migrate
 *  will be PDF printout system
 *  Created by Peter Beerli on 7/25/05.
 *  (c) 2005 Peter Beerli. All rights reserved.
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
 *
 */
#ifdef PRETTY
#include "migration.h"
#include "sighandler.h"
#include "tools.h"
#ifdef BEAGLE
#include "calculator.h"
#endif
#undef USE_ENCRYPTION
#define HAVE_BOOLEAN 
#include <hpdf.h>
//#include "haru/libharuc.h"
#define SEARCHUP 0
#define SEARCHDOWN 1
// GENERAL
// initialize PDF structures
extern void pdf_master_init(world_fmt *world, option_fmt *options, data_fmt *data);
extern int pdf_init(void);
// print title on every new page
extern int pdf_new_page(char *title);
// print master title on first page 
extern int pdf_master_title(char *title);

extern void pdf_print_end_time(double *page_height);
// OPTIONS
// print options onto the first page
extern void pdf_print_options(world_fmt * world, option_fmt *options, data_fmt * data);
#ifdef BEAGLE
extern void print_beagle_resources_pretty(BeagleInstanceDetails instDetails);
#endif
extern void pdf_print_mutationrate_weights(MYREAL *murates, long *segregs, MYREAL *wattersons, long loci);
// DATA
// print data summary
extern void pdf_print_data_summary(world_fmt * world, option_fmt *options, data_fmt * data, 
                              double *orig_page_height, double *orig_left_margin);
extern void pdf_print_random_subset(data_fmt * data, option_fmt *options);
// print data
extern void pdf_print_data (world_fmt * world, option_fmt * options, data_fmt * data);
extern void pdf_print_spectra(world_fmt *world, data_fmt *data, option_fmt *options, MYREAL ***freq, MYREAL ** total, MYREAL *grandtotal, MYREAL *avg_het, MYREAL avghet1, long * maxalleles);
extern void pdf_print_averageheat(world_fmt **universe, option_fmt *options);

// BAYES OUTPUTS
// print bayes table
extern void pdf_print_bayestable(world_fmt *world);
// print histogram at location lx ly with widht and height
extern void pdf_histogram(double *binvals, char *set50, char *set95, long bins, double bindelta, double binmin, double binmax, double lx, double ly, double width, double height, boolean nofreq, MYREAL *priors);
extern void pdf_histogram_plus(double *binvals, MYREAL *std, char *set50, char *set95, long bins, double bindelta, double binmin, double binmax, double lx, double ly, double width, double height, MYREAL scaler, boolean nofreq, world_fmt * world, double *confidence, long topop);
// print histogram for each locus and overall loci
extern double pdf_locus_histogram(world_fmt *world, long locus);
extern double pdf_loci_histogram(world_fmt *world);
// print acceptance ratios
extern void pdf_bayes_print_accept(world_fmt *world);
extern void pdf_bayes_print_hyperpriors(world_fmt *world);
// print ESS and autocorrelation
extern void pdf_bayes_print_ess(world_fmt *world);
extern void pdf_bayes_factor_header(world_fmt *world, option_fmt *options);
extern void pdf_bayes_factor_rawscores_header(world_fmt *world, option_fmt *options);
extern void pdf_bayes_factor_rawscores(long locus, MYREAL rawtermo, MYREAL beziertermo, MYREAL ss, MYREAL harmo);
extern void pdf_bayes_factor_rawscores_minmax(int minmax, MYREAL rawtermo, MYREAL beziertermo, MYREAL harmo);
void pdf_bayes_factor_rawscores_harmo(long locus, MYREAL harmo);
extern void pdf_bayes_factor_comment(world_fmt *world, MYREAL scaling_factor);
extern void pdf_burnin_stops(world_fmt *world, long maxreplicate);
// MIGHIST output
//extern void 
//pdf_mig_histogram(histogram_fmt ** histogram,
//                  plotfield_fmt ** plotfield, long loci, long nmigs,
//                  long bins, long *sum, MYREAL ***migtable,  boolean precalc, world_fmt *world);
//extern void pdf_print_mighist_table (world_fmt * world, MYREAL ***migtable,
extern void pdf_print_event_table (world_fmt * world, float *migtable);
extern void pdf_print_eventtime_table(world_fmt *world);
extern void pdf_print_time_table(world_fmt *world, 
				 double *meantime, double *mediantime, double *stdtime, double *freq,
				 boolean mrca);
extern void pdf_event_histogram(long loci, long numparams,  world_fmt *world);
// SKYLINE plot
extern void pdf_skyline_histogram(long loci, long numparams,  world_fmt *world, boolean enlarged);

// ML OUPUTS
// print mcmc table
//extern void pdf_print_mcmctable(float *page_height, float page_width, world_fmt *world, data_fmt *data, option_fmt *options);
extern void pdf_print_results (world_fmt ** universe, option_fmt * options, data_fmt * data);
extern void pdf_print_profile_title(world_fmt *world);
extern void pdf_print_profile_table(float left_margin, float * page_height, char method, char *buffer, world_fmt *world);
extern long prepare_profile_table(char method, long which, MYREAL *likes, MYREAL **param, world_fmt *world, long bufsize);
extern void pdf_print_profile_percentile (world_fmt * world);

extern void pdf_histogram_legend(void);
extern void pdf_print_lrt_header(MYREAL aicfull, long aicfullparamnum,  float *page_width,  float *page_height);
extern void pdf_print_LRT_box(world_fmt* world,char **box,long size,
			      MYREAL like0, MYREAL like1,MYREAL lrt, MYREAL chiprob, MYREAL chiprob2, MYREAL aic, long df, 
			      long aicparamnum, float *page_height, float *page_width);

// haplotyping 
extern void pdf_print_haplotype_title(void);
extern void pdf_print_haplotypes(world_fmt *world);
extern void pdf_print_stored_warnings(world_fmt *world);
extern void pdf_print_haplotypes2(char *buffer, long buflen);

// assignment
extern void pdf_report_unassigned(world_fmt *world);

// write PDF content to file
extern int pdf_write_file(option_fmt *options);

extern  HPDF_Doc doc;
extern  HPDF_Page page;
extern  HPDF_Page canvas;
extern  HPDF_Font helv;
extern  HPDF_Font helvob;
extern  HPDF_Font symb;
extern HPDF_Font cour;
extern  int page_counter;
extern  char pdf_pagetitle[LINESIZE+1];
extern  char pdf_time[LINESIZE+1];
extern  double page_height;
extern  double left_margin;
extern int numcpu;
extern HPDF_Page firstcanvas;
extern  double timestampy;


#endif
