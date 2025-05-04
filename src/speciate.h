/* speciate.h */
// contains declarations for speciate.c
//
#ifndef SPECIATE_H
#define SPECIATE_H
/*
  Copyright (c) 2015 Peter Beerli

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  THE SOFTWARE.
*/


extern MYREAL log_prob_wait_speciate_weibull(MYREAL t0, MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt *s);
extern MYREAL log_prob_wait_speciate_normal(MYREAL t0, MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt *s);
extern MYREAL log_prob_wait_speciate_normalorig(MYREAL t0, MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt *s);
extern MYREAL log_prob_wait_speciate_exp(MYREAL t0, MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt *s);
extern MYREAL log_point_prob_speciate_weibull(MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt * s);
extern MYREAL log_point_prob_speciate_normal(MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt * s);
extern MYREAL log_point_prob_speciate_normalorig(MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt * s);
extern MYREAL log_point_prob_speciate_exp(MYREAL t1, MYREAL mu, MYREAL sigma, species_fmt * s);
extern MYREAL eventtime_single(proposal_fmt *proposal, world_fmt *world, long pop, long timeslice, long *lineages, double age, char * event, long *to, long *from);
extern node * set_type2(world_fmt *world, node *p, node *q, char *custm2);
extern char set_type(world_fmt *world, long topop, long frompop, char *custm2, long numpop);
extern long newtree_update (world_fmt * world, long g, boolean assign);
extern long speciation_from(long to, proposal_fmt * proposal);
extern void loopcleanup(boolean assign, world_fmt * world, proposal_fmt *proposal, long oldpop, long newpop, timelist_fmt *timevector);
extern long init_speciesvector(world_fmt * world, option_fmt *options);
extern void fill_speciesvector(world_fmt*world, option_fmt *options);
extern long propose_new_spec_mu_old(world_fmt * world, MYREAL *oldmu, MYREAL *oldsigma, MYREAL *newmu, MYREAL *newsigma, long which, boolean *is_mu);
extern void adjust_averagediv(world_fmt * world, long which, MYREAL *newparam);
extern long propose_new_spec_mu(world_fmt * world, long which, boolean *is_mu, MYREAL *newparam);
extern MYREAL wait_event_species(world_fmt *world, vtlist *tli, MYREAL t0, MYREAL t1, boolean waitonly, MYREAL *eventprob);
extern MYREAL wait_D(long pop, MYREAL t0, MYREAL t1, long *lineages, world_fmt *world);
extern void set_first_speciestree(node *mrca, world_fmt *world);
extern void species_datarecorder(world_fmt *world);
extern species_fmt * get_which_species_model( long which, species_fmt * s ,  long ssize);
extern species_fmt * get_fixed_species_model( long from, long to, species_fmt * s ,  long ssize);
extern species_fmt * get_species_model( long to, species_fmt * s ,  long ssize);
extern void species_datarecorder(world_fmt *world);
extern void construct_locusspecies_histogram(world_fmt *world, long locus, MYREAL *mini, MYREAL *maxi, float **results);
extern void read_species_record(long masternumbinall, long locus, MYREAL *params, long *n, MYREAL *oldmeans, MYREAL *lowerbound, MYREAL *upperbound, MYREAL *delta, world_fmt *world, char **inptr);

extern void record_parameters(world_fmt *world);

extern log_prob_wait_speciate_func log_prob_wait_speciate;
extern log_point_prob_speciate_func log_point_prob_speciate;
extern time_to_speciate_func time_to_speciate;



#if defined(MPI) && !defined(PARALIO) /* */
void print_species_record(float *temp, long *z, world_fmt *world);
#else
void  print_species_record(char *temp, long *c, world_fmt * world);
#endif
void set_speciate_functions(world_fmt *world);
#endif
