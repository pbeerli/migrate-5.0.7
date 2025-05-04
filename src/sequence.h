#ifndef SEQUENCE_INCLUDE
#define SEQUENCE_INCLUDE
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 S E Q U E N C E S   R O U T I N E S 
 
 
 Peter Beerli 1996, Seattle
 beerli@fsu.edu
 
(c) 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
(c) 2003-2004 Peter Beerli, Tallahassee FL
 
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

$Id: sequence.h 2157 2013-04-16 22:27:40Z beerli $
-------------------------------------------------------*/

void initratio (option_fmt * options);
void initfreqs (MYREAL *freqa, MYREAL *freqc, MYREAL *freqg,
                    MYREAL *freqt);
void initcatn (long *categs);
boolean initcategs (long categs, MYREAL *rate, MYREAL *probcat);
void initprobcat (long categs, MYREAL *probsum, MYREAL *probcat);
void constrain_rates(long categs, MYREAL *rate, MYREAL *probcat);
void initlambda (option_fmt * options);

void init_sequences (world_fmt * world, option_fmt * options,
			    data_fmt * data, long locus);
void init_sequences2 (world_fmt * world, seqmodel_fmt * seq,
			     long locus);
void init_tbl (world_fmt * world, long locus);
void empiricalfreqs (world_fmt * world, option_fmt * options, mutationmodel_fmt * s, long sublocus);

void print_weights (FILE * outfile, world_fmt * world,
			   option_fmt * options, long locus);
void print_seqfreqs (FILE * outfile, world_fmt * world,
			    option_fmt * options);
MYREAL treelike_seq (mutationmodel_fmt *s, long sublocus, world_fmt * world, long locus);
MYREAL treelike_snp (mutationmodel_fmt *s, long sublocus, world_fmt * world, long locus);
void snp_invariants (mutationmodel_fmt *s, long sublocus, contribarr invariants, world_fmt *world, long locus,  phenotype x1, MYREAL *scale);
//void make_invarsites (world_fmt * world, data_fmt * data, long locus);
void make_invarsites_unlinked (world_fmt * world, data_fmt * data,
				      long locus);
void make_snp (world_fmt * world, option_fmt * options,
		      data_fmt * data, long locus);

MYREAL treelike_snp_unlinked (mutationmodel_fmt *s, long sublocus, world_fmt * world, long locus);
void copy_seq (world_fmt * original, world_fmt * kopie);
void init_sequences_aliases (world_fmt * world, option_fmt * options,
                                 data_fmt * data, long locus);
extern void find_rates_fromdata(data_fmt * data, option_fmt * options, world_fmt *world);
extern void find_rates_fromdata_alleles(data_fmt * data, option_fmt * options, world_fmt *world, MYREAL mean);
void free_seq(seqmodel_fmt **seq, long seqnum);
void set_nucleotide(MYREAL *treedata, const char nucleotide, const MYREAL * seqerr);
void check_basefreq (option_fmt * options);

#endif /*SEQUENCE_INCLUDE */
