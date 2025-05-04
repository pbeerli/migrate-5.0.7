#ifndef HAPLOTYPE_INCLUDE
#define HAPLOTYPE_INCLUDE
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 D A T A   R O U T I N E S 
 
 creates and manipulates data structures for haplotyping.
 
 Peter Beerli
 beerli@fsu.edu
 

Copyright 2010 Peter Beerli, Tallahassee FL
 
 This software is distributed free of charge for non-commercial use
 and is copyrighted. Of course, we do not guarantee that the software
 works and are not responsible for any damage you may cause or have.
 
*/
#include "migration.h"
#include "hash.h"

void init_individuals(data_fmt **data, long loci);
void insert_individual_inlist(char *input, data_fmt *data, option_fmt *options);
void link_individual_node(char *rawname, node *thenode, long region, world_fmt *world);
void reset_ID_nodelist(long locus, world_fmt * world);
void copy_individuals_from_data(data_fmt * data, worlddata_fmt * wdata, long locus);
void copy_individuals_differences(MYREAL heat, worlddata_fmt * ndata, worlddata_fmt * odata);
void add_state_counter(hash_fmt **hash, size_t *numhash, float *total, int *states, long numstates, int count);
void store_indname(char *name, long pop, long ind, long region, data_fmt *data);
//long calc_individual_checksum(int *difference, long allsubloci);
long calc_individual_checksum(int *difference, world_fmt * world, long locus);
void check_individual_nodes(long region, world_fmt *world);
void report_states(char *buffer, long *bufsize, world_fmt *world, long locus, boolean final, long numreplicates);
long swap_haplotypes(world_fmt *world);

void read_individuals_dictionary(data_fmt * data, option_fmt *options);
void  set_individuals_request_haplotyping(world_fmt * world, data_fmt * data, long locus);
void cleanup_individual_nodes(long region, world_fmt *world);
void print_haplotypes2(char **buffer, long *allocbufsize, long *bufsize, world_fmt *world, long locus, long pop, long ind, boolean fraction);
void print_haplotypes(char **buffer, long *allocbufsize, long *bufsize, world_fmt* world, data_fmt *data);
void print_haplotype_stat(world_fmt *world, data_fmt *data);
void get_haplotypes (world_fmt * world, option_fmt * options);
void save_haplotypes(world_fmt* world, long locus);
void reset_haplotypes(world_fmt* world, long locus);
long set_haplotype_states(individualDB_fmt *winner, char *word);
long find_inDB(char * thisname, long region, individualDB_fmt *individuals, long numindividuals);
void free_haplotypes(world_fmt *world, long locus);
#endif
