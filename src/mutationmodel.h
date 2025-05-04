#ifndef MUTATIONMODEL_H
#define MUTATIONMODEL_H
#include "migration.h"
#include "sighandler.h"
/*
 
 (c) Peter Beerli 2013 Tallahassee FL
 
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
*/

///
/// initialize the mutation model structure
extern void init_mutationmodel_second(world_fmt *world, data_fmt *data);
extern void init_mutationmodel_first(world_fmt *world, data_fmt *data);
extern void init_mutationmodel_readsites(mutationmodel_fmt *mumod, char datatype, char *  sites);
extern void init_mutationmodel_readsites2(mutationmodel_fmt *mumod, char datatype, long numsites);
extern void init_mutationmodel_readsites3(mutationmodel_fmt *mumod, char datatype, long numsites); //long sites)
extern void finish_mutationmodel(world_fmt *world, data_fmt *data, option_fmt *options, long locus);
extern void print_mutationrate_weights(FILE *file, MYREAL *murates, long *segregs, MYREAL *wattersons, long loci);
// klone for heating 
extern void klone_mutationmodel(world_fmt *newcopy, world_fmt *original, data_fmt *data, long locus);
void read_mutationmodel_comments(char *input, data_fmt * data, world_fmt * world);
extern int get_states(mutationmodel_fmt *s, data_fmt *data, long locus);
//extern void set_subloci_basefrequencies(mutationmodel_fmt *mumod, world_fmt *world, option_fmt *options, data_fmt *data, long sublocus);
extern void set_subloci_frequencies_alleles(world_fmt *world, option_fmt *options, data_fmt *data, long locus);
extern void set_subloci_frequencies(world_fmt *world, option_fmt *options, data_fmt *data, long locus);
extern void set_subloci_basedefaults(mutationmodel_fmt *s, world_fmt *world, option_fmt *options, data_fmt *data, long sublocus);
extern void get_mutationmodel_nameparam(char *modelname, char *modelparams, mutationmodel_fmt *s);
extern void set_siterates(long z, world_fmt *world, option_fmt *options);
extern void force_basefreqs(MYREAL ** basefreqs, MYREAL pA, MYREAL pC, MYREAL pG);
extern void nuview_tn93 (mutationmodel_fmt *s, long sublocus, long xs, node * mother, world_fmt * world, const long locus);
extern void pseudonu_tn93 (mutationmodel_fmt *s, proposal_fmt *proposal, xarray_fmt *xxx1, MYREAL *lx1, MYREAL v1, xarray_fmt *xxx2, MYREAL *lx2, MYREAL v2, long xs);
extern int get_mutationmodel(char x);
extern void destroy_mutationmodel(world_fmt* world);
#endif
