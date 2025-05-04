#ifndef __CALCULATOR__
#define __CALCULATOR__
/*
 
(c) Peter Beerli Tallahassee 2013

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
#include "beagle.h"
#include "migration.h"
void print_beagle_available_resources(world_fmt *world);
void reset_beagle(beagle_fmt *beagle);
void init_beagle(world_fmt *world, long locus);
void reinit_beagle(world_fmt *world, long locus);
void change_beagle_scalingIndex(beagle_fmt *beagle);
void change_beagle(node *theNode, beagle_fmt *beagle, long sumtips);
void  set_beagle_instances(world_fmt *world, long locus);
void prepare_beagle_instances(node *theNode, node * left, node *right, beagle_fmt *beagle);
void prepare_beagle_instances_proposal(proposal_fmt *proposal, long trueparentid, 
				       long leftid, long leftbid, double leftbranch, 
				       long rightid, long rightbid, double rightbranch, beagle_fmt *beagle);
void fill_beagle_instances(world_fmt *world, long locus);
void evaluate_beagle_instances_proposal (proposal_fmt * proposal,
				    node * mother,  
				    node * newdaughter, long newdaughter_id, long newdaughter_bid, 
					 MYREAL v);
double calcLnL(world_fmt *world, int scalingFactorIndex);
double force_beagle_recalculate(world_fmt *world, long locus);
void beagle_stop(world_fmt **universe, long usize);
long set_branch_index (node * p,  long *bid);
long new_id(long id, long sumtips);
#endif