#ifndef EVENTUPDATE_H
#define EVENTUPDATE_H
/*------------------------------------------------------
Maximum likelihood estimation
of migration rate  and effectice population size
using a Metropolis-Hastings Monte Carlo algorithm
-------------------------------------------------------
    gathering migration and coalescence events  routines   R O U T I N E S

    Peter Beerli 2006, Tallahassee
    beerli@fsu.edu

    Copyright 2006 Peter Beerli, Tallahassee

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


$Id$
  */
#include "migration.h"
extern void increase_mighist (mighistloci_fmt * mighistlocus);
extern void setup_mighist (world_fmt * world, option_fmt * options);
extern void destroy_mighist (world_fmt * world);
extern void calculate_event_values(duo **eventbins, long *eventbinnum, 
				   MYREAL eventinterval, MYREAL interval, 
				   MYREAL age,  char type, long from, long to, 
				   long * lineages, long numpop, boolean is_last, world_fmt *world);
extern void setup_event_events (world_fmt * world, option_fmt * options);
extern void destroy_event_events (world_fmt * world);
extern void print_mighist_output (FILE * out, world_fmt * world, MYREAL *sums, boolean mrca);
extern void print_event_values(world_fmt * world);
extern void store_events (world_fmt * world, timelist_fmt * ltl, long np, long rep);

#endif /*_EVENTUPDATE_*/
