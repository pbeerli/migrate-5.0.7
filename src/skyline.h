#ifndef SKYLINEUPDATE_H
#define SKYLINEUPDATE_H
/*------------------------------------------------------
    variation over time routines   R O U T I N E S

    Peter Beerli 2006, Tallahassee
    beerli@fsu.edu

    (c) 2006 Peter Beerli, Tallahassee

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
extern void calculate_expected_values(tetra **eventbins, long *eventbinnum, MYREAL eventinterval, MYREAL interval, MYREAL age, char type, long from, long to, long * lineages, long numpop, world_fmt *world);
extern void setup_expected_events (world_fmt * world, option_fmt * options);
extern void destroy_expected_events (world_fmt * world);
extern void print_expected_values(world_fmt * world, option_fmt *options);
extern void debug_skyline(world_fmt *world, char text[]);
#endif /*_SKYLINEUPDATE_*/
