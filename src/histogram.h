// histogram.h
#ifndef BAYESREREAD_H
#define BAYESREREAD_H
/*
 started October 2007
 (c) Peter Beerli Tallahassee 2007
 $Id:$

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

#include "migration.h"
#ifdef ZNZ
extern void read_from_bayesmdim_minimal_info(znzFile mdimfile, world_fmt *world,option_fmt *options, data_fmt *data);
extern void read_bayes_fromfile(znzFile mdimfile, world_fmt *world, option_fmt *options, char **files, long fnum);
#else
extern void read_from_bayesmdim_minimal_info(FILE *mdimfile, world_fmt *world,option_fmt *options, data_fmt *data);
extern void read_bayes_fromfile(FILE *mdimfile, world_fmt *world, option_fmt *options, char **files, long fnum);
#endif
extern boolean checking_bayesallfile(world_fmt *world, option_fmt *options, data_fmt *data, long ***unfinished);

#endif
