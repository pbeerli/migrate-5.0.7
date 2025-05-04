#ifndef CORRELATION_
#define CORRELATION_
/* Correlation table of all parameters for Bayesian inference
 and using the profile likelihood tables
 (1) Bayesian inference: read parameter list and calculate correlation table
 (2) ML inference: read profile table and construct correlation table
     [uses print_cov .... already present, but needs testing]
 

 (c) Peter Beerli 2012

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
#include "definitions.h"
#include "migration.h"

extern void covariance_bayes(world_fmt *world, long locus);
//extern void covariance_bayes2(world_fmt *world, long locus, MYREAL *params);
extern void covariance_summary(world_fmt *world);
extern void adjust_covariance(MYREAL **cov, long size, long n);
#endif
