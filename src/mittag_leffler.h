#ifndef MITTAGLEFFLER
#define MITTAGLEFFLER
/*-----------------------------------------------------------------
  Bayesian inference of population genetic forces: drift, migraiton, divergence
  allowing for the n-coalescent, the f-coalescent, and the BSC-coalescent
 
  Peter Beerli
  Department of Scientific Computing
  Florida State University
  Tallahassee FL 32306-4120
  beerli@fsu.edu
 
  Copyright 2017 Peter Beerli, Tallahassee FL

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
*-----------------------------------------------------------------
*/

#include <complex.h>
#include "definitions.h"
#include "migration.h"

extern MYCOMPLEX mittag_leffler(double alpha, double beta, MYCOMPLEX z);
extern double interval_mittag_leffler(double r, double alpha, double lambda, double tmin, double tmax);
extern MYCOMPLEX  ML(MYCOMPLEX z, double alpha, double beta, double Q, double X0);
extern void set_mittag_leffler(option_fmt * options);
extern double interval_mittag_leffler_func(double r, double alpha, double t0, double mu, double sigma, species_fmt *s, double tmin, double tmax);
extern double propose_new_mlftime(double lambda, double alpha, double r1, double r2);
extern void change_mittag_leffler(world_fmt* world);
#endif

