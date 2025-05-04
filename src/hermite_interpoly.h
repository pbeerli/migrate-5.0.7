#ifndef HERMITE_INTERPOL_H
#define HERMITE_INTERPOL_H
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

// Piece-wise interpolation
#include <complex.h>
#include "definitions.h"
extern MYCOMPLEX MLinterpol(MYCOMPLEX z, double alpha, double beta, double Q, double X0);
#endif
