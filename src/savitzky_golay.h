/*------------------------------------------------------
 Bayesian inference of
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm
 -------------------------------------------------------
 smoothing routines

 The Savitzky-Golay smoother implementation used the resources of
 a pascal code Peter A. Gorry: 
 General least-squares smoothing and differentiation by the convolution (Savitzky-Golay) method.
 Analytical Chemistry 1990, 62, 6, 570–573
 and its extension by the machine learning java script repo: 
 https://github.com/mljs/savitzky-golay-generalized

 PB translated the Pascal code and and the JS code to python (see savitzky_golay.py) and that to C.

 Send questions concerning this software to:
 Peter Beerli
 beerli@fsu.edu
 
 Copyright 2020 Peter Beerli
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 of the Software, and to permit persons to whom the Software is furnished to do
 so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all copies
 or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 
 */

extern double savitzky_golay(double *x, long xelem, double *xvals, long windowsize, long derivative, long polynomial);
