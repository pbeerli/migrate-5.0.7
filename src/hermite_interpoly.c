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

// interpolation algorithm is based on the algorithm 
// using http://paulbourke.net/miscellaneous/interpolation/
// His Copyright : 	
// The contents of this web site are © Copyright Paul Bourke or a
// third party contributor where indicated. You may print or save an
// electronic copy of parts of this web site for your own personal use.
// Permission must be sought for any other use. Any source code found
// here may be freely used provided credits are given to the author.
// Purchase of credit free licenses of material found on this site can be
// negotiated with the author. The author can also quote for unique
// variations and/or higher resolution versions of the images found on
// this site.
//
#include "definitions.h"
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include "hermite_interpoly.h"
#include "mittag_leffler.h"
#include "mittag_leffler_interpol_data.h"


long find_table_value(double mu, double *values, long n);
double linear_interpolate(double y1, double y2, double mu);
double interpoly(double x, double alpha, double beta, double z[NCOLS], double tables[NCOLS], long ncols, double Q, double X0);

/*
   Tension: 1 is high, 0 normal, -1 is low
   Bias: 0 is even,
         positive is towards first segment,
         negative towards the other
*/
#include "hermite_interpoly.h"
// mu is the value to interpolate
// y0 to y4 are values from the table
// mu is between y1 and y2

MYCOMPLEX MLinterpol(MYCOMPLEX z, double alpha, double beta, double Q, double X0);

double hermite_interpolate(
   double y0,double y1,
   double y2,double y3,
   double x0,double x1,
   double x2,double x3,
   double mu,
   double tension,
   double bias);


double hermite_interpolate(
   double y0,double y1,
   double y2,double y3,
   double x0,double x1,
   double x2,double x3,
   double mu,
   double tension,
   double bias)
{
   double m0,m1,mu2,mu3;
   double a0,a1,a2,a3;
   double plus = (1+bias)*(1-tension)/2;
   double minus = (1-bias)*(1-tension)/2;
   double s0 = (x2 - x1)/(x1 - x0);
   double s3 = (x2 - x1)/(x3 - x2);

   mu2 = mu * mu;
   mu3 = mu2 * mu;
   m0  = s0*(y1-y0)*plus;
   m0 += (y2-y1)*minus;
   m1  = (y2-y1)*plus;
   m1 += s3*(y3-y2)*minus;
   a0 =  2*mu3 - 3*mu2 + 1;
   a1 =    mu3 - 2*mu2 + mu;
   a2 =    mu3 -   mu2;
   a3 = -2*mu3 + 3*mu2;

   return(a0*y1+a1*m0+a2*m1+a3*y2);
}

double linear_interpolate(double y1, double y2, double mu)
{
   return(y1 * (1-mu) + y2*mu);
}


long find_table_value(double mu, double *values, long n)
{
  long first = 0;
  long last  = n-1;
  long index = ( last + first ) /2;
  while (last - first > 1)
    {
      if (mu < values[index])
	last = index;
      else
	first = index;
      if (fabs(mu - first) < PRECISION)
	return index;
      index = (last + first) / 2;
    }
  return index;
}


double interpoly(double x, double alpha, double beta, double z[NCOLS], double tables[NCOLS], long ncols, double Q, double X0)
{
  double y0, y1, y2, y3;
    double x0, x1, x2, x3;
  // we assume that nrows is always 99
  // tables has 99 rows for alpha 0.01 ... 0.99
  
  if (x < z[0])
    {
      /* far tail (|z| > 2e6): asymptotic series
	 E_{a,b}(z) ~ -sum_k z^{-k} / Gamma(b - a k); terms with
	 b - a k at a pole of Gamma vanish. The former fallback called ML()
	 with Q and X0 swapped and returned log(0) here. */
      double sum = 0.0, g;
      long k;
      for (k = 1; k <= 4; k++)
	{
	  g = beta - alpha * (double) k;
	  if (g <= 0.0 && fabs(g - floor(g + 0.5)) < 1e-12)
	    continue;
	  sum -= pow(x, (double) -k) / tgamma(g);
	}
      return (sum > 0.0) ? log(sum) : -(double) HUGE;
    }
  if (x>z[ncols-1])
    {
#ifdef WINDOWS
      MYCOMPLEX xx = {x, 0.0};
      MYCOMPLEX b = ML(xx, alpha, beta, Q, X0);
#else
      MYCOMPLEX b = ML(x, alpha, beta, Q, X0);
#endif
      if (creal(b) <= DBL_MIN)
	return -(double) HUGE;
      else
	return log(creal(b));
    }
  
  long index = find_table_value(x,z,ncols);
  if (fabs(x-z[index])<PRECISION)
    return tables[index];
  if (index == 0)
    {
      y0 = tables[index];
      y1 = tables[index+1];
      return linear_interpolate(y0, y1, (x-z[index])/(z[index+1]-z[index]));
    }
  if (index == ncols-2)
    {
      y0 = tables[index];
      y1 = tables[index+1];
      return linear_interpolate(y0, y1, (x-z[index])/(z[index+1]-z[index]));
    }
  y0 = tables[index-1];
  y1 = tables[index];
  y2 = tables[index+1];
  y3 = tables[index+2];
  x0 = z[index-1];
  x1 = z[index];
  x2 = z[index+1];
  x3 = z[index+2];
  return hermite_interpolate(y0, y1, y2, y3, x0, x1, x2, x3, (x-z[index])/(z[index+1]-z[index]), 0.0, 0.0);
}

// mlf_lookup and mlf_lookup1 and z_lookup are global static tables
/* log E from table row r (beta = alpha of that row, or beta = 1) */
static double mlrow(MYCOMPLEX z, long r, boolean aa, double Q, double X0)
{
  if (aa)
    return interpoly(creal(z), alpha_lookup[r], alpha_lookup[r], z_lookupaa[r], mlf_lookupaa[r], nz_lookupaa[r], Q, X0);
  return interpoly(creal(z), alpha_lookup[r], 1.0, z_lookupa1[r], mlf_lookupa1[r], nz_lookupa1[r], Q, X0);
}

MYCOMPLEX MLinterpol(MYCOMPLEX z, double alpha, double beta, double Q, double X0)
{
  /* log E_{alpha,beta}(z) for beta==alpha or beta==1 from the tables; the
     table rows have their own alpha values (alpha_lookup, finer near 1) and
     values in between are interpolated linearly in alpha (formerly alpha was
     truncated to two digits) */
  double b, w;
  long lo, hi, mid;
  const boolean aa = (fabs(alpha-beta)<PRECISION);
  if (fabs(alpha-1.0)<PRECISION && fabs(beta-1.0)<PRECISION)
    {
      return z;//we return the log of exp(z) (for alpha=1 beta=1)
    }
  if (alpha <= alpha_lookup[0])
    lo = hi = 0;
  else if (alpha >= alpha_lookup[NROWS-1])
    lo = hi = NROWS-1;
  else
    {
      lo = 0; hi = NROWS-1;
      while (hi - lo > 1)
	{
	  mid = (lo + hi) / 2;
	  if (alpha_lookup[mid] <= alpha) lo = mid; else hi = mid;
	}
    }
  /* interpolate in u = log(1 - alpha): the tail depends on alpha through
     log Gamma(1-alpha) ~ -log(1-alpha), which is nearly linear in u but
     bends sharply in alpha near 1 (interpolating in alpha left errors of
     1e-2 .. 6e-2 in log E there) */
  if (hi == lo)
    b = mlrow(z, lo, aa, Q, X0);
  else if (lo >= 1 && hi <= NROWS-2)
    {
      long r[4] = {lo-1, lo, hi, hi+1};
      long i, j;
      double lw, u = log(1.0 - alpha), ur[4];
      for (i = 0; i < 4; i++)
	ur[i] = log(1.0 - alpha_lookup[r[i]]);
      b = 0.0;
      for (i = 0; i < 4; i++)
	{
	  lw = 1.0;
	  for (j = 0; j < 4; j++)
	    if (j != i)
	      lw *= (u - ur[j]) / (ur[i] - ur[j]);
	  b += lw * mlrow(z, r[i], aa, Q, X0);
	}
    }
  else
    {
      w = (log(1.0 - alpha) - log(1.0 - alpha_lookup[lo]))
	/ (log(1.0 - alpha_lookup[hi]) - log(1.0 - alpha_lookup[lo]));
      b = (1.0 - w) * mlrow(z, lo, aa, Q, X0) + w * mlrow(z, hi, aa, Q, X0);
    }
#ifdef WINDOWS
#ifdef NMAKE
  return {b, 0.0};
#else
  return (_Dcomplex) {b, 0.0};
#endif
#else
  return (MYCOMPLEX) b;
#endif
}


#ifdef INTERPOLYTESTER
int main(int argc, char **argv)
{
  double x = atof(argv[1]);
  double alpha = atof(argv[2]);
  double beta = alpha;
  double logx = interpoly(x, alpha, beta, z_lookup, mlf_lookup, NCOLS);
  double xx = exp(logx);
  printf("log(mlf(alpha=%f, beta==%f x=%f)) = %f --> %f\n",alpha, beta,x,
	 logx, xx);

}
#endif
