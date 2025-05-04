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

#include <assert.h>
#include <stdlib.h>
#include "definitions.h"
#include "migration.h"
#include "options.h"
#include "bayes.h"
#include "random.h"
#include "tools.h"
#include "sighandler.h"
#include "world.h"
#include "tree.h"
#include "slice.h"
#include "reporter.h"
#include "haplotype.h"
#include "mcmc.h"
#include "skyparam.h"
#include "speciate.h"
#include "correlation.h"
#include "autotune.h"
#include "priors.h"

#include "savitzky_golay.h"

double savitzky_golay(double *x, long xelem, double *xvals, long windowsize, long derivative, long polynomial);
void fullweights(double ***weights, long m, long n, long s);
double sg_weight(long i, long t, long m, long n, long s);
double genfact(long a, long b);
double grampoly(long i,long m, long k,long s);
double getHs(double *xvals, long numxvals, long center, long half, long derivative);




// savitzky_golay(testdata,x,windowsize=150,derivative=0,polynomial=5)
double savitzky_golay(double *x, long xelem, double *xvals,  long windowsize, long derivative, long polynomial)
{
  double mse;
  long np1 = xelem;
  double *xr = (double *) mycalloc(xelem,sizeof(double));

  if (windowsize %2 == 0)
    windowsize += 1;

  if (windowsize<5)
    windowsize = 5;
  
  if (windowsize>xelem)
    {
      printf("windowsize is too big, it is bigger than number of data elements");
      exit(-1);
    }
  if (polynomial>=6)
    printf("Savitzky-Golay smooth: You should not use polynomial grade higher than 5");

  long half = (long) windowsize / 2;
  double **weights;
  doublevec2d(&weights,windowsize, windowsize);
  fullweights(&weights, windowsize, polynomial, derivative);
  double hs=0;
  double h=1.0;
  // currently we assume that all histogram cells are filled (does *x has no empty cell
  // this is incorrect but close for long runs
  boolean constantH = TRUE;

  if (xvals[2] > (double) -HUGE)
    constantH = FALSE;
  
  if (constantH)
    {
      h = xvals[1]-xvals[0]; //if we only look at smoothing and
      //not derivatives h does not matter h**derivative=h^0=1
      hs = pow(h,derivative);
    }
  // for the border points
  long i;
  long l;
  double d1, d2;
  for (i=0; i < half; i++)
    {
      //double *wg1 = weights[half - i - 1];
      //double *wg2 = weights[half + i + 1];
      d1 = 0;
      d2 = 0;
      for (l=0;l<windowsize;l++)
	{
	  //d1 += wg1[l] * x[l];
	  //d2 += wg2[l] * x[xelem - windowsize + l];
	  d1 += weights[half-i-1][l] * x[l];
	  d2 += weights[half+i+1][l] * x[xelem - windowsize + l];
	  if (constantH)
	    {
	      xr[half - i - 1] = d1 / hs;
	      xr[np1 - half + i] = d2 / hs;
	    }
	  else
	    {
	      hs = getHs(xvals, xelem, half - i - 1, half, derivative);
	      xr[half - i - 1] = d1 / hs;
	      hs = getHs(xvals, xelem, np1 - half + i, half, derivative);
	      xr[np1 - half + i] = d2 / hs;
	    }
	}
    }
  // internal points
  //double *wg = weights[half];

  for (i=windowsize;i<xelem+1;i++)
    {
      double d = 0.0;
      //..printf ("%f ",wg[l]);fflush(stdout);
      for (l=0;l<windowsize;l++)
	{
	  //d += wg[l] * x[l + i - windowsize];
	  d += weights[half][l] * x[l + i - windowsize];
	}
      if (!constantH)
	{
	  hs = getHs(xvals,xelem,i-half-1,half,derivative);
	}
      xr[i-half-1] = d / hs;
    }
  mse = 0;
  for (i=0;i<xelem;i++)
    {
      double diff = x[i]-xr[i];
      mse += diff * diff;
	
      if (xr[i]<0.0)
	{
	  x[i]=0.0;
	}
      else
	{
	  x[i] = xr[i];
	}
    }
  myfree(xr);
  return sqrt(mse);
}

void fullweights(double ***weights, long m, long n, long s)
{
  long t,j;
  long np1 = (long) m / 2;
  for (t = -np1; t < (np1+1); t++)
    {
      for (j = -np1; j < (np1+1); j++)
	(*weights)[t + np1][j + np1] = sg_weight(j, t, np1, n, s);
    }
}

double sg_weight(long i, long t, long m, long n, long s)
{
  // calculates the weith of the ith data point for the tth least square
  // point of the sth derivative, over 2m+1 points, order n
  long k;
  double sum = 0.0;
  
  for (k=0; k < n + 1; k++)
    {
      sum += (2.0*k+1.0)*(genfact(2*m,k)/genfact(2*m+k+1,k+1))*grampoly(i,m,k,0)*grampoly(t,m,k,s);
    }
  //printf("%f\n",sum);
  return sum;
}

double genfact(long a, long b)
{
  long c = 1;
  long j = a-b+1;
  if (a>=b)
    {
      while(j<=a)
	{
	  c *= j;
	  j++;
	}
    }
  return (double) c;
}

double grampoly(long i,long m, long k,long s)
{
  // calculates the chebyshev/gram polynomial (s=0) for its sth
  // derivative evaluated at i, order k, over 2m+1 points
  double g = 0;
  if (k>0)
    {
      g = (4.0*k-2.0)/(k*(2.0*m-k+1.0))*(i * grampoly(i,m,k-1,s)+ \
				 s*grampoly(i,m,k-1,s-1)) - \
	(((k - 1.0) * (2.0 * m + k)) / (k * (2.0 * m - k + 1.0)))*grampoly(i,m,k-2,s);
    }
  else
    {
      if ((k==0) && (s==0)) 
	g = 1.0;
      else
	g = 0.0;
    }
  return g;
}

double getHs(double *xvals, long numxvals, long center, long half, long derivative)
{
  //assumes h is list
  double hs = 0.0;
  
  long i;
  long count = 0;
  for (i=center - half; i <  center + half; i++)
    {
      if ((i >= 0) && (i < numxvals - 1))
	{
            hs += xvals[i + 1] - xvals[i];
            count+=1;
	}
    }
  return pow((hs / count), (double) derivative);
}



#ifdef  STANDALONE
// compile as gcc -o savgo savitzky_golay.c -DSTANDALONE
#include "savitzky-golay-helper.h"
void doublevec2d (MYREAL ***v, long size1, long size2);

/// allocate a 2D array with MYREALS
void
doublevec2d (MYREAL ***v, long size1, long size2)
{
    long i;
    *v = (MYREAL **) mycalloc (size1, sizeof (MYREAL *));
    (*v)[0] = (MYREAL *) mycalloc ((size1 * size2), sizeof (MYREAL));
    for (i = 1; i < size1; i++)
    {
        (*v)[i] = (*v)[0] + size2 * i;
    }
}


int main(int argc, char **argv)
{
  double mse;
  long window=150;
  if (argc>1)
    window = atol(argv[1]);
  
  // histogram is in range 0.0 to 0.1 --> delta is 0.1/1500
  double x[]=THETAS1500;
  double y[1500];
  long  i;
  double delta = 1000./1500.;
  //double delta = 0.1/1500.;
  double xi = 0.0;
  double xval[1500];
  xval[0]=0.0;
  for (i=1; i < 1500; i++)
  {
    xval[i] = xval[i-1] + delta;
  }
  long w = 10;
  double bestmse = 1e100;
  if (window==0)
    {
      while(w<1000)
	{
	  memcpy(y,x,sizeof(double)*1500);
	  mse = savitzky_golay(&y[0],1500, xval, w, 0, 3);
	  mse = mse / sqrt(w);
	  fprintf(stderr,"MSE=%f (%f) el=%li\n",mse, bestmse,w);
	  if (mse < bestmse)
	    {
	      bestmse = mse;
	      window = w;
	    }
	  w = (long) (1.2 * w);
	}
      mse = savitzky_golay(&x[0],1500, xval, window, 0, 3);
      fprintf(stderr,"using MSE=%f (%f) el=%li\n",mse/sqrt(w), bestmse,window);
    }
  else
    {      
      mse = savitzky_golay(&x[0],1500, xval, window, 0, 3);
      fprintf(stderr,"using MSE=%f el=%li\n",mse/sqrt(w), window);
    }
  for (i=0; i < 1500; i++)
  {
    printf("%f %f\n",xval[i],x[i]);
  }
  return 1;
}
#endif
