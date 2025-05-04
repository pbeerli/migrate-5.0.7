/*------------------------------------------------------
 Bayesian inference of
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm
 -------------------------------------------------------
 kernel density smoothing using an epanetchnikov kernel


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

#include "definitions.h"
#include "migration.h"

void kernel_smooth(double *x, long xelem, double *fx, long fxelem, long el, MYREAL mini, MYREAL maxi);


//
// kernel smoother uses params not histogram cells
// a window of size 2*el + 1 is averaged and replaces the central value
// the array with the values is prepended with el zeros and also appended with el zeros,
// and the output replaces the input array,
// this methods uses an Epanechnikov kernel [k(x) = 3/4 (1-x^2) for -1<=x<=1
//
// use this to generate the histograms from the parameter list [2020 June 24]
//
void kernel_smooth(double *x, long xelem, double *fx, long fxelem, long el, MYREAL mini, MYREAL maxi)
{
    double d;
    double dd;
    double *weight;
    double total;
    long i, j, jj;
    long el2 = 2 * el + 1;
    weight = (double *) calloc((size_t) el2, sizeof(double));
    if(xelem==0)
        return;
    weight[el] = 0.75;
    total = 1.;
    d =  (double) (2./(el2-1.));
    for(j = el+1; j < el2; j++)
      {
	dd = (-1.0 + d*j);
        
        dd *= dd;
        weight[j] =  0.75 * (1.0 - dd);
        total += 2 * weight[j];
      }
    //    printf("weight=%f\n",total);
    weight[el] /= total;
    
    for(j=el+1, jj = el-1; j < el2; j++, jj--)
      {
        weight[j] /= total;
        weight[jj] = weight[j];
      }
    total = 0.0;
    d = (maxi - mini)/fxelem;
    for(i=0; i < xelem; i++)
      {
	//find fx bin that contains xx[i]
	long fi = (long) (0.5 + (x[i] - mini)/d);
	fx[fi] += weight[el];
	for(j=el+1, jj = el-1; j < el2; j++, jj--)
	  {
	    if ((fi+j-el) < fxelem)
	      {
		double w1 = weight[j];
		total += w1;
		fx[fi+j-el]  += w1;
	      }
	    if ((fi-jj-el)>=0)
	      {
		double w2 = weight[jj];
		total += w2;
		fx[fi-jj-el] += w2; 
	      }
	  }
      }
    //    total = 0.0;
    //for(i=0; i < fxelem; i++)
    //total += fx[i];
    total = 1.0/total;
    for(i=0; i < fxelem; i++)
      {
	fx[i] *=  total;
      }
    free(weight);
}


#ifdef STANDALONE
#include "kernel-density-helper.h"
int main()
{
  double data[] = THETAS;
  long dnum= 9999;
  double fx[1500] = { 0 };
  long fxnum = 1500;
  double mini = 0.0;
  double maxi = 0.1;
  double delta = (maxi-mini)/fxnum;
  long i;
  //window = 2*el + 1;
  double el = 20;
  kernel_smooth(&data[0], dnum, &fx[0], fxnum, el, mini, maxi);
  for (i=0;i<fxnum;i++)
    {
      printf("%f %f %li\n",fx[i],mini+i*delta,i); 
    }
  return 1;
}

#endif
