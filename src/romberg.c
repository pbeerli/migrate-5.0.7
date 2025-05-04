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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "definitions.h"
#include "sighandler.h"
#include "romberg.h"

// ROMBERG code is from Wikipedia (downloaded 2017)
// https://en.wikipedia.org/wiki/Romberg%27s_method
// Code from wikipedia
void dump_row(size_t i, double *R);
double romberg(double (*f)(double, MYCOMPLEX *, long), MYCOMPLEX *args, long arglen, double a, double b, size_t max_steps, double acc);
double check_bounds(double (*f)(double, MYCOMPLEX *, long), MYCOMPLEX *args, long arglen, double a, double b, size_t max_steps, double acc);

double check_bounds(double (*f)(double, MYCOMPLEX *, long), MYCOMPLEX *args, long arglen, double a, double b, size_t max_steps, double acc)
{
  double x = 0.0;
  //double bb = b;
  double maxu = b;
  double minu = a;
  size_t maxcount=0;
  //boolean done=FALSE;
  //boolean second=FALSE;
  double bb = (maxu-minu)/2.0;
  while(maxu-minu>acc && maxcount++ < max_steps)
    {
      x = f(bb,args, arglen);
      if(!isfinite(x))
	{
	  maxu = bb;
	}
      else
	{
	  minu = bb;
	}
      bb = (maxu-minu)/2.0;
    }
      //do
      //{
      //x = f(bb,args, arglen);
      //bb = bb/2.;
      //}
      //while(!isfinite(x));
  printf("maxcount: %li upper bound: %f with value %f\n", maxcount, bb, x);
  return bb;
}

double romberg(double (*f)(double, MYCOMPLEX *, long), MYCOMPLEX *args, long arglen,
	       double a, double b, size_t max_steps, double acc)
{
  //double R1[max_steps], R2[max_steps]; //buffers
  //double *Rp = &R1[0], *Rc = &R2[0]; //Rp is previous row, Rc is current row
  double *R1;
  double *R2;
  R1 = (double *) mycalloc(max_steps,sizeof(double));
  R2 = (double *) mycalloc(max_steps,sizeof(double));
  double *Rp = R1;
  double *Rc = R2;
  double h = (b-a); //step size
  Rp[0] = (f(a,args, arglen) + f(b, args, arglen))*h*.5; //first trapezoidal step
   size_t i;
   size_t j;
   if (max_steps > 15)
     max_steps = 15;
   for(i = 1; i < max_steps; ++i){
      h /= 2.;
      double c = 0;
      size_t ep = 1 << (i-1); //2^(n-1)
      for(j = 1; j <= ep; ++j){
	c += f(a+(2*j-1)*h, args, arglen);
      }
      Rc[0] = h*c + .5*Rp[0]; //R(i,0)

      for(j = 1; j <= i; ++j){
	unsigned long n_k = ( 1 << (2*j));
         Rc[j] = (n_k*Rc[j-1] - Rp[j-1])/(n_k-1); //compute R(i,j)
      }

#ifdef ROMBERGERTESTER
      //Dump ith column of R, R[i,i] is the best estimate so far
      //dump_row(i, Rc);
#endif
      if(isinf(Rc[i]))
	{
	  //return Rc[i-1];
	  double rcval = Rc[i-1];
	  myfree(Rc);
	  myfree(Rp);
	  return rcval;
	}
      if(i > 1 && fabs(Rp[i-1]-Rc[i]) < acc){
	{
	  //return Rc[i-1];
	  double rcval = Rc[i-1];
	  myfree(Rc);
	  myfree(Rp);
	  return rcval;
	}
      }
      //swap Rn and Rc as we only need the last row
      double *rt = Rp;
      Rp = Rc;
      Rc = rt;
   }
   //return Rp[max_steps-1]; //return our best guess
   double rpval = Rp[max_steps-1];
   myfree(Rc);
   myfree(Rp);
   return rpval;
}

#ifdef ROMBERGERTESTER
double f(double x, _Complex double *args, long arglen)
{
  double v = args[0] * exp(-(x*x) * args[0]);
  return v;
}

double func_probg_div(double x, _Complex double *args, long arglen)
{
  double t0=0.0;
  double mu=0.0025;
  double sigma=0.0025;
  double theta = 0.01;
  double M = 100;
  long k = 2;
  if (arglen>=6)
    {
      t0 = creal(args[0]);
      theta = creal(args[1]);
      M = creal(args[2]);
      k = (long) creal(args[3]);
      mu = creal(args[4]);
      sigma = creal(args[5]);
    }
  double t=t0+x;
  double lambda2 = k * (k-1)/theta;
  double lambda3 = k * M;
  double mut = mu-t;
  double mut0 = mu-t0;
  double twosigma2 = 2.0*sigma*sigma;
  double sqrt2sigma = sqrt(2.0) * sigma;
  
  double explambda1lambda1;
  double first = -mut*mut/twosigma2 + log(k) +
    log(sqrt(2.0/3.1415926));
  double v = erf(mut/sqrt2sigma);
  double second;
  if (v > -1)
    second = log(1.0+v) * k;
  else
    second = -100000.;
  double third = log(1.0 + erf(mut0/sqrt2sigma)) * (-k);
  double fourth = log((sigma * (1.0-erf((-mu + t)/sqrt2sigma))));  
  // Exp[-lambda1[t]]*Exp[-(t-to) (lambda2 + lambda3)
  double myf = exp(-(t-t0) * (lambda2 + lambda3)) * exp(first+second+third - fourth);
  //printf ("t=%f f(t)=%f (exp[%f+%f+%f-%f])\n",t,myf,first,second,third,fourth);
  return myf;
}
double func_probg_theta(double x, _Complex double *args, long arglen)
{
  double t0=0.0;
  double mu=0.0025;
  double sigma=0.0025;
  double theta = 0.01;
  double M = 100;
  long k = 2;
  if (arglen>=6)
    {
      t0 = creal(args[0]);
      theta = creal(args[1]);
      M = creal(args[2]);
      k = (long) creal(args[3]);
      mu = creal(args[4]);
      sigma = creal(args[5]);
    }
  double t=t0+x;
  double lambda2 = k * (k-1)/theta;
  double lambda3 = k * M;
  double mut = mu-t;
  double mut0 = mu-t0;
  double twosigma2 = 2.0*sigma*sigma;
  double sqrt2sigma = sqrt(2.0) * sigma;
  
  double explambda1lambda1;
  //double first = -mut*mut/twosigma2 + log(k) +
  //  log(sqrt(2.0/3.1415926));
  double v = erf(mut/sqrt2sigma);
  double second;
  if (v > -1)
    second = log(1.0+v) * k;
  else
    second = -100000.;
  double third = log(1.0 + erf(mut0/sqrt2sigma)) * (-k);
  //double fourth = log((sigma * (1.0-erf((-mu + t)/sqrt2sigma))));  
  // Exp[-lambda1[t]]*Exp[-(t-to) (lambda2 + lambda3)
  double myf = lambda2 * exp(-(t-t0) * (lambda2 + lambda3)) * exp(second+third);
  //printf ("t=%f f(t)=%f (exp[%f+%f+%f-%f])\n",t,myf,first,second,third,fourth);
  return myf;
}
double func_probg_M(double x, _Complex double *args, long arglen)
{
  double t0=0.0;
  double mu=0.0025;
  double sigma=0.0025;
  double theta = 0.01;
  double M = 100;
  long k = 2;
  if (arglen>=6)
    {
      t0 = creal(args[0]);
      theta = creal(args[1]);
      M = creal(args[2]);
      k = (long) creal(args[3]);
      mu = creal(args[4]);
      sigma = creal(args[5]);
    }
  double t=t0+x;
  double lambda2 = k * (k-1)/theta;
  double lambda3 = k * M;
  double mut = mu-t;
  double mut0 = mu-t0;
  double twosigma2 = 2.0*sigma*sigma;
  double sqrt2sigma = sqrt(2.0) * sigma;
  
  double explambda1lambda1;
  //double first = -mut*mut/twosigma2 + log(k) +
  //  log(sqrt(2.0/3.1415926));
  double v = erf(mut/sqrt2sigma);
  double second;
  if (v > -1)
    second = log(1.0+v) * k;
  else
    second = -100000.;
  double third = log(1.0 + erf(mut0/sqrt2sigma)) * (-k);
  //double fourth = log((sigma * (1.0-erf((-mu + t)/sqrt2sigma))));  
  // Exp[-lambda1[t]]*Exp[-(t-to) (lambda2 + lambda3)
  double myf = lambda3 * exp(-(t-t0) * (lambda2 + lambda3)) * exp(second+third);
  //printf ("t=%f f(t)=%f (exp[%f+%f+%f-%f])\n",t,myf,first,second,third,fourth);
  return myf;
}


void dump_row(size_t i, double *R){
   printf("R[%2zu] = ", i);
   for(size_t j = 0; j <= i; ++j){
      printf("%f ", R[j]);
   }
   printf("\n");
}



int main(int argc,char **argv)
{
  _Complex double *args;
  long order = atol(argv[1]);
  long arglen = 1;
  if(argc>=7)
    args = calloc(6,sizeof(_Complex double));
  else
    args = calloc(arglen,sizeof(_Complex double));
  if (argc>=8)
    {
      args[0] = atof(argv[2]) + 0.0*I; //t0
      args[1] = atof(argv[3]) + 0.0*I; //theta
      args[2] = atof(argv[4]) + 0.0*I; // M
      args[3] = atof(argv[5]) + 0.0*I; // k
      args[4] = atof(argv[6]) + 0.0*I; // mu
      args[5] = atof(argv[7]) + 0.0*I; // sigma
      arglen=6;
    }
  else
    {
      args[0] = 1.0;
      arglen = 1;
    }
  struct timespec stop;
  struct timespec start;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
  //double a = romberg(&f,args,3, 0.0,1.0,order,0.0000001);
  double b = check_bounds(&func_probg_div,args,arglen, 0.0,100.0,order,0.0001);
  double a = romberg(&func_probg_div,args,arglen, 0.0,b,order,0.0001);
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
  double elap = (stop.tv_sec - start.tv_sec) * 1e6 + (stop.tv_nsec - start.tv_nsec) / 1e3;
  printf("Divergence Elapsed: %f micro-seconds: Romberg(l*exp(-u*l),a,b,3)=%f\n",elap, a);
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
  //double a = romberg(&f,args,3, 0.0,1.0,order,0.0000001);
  b = check_bounds(&func_probg_theta,args,arglen, 0.0,1.0,order,0.0001);
  a = romberg(&func_probg_theta,args,arglen, 0.0,b,order,0.0001);
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
  elap = (stop.tv_sec - start.tv_sec) * 1e6 + (stop.tv_nsec - start.tv_nsec) / 1e3;
  printf("Theta Elapsed: %f micro-seconds: Romberg(l*exp(-u*l),a,b,3)=%f\n",elap, a);
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
  //double a = romberg(&f,args,3, 0.0,1.0,order,0.0000001);
  b = check_bounds(&func_probg_M,args,arglen, 0.0,1.0,order,0.0001);
  a = romberg(&func_probg_M,args,arglen, 0.0,b,order,0.0001);
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
  elap = (stop.tv_sec - start.tv_sec) * 1e6 + (stop.tv_nsec - start.tv_nsec) / 1e3;
  printf("Migration Elapsed: %f micro-seconds: Romberg(l*exp(-u*l),a,b,3)=%f\n",elap, a);

  return 0;
}
 
#endif
