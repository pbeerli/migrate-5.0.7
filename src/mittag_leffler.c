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

// to test standalone use
//gcc -g -DSTANDALONEMITTAGLEFFLER mittag_leffler.c hermite_interpoly.c romberg.c mittag_leffler_interpol_data.c -o mlf

#include "mittag_leffler.h"
#include <math.h>
#include "bayes.h"
#include "migration.h"
#include "random.h"
#include "romberg.h"
#include "hermite_interpoly.h"
#include "priors.h"
#include "sighandler.h"
#include "speciate.h"

//!Mittage-Leffler function in the real case
//!Somayeh Mashayekhi March 2017
// translated from Fortran to C Peter Beerli March 2017
#include <complex.h>    // Standard Library of Complex Numbers
//#undef I
#ifdef WINDOWS
#include "mycomplex.h"
#ifdef NMAKE
#define J {0.0, 1.0}
#else
#define J ((_Dcomplex) {0.0, 1.0})
#endif
#else
#define J ((complex double) I)
#endif
#ifndef Pi
#define Pi 3.14159265358979323846264338328
#endif
#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif
#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif

MYCOMPLEX mittag_leffler(double alpha, double beta, MYCOMPLEX z);
MYCOMPLEX  K(double alpha, double beta, double x, MYCOMPLEX z);
MYCOMPLEX P(double alpha, double beta, double eps, double phi, MYCOMPLEX z); // returns , MYCOMPLEX P1
MYCOMPLEX  ML(MYCOMPLEX z, double alpha, double beta, double Q, double X0);
double KK(double x, MYCOMPLEX *args, long arglen);
double PP(double phi, MYCOMPLEX *args, long arglen);
double interval_mittag_leffler(double r, double alpha, double lambda, double tmin, double tmax);
void set_mittag_leffler(option_fmt * options);
double interval_mittag_leffler_func(double r, double alpha, double t0, double mu, double sigma, species_fmt *s, double tmin, double tmax);
double propose_new_mlftime(double lambda, double alpha, double r1, double r2);

void change_mittag_leffler(world_fmt* world);

// calculates internals of the generalized mittag-leffler function for z with alpha and beta
// with precision Q and X0
MYCOMPLEX  ML(MYCOMPLEX z, double alpha, double beta, double Q, double X0) 
{
  long K0;
#ifdef WINDOWS
  MYCOMPLEX b = {0.0, 0.0};
  MYCOMPLEX aa = {alpha, 0.0};
  MYCOMPLEX bb = {beta, 0.0};
  MYCOMPLEX one = {1.0,0.0};
  MYCOMPLEX kargs[3] = {aa, bb, z};
  MYCOMPLEX pargs[4] = {aa, bb, one, z};
#else
  MYCOMPLEX b = 0.0;
  MYCOMPLEX kargs[3] = {alpha, beta, z};
  MYCOMPLEX pargs[4] = {alpha, beta, 1.0, z};
#endif
  if (1.0 < alpha) //alpha > 1;  beta in Real;  z in Complex; 
    {
      K0 = (long) (floor(alpha) + 1.0);
      int i;
#ifdef WINDOWS
#ifdef NMAKE
      MYCOMPLEX invK0 = {1.0/K0, 0.0};
#else
      MYCOMPLEX invK0 = (_Dcomplex) {1.0/K0, 0.0};
#endif
      for (i=0; i<K0; ++i)
	{
	  MYCOMPLEX k0temp = _Cmulcr(J, 2.0*Pi*i/K0);
	  MYCOMPLEX ctemp = _Cmulcc(cpow(z,invK0),cexp(k0temp));
	  b = myadd(b,ctemp);
	}
#else
      double invK0 = 1.0/K0;
      for (i=0; i<K0; ++i)
	{
	  b += cpow(z,invK0) * cexp(2.0*Pi*J*i/K0);
	}
#endif
      //printf("1 z= %.5f alpha > 1.0 \n", creal(z));
      return b;
    }
  else if (alpha == 1.0  && beta == 1.0)
    {
      b = cexp(z);
      return b;
    }
#ifdef WINDOWS
  else if (creal(z) == 0.0)
#else
  else if (z==0.0) // z==0
#endif
    {
      //printf("2 (z=%.5f)==0 alpha < 1.0 \n", creal(z));
#ifdef WINDOWS
#ifdef NMAKE
      MYCOMPLEX b2 = { 1.0/tgamma(beta) , 0.0};
#else
      MYCOMPLEX b2 = (_Dcomplex) { 1.0/tgamma(beta) , 0.0};
#endif
      return b2;
#else
      b = 1.0 /tgamma(beta);
      return b;
#endif
    }
  else if (cabs(z)<1.0) // |z| < 1.0 z:{-1,1}
    {
      //K0 = (long) max(creal((ceil ((1.0-beta)/alpha))),ceil(creal(clog(Q*(1.0-cabs(z)))/clog(cabs(z)))));
#ifdef WINDOWS
      double tempa = ceil ((1.0-beta)/alpha);
      double invdenom = 1.0/log(cabs(z));
      MYCOMPLEX tmpc = {Q*(1.0 - cabs(z)),0.0};
      MYCOMPLEX num = clog(tmpc);
      MYCOMPLEX tempd = _Cmulcr(num, invdenom);
      double tempb = ceil(creal(tempd));
      K0 = (long) max(tempa,tempb);
      MYCOMPLEX b = {0.0, 0.0};
#else   
      K0 = (long) max((ceil ((1.0-beta)/alpha)),ceil(creal(clog(Q*(1.0-cabs(z)))/clog(cabs(z)))));
      b = 0.0;
#endif
      K0 = (long) min(200,K0);
      int i;
      for (i=0;i<K0+1;i++)
	{
#ifdef WINDOWS
	  double x = 1.0 / tgamma(beta+alpha*i);
	  MYCOMPLEX ic = {(double) i, 0.0};
	  MYCOMPLEX temp0 = cpow(z,ic);
	  MYCOMPLEX temp = _Cmulcr(temp0, x);
	  b = myadd(b, temp);
#else
	  b += (cpow(z,i))/ tgamma(beta+alpha*i);
#endif
	}
      return b;
    }
  else if (cabs(z)>floor(10.0+5.0*alpha))
    {
#ifdef WINDOWS
       MYCOMPLEX tempQ = {Q,0.0};
       MYCOMPLEX temp1 = {cabs(z),0.0};
       MYCOMPLEX ltempQ = clog(tempQ);
       MYCOMPLEX ltemp1 = clog(temp1);
       MYCOMPLEX div =  mydiv(ltempQ,ltemp1);
       K0= (long) (floor(-creal(div)));
       MYCOMPLEX c = {0.0, 0.0};
#else
       K0= (long) (floor(creal(-clog(Q)/clog(cabs(z)))));
       double complex c = 0.0;
#endif
      K0 = min(200,K0);
      int i;
      for(i=1; i<K0+1; i++)
	{
#ifdef WINDOWS
	  double x = 1.0 / tgamma(beta-alpha*i);
	  MYCOMPLEX negi = {(double) -i, 0.0};
	  MYCOMPLEX temp0 = cpow(z,negi);
	  c = myadd(c, _Cmulcr(temp0,x));
#else
	  c += cpow(z,(-i)) / tgamma(beta-alpha*i);
#endif
	}
      if(fabs(carg(z)) < (alpha*Pi/4.0 + 0.5 * min(Pi, alpha*Pi)))
	{
#ifdef WINDOWS
	  MYCOMPLEX tempd= {(1.0/alpha), 0.0};
	  MYCOMPLEX tempe = {(1-beta)/alpha,0.0};
	  MYCOMPLEX temp0 = cpow(z,tempe);
	  MYCOMPLEX temp1 =  cexp(cpow(z,tempd));
	  MYCOMPLEX temp2 = _Cmulcc(temp0,temp1);
	  b = _Cmulcc(temp2,tempd);
	  b = mysub(b,c);
#else
	  b=(1.0/alpha)*cpow(z,((1-beta)/alpha))*cexp(cpow(z,(1.0/alpha))) - c;
#endif
	}
      else 
	{
#ifdef WINDOWS
	  b = neg(c);
#else
	  b = -c;
#endif
	}
      return b;
    }
  else
    {
      // compute X0 (done outside)
      if(fabs(carg(z)) > alpha * Pi)
	{
#ifdef WINDOWS
	  if (beta <= 1.0)
	    {
	      double romres = romberg(&KK, kargs, 3, 0.0, X0, 100, Q);
#ifdef NMAKE
	      b = {romres, 0.0};
#else
	  b = (_Dcomplex) {romres, 0.0};
#endif
	    }
	  else
	    {
	      double romres = romberg(&KK, kargs, 3, 1.0, X0, 100, Q);
	      romres += romberg(&PP, pargs, 3, -alpha*Pi, alpha*Pi, 100, Q);
#ifdef NMAKE
	      b =  {romres, 0.0};
#else
	      b =  (_Dcomplex) {romres, 0.0};
#endif
	    }
#else
	  if (beta <= 1.0)
	    {
	      b = romberg(&KK, kargs, 3, 0.0, X0, 100, Q);
	    }
	  else
	    {
	      b = romberg(&KK, kargs, 3, 1.0, X0, 100, Q);
	      b += romberg(&PP, pargs, 3, -alpha*Pi, alpha*Pi, 100, Q);
	    }
#endif
	  return b;
	}
      else if(fabs(carg(z)) < alpha * Pi)
	{
	  if (beta <= 1)
	    {
#ifdef WINDOWS
	      double romres = romberg(&KK, kargs, 3, 0.0, X0, 100, Q);
#ifdef NMAKE
	      b = {romres, 0.0};	    
	      MYCOMPLEX tempd= {(1.0/alpha), 0.0};
	      MYCOMPLEX tempe = {(1-beta)/alpha,0.0};
#else
	      b = (_Dcomplex) {romres, 0.0};	    
	      MYCOMPLEX tempd= (_Dcomplex) {(1.0/alpha), 0.0};
	      MYCOMPLEX tempe = (_Dcomplex) {(1-beta)/alpha,0.0};
#endif
	      MYCOMPLEX temp0 = cpow(z,tempe);
	      MYCOMPLEX temp1 =  cexp(cpow(z,tempd));
	      MYCOMPLEX temp2 = _Cmulcc(temp0,temp1);
	      b = myadd(b,_Cmulcc(temp2,tempd));
#else
	      b = romberg(&KK, kargs, 3, 0.0, X0, 100, Q);
	      b += (1.0/alpha)*cpow(z,((1.0-beta)/alpha)) * cexp(cpow(z,(1.0/alpha)));
#endif
	    }
	  else
	    {
#ifdef WINDOWS
	      double tempabs = cabs(z)/2.0;
#ifdef NMAKE
	      pargs[2]= {tempabs,0.0};
#else
	      pargs[2]= (_Dcomplex) {tempabs,0.0};
#endif
	      double romres = romberg(&KK, kargs, 3, creal(pargs[2]) , X0, 100, Q);
	      romres += romberg(&PP, pargs, 3, -alpha*Pi, alpha*Pi, 100, Q);
#ifdef NMAKE
	      b = {romres, 0.0};
#else
	      b = (_Dcomplex) {romres, 0.0};
#endif
	      MYCOMPLEX tempd= {(1.0/alpha), 0.0};
	      MYCOMPLEX tempe = {(1.0-beta)/alpha,0.0};
	      MYCOMPLEX temp0 = cpow(z,tempe);
	      MYCOMPLEX temp1 =  cexp(cpow(z,tempd));
	      MYCOMPLEX temp2 = _Cmulcc(temp0,temp1);
	      b = myadd(b,_Cmulcc(temp2,tempd));
#else
	      pargs[2]= cabs(z)/2.0;
	      b = romberg(&KK, kargs, 3, creal(pargs[2]) , X0, 100, Q);
	      b += romberg(&PP, pargs, 3, -alpha*Pi, alpha*Pi, 100, Q);
	      b += (1.0/alpha)*cpow(z,((1.0-beta)/alpha)) * cexp(cpow(z,(1.0/alpha)));
#endif
	    }
	  return b;
	}
      else
	{
	  //double lowbound = creal((cabs(z) + 1.0)/ 2.0);
	  double lowbound = (cabs(z) + 1.0)/ 2.0;
#ifdef WINDOWS
#ifdef NMAKE
	  pargs[2]= {lowbound, 0.0};
#else
	  pargs[2]= (_Dcomplex) {lowbound, 0.0};
#endif
	  double romres = romberg(&KK, kargs, 3, creal(pargs[2]) , X0, 100, Q);
	  romres += romberg(&PP, pargs, 3, -alpha*Pi, alpha*Pi, 100, Q);
#ifdef NMAKE
	  b = {romres, 0.0};
#else
	  b = (_Dcomplex) {romres, 0.0};
#endif
		  
#else
	  pargs[2]= lowbound;
	  b = romberg(&KK, kargs, 3, creal(pargs[2]) , X0, 100, Q);
	  b += romberg(&PP, pargs, 3, -alpha*Pi, alpha*Pi, 100, Q);
#endif
	}
      return b;
    }
}

MYCOMPLEX  K(double alpha, double beta, double x, MYCOMPLEX z) 
{
#ifdef WINDOWS
  if(x==0.0)
#ifdef NMAKE
    return {0.0, 0.0};
#else
  return (_Dcomplex) {0.0, 0.0};
#endif
#else
  if(x==0.0)
    return (double complex) 0.0;
#endif
  double a = (1.0-beta)/alpha;
  double b = alpha * Pi;
  double a0 = 1/b;
#ifdef WINDOWS
#ifdef NMAKE
  MYCOMPLEX xx = {x,0.0};
  MYCOMPLEX aa = {a,0.0};
#else
  MYCOMPLEX xx = (_Dcomplex) {x,0.0};
  MYCOMPLEX aa = (_Dcomplex) {a,0.0};
#endif  
  MYCOMPLEX temp0 = cpow(xx,aa);
  MYCOMPLEX inva = {1.0/alpha,0.0};
  MYCOMPLEX temp1 = cexp(neg(cpow(xx,inva)));
  MYCOMPLEX a1 = {Pi*(1.0-beta),0.0};
  MYCOMPLEX a2 = {Pi*(1.0-beta+alpha),0.0};
  MYCOMPLEX a3 = {b,0};
  MYCOMPLEX a4 = _Cmulcr(csin(a1),x);
  MYCOMPLEX a5 = _Cmulcc(csin(a2),z);
  MYCOMPLEX K0 = mysub(a4,a5);
  MYCOMPLEX a7 = _Cmulcr(z, 2.0*x);
  MYCOMPLEX a6 = _Cmulcc(ccos(a3),a7);
  MYCOMPLEX K1 = _Cmulcr(temp0,a0);
  K1 = _Cmulcc(K1,temp1);
#ifdef NMAKE
  MYCOMPLEX xc = {x*x,0.0};
#else
    MYCOMPLEX xc = (_Dcomplex) {x*x,0.0};
#endif
  MYCOMPLEX denom = myadd(xc,a6);
  denom = myadd(denom, _Cmulcc(z,z));
  K1 = _Cmulcc(K1,K0);
  K1 = mydiv(K1,denom);  // if this does not fail than lots of things can be simplified
#else
  //  ((x*csin(Pi*(1.0-beta))-z*csin(Pi*(1.0-beta+alpha)))/(x*x-2.0*x*z*ccos(alpha*Pi)+z*z))
  MYCOMPLEX K1=a0*cpow(x,a)*
    cexp(-cpow(x,(1.0/alpha)))*
    ((x*csin(Pi*(1.0-beta))-z*csin(Pi*a))/(x*x-2.0*x*z*ccos(b)+z*z));
#endif
    return K1;
}

MYCOMPLEX P(double alpha, double beta, double eps, double phi, MYCOMPLEX z) // returns , MYCOMPLEX P1
{
#ifdef WINDOWS
  double phistar = phi * (1.0+ (1.0-beta)/alpha);
  MYCOMPLEX omega = {phistar + pow(eps,(1.0/alpha)) * sin(phi/alpha), 0.0};    
  double P1star = (1.0/(2.0*alpha*Pi)) * pow(eps,(1.0+(1.0-beta)/alpha));
  MYCOMPLEX temp1 = {pow(eps,(1.0/alpha)), 0.0};
  MYCOMPLEX P1 = _Cmulcr(cexp(temp1),P1star);
  double costemp = cos(phi/alpha);
  P1 = _Cmulcr(P1,costemp);
  MYCOMPLEX temp2 = myadd(ccos(omega),_Cmulcc(J,csin(omega)));
  P1 = _Cmulcc(P1, temp2);
  MYCOMPLEX temp3 = _Cmulcr(cexp(_Cmulcr(J,phi)),eps);
  MYCOMPLEX temp4 = mysub(temp3,z);
  P1 = mydiv(P1,temp4);
#else
  MYCOMPLEX omega = phi * (1.0+ (1.0-beta)/alpha) + cpow(eps,(1.0/alpha)) * csin(phi/alpha);    
  MYCOMPLEX P1 = (1.0/(2.0*alpha*Pi)) * cpow(eps,(1.0+(1.0-beta)/alpha)) * cexp(cpow(eps,(1.0/alpha)) * ccos(phi/alpha)) * (ccos(omega)+J*csin(omega))/(eps*cexp(J*phi)-z);
#endif
  return P1;
}

//double mlf(double alpha, double beta, MYCOMPLEX z)
MYCOMPLEX mittag_leffler(double alpha, double beta, MYCOMPLEX z)
{
  
  double Q=PRECISION;  // set in definitions.h
  double X0;
  double dz = creal(z);
  
  if(beta >= 0.0)
    {
      //            X0 = max(
      //       max(creal(1.0), creal(2.0*creal(z))),
      //       creal(cpow((-clog(Pi*Q/creal(6.0))),(alpha))));
      X0 = max(
	       max(1.0, 2.0*dz),
	       pow((-log(Pi*Q/6.0)),alpha));
    }
  else
    {
      //X0 = max(
      //       max(creal(cpow((fabs(beta)+1.0),(alpha))), creal(2.0*creal(z))),
      //       creal(cpow((-2.0*clog(Pi*Q/creal(6*(fabs(beta)+2.0)* cpow((2.0*fabs(beta)),fabs(beta))))),(alpha))));
      double fb = fabs(beta);
      X0 = max(
	       max(pow((fb+1.0), alpha),
		   2.0*dz),
	       pow(-2.0 * log(Pi*Q/(6*(fb+2.0)* pow(2.0*fb,fb))), alpha));
    }
#ifndef MLF_SLOW
  MYCOMPLEX b;
  if (beta==1.0 || fabs(alpha-beta)>EPSILON)
    b = MLinterpol(z,alpha, beta, Q, X0);
  else
    b = clog(ML(z,alpha, beta, Q, X0));
#ifdef MLFCHECK
  MYCOMPLEX c = clog(ML(z,alpha, beta, Q, X0));
  printf("@CHECK@ interpol= %f mlf= %f z= %f a= %f b= %f\n",creal(b),creal(c),creal(z),alpha,beta);
#endif
#else
  MYCOMPLEX b = clog(ML(z,alpha, beta, Q, X0));
#endif
  return b;
}

double KK(double x, MYCOMPLEX *args, long arglen)
{
  (void) arglen;
  double alpha = creal(args[0]);
  double beta  = creal(args[1]);
  MYCOMPLEX z = args[2];
  MYCOMPLEX b = K(alpha,beta,x,z);
  return creal(b);
}

double PP(double phi, MYCOMPLEX *args, long arglen)
{
  (void) arglen;
  double alpha = creal(args[0]);
  double beta  = creal(args[1]);
  double eps   = creal(args[2]);
  MYCOMPLEX z = args[3];
  MYCOMPLEX b = P(alpha,beta,eps,phi,z);
  return creal(b);
}
//1) For drawing time you need to solve
//   r=E_alpha(-lambda *t^alpha).
//   Before you had r=E_alpha(-lambda *t).

//2) For probability you need to write
//   t^(alpha-1)*E_alpha,alpha(-lambda *t^alpha).
//   Before you       had E_alpha,alpha(-lambda *t).


double interval_mittag_leffler(double r, double alpha, double lambda, double tmin, double tmax)
{
  double minu=tmin;
  double maxu=tmax;
  double u=(maxu+minu)/2.0;
  double prob;
  long maxcount=200;
  double beta = 1.0;
  double logr = log(r);
  while(maxu-minu>SMALLEPSILON && maxcount-- > 0)
    {
#ifdef WINDOWS
      MYCOMPLEX  arguments = { (-(pow(u,alpha)) * lambda), 0.0};
#else
      MYCOMPLEX arguments = (MYCOMPLEX) (-(pow(u,alpha)) * lambda);
#endif
      prob = creal(mittag_leffler(alpha,beta, arguments));
      //printf("@ %f %f %f prob=%f r=%f \n",minu,u,maxu, prob, r);
      if (logr>prob)
	{
	  maxu = u;
	}
      else
	{
	  minu = u;
	}
      u = (maxu+minu)/2.0;
    }
  if (maxcount<=0)
    {
      fprintf(stderr,"WARNING: drawing a random time for MLF failed: minu=%f u=%f maxu=%f\n",minu,u,maxu);
    }
  return u;
}

#ifndef STANDALONEMITTAGLEFFLER
double interval_mittag_leffler_func(double r, double alpha, double t0, double mu, double sigma, species_fmt *s, double tmin, double tmax)
{
  double minu=tmin;
  double maxu=tmax;
  double u=(maxu+minu)/2.0;
  double prob;
  long maxcount=200;
  double beta = 1.0;
  double logr = log(r);
  double lambda;
  while(maxu-minu>EPSILON && maxcount-- > 0)
    {
      double t1 = t0 + u;
      lambda = (*log_prob_wait_speciate)(t0,t1,mu,sigma,s);
#ifdef WINDOWS
      MYCOMPLEX  arguments = { (-(pow(u,alpha)) * lambda), 0.0};
#else
      MYCOMPLEX arguments = (MYCOMPLEX) (-(pow(u,alpha)) * lambda);
#endif
      prob = creal(mittag_leffler(alpha,beta, arguments));
      //printf("@ %f %f %f prob=%f r=%f \n",minu,u,maxu, prob, r);
      if (logr>prob)
	{
	  maxu = u;
	}
      else
	{
	  minu = u;
	}
      u = (maxu+minu)/2.0;
    }
  if (maxcount<=0)
    warning("drawing a random time for MLF failed: minu=%f u=%f maxu=%f\n",minu,u,maxu);
  return u;
}

// Somayeh Mashayekhi 
double propose_new_mlftime(double lambda, double alpha, double r1, double r2)
{
  double pia = PI * alpha;
  double denoma = 1.0 / alpha;
  double denomlambda = 1.0 / lambda;
  //return -pow(denomlambda,denoma) * pow(sin(pia)/(tan(pia*(1-r1))) - cos(pia),denoma) * log(r2);
  return -pow(denomlambda * (sin(pia)/(tan(pia*(1-r1))) - cos(pia)),denoma) * log(r2);
}

void change_mittag_leffler(world_fmt * world)
{
  //CAREFUL this still has hardcoded elements for the prior
  double oldalpha = world->mlalpha;
  MYREAL r = UNIF_RANDUM();
  MYREAL delta = 0.05;
  MYREAL np = propose_uniform(oldalpha, 0.4, 1.0, &r, delta);
  MYREAL oldval = probg_treetimes(world);
  world->mlalpha=np;
  MYREAL newval = probg_treetimes(world);
  boolean success = bayes_accept(newval, oldval,world->heat, 0.0);
  if (!success)
    world->mlalpha=oldalpha;
}

void set_mittag_leffler(option_fmt * options)
{
  char *input;
  input = (char *) mycalloc(LINESIZE,sizeof(char));
  printf("Mittag-Leffler distribution\n");
  printf("----------------------------------------------------------------\n");
  printf("The standard Coalescent is based on the Exponential distribution\n");
  printf("Mittag-Leffler is a generalization with an additional parameter\n");
  printf("alpha, ranging between 0 and 1; With an alpha=1.0 we recover the\n");
  printf("standard Exponential; with alpha < 1 the distribution has\n");
  printf("shorter coalescences near today, values are allowed\n");
  printf("between 0.01 and 1.00, using 2-digit precision.\n");
  printf("[Default: alpha=1.00]\n");
  printf("> ");fflush(stdout);
  FGETS(input,LINESIZE,stdin);
  if(input[0]!='\0')
    {
      options->mlalpha = atof(input);
    }
  else
    {
      options->mlalpha = 1.0;
    }
}
#endif



#ifdef STANDALONEMITTAGLEFFLER
void run_grid(double alpha,double beta)
{
  MYCOMPLEX zmax = 100;
  MYCOMPLEX zmin = -100;
  MYCOMPLEX z;
  for (z=zmin; creal(z) < creal(zmax); z += 1.0)
    {
      MYCOMPLEX b = mittag_leffler(alpha, beta, z);
      printf("mlf(z=%f,alpha=%f,beta=%f) =  %.8g %+.8gi\n", creal(z), alpha, beta, creal(b), cimag(b));
    }
}


void help(int argc, char **argv);

void help(int argc, char **argv)
{
  if ((argc==1) || (argc>1 && argv[1][0]=='-' && argv[1][1]=='h'))
    {
      printf("Syntax: mlf g           #prints mlf(lambda=[-100..99],alpha=0.6,beta=0.6)\n");
      printf("Syntax: mlf z           #prints mlf(lambda=z,alpha=0.6,beta=0.6)\n");
      printf("Syntax: mlf z a         #prints mlf(lambda=z,alpha=a,beta=a)\n");
      printf("Syntax: mlf z a b       #prints mlf(lambda=z,alpha=a,beta=b)\n");
      exit(-1);
    }
}

int main(int argc, char ** argv)
{
  double alpha=0.6;
  double beta = alpha;
  //double eps = 1e-8;
  MYCOMPLEX z = 20.0 + 0.0 * I;
  help(argc,argv);
  switch (argc)
    {
    case 2:
      if (argv[1][0]=='t')
	{
	  long i;
	  for(i=0;i<1000000;i++)
	    {
	      mittag_leffler(alpha, beta, z);
	    }
	  return 0;
	}
      if (argv[1][0]=='g')
	{
	  run_grid(alpha,beta);
	  return 0;
	}
      else
      z = atof(argv[1]) + 0.0 * I;
      break;
    case 3:
      z = atof(argv[1]) + 0.0 * I;
      alpha = atof(argv[2]);
      beta = alpha;
      break;
    case 4:
      z = atof(argv[1]) + 0.0 * I;
      alpha = atof(argv[2]);
      beta = atof(argv[3]);
      break;
    default:
      break;
    }
  MYCOMPLEX b = mittag_leffler(alpha, beta, z);
  printf("log(mlf(z=%f,alpha=%f,beta=%f)) =  %.8g %+.8gi\n", creal(z), alpha, beta, creal(b), cimag(b));
  return 0;
}
  
#endif /*end standalonemittagleffler*/




