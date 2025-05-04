/* Correlation table of all parameters for Bayesian inference
   and using the profile likelihood tables
   (1) Bayesian inference: read parameter list and calculate correlation table
   (2) ML inference: read profile table and construct correlation table 
      [uses print_cov .... already present, but needs testing]
  

  (c) Peter Beerli 2012

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
#include "bayes.h"
#include "tools.h"
#include "world.h"
#include "sighandler.h"
#include "reporter.h"
#include "correlation.h"

// functions
void covarianceBayes(world_fmt *world, long T, MYREAL *params, long offset, long locus, MYREAL ***cov);
void correlationBayes(world_fmt *world, long locus, MYREAL **cov, MYREAL ***corr);
//void covariance_bayes(world_fmt *world, long locus);
//void covariance_summary(world_fmt *world);
//void adjust_covariance(MYREAL **cov, long size, long n);

// unit test compile using clang correlation.c -o corrtest
#ifdef FUNCTIONTEST
#define NOJPEG 
#define NOPNG
#define NOZLIB
#endif
#include "pretty.h"
#include <errno.h>
extern int myID;
extern int numcpu;
#ifdef FUNCTIONTEST
void
doublevec2d (MYREAL ***v, long size1, long size2)
{
    long i;
    *v = (MYREAL **) mycalloc (size1, sizeof (MYREAL *));
    (*v)[0] = (MYREAL *) mycalloc (size1 * size2, sizeof (MYREAL));
    for (i = 1; i < size1; i++)
    {
        (*v)[i] = (*v)[0] + size2 * i;
    }
}
/// free 2D vector of MYREAL
void
free_doublevec2d (MYREAL **v)
{
    myfree(v[0]);
    myfree(v);
}

#endif



// function implementation: internal functions 
void covarianceBayes(world_fmt *world, long T, MYREAL *params, long offset, long locus, MYREAL ***cov)
{
  (void) locus;
  const long o=offset;
  const long nn = world->numpop2+(long)(world->bayes->mu) + world->species_model_size * 2 + world->grownum;
  const long nno = nn + o;
  MYREAL nk = 0;
  MYREAL *x;
  MYREAL *xn;
  MYREAL **cn;
  MYREAL *n;
  long i,k0,k,l0,l;
  doublevec1d(&xn,nn);
  doublevec2d(cov,nn,nn);
  cn = *cov;
  n = (MYREAL *) mycalloc(nn,sizeof(MYREAL));
  for (i=0;i<T;i++)
    {
      x = &params[i*nno];
      for(k0=0;k0<nn; k0++)
	{
	  if(shortcut(k0,world,&k))
	    continue;
	  else
	    {
	      n[k] += 1;
	      nk = 1./n[k];
	      xn[k] += (x[k] - xn[k])*nk;
	    }
	}
    }
  for (i=0;i<T;i++)
    {
      x = &params[i*nno+o];
      for(k0=0;k0<nn; k0++)
	{
	  if(shortcut(k0,world,&k))
            continue;
          else
            {
	      for(l0=0;l0<nn;l0++)
		{
		  if(shortcut(l0,world,&l))
		    continue;
		  else
		    {
		      cn[k][l] += (x[k] - xn[k]) * (x[l] - xn[l])*nk;//nk is correct after
		      // loop to calculate the means
		    }
		}
	    }
	}
    }
  myfree(xn);
  myfree(n);
}

void correlationBayes(world_fmt *world, long locus, MYREAL **cov, MYREAL ***corr)
{
  (void) locus;
  const long nn = world->numpop2;
  long k,l;
  MYREAL sqrtcovkk;
  doublevec2d(corr,nn,nn);
  for(k=0;k<nn; k++)
    {
      sqrtcovkk = sqrt(cov[k][k]);
      for(l=0;l<nn;l++)
	{
	  (*corr)[k][l] = cov[k][l]/(sqrt(cov[l][l])*sqrtcovkk);
	}
    }
}

// function implementation: this is the API 
void covariance_bayes(world_fmt *world, long locus)
{
  long offset=2;
  world->bayes->histogram[locus].n = world->bayes->numparams;
  covarianceBayes(world,world->bayes->numparams,world->bayes->params, offset, locus,&world->bayes->histogram[locus].covariance);
}

//void covariance_bayes2(world_fmt *world, long locus, MYREAL *params)/
//{
//  covarianceBayes(world,world->bayes->numparams,params,locus,&world->bayes->histogram[locus].covariance);
//}

void covariance_summary(world_fmt *world)
{
  const long nn = world->numpop2 + world->bayes->mu + world->species_model_size * 2 + world->grownum;
  bayeshistogram_fmt *target = &world->bayes->histogram[world->loci];
  const MYREAL invn = 1.0 / (world->loci);
  long i,j,locus;
  MYREAL **cov;
  if(world->bayes->histogram[0].covariance==NULL)
    return;
  if (target->covariance==NULL)
    {
      doublevec2d(&target->covariance,nn,nn);
    }
  for (locus=0;locus<world->loci;locus++)
    {
      if (!world->data->skiploci[locus])
	{
	  cov = world->bayes->histogram[locus].covariance;
	  for(i=0;i<nn;i++)
	    {
	      for(j=0;j<nn;j++)
		{
		  target->covariance[i][j] += cov[i][j]*invn;
		}
	    }
	}
    }
}

void adjust_covariance(MYREAL **cov, long size, long n)
{
  const long nn = size;
  const MYREAL invn = 1.0 / (n-1);
  long i,j;
  for(i=0;i<nn;i++)
    {
      for(j=0;j<nn;j++)
	{
	  cov[i][j] *= invn;
	}
    }
}


// testing of the function
//
//
#ifdef FUNCTIONTEST
int myID;

int main(int argc, char **argv)
{
  long i,j;
  long locus = 0;
  MYREAL **cov=NULL;
  MYREAL **corr=NULL;
  world_fmt *world;
  world = (world_fmt *) calloc(1,sizeof(world_fmt));
  MYREAL *params = calloc(3*4, sizeof(MYREAL));
  for (i=0;i<3;i++)
    {
      for(j=0;j<4;j++)
	{
	  params[i*3+j] = (1+i)*(1+j);
	  printf("%f ",params[i*3+j]);
	}
      printf("\n");
    }
      printf("\n");
  world->numpop2 = 4;
  
  covarianceBayes(world, 3, params, 0, locus, &cov);
  correlationBayes(world, locus, cov, &corr);
  printf("Covariance\n");
  for (i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  printf("%f ",cov[i][j]);
	}
      printf("\n");
    }
  printf("\nCorrelation\n");
  for (i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  printf("%f ",corr[i][j]);
	}
      printf("\n");
    }
  printf("\n");
  return 0;
}

#endif
