/*------------------------------------------------------
Maximum likelihood estimation
of migration rate  and effectice population size
using a Metropolis-Hastings Monte Carlo algorithm
-------------------------------------------------------
    slice sampler routines   R O U T I N E S

    Peter Beerli 2006, Tallahassee
    beerli@fsu.edu

    Copyright 2006-2020 Peter Beerli, Tallahassee

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


 $Id$
    */
/*! \file slice.c 
this file implements the slice sampler that can replace the proposal distribution
for parameters as an alternative to the more standard proposal. The code is based on this
mathematica function (created after talking to Paul Lewis):

slice[start_, stick_] := Module[{x = start, r, y, z, z2},
    r = f[x];
y = Random[Real, {0, r}];
r2 = Random[];
z = {x - (1 - r2)*stick, x + r2 * stick};
z2 = Map[f, z];
While[y < z2 [[1]] ,
      z2[[1]] = f[z[[1]] -= stick]];
While[z2 [[2]] >  y, z2[[2]] = f[z[[2]] += stick]];
    x = Random[Real, z];
While[(f[x]) < y && count++ < 100,
      If[(x - z[[1]]) < z[[2]] - x, z[[1]] = x, z[[2]] = x];
      Print[{x, f[x]}];
      x = Random[Real, z];
      ];
      {x, y}]

stick is arbitrary in size, too small to much work, too large to much work
*/

#include <stdlib.h>
#include "migration.h"
#include "bayes.h"
#include "slice.h"
#include "random.h"
#include "tools.h"
#include "sighandler.h"
#include "world.h"
#include "speciate.h"
#ifdef PRETTY
#include "pretty.h"
#endif
#ifdef MPI
#include "migrate_mpi.h"
#else
extern int myID;
#endif
extern MYREAL norm (MYREAL *d, long size);
// functions
void set_slice_param(MYREAL newparam, long which, world_fmt * world);
void set_slice_param_timek(MYREAL newparam, long which, world_fmt * world);
void set_multi_slice_param(MYREAL *newparam, long numparam, world_fmt * world);
void random_direction(MYREAL *coord,long dim);
void sphere_random(MYREAL *lower, MYREAL *values, MYREAL *upper, 
		   MYREAL *r, MYREAL stick, long numparam, MYREAL *minparam, MYREAL *maxparam);
void sphere_x_random(MYREAL *lower, MYREAL *values, MYREAL *upper, 
		     MYREAL *r, long numparam, MYREAL * minparam, MYREAL * maxparam);
void add_stick(MYREAL *values, MYREAL stick, MYREAL *r, long numparam, MYREAL * minparam, MYREAL * maxparam);
//##
// functions implementation
void set_slice_param(MYREAL newparam, long which, world_fmt * world)
{
  MYREAL oldparam;
  const long numpop = world->numpop;  
  const long numpop2 = world->numpop2;  
 
  if(newparam <= world->bayes->minparam[which])
    {
      newparam = world->bayes->minparam[which];
    }
  if(newparam > world->bayes->maxparam[which])
    newparam = world->bayes->maxparam[which];
  
  if(which<numpop)
    {
      if (newparam <=0.0)
	newparam = SMALL_VALUE;
    }
  //if (which > numpop2 + world->bayes->mu)
  //  {
  //    adjust_averagediv(world, which, &newparam);
  //  }
  //else
    bayes_set_param(world->param0,newparam,which,world->options->custm2, world->numpop);
  reprecalc_world(world, which);

  if (which==numpop2 && world->bayes->mu)
    {
      oldparam = world->options->mu_rates[world->locus];
      world->options->mu_rates[world->locus] = newparam;
      //printf("@@@ rate = %f\n",newparam);
      recalc_timelist(world, world->options->mu_rates[world->locus] , oldparam);
      return;
    }
}

void set_slice_param_timek(MYREAL newparam, long which, world_fmt * world)
{
  //MYREAL oldparam;
  //long numpop2 = world->numpop2;  
 
  if(newparam < world->bayes->minparam[which])
    newparam = world->bayes->minparam[which];
  if(newparam > world->bayes->maxparam[which])
    newparam = world->bayes->maxparam[which];

  world->timek[which]=newparam;
  //  reprecalc_world(world, which);
  //      oldparam = world->options->mu_rates[world->locus];
  //  world->options->mu_rates[world->locus] = newparam;
  //  recalc_timelist(world, world->options->mu_rates[world->locus] , oldparam);
}


void set_multi_slice_param(MYREAL *newparam, long numparam, world_fmt * world)
{
  (void) numparam;
  MYREAL oldparam;
  long numpop2 = world->numpop2;  
  // set new parameters and allow for custom migration matrix
  memcpy(world->param0,newparam,sizeof(MYREAL)*(size_t) numpop2);
  precalc_world(world);
  if(world->bayes->mu)
    {
	oldparam = world->options->mu_rates[world->locus];
	world->options->mu_rates[world->locus] = newparam[numpop2];
	recalc_timelist(world, world->options->mu_rates[world->locus] , oldparam);
    }
}


///
/// slice() calculates a new parameter value that comes from the proposal distribution
/// returns log(prob) and replaces startval, slice uses logs internally and so changes the
/// prior distribution 
MYREAL slice (MYREAL *startval, long which, world_fmt * world, MYREAL  (*func) (MYREAL, MYREAL, bayes_fmt *, long))
{
  bayes_fmt * bayes = world-> bayes;
  MYREAL *minparam = bayes->minparam;
  MYREAL *maxparam = bayes->maxparam;
  MYREAL newstartval = *startval;
  MYREAL log_newstartval;
  MYREAL ra;
  MYREAL x;
  MYREAL fx;
  MYREAL y;
  MYREAL r;
  MYREAL zl;
  MYREAL zr;
  MYREAL ezl;
  MYREAL ezr;
  MYREAL ex;
  MYREAL fzl;
  MYREAL fzr;
  MYREAL priorratio;
  register MYREAL stick = world->options->slice_sticksizes[which];
  const MYREAL stick_expand = 1.02;
  const MYREAL stick_reduce = 1.0/1.02;
  register long w;
  long count;
  // find parameter values to honor the custm-migration matrix
  ///if(which < world->numpop2)
  //  {
      if(bayes->map[which][1] == INVALID)
	return -MYREAL_MAX;
      else
	w = bayes->map[which][1];
      // }
      //else
      //{
      //w = which;
      //}
      // if user chose prior bounds that essentially fix value simple return old log(prob)
  if(maxparam[w] -  minparam[w] < SMALLEPSILON)
    {
       return world->bayes->oldval;
    }

  // normal execution pattern
  if((newstartval > maxparam[w]) || (newstartval < minparam[w]))
    {
      fprintf(stdout,"Adjusting startvalue in slice sampler oldval=%f",newstartval);
      newstartval = minparam[w] + RANDUM() * (maxparam[w] - minparam[w]);
      fprintf(stdout," newval=%f\n", newstartval);
    }
  log_newstartval = log(newstartval);
  x = log_newstartval;
  priorratio = (MYREAL) func(newstartval, -1. , bayes, w);
  set_slice_param(newstartval,which, world);
  // calculate the function
  fx = probg_treetimes(world)+ priorratio;
  // find a random point between zero and fx
  ra = log(RANDUM());
  y  = ra + fx;
  // position horizontal stick on coordinates (x,y)
  ra = RANDUM();
  r  = ra * stick;
  zl = x - r;
  zr = zl + stick;
  // extend the stick until it crosses the function
  ezl = EXP(zl);
  ezr = EXP(zr);
  priorratio = func(ezl, -1. , bayes, w);
  set_slice_param(ezl,which, world);
  fzl= probg_treetimes(world) + priorratio;
  priorratio = func(ezr, -1. , bayes, w);
  set_slice_param(ezr,which, world);
  fzr= probg_treetimes(world) + priorratio;
  while(y < fzl && y > -DBL_MAX)
    {
      stick *= stick_expand;
      zl -= stick;
      ezl = EXP(zl);
      priorratio = func(ezl, -1. , bayes, w);
      set_slice_param(ezl,which, world);
      fzl= probg_treetimes(world) + priorratio;
    }
  while(y < fzr && y > -DBL_MAX)
    {
      stick *= stick_expand;
      zr += stick;
      ezr = EXP(zr);
      priorratio = func(ezr, -1. ,bayes, w);
      set_slice_param(ezr,which, world);
      fzr= probg_treetimes(world) + priorratio;
    }
  // pick new value at random
  ra = RANDUM();
  x  = ra * (zr - zl) + zl;
  ex = EXP(x);
  priorratio = func(ex, -1. , bayes, w);
  set_slice_param(ex,which, world);
  fx = probg_treetimes(world) + priorratio;
  count = 0;
  while(fx < y && ((zr-zl) > EPSILON))
    {
      stick *= stick_reduce;
      count++;
      //      if((x - zl) < (zr - x))
      if(x < log_newstartval)
	{
	  zl = x;
	}
      else
	{
	  zr = x;
	}
      ra = RANDUM();
      x  = ra * (zr-zl)+ zl;
      priorratio = func(EXP(x), -1. , bayes, w);
      set_slice_param(EXP(x),which, world);
      fx = probg_treetimes(world) + priorratio;
    }
  //fprintf(stdout,"{{%f,%f},{%f,%f},{%f,%f}} (* y=%f, c=%li, %li, lv=%f *)\n",zl,fzl,x,fx,zr,fzr,y, count, which,startval);
  *startval = EXP(x);
  world->options->slice_sticksizes[which] = stick;
  //fprintf(stdout,"### %li %f\n",which, stick);
  return fx - priorratio;
}


/// slicing for rates of parameters using the times
/// slice() calculates a new parameter value that comes from the proposal distribution
/// returns log(prob) and replaces startval, slice uses logs internally and so changes the
/// prior distribution 
MYREAL sliceRatio (MYREAL *startval, long which, world_fmt * world, MYREAL  (*func) (MYREAL, MYREAL, bayes_fmt *, long))
{
  bayes_fmt * bayes = world-> bayes;
  MYREAL minparam = -4;
  MYREAL maxparam = 4;
  MYREAL newstartval = *startval;
  MYREAL log_newstartval;
  MYREAL ra;
  MYREAL x;
  MYREAL fx;
  MYREAL y;
  MYREAL r;
  MYREAL zl;
  MYREAL zr;
  MYREAL ezl;
  MYREAL ezr;
  MYREAL ex;
  MYREAL fzl;
  MYREAL fzr;
  MYREAL priorratio;
  register MYREAL stick = 2.0;
  //  const MYREAL stick_expand = 1.02;
  //const MYREAL stick_reduce = 1.0/1.02;
  register long w;
  long count;
  long wbase;
  // find parameter values to honor the custm-migration matrix
  wbase = (long) (which / (world->timeelements * (world->numpop2+1)));
  if (wbase == 0)
    {
      wbase = 0;
    }
  else
    {
      wbase = wbase * (world->timeelements * (world->numpop2+1));
      which = which - wbase;
    }
  /*  if(which < world->numpop2)
    {
      if(bayes->map[which][1] == INVALID)
	return -MYREAL_MAX;
      else
	w = bayes->map[which][1];
    }
  else
    {
      w = which;
      }*/
  w=which;
      // if user chose prior bounds that essentially fix value simple return old log(prob)
  if(fabs(maxparam - minparam) > (double) FLT_EPSILON)
    {
       return world->bayes->oldval;
    }

  // normal execution pattern
  log_newstartval = log(newstartval);
  if((log_newstartval > maxparam) || (log_newstartval < minparam))
    {
      fprintf(stdout,"Adjusting startvalue in slice sampler oldval=%f",newstartval);
      log_newstartval = minparam + RANDUM() * (maxparam - minparam);
      fprintf(stdout," log_newval=%f\n", log_newstartval);
      newstartval = exp(log_newstartval);
    }
  x = log_newstartval;
  priorratio = 1.0;//(MYREAL) func(newstartval, -1. , bayes, w);
  set_slice_param_timek(newstartval,wbase + which, world);
  // calculate the function
  fx = probg_treetimes(world)+ priorratio;
  // find a random point between zero and fx
  ra = log(RANDUM());
  y  = ra + fx;
  // position horizontal stick on coordinates (x,y)
  ra = RANDUM();
  r  = ra * stick;
  zl = x - r;
  zr = zl + stick;
  // extend the stick until it crosses the function
  ezl = EXP(zl);
  ezr = EXP(zr);
  priorratio = 1.0;//func(ezl, -1. , bayes, w);
  set_slice_param_timek(ezl,wbase+which, world);
  fzl= probg_treetimes(world) + priorratio;
  priorratio = func(ezr, -1. , bayes, w);
  set_slice_param_timek(ezr,wbase+which, world);
  fzr= probg_treetimes(world) + priorratio;
  while(y < fzl && y > -DBL_MAX)
    {
      //stick *= stick_expand;
      zl -= stick;
      ezl = EXP(zl);
      priorratio = 1.0;//func(ezl, -1. , bayes, w);
      set_slice_param_timek(ezl,wbase+which, world);
      fzl= probg_treetimes(world) + priorratio;
    }
  while(y < fzr && y > -DBL_MAX)
    {
      //stick *= stick_expand;
      zr += stick;
      ezr = EXP(zr);
      priorratio = 1.0;//func(ezr, -1. ,bayes, w);
      set_slice_param_timek(ezr,wbase+which, world);
      fzr= probg_treetimes(world) + priorratio;
    }
  // pick new value at random
  ra = RANDUM();
  x  = ra * (zr - zl) + zl;
  ex = EXP(x);
  priorratio = 1.0;//func(ex, -1. , bayes, w);
  set_slice_param_timek(ex,wbase+which, world);
  fx = probg_treetimes(world) + priorratio;
  count = 0;
  while(fx < y && ((zr-zl) > EPSILON))
    {
      //      stick *= stick_reduce;
      count++;
      //      if((x - zl) < (zr - x))
      if(x < log_newstartval)
	{
	  zl = x;
	}
      else
	{
	  zr = x;
	}
      ra = RANDUM();
      x  = ra * (zr-zl)+ zl;
      priorratio = 1.0;//func(EXP(x), -1. , bayes, w);
      set_slice_param_timek(EXP(x),wbase+which, world);
      fx = probg_treetimes(world) + priorratio;
    }
  //fprintf(stdout,"{{%f,%f},{%f,%f},{%f,%f}} (* y=%f, c=%li, %li, lv=%f *)\n",zl,fzl,x,fx,zr,fzr,y, count, which,startval);
  *startval = EXP(x);
  //  world->options->slice_sticksizes[wbase+which] = stick;
  //fprintf(stdout,"### %li %f\n",which, stick);
  return fx - priorratio;
}

///
/// slice() calculates a new parameter value  that comes from the proposal distribution
/// returns the log(probability) and changes startval]
MYREAL expslice (MYREAL *startval, long which, world_fmt * world, MYREAL  (*func) (MYREAL, MYREAL, bayes_fmt *, long))
{
  register long w;
  bayes_fmt * bayes = world-> bayes;
  const MYREAL *minparam = bayes->minparam;
  const MYREAL *maxparam = bayes->maxparam;
  MYREAL newstartval = *startval;
  register MYREAL ra;
  register MYREAL x;
  register MYREAL fx;
  register MYREAL y;
  register MYREAL r;
  register MYREAL zl;
  register MYREAL zr;
  register MYREAL fzl;
  register MYREAL fzr;
  MYREAL priorratio;
  long count = 3L;
  register MYREAL   stick = world->options->slice_sticksizes[which];
  const MYREAL stick_expand = 1.02;
  const MYREAL stick_reduce = 1.0/1.02;


  // find parameter values to honor the custm-migration matrix
  //  if(which < world->numpop2)
  //  {
      if(bayes->map[which][1] == INVALID)
	return -MYREAL_MAX;
      else
	w = bayes->map[which][1];	  
      //  }
      //else
      //{
      //w = which;
      //}
  // if user chose prior bounds that essentially fix value simple return old log(prob)
  if(maxparam[w] -  minparam[w] < SMALLEPSILON)
    {
       return world->bayes->oldval;
    }
  // normal execution pattern
  // replace the parameter value at which with startval
  if((newstartval > maxparam[w]) || (newstartval < minparam[w]))
  {
    // in principle we never should come here, but the combine parameters particularly xNm 
    // may lead to here
    fprintf(stdout,"Adjusting startvalue in slice sampler:\n   Out of bound  start value =%f\n",newstartval);
    newstartval = minparam[w] + RANDUM() * (maxparam[w] - minparam[w]);
    fprintf(stdout,"    In bound new value =%f\n", newstartval);
  }
  x = newstartval;
  priorratio = (MYREAL) func(newstartval, -1. , bayes, w);
  set_slice_param(newstartval,which, world);
  // calculate the function
  fx = probg_treetimes(world)+ priorratio;
  // find a random point between zero and fx
  ra = log(RANDUM());
  y  = ra + fx;
  // position horizontal stick on coordinates (x,y) --> |-------x--|
  ra = RANDUM();
  r  = ra * stick;
  zl = x - r;
  zr = zl + stick;
  // extend the stick until it crosses the function
  priorratio = func(zl, -1. , bayes, w);
  set_slice_param(zl,which, world);
  fzl= probg_treetimes(world) + priorratio;
  priorratio = func(zr, -1. , bayes, w);
  set_slice_param(zr,w, world);
  fzr= probg_treetimes(world) + priorratio;
  // 
  while(y < fzl && zl>minparam[w])
    {
      stick *= stick_expand;
      zl -= stick;
      //      ezl = EXP(zl);
      priorratio = func(zl, -1. , bayes, w);
      set_slice_param(zl,which, world);
      fzl= probg_treetimes(world) + priorratio;
      count++;
#ifdef SLICEREPORTER
  fprintf(stdout,"#@ %f %f %s %li\n", zl, fzl,"lower", count);
#endif
    }
  while(y < fzr && zr< maxparam[w])
    {
      stick *= stick_expand;
      zr += stick;
      //      ezr = EXP(zr);
      priorratio = func(zr, -1. ,bayes, w);
      set_slice_param(zr,which, world);
      fzr= probg_treetimes(world) + priorratio;
      count++;
#ifdef SLICEREPORTER
  fprintf(stdout,"#@ %f %f %s %li\n", zr, fzr,"upper", count);
#endif
    }
  // pick new value at random
  ra = RANDUM();
  x  = ra * (zr - zl) + zl;
  //  ex = EXP(x);
  priorratio = func(x, -1. , bayes, w);
  set_slice_param(x,which, world);
  fx = probg_treetimes(world) + priorratio;
  count++;
#ifdef SLICEREPORTER
  fprintf(stdout,"#@ %f %f %s %li\n", x, fx,"new_x", count);
#endif
  while(fx < y && ((zr-zl) > EPSILON))
    {
      stick *= stick_reduce;
      count++;
      if(x < newstartval)
	{
	  zl = x;
#ifdef SLICEREPORTER
  fprintf(stdout,"#@ %f %f %s %li\n", x, fx,"new_lower", count);
#endif
	}
      else
	{
	  zr = x;
#ifdef SLICEREPORTER
  fprintf(stdout,"#@ %f %f %s %li\n", x, fx,"new_upper", count);
#endif
	}
      ra = RANDUM();
      x  = ra * (zr-zl)+ zl;
      priorratio = func(x, -1. , bayes, w);
      set_slice_param(x,which, world);
      fx = probg_treetimes(world) + priorratio;
#ifdef SLICEREPORTER
      fprintf(stdout,"#@ %f %f %s %li\n", x, fx,"new_x", count);
#endif
    }
  //fprintf(stdout,"{{%f,%f},{%f,%f},{%f,%f}} (* y=%f, c=%li, %li, lv=%f *)\n",zl,fzl,x,fx,zr,fzr,y, count, which,startval);
  //fprintf(stdout,"%3li: %f %f %li\n",which, x, fx, count);
#ifdef SLICEREPORTER
  fprintf(stdout,"#@ %f %f %s %li\n#----------------------\n", x, fx,"return_x", count);
#endif
  if((x < maxparam[w]) && (x > minparam[w]))
    {
      *startval = x;
      world->options->slice_sticksizes[which] = stick;
      return fx - priorratio;
    }
  else
    {
      if(x < minparam[w])
	{
	  *startval = minparam[w];
	  world->options->slice_sticksizes[which] = stick;
	  return fx - priorratio;
	}
      else
	{
	  *startval = maxparam[w];
	  world->options->slice_sticksizes[which] = stick;
	  return fx - priorratio;
	}
    }
}

//#ifdef TESTING2
// log_prior_ratioall needs a body in bayes.c
void random_direction(MYREAL *coord,long dim)
{
  long i;
  MYREAL ci;
  MYREAL norm=0.0;
  for(i=0;i<dim;i++)
    {
      ci = rannor(0.0,1.0);
      coord[i] = ci;
      norm += ci * ci;
    }
  norm = sqrt(norm);
  for(i=0;i<dim;i++)
    {
      coord[i] /= norm;
    }
}

void sphere_random(MYREAL *lower, MYREAL *values, MYREAL *upper, 
		   MYREAL *r, MYREAL stick, long numparam, MYREAL *minparam, MYREAL *maxparam)
{
  long i;
  random_direction(r,numparam);
  for(i=0;i<numparam;i++)
    {
      lower[i] = values[i] - (r[i] * stick);
      if(lower[i] < minparam[i])
	lower[i]=minparam[i];
      upper[i] = lower[i] + stick;
      if(upper[i]>maxparam[i])
	upper[i] = maxparam[i];
    }
}

void sphere_x_random(MYREAL *lower, MYREAL *values, MYREAL *upper, 
		   MYREAL *r, long numparam, MYREAL * minparam, MYREAL * maxparam)
{
  long i;
  random_direction(r,numparam);
  for(i=0;i<numparam;i++)
    {
      values[i] = (RANDUM() * (upper[i]-lower[i]) + lower[i]);
      if(values[i] < minparam[i])
	values[i]=minparam[i];
      else
	{
	  if(values[i] > maxparam[i])
	    values[i] = maxparam[i];
	}
    }
}

void add_stick(MYREAL *values, MYREAL stick, MYREAL *r, long numparam, MYREAL * minparam, MYREAL * maxparam)
{
  long i;
  for(i=0;i<numparam;i++)
    {
      values[i] = values[i] + (r[i] * stick);
    }
  if(values[i] < minparam[i])
    values[i]=minparam[i];
  else
    {
      if(values[i] > maxparam[i])
	values[i] = maxparam[i];
    }
}

/*
///
/// slice() calculates a new point that comes from the proposal distribution changing all parameters
/// at once using Michal's sphere approach 
void expallslice (MYREAL * values, MYREAL *ratevalue, long which, world_fmt * world)
{
  long i;
  long numpop2 = world->numpop2;
  long numparam = numpop2 + world->bayes->mu * 1;// we use slice per locus so we do not need
  // to keep track of all loci at once
  MYREAL norm_newstartval;
  MYREAL *minparam = world->bayes->minparam;
  MYREAL *maxparam = world->bayes->maxparam;
  MYREAL *newstartval;
  //  MYREAL log_newstartval;
  MYREAL ra;
  MYREAL *x;
  MYREAL fx;
  MYREAL y;
  MYREAL *r;
  MYREAL *zl;
  MYREAL *zr;
  //  MYREAL ezl;
  //  MYREAL ezr;
  //  MYREAL ex;
  MYREAL *v;
  MYREAL fzl;
  MYREAL fzr;
  MYREAL priorratio=1.0;
  MYREAL stick = 1.0; //this needs a remedy! (maxparam[which] - minparam[which]) / 20.;
  long count = 3L;

  v = (MYREAL *) mycalloc(numparam * 5, sizeof(MYREAL));
  r = v;
  x = r + numparam;
  zl = x + numparam;
  zr = zl + numparam;
  newstartval = zr + numparam;
  memcpy(newstartval,world->param0,sizeof(MYREAL)*world->numpop2);
  if(world->bayes->mu)
    newstartval[world->numpop2] = world->options->mu_rates[world->locus];
  for(i=0; i < numparam; i++)
    {
      // safeguard against slicing towards zero, eplxore whether zero is harmful, should not be!
      //      if(minparam[i] < SMALL_VALUE)
      //	minparam[i] = SMALL_VALUE;
      // replace the parameter value at which with startval
      //      if((newstartval[i] > maxparam[i]) || (newstartval[i] < minparam[i]))
      //{
      //
      //  newstartval[i] = minparam[i] + RANDUM() * (maxparam[i] - minparam[i]);
      //  fprintf(stdout,"Adjusting startvalue in slice sampler oldval[%li]=%f",i,startval[i]);
      //  fprintf(stdout," newval[%li]=%f\n", i, newstartval[i]);
      //}
      //else
      //	{
      //  newstartval[i] = startval[i];
	  //	}
      x[i] = newstartval[i];
    }
  //  norm_newstartval = norm(newstartval,numparam, minparam, maxparam);
  norm_newstartval = norm(newstartval,numparam);
  //  priorratio = (MYREAL) log_prior_ratio_all(world,x);
  error("needs implementation for log_prior_ratio_all");
  set_multi_slice_param(newstartval,numparam, world);
  // calculate the function
  fx = probg_treetimes(world)+ priorratio;
  // find a random point between zero and fx
  ra = log(RANDUM());
  y  = ra + fx;
  // position horizontal stick on coordinates (x1,x2,...,xn,y)
  sphere_random(zl,x,zr,r,stick, numparam, minparam, maxparam);//zl = x - (r=sphererand()) * stick, zr = zl + stick
  //  zl = x - r;
  //zr = zl + stick;
  // extend the stick until it crosses the function
  //  ezl = EXP(zl);
  //ezr = EXP(zr);
  priorratio = log_prior_ratio_all(world,zl);
  set_multi_slice_param(zl, numparam, world);
  fzl= probg_treetimes(world) + priorratio;
  priorratio = log_prior_ratio_all(world,zr);
  set_multi_slice_param(zr, numparam, world);
  fzr= probg_treetimes(world) + priorratio;
  while(y < fzl)
    {
      add_stick(zl, -stick, r, numparam, minparam, maxparam);//zl -= stick;
      //      ezl = EXP(zl);
      priorratio = log_prior_ratio_all(world,zl);
      set_multi_slice_param(zl,numparam, world);
      fzl= probg_treetimes(world) + priorratio;
      count++;
    }
  while(y < fzr)
    {
      add_stick(zr, stick, r, numparam, minparam, maxparam); //zr += stick;
      //      ezr = EXP(zr);
      priorratio = log_prior_ratio_all(world,zr);
      set_multi_slice_param(zr,which, world);
      fzr= probg_treetimes(world) + priorratio;
      count++;
    }
  // pick new value at random
  sphere_x_random(zl,x,zr,r, numparam, minparam, maxparam);
  //ra = RANDUM();
  //x  = ra * (zr - zl) + zl;
  //  ex = EXP(x);
  priorratio = log_prior_ratio_all(world,x);
  set_multi_slice_param(x,which, world);
  fx = probg_treetimes(world) + priorratio;
  count++;
  while(fx < y && ((norm(zr,numparam)-norm(zl,numparam)) > EPSILON))
    {
      count++;
      //      if((x - zl) < (zr - x))
      if(norm(x,numparam)<norm_newstartval)//vector comparison?
	{
	  memcpy(zl,x,sizeof(MYREAL)*numpop2);
	}
      else
	{
	  memcpy(zr,x,sizeof(MYREAL)*numpop2);
	}
      sphere_x_random(zl,x, zr,r, numparam, minparam, maxparam);
      //      ra = RANDUM();
      //x  = ra * (zr-zl)+ zl;
      priorratio = log_prior_ratio_all(world,x);
      set_multi_slice_param(x,which, world);
      fx = probg_treetimes(world) + priorratio;
    }
  //fprintf(stdout,"{{%f,%f},{%f,%f},{%f,%f}} (* y=%f, c=%li, %li, lv=%f *)\n",zl,fzl,x,fx,zr,fzr,y, count, which,startval);
  //fprintf(stdout,"%3li: %f %f %li\n",which, x, fx, count);
  memcpy(values,x,sizeof(MYREAL)*world->numpop2);
  if(world->bayes->mu)
    *ratevalue = x[world->numpop2];
  myfree(v);
}
//#endif

*/
