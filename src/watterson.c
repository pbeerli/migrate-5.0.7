/*
 *  watterson.c
 *  migrate-n
 *
 *  Created by Peter Beerli on June 26 2006.
 *  Copyright 2006 Peter Beerli. All rights reserved.
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
 *
 */

#include "watterson.h"
double watterson_a(long n);
double watterson_b(long n);


double watterson_a(long n)
{
  long k;
  double sum=0.;
  for(k=1; k<n; k++)
    {
      sum += 1. / k;
    }
  return sum;
}
double watterson_b(long n)
{
  long k;
  double sum=0.;
  for(k=1; k<n; k++)
    {
      sum += 1. / (k*k);
    }
  return sum;
}

MYREAL watterson(long segreg, long n)
{
  return (MYREAL) (segreg / watterson_a(n));
}

//wattersonvar[s_, n_] :=
//((an = (Sum[1./k, {k, 1,
// n}]))(thetaw = watterson[s, n]) +  (bn = (Sum[
// 1./(k^2), {k, 1, n}]))thetaw^2 )/(an^2)
MYREAL wattersonvar(double  thetaw, long n)
{
  double an = watterson_a(n);
  double bn = watterson_b(n);
  return (MYREAL) ((an * thetaw + bn * thetaw * thetaw) / (an*an));
}

MYREAL wattersonstd(double  thetaw, long n)
{
  return (MYREAL)(sqrt(wattersonvar(thetaw,n)));
}


