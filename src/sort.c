/* ----------------------------------------------------- */
/* sort.c                       */
/* comparison routines for various qsorts and bsearch's */
/* ----------------------------------------------------- */

/*                              */
/* P. Beerli                        */
/*
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2004 Peter Beerli, Tallahassee FL
 
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

 
$Id: sort.c 2158 2013-04-29 01:56:20Z beerli $ */
/* ----------------------------------------------------- */
/*! \file sort.c */
#pragma clang diagnostic ignored "-Wcast-qual"

/* include files */
#include "migration.h"
#include "sort.h"
/* private functions */

/* public functions */
int
charcmp (const void *v1, const void *v2)
{
    if (*(const char *) v1 < *(const char *) v2)
    {
        return -1;
    }
    else
    {
        if (*(const char *) v1 > *(const char *) v2)
        {
            return 1;
        }
        else
            return 0;
    }
}

int
stringcmp (const void *v1, const void *v2)
{
    if (strcmp ((const char *) v1, (const char *) v2) < 0)
    {
        return -1;
    }
    else
    {
        if (strcmp ((const char *) v1, (const char *) v2) > 1)
            return 1;
        else
            return 0;
    }
}

int
numcmp (const void *v1, const void *v2)
{
    if (*(const MYREAL *) v1 < *(const MYREAL *) v2)
    {
        return -1;
    }
    else
    {
        if (*(const MYREAL *) v1 > *(const MYREAL *) v2)
        {
            return 1;
        }
        else
            return 0;
    }
}

int
paircmp (const void *v1, const void *v2)
{
    if (((const MYREAL *) v1)[1] < ((const MYREAL *) v2)[1])
    {
        return -1;
    }
    else
    {
        if (((const MYREAL *) v1)[1] > ((const MYREAL *) v2)[1])
        {
            return 1;
        }
        else
            return 0;
    }
}

int
paircmp_first (const void *v1, const void *v2)
{
    if (((const MYREAL *) v1)[0] < ((const MYREAL *) v2)[0])
    {
        return -1;
    }
    else
    {
        if (((const MYREAL *) v1)[0] > ((const MYREAL *) v2)[0])
        {
            return 1;
        }
        else
            return 0;
    }
}



/* by the second element, ties by the first */
static int
paircmp_second_first (const void *v1, const void *v2)
{
    const int c = paircmp (v1, v2);
    return c != 0 ? c : paircmp_first (v1, v2);
}

/// sort pairs by their second element (the histogram height in
/// calc_hpd_credibility()), ties by the first (the bin). The former second
/// pass meant to order ties re-sorted, at every tie, the run of DIFFERENT
/// heights before it by the first element (and its tie test was an absolute
/// DBL_EPSILON), so the order and with it the HPD sets depended on the scale
/// of the heights.
void paired_qsort2(pair *x, long xelem)
{
    qsort(x,(size_t) xelem,sizeof(pair),paircmp_second_first);
}


int
floatcmp (const void *v1, const void *v2)
{
    if (*(const float *) v1 < *(const float *) v2)
    {
        return -1;
    }
    else
    {
        if (*(const float *) v1 > *(const float *) v2)
        {
            return 1;
        }
        else
            return 0;
    }
}

int
doublecmp (const void *v1, const void *v2)
{
    if (*(const double *) v1 < *(const double *) v2)
    {
        return -1;
    }
    else
    {
        if (*(const double *) v1 > *(const double *) v2)
        {
            return 1;
        }
        else
            return 0;
    }
}


int
longcmp (const void *v1, const void *v2)
{
    if (*(const long *) v1 < *(const long *) v2)
    {
        return -1;
    }
    else
    {
        if (*(const long *) v1 > *(const long *) v2)
        {
            return 1;
        }
        else
            return 0;
    }
}

int
intcmp (const void *v1, const void *v2)
{
    if (*(const int *) v1 < *(const int *) v2)
    {
        return -1;
    }
    else
    {
        if (*(const int *) v1 > *(const int *) v2)
        {
            return 1;
        }
        else
            return 0;
    }
}
 
int
agecmp (const void *x, const void *y)
{
  const vtlist * xx = (const vtlist *) x;
  const vtlist * yy = (const vtlist *) y;
  MYREAL xa = xx->age;
  MYREAL ya = yy->age;

  if (xa < ya)
    {
        return -1;
    }
    else
    {
      if (xa > ya)
        {
            return 1;
        }
        else
            return 0;
    }
}

int
delcmp (const void *x, const void *y)
{
    if ((*((const node **) x))->id < (*((const node **) y))->id)
    {
        return -1;
    }
    else
    {
        if ((*((const node **) x))->id > (*((const node **) y))->id)
        {
            return 1;
        }
        else
            return 0;
    }
}

int
migr_time_cmp (const void *x, const void *y)
{
    if (((const migr_table_fmt *) x)->time < ((const migr_table_fmt *) y)->time)
    {
        return -1;
    }
    else
    {
      if (((const migr_table_fmt *) x)->time > ((const migr_table_fmt *) y)->time)
        {
            return 1;
        }
        else
            return 0;
    }
}


int
searchagecmp (const void *x, const void *y)
{
    const MYREAL xx = (const MYREAL) *((const MYREAL *) x);
    const MYREAL age = (const MYREAL) ((const vtlist *) y)->age;

    if (xx < age)
    {
        return -1;
    }
    else
    {
      if (xx > age)
        {
            return 0;
        }
        else
            return 1;
    }

}


int
aiccmp (const void *x, const void *y)
{
    if (((const aic_fmt *) x)->aic < ((const aic_fmt *) y)->aic)
    {
        return -1;
    }
    else
    {
      if (((const aic_fmt *) x)->aic > ((const aic_fmt *) y)->aic)
        {
            return 1;
        }
        else
            return 0;
    }
}
