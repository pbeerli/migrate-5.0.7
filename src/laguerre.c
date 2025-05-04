/*
   (mutation) rates following a Gamma distribution
   using orthogonal polynomials for finding rates and
   LOG(probabilities)
   [based on initgammacat of Joe Felsenstein]
 
   - Generalized Laguerre (routine by Joe Felsenstein 2000)
     defining points for a Gamma distribution with 
     shape parameter alpha and location parameter beta=1/alpha
     [mean=1, std = 1/alpha^2]
   - Hermite (approximates a normal and is activated when 
     the shape parameter alpha is > 100.)
 
   Part of Migrate
   http://popgen.csit.fsu.edu/migrate.html
 
   Peter Beerli, Seattle 2001
   
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2005 Peter Beerli, Tallahassee FL
 
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
 
$Id: laguerre.c 2067 2012-07-27 20:59:32Z beerli $
*/
/* \file laguerre.c
Calculates the Laguerre Quadrature points, used for gamma deviated site rate variation

*/
#include "migration.h"
#include "laguerre.h"
#include "tools.h"
#include "sighandler.h"

#define SQRTPI 1.7724538509055160273
#define SQRT2  1.4142135623730950488


/*this triggers the test  main()
and is called with
gcc -DLAGUERRE_TEST -g laguerre.c -o laguerre -lm*/

#ifdef LAGUERRE_TEST
/* at the end is a test main to help test if the root/weights finding
   is OK*/

/* if machine has lgamma() use it otherwise use lgamma from
   tools.h*/
#undef LGAMMA
#define LGAMMA lgamma


/* for migrate this is defined in tools.h */
MYREAL
logfac (long n)
{
    /* log(n!) values were calculated with Mathematica
       with a precision of 30 digits */
    switch (n)
    {
    case 0:
        return 0.;
    case 1:
        return 0.;
    case 2:
        return 0.693147180559945309417232121458;
    case 3:
        return 1.791759469228055000812477358381;
    case 4:
        return 3.1780538303479456196469416013;
    case 5:
        return 4.78749174278204599424770093452;
    case 6:
        return 6.5792512120101009950601782929;
    case 7:
        return 8.52516136106541430016553103635;
    case 8:
        return 10.60460290274525022841722740072;
    case 9:
        return 12.80182748008146961120771787457;
    case 10:
        return 15.10441257307551529522570932925;
    case 11:
        return 17.50230784587388583928765290722;
    case 12:
        return 19.98721449566188614951736238706;
    default:
        return LGAMMA (n + 1.);
    }
}
#endif

/* prototypes */
MYREAL hermite (long n, MYREAL x);
void root_hermite (long n, MYREAL *hroot);
MYREAL halfroot (MYREAL (*func) (long m, MYREAL x),
                 long n, MYREAL startx, MYREAL delta);
void hermite_weight (long n, MYREAL *hroot, MYREAL *weights);
void inithermitcat (long categs, MYREAL alpha, MYREAL theta1,
                    MYREAL *rate, MYREAL *probcat);

MYREAL glaguerre (long m, MYREAL b, MYREAL x);
void initlaguerrecat (long categs, MYREAL alpha, MYREAL theta1,
                      MYREAL *rate, MYREAL *probcat);
void roots_laguerre (long m, MYREAL b, MYREAL **lgroot);

void initgammacat (long categs, MYREAL alpha, MYREAL theta1,
                   MYREAL *rate, MYREAL *probcat);

void integrate_laguerre (long categs, MYREAL *rate,
                         MYREAL *probcat,
                         MYREAL (*func) (MYREAL theta, helper_fmt * b),
                         helper_fmt * helper, MYREAL *result, MYREAL *rmax);


/*------------------------------------------------------
  Generalized Laguerre polynomial computed recursively.
  For use by initgammacat 
*/
MYREAL
glaguerre (long m, MYREAL b, MYREAL x)
{
    long i;
    MYREAL gln, glnm1, glnp1; /* L_n, L_(n-1), L_(n+1) */

    if (m == 0)
        return 1.0;
    else
    {
        if (m == 1)
            return 1.0 + b - x;
        else
        {
            gln = 1.0 + b - x;
            glnm1 = 1.0;
            for (i = 2; i <= m; i++)
            {
                glnp1 =
                    ((2 * (i - 1) + b + 1.0 - x) * gln - (i - 1 + b) * glnm1) / i;
                glnm1 = gln;
                gln = glnp1;
            }
            return gln;
        }
    }
}    /* glaguerre */


/* calculates hermite polynomial with degree n and parameter x */
/* seems to be unprecise for n>13 -> root finder does not converge*/
MYREAL
hermite (long n, MYREAL x)
{
    MYREAL h1 = 1.;
    MYREAL h2 = 2. * x;
    MYREAL xx = 2. * x;
    long i;
    for (i = 1; i < n; i++)
    {
        xx = 2. * x * h2 - 2. * (i) * h1;
        h1 = h2;
        h2 = xx;
    }
    return xx;
}

void
root_hermite (long n, MYREAL *hroot)
{
    long z = 0;
    long ii;
    long start;
    if (n % 2 == 0)
    {
        start = n / 2;
        z = 1;
    }
    else
    {
        start = n / 2 + 1;
        z = 2;
        hroot[start - 1] = 0.0;
    }
    for (ii = start; ii < n; ii++)
    {
        /* search only upwards */
        hroot[ii] = halfroot (hermite, n, hroot[ii - 1] + EPSILON, 1. / n);
        hroot[start - z] = -hroot[ii];
        z++;
    }
}

/*searches from the bound (startx) only in one direction
  (by positive or negative delta, which results in 
  other-bound=startx+delta)
  delta should be small.
  (*func) is a function with two arguments
*/
MYREAL
halfroot (MYREAL (*func) (long m, MYREAL x),
          long n, MYREAL startx, MYREAL delta)
{
    MYREAL xl;
    MYREAL xu;
    MYREAL xm = 0.;
    MYREAL fu;
    MYREAL fl;
    MYREAL fm = 100000.;
    MYREAL gradient;
    boolean down = FALSE;
    /* decide if we search above or below startx and escapes to trace back
       to the starting point that most often will be 
       the root from the previous calculation */
    if (delta < 0)
    {
        xu = startx;
        xl = xu + delta;
    }
    else
    {
        xl = startx;
        xu = xl + delta;
    }
    delta = fabs (delta);
    fu = (*func) (n, xu);
    fl = (*func) (n, xl);
    gradient = (fl - fu) / (xl - xu);

    while (fabs (fm) > EPSILON)
    {
        /* is root outside of our bracket? */
        if ((fu < 0.0 && fl < 0.0) || (fu > 0.0 && fl > 0.0))
        {
            xu += delta;
            fu = (*func) (n, xu);
            fl = (*func) (n, xl);
            gradient = (fl - fu) / (xl - xu);
            down = gradient < 0 ? TRUE : FALSE;
        }
        else
        {
            xm = xl - fl / gradient;
            fm = (*func) (n, xm);
            if (down)
            {
                if (fm > 0.)
                {
                    xl = xm;
                    fl = fm;
                }
                else
                {
                    xu = xm;
                    fu = fm;
                }
            }
            else
            {
                if (fm > 0.)
                {
                    xu = xm;
                    fu = fm;
                }
                else
                {
                    xl = xm;
                    fl = fm;
                }
            }
            gradient = (fl - fu) / (xl - xu);
        }
    }
    return xm;
}


// calculate the weights for the hermite polynomial
// at the roots
// using formula Abramowitz and Stegun chapter 25.4.46 p.890
void
hermite_weight (long n, MYREAL *hroot, MYREAL *weights)
{
    long i;
    MYREAL hr2;
    MYREAL nominator = EXP (LOG2 * (n - 1.) + logfac (n)) * SQRTPI / (n * n);
    for (i = 0; i < n; i++)
    {
        hr2 = hermite (n - 1, hroot[i]);
        weights[i] = nominator / (hr2 * hr2);
    }
}

/* calculates rates and LOG(probabilities) */
void
inithermitcat (long categs, MYREAL alpha, MYREAL theta1,
               MYREAL *rate, MYREAL *probcat)
{
    long i;
    MYREAL *hroot;
    MYREAL std = SQRT2 * theta1 / sqrt (alpha);
    hroot = (MYREAL *) mycalloc (categs + 1, sizeof (MYREAL));
    root_hermite (categs, hroot); // calculate roots
    hermite_weight (categs, hroot, probcat); // set weights
    for (i = 0; i < categs; i++) // set rates
    {
        rate[i] = theta1 + std * hroot[i];
        probcat[i] = LOG (probcat[i]);
    }
    myfree(hroot);
}

///
/// For use by initgammacat().
/// Get roots of m-th Generalized Laguerre polynomial, given roots 
/// of (m-1)-th, these are to be stored in lgroot[m][]
void
roots_laguerre (long m, MYREAL b, MYREAL **lgroot)
{
    long i;
    long count=0;
    MYREAL upperl, lower, x, y;
    boolean dwn = FALSE;
    //MYREAL tmp;
    /* is function declining in this interval? */
    if (m == 1)
    {
        lgroot[1][1] = 1.0 + b;
    }
    else
    {
        dwn = TRUE;
        for (i = 1; i <= m; i++)
        {
            if (i < m)
            {
                if (i == 1)
                    lower = 0.0;
                else
                    lower = lgroot[m - 1][i - 1];
                upperl = lgroot[m - 1][i];
            }
            else
            {   /* i == m, must search above */
                lower = lgroot[m - 1][i - 1];
                x = lgroot[m - 1][m - 1];
                do
                {
                    x = 2.0 * x;
                    y = glaguerre (m, b, x);
                }
                while ((dwn && (y > 0.0)) || ((!dwn) && (y < 0.0)));
                upperl = x;
            }
            count = 0;
            while (upperl - lower > 0.000000001 && count++  < 1000)
            {
                x = (upperl + lower) / 2.0;
                if (glaguerre (m, b, x) > 0.0)
                {
                    if (dwn)
                        lower = x;
                    else
                        upperl = x;
                }
                else
                {
                    if (dwn)
                        upperl = x;
                    else
                        lower = x;
                }
            }
            lgroot[m][i] = (lower + upperl) / 2.0;
            dwn = !dwn;  /* switch for next one */
        }
    }
}    /* root_laguerre */


void
initgammacat (long categs, MYREAL alpha, MYREAL theta1,
              MYREAL *rate, MYREAL *probcat)
{
    /* calculate rates and probabilities to approximate Gamma distribution
       of rates with "categs" categories and shape parameter "alpha" using
       rates and weights from Generalized Laguerre quadrature */

    if (alpha >= 100.)
    {
        inithermitcat (categs, alpha, theta1, rate, probcat);
    }
    else
    {
        initlaguerrecat (categs, alpha, theta1, rate, probcat);
    }
}

void
initlaguerrecat (long categs, MYREAL alpha, MYREAL theta1, MYREAL *rate,
                 MYREAL *probcat)
{
    long i;
    MYREAL **lgroot;  /* roots of GLaguerre polynomials */
    MYREAL f, x, xi, y;

    lgroot = (MYREAL **) mycalloc (categs + 1, sizeof (MYREAL *));
    lgroot[0] =
      (MYREAL *) mycalloc (((categs + 1) * (categs + 1)), sizeof (MYREAL));
    for (i = 1; i < categs + 1; i++)
    {
        lgroot[i] = lgroot[0] + i * (categs + 1);
    }
    lgroot[1][1] = 1.0 + alpha;
    for (i = 2; i <= categs; i++)
        roots_laguerre (i, alpha, lgroot); /* get roots for L^(a)_n */
    /* here get weights */
    /* Gamma weights are
       (1+a)(1+a/2) ... (1+a/n)*x_i/((n+1)^2 [L_{n+1}^a(x_i)]^2)  */
    f = 1;
    for (i = 1; i <= categs; i++)
        f *= (1.0 + alpha / i);
    for (i = 1; i <= categs; i++)
    {
        xi = lgroot[categs][i];
        y = glaguerre (categs + 1, alpha, xi);
        x = f * xi / ((categs + 1) * (categs + 1) * y * y);
        rate[i - 1] = xi / (1.0 + alpha);
        probcat[i - 1] = x;
    }
    for (i = 0; i < categs; i++)
    {
        probcat[i] = LOG (probcat[i]);
        rate[i] *= theta1;
    }
    myfree(lgroot[0]);
    myfree(lgroot);
}    /* initgammacat */


void
integrate_laguerre (long categs, MYREAL *rate,
                    MYREAL *probcat,
                    MYREAL (*func) (MYREAL theta, helper_fmt * b),
                    helper_fmt * helper, MYREAL *result, MYREAL *rmax)
{
    MYREAL summ = 0.;
    long i;
    MYREAL *temp;
    int *stemp;
    *rmax = -MYREAL_MAX;

    temp = (MYREAL *) mycalloc (categs, sizeof (MYREAL));
    stemp = (int *) mycalloc (categs, sizeof (int));
    for (i = 0; i < categs; i++)
    {
        temp[i] = (*func) (rate[i], /*(void *)*/ helper);
        stemp[i] = (int) helper->sign;
        if (temp[i] > *rmax)
            *rmax = temp[i];
    }
    for (i = 0; i < categs; i++)
    {
        summ += ((MYREAL) stemp[i]) * probcat[i] * EXP (temp[i] - *rmax);
    }
    myfree(temp);
    myfree(stemp);
    *result = summ;
}


/* initgammacat test function*/
#ifdef LAGUERRE_TEST
int
main ()
{
    long categs = 10;
    MYREAL alpha;
    MYREAL theta1;
    MYREAL *rate;
    MYREAL *probcat;
    long i;
    long retval;
    printf ("Enter alpha, theta1,  and categs\n");
    fscanf (stdin, "%lf%lf%li", &alpha, &theta1, &categs);
    rate = (MYREAL *) mycalloc (categs + 1, sizeof (MYREAL));
    probcat = (MYREAL *) mycalloc (categs + 1, sizeof (MYREAL));
    initgammacat (categs, alpha, theta1, rate, probcat);
    printf ("Rate                   Log(prob)              Prob\n");
    for (i = 0; i < categs; i++)
    {
      printf ("%20.20f %20.20f %20.20f\n", rate[i], probcat[i], exp(probcat[i]));
    }
    return 0;
}
#endif



