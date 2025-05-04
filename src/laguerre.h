#ifndef LAGUERRE_H_
#define LAGUERRE_H_
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 Laguerre integration  R O U T I N E S 
 
 Peter Beerli 2000, Seattle
 beerli@fsu.edu
 
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

$Id: laguerre.h 2157 2013-04-16 22:27:40Z beerli $
-------------------------------------------------------*/

extern void integrate_laguerre (long categs, MYREAL *rate,
                                    MYREAL *probcat,
                                    MYREAL (*func) (MYREAL, helper_fmt *),
                                    helper_fmt * helper, MYREAL *result,
                                    MYREAL *rmax);
extern void initgammacat (long categs, MYREAL alpha, MYREAL theta1,
                              MYREAL *rate, MYREAL *probcat);


#endif
