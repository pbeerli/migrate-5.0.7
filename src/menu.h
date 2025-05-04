#ifndef MENU_INCLUDE
#define MENU_INCLUDE
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 M E N U   R O U T I N E S 
 
 presents the menu and its submenus.                                                                                                               
 Peter Beerli 1996, Seattle
      updated 2009, Tallahassee
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2009 Peter Beerli, Tallahassee FL
 
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

$Id: menu.h 2157 2013-04-16 22:27:40Z beerli $
-------------------------------------------------------*/

#include "migration.h"
///
/// print the title of the analysis to the log and the screen and the text outputfile
long  print_menu_title (FILE * file, option_fmt * options, world_fmt *world);
///
/// print the acceptance ratio to the screen
/// \param a {likelihood of the newly proposed MCMC state}
/// \param b {likelihood of the old MCMC state
/// \param world {holds all runtime related material including the data and run-options structures}
void print_menu_accratio (long a, long b, world_fmt *world);
///
/// print the title to the outfile
/// \param world {holds all runtime related material including the data and run-options structures}
/// \param options  {is the option structure that holds all general option data}
/// \returns {the position in the output file}
long print_title (world_fmt * world, option_fmt * options);
///
/// presents the menu to the user
/// \param options {is the option structure that holds all general option data}
/// \param world {holds all runtime related material including the data and run-options structures, 
///               but this is not really filled yet because all data related material is missing}
/// \param data  {holds all data related parts}
void get_menu (option_fmt * options, world_fmt *world, data_fmt *data);
///
/// presents title of the prior method
/// \returns {a string with the name of the proposal method}
extern long  is_priortype(prior_fmt *p, long pnum, int priortype);
///
/// presents title of the posterior method
/// \param proposalset {is a boolean to specify whether this is SPLICE or MH sampling scheme}
/// \returns {a string with the name of the proposal method}
char * is_proposaltype(boolean proposalset);
#endif /*MENU_INCLUDE */
