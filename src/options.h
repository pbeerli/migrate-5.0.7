/*! \file options.h */
#ifndef OPTIONS_INCLUDE
#define OPTIONS_INCLUDE
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 O P T I O N S   R O U T I N E S 
 
 creates options structures,
 reads options from parmfile if present
 set first parameter for mcmc run,
 prints options,
 and finally helps to destroy itself.
                                                                                                               
 Peter Beerli 1996, Seattle
    updated   2009, Tallahassee
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

$Id: options.h 2157 2013-04-16 22:27:40Z beerli $
-------------------------------------------------------*/

#include "migration.h"
//
// initializes the memory for the option structure
// \param ** options {pointer to the optionstructure pointer}
extern void create_options (option_fmt ** options);
//
// populates the content of the option structure with default values 
// \param options  {is the option structure that holds all general option data}
// \returns None
extern void init_options (option_fmt * options);
extern void destroy_options (option_fmt * options);

//
// reads the options from the parmfile, the parmfile name default is "parmfile"
// but this can be overidden using the arguments to the main(). The parmfile is already
// open and has a pointer to options->parmfile.
// \param options  {is the option structure that holds all general option data}
// \returns None
extern void read_options_master (option_fmt * options);
//
// reads the options from the buffer that was filled by the master in save_options_buffer().
// \param buffer { pointer to the buffer string that contains all parameter values 
//                 this buffer is filled by the function save_options_buffer() } 
// \param options  {is the option structure that holds all general option data}
// \returns None
extern void read_options_worker (char **buffer, option_fmt * options);
//
// reads the custom migration matrix from the parmfile
// \param file { pointer to the parmfile}  
// \param options  {is the option structure that holds all general option data}
// \param value {value string holding }
// \param customnumpop {number of populations specified in the custom-migration string read from the parmfile}
// \returns None
extern void read_custom_migration (FILE * file, option_fmt * options, char *value, long customnumpop,long pos);

extern long scan_connect (char *custm2, long start, long stop, int check);
extern void fillup_custm (long len, world_fmt * world, option_fmt * options);

extern void set_param (world_fmt * world, data_fmt * data, option_fmt * options, long locus);
extern void set_profile_options (option_fmt * options);
extern void decide_plot (worldoption_fmt * options, long chain, long chains, char type);
extern void set_plot (option_fmt * options);

extern long save_options_buffer (char **buffer, long *allocbufsize, option_fmt * options, data_fmt *data);
extern long save_mu_rates_buffer (char **buffer, long *allocbufsize, option_fmt * options);
extern long save_parmfile (option_fmt * options, world_fmt * world, data_fmt *data);
extern void print_menu_options (world_fmt * world, option_fmt * options, data_fmt * data);
extern void print_options (FILE * file, world_fmt * world, option_fmt * options, data_fmt * data);

//
// prints the value of minimum value of the prior distribution
// \param tmp {string to hold the minimum value of the prior distribution}
// \param prior {pointer to the prior structure}
// \param priortype {indicator of the priortype}
// \returns {pointer to tmp}
extern char * show_priormin(char *tmp, prior_fmt *prior);
//
// prints the value of the mean of the prior distribution
// \param tmp {string to hold the mean value of the prior distribution}
// \param prior {pointer to the prior structure}/// \param priortype {indicator of the priortype}
// \returns {pointer to tmp}
extern char * show_priormean(char *tmp, prior_fmt *prior);
//
// prints the value of maximum value of the prior distribution
// \param tmp {string to hold the maximum value of the prior distribution}
// \param prior {pointer to the prior structure}
// \param priortype {indicator of the priortype}
// \returns {pointer to tmp}
extern char * show_priormax(char *tmp, prior_fmt *prior);
//
// prints the value of delta used for the prior distribution
// \param tmp {string to hold the delta value used for the prior distribution}
// \param prior {pointer to the prior structure}
// \param priortype {indicator of the priortype}
// \returns {pointer to tmp}
extern char * show_priordelta(char *tmp, prior_fmt *prior);
//
// prints the value of the number of bins of the posterior and prior distribution
// \param tmp {string to hold the number of bins for recording of the posterior}
// \param prior {pointer to the prior structure}
// \param priortype {indicator of the priortype}
// \returns {pointer to tmp}
extern char * show_priorbins(char *tmp, prior_fmt *prior);

extern char * show_priorupdatefreq(char *tmp, prior_fmt *prior);
//
// adjust the parameter values according to the custom migration matrix
// to set the parameters to specific values
// that allow for parameter reduction by using only one average population size
// or one average migration rate, or symmetric migration rates [symmetric m/mu
// or symmetric Nm values] or constant parameter values that do not change during the run.
// \param world {holds all runtime related material including the data and run-options structures}
// \param options  {is the option structure that holds all general option data}
// \returns {None}
extern void synchronize_param (world_fmt * world, option_fmt * options);
//
// resynchronizes the parameters with the custm migraiton matrix but only recalculates the material that was 
// changed, is supposed to do less work than the synchronize_param() function.
// \param world {holds all runtime related material including the data and run-options structures}
// \returns {None}
extern void resynchronize_param (world_fmt * world);
//
// Reorders and relabels the locations(populations), so that the datafile may contain n number of locations that
// can be combined into m number of new locations that are then used for the anlysis.
// Examples:\n
// \verbatim
// (1,2,3,4,5) --> (1,1,1,2,2)
// (1,2,3,4,5) --> (2,3,1,1,1)
// \endverbatim
// The reordering parameters are held in options->newposts and its size is data->numpop 
// \param world {holds all runtime related material including the data and run-options structures}
// \param options  {is the option structure that holds all general option data}
// \param data  {holds all data related parts}
// \returns {None}
extern void reorder_populations(world_fmt *world, option_fmt *options, data_fmt *data);
// 
// fills the option->newposts array with localities to reorder 
extern void set_localities(char **value, char **tmp, option_fmt *options);
//
// sets the average of the mutation rates among loci and updates option structures
// in world->options and options.
// \param wopt  {is the lean local copy of options in world}
// \param options  {is the option structure that holds all general option data}
// \param loci {number of loci}
// \returns {None}
extern void set_meanmu(worldoption_fmt * wopt, option_fmt * options, long loci);
//
// sets the vector of cumumulative probabilities to choose one of the avaiable
// updates: tree, parameter, haplotype in world->options.
// \param choices  {vector of comumulative probabilities
// \param options  {is the option structure that holds all general option data}
// \param flag {holds particular option: STANDARD allows input for all updates, NOPARAMETER
//        sets parameter update to zero, NOTREE sets tree updates to zero.}
// \returns {None}
extern void set_updating_choices(double *choices, option_fmt * options, int flag);
//
// initialize or reinitialize the start parameters, so that all parameter values can have a setting
// that either comes from the PRIOR (using a value at a percentage of PDF or a random value
// from the prior distribution.
extern void allocate_startparam(option_fmt *options, long numpop);
extern void realloc_startparam(option_fmt *options, long numpop);
//
// fill the startparam, this is called in the menuParam() function
extern void fill_startparam(option_fmt *options, short key, long index, float value);
extern void fill_printvar_startparam(option_fmt *options, char *paramtgen[], char *parammgen[]);
// delivers the keys values for the startparameters
extern long startguessvalue(option_fmt *options, short key);
extern void show_allstartparam(FILE *file, option_fmt *options);
// shows values for each of the keys
extern long show_thetaownparam(FILE *file, option_fmt *options, char **temp);
extern long show_migownparam(FILE *file, option_fmt *options, char **temp);
extern long show_rateownparam(FILE *file, option_fmt *options, char **temp);
extern long show_splitownparam(FILE *file, option_fmt *options, char **temp);
extern long show_splitstdownparam(FILE *file, option_fmt *options, char **temp);
// deliver the correct type of startparam for the options in the outfile
extern void show_startparamtype(twin_fmt *guess, long index, char **temp);
// show the parameter start values for each key
extern void show_startparam(FILE *file, option_fmt *options, short key, boolean verbose);
#endif /*OPTIONS_INCLUDE */
