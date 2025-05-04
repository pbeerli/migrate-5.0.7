#ifndef HASH_INCLUDE 
#define HASH_INCLUDE
/// hash functionality
/// for migrate, currently uses string with 0 1 characters of 'any' length
/// and adding a count to the values, this is geared towards counting patterns
/// but may be easily expanded to more complex storage devices
/// currently the keys are partitioned into 1<<9 blocks and subkeys are link hashtables
/// I guess they call this a perfect hash, but this one is sparse and not all subtrees are completely
/// allocated.
///
/*
  Copyright (c) 2015 Peter Beerli

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  THE SOFTWARE.
*/
#include "migration.h"

extern boolean insert_hash(char *key, int *offset, size_t numstates, size_t realnumstates, long ptr, hash_fmt **hash, size_t *numhash);
extern long find_hash(char *key, size_t numstates, hash_fmt *hash, size_t *numhash);
extern void hash_elements(hash_fmt *hash, size_t numhash, long *total);
extern void hash_leaves(hash_fmt *hash, size_t numhash, long *leaves);
extern void show_hash(hash_fmt *hash, size_t numhash, size_t numstates, char ***showkeys, int **showvalues, long *z);
extern void free_hash(hash_fmt *hash, size_t numhash, size_t numstates);
#endif
