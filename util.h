 /***********************************************************
 * Created: Mon 04 Feb 2013 03:35:22 PM BRST
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************/

 
#ifndef __UTIL_H
#define __UTIL_H

#include <stdio.h>
#include "permutations.h"

#define INC 1
#define DEC 0

struct information {
    int nmoves;
    long double alpha;
    long double weight;
    int nof_pr, nof_sr, nof_pt, nof_st, nof_r, nof_t;
};
typedef struct information info_t;


void* Malloc(int size);
void init_info(info_t *info);
void print_verbose(permutation_t *p, int verbose, char r, int i, int j, int k);
long double pot(long double x, long double y);
int inc_mod(int i, int n);
int dec_mod(int i, int n);
long double pot(long double x, long double y);
void exchange(int *a, int *b);
 
#endif /* __UTIL_H */

