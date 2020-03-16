/***********************************************************
 * Created: Mon 04 Feb 2013 03:35:22 PM BRST
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************/
 
 
#ifndef __WEIGHTED_PREFIX_TRANSPOSITION_H
#define __WEIGHTED_PREFIX_TRANSPOSITION_H

#include "../util.h"
#include "../permutations.h"

void alg_WPTg(permutation_t *p, int verbose, info_t* info);


long double WPT(permutation_t *p, int end, int verbose, info_t *info);
long double partitionWPT(permutation_t *p, int fim, int divisor, int verbose, info_t *info);
void alg_WPT(permutation_t *p, int verbose, info_t* info);
 
#endif /* __WEIGHTED_PREFIX_TRANSPOSITION_H */

