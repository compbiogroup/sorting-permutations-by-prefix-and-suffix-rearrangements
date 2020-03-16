/**********************************************************
 * Created: Seg 17 Mar 2014 02:52:06 BRT
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 **********************************************************/
 
 
#ifndef __WPRPT_H
#define __WPRPT_H
 
#include "../util.h"
#include "../permutations.h"

void alg_WPRTg(permutation_t *p, int verbose, info_t *info);
long double partitionWPRT(permutation_t *p, int end, int type, int divisor, int verbose, info_t *info);
long double WPRT(permutation_t *p, int end, int type, int verbose, info_t *info);
void alg_WPRT(permutation_t *p, int verbose, info_t *info);



void alg_WSPRTg(permutation_t *p, int verbose, info_t *info);
long double partitionWSPRT(permutation_t *p, int end, int type, int divisor, int verbose, info_t *info);
long double WSPRT(permutation_t *p, int end, int type, int verbose, info_t *info);
void alg_WSPRT(permutation_t *p, int verbose, info_t *info);

#endif /* __WPRPT_H */

