/**********************************************************
 * Created: Seg 17 Mar 2014 02:52:06 BRT
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 **********************************************************/
 
 
#ifndef __WPRPTSRST_H
#define __WPRPTSRST_H
 
#include "../util.h"
#include "../permutations.h"


void alg_WPSRTg(permutation_t *p, int verbose, info_t *info);
long double partitionWPSRT(permutation_t *p, int verbose, info_t *info);
long double WPSRT(permutation_t *p, int verbose, info_t *info);
void alg_WPSRT(permutation_t *p, int verbose, info_t *info);




void alg_WSPSRTg(permutation_t *p, int verbose, info_t *info);
long double partitionWSPSRT(permutation_t *p, int verbose, info_t *info);
long double WSPSRT(permutation_t *p, int verbose, info_t *info);
void alg_WSPSRT(permutation_t *p, int verbose, info_t *info);

#endif /* __WPRPTSRST_H */

