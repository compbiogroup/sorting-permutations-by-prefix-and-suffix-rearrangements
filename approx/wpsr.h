/**********************************************************
 * Created: Seg 17 Mar 2014 02:52:06 BRT
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 **********************************************************/
 
 
#ifndef __WPRSR_H
#define __WPRSR_H
 
#include "../util.h"
#include "../permutations.h"

void wpsr_star_perm(permutation_t *p, int verbose, info_t *info);
void alg_WPSRg(permutation_t *p, int verbose, info_t *info);



long double partitionWPSR(permutation_t *p, int verbose, info_t *info);
long double WPSR(permutation_t *p, int verbose, info_t *info);
void alg_WPSR(permutation_t *p, int verbose, info_t *info);







void wspsr_star_perm(permutation_t *p, int verbose, info_t *info);
void alg_WSPSRg(permutation_t *p, int verbose, info_t *info);



long double partitionWSPSR(permutation_t *p, permutation_t *sp, int verbose, info_t *info);
long double WSPSR(permutation_t *p, permutation_t *sp, int verbose, info_t *info);
void alg_WSPSR(permutation_t *p, int verbose, info_t *info);


#endif /* __WPRSR_H */

