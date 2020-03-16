/**********************************************************
 * Created: Seg 17 Mar 2014 02:52:06 BRT
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 **********************************************************/
 
 
#ifndef __WPR_H
#define __WPR_H
 
#include "../util.h"
#include "../permutations.h"

void alg_WPRg(permutation_t *p, int verbose, info_t *info);
void alg_WPRm(permutation_t *p, int verbose, info_t *info);

long double partitionWPR(permutation_t *p, int end, int type, int divisor, int verbose, info_t *info);
long double WPR(permutation_t *p, int end, int type, int verbose, info_t *info);
void alg_WPR(permutation_t *p, int verbose, info_t *info);



void alg_WSPRg(permutation_t *p, int verbose, info_t *info);

int partitionWSPR(permutation_t *p, permutation_t *sp, int end, int type, int divisor, int verbose, info_t *info);
int WSPR(permutation_t *p, permutation_t *sp, int end, int type, int verbose, info_t *info);
void alg_WSPR(permutation_t *p, int verbose, info_t *info);

#endif /* __WPR_H */

