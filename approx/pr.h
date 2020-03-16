/***********************************************************
 * Created: Mon 04 Feb 2013 03:35:22 PM BRST
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************/
 
 
#ifndef __PREFIX_REVERSAL_H
#define __PREFIX_REVERSAL_H

#include "../util.h"
#include "../permutations.h"

void pr_star_perm(permutation_t *p, int verbose, info_t *info);
int pr_pr_edge_type1(int *edge, permutation_t *p);
int pr_pr_edge_type2(int *edge, permutation_t *p);
int pr_pr_edge_type3A(int *edge, permutation_t *p);
int pr_pr_edge_type3B(int *edge, permutation_t *p);
void alg_2PR(permutation_t *p, int verbose, info_t *info);


void spr_star_perm(permutation_t *p, int verbose, info_t *info);
void alg_2SPR(permutation_t *p, int verbose, info_t *info);
 
#endif 

