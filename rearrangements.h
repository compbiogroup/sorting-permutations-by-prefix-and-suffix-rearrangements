/**********************************************************
 * Created: Sáb 19 Jan 2013 11:39:18 BRST
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************/
 
 
#ifndef __REARRANGEMENTS_H
#define __REARRANGEMENTS_H

#include "permutations.h"

void op_reversal(permutation_t *p, int i, int j);
void op_prefix_reversal(permutation_t *p, int j);
void op_suffix_reversal(permutation_t *p, int j);
void op_transposition(permutation_t *p, int i, int j, int k);
void op_prefix_transposition(permutation_t *p, int j, int k);
void op_suffix_transposition(permutation_t *p, int i, int j);

double lower_bound_pr(permutation_t *p);
double lower_bound_psr(permutation_t *p);
double lower_bound_pt(permutation_t *p);
double lower_bound_pst(permutation_t *p);
double lower_bound_prt(permutation_t *p);
double lower_bound_psrt(permutation_t *p);

double lower_bound_spr(permutation_t *p);
double lower_bound_spsr(permutation_t *p);
double lower_bound_sprt(permutation_t *p);
double lower_bound_spsrt(permutation_t *p);

double lower_bound_wpr(permutation_t *p);
double lower_bound_wpt(permutation_t *p);
double lower_bound_wprt(permutation_t *p);
/* It only works for alpha = 1*/
double lower_bound_wpsr(permutation_t *p);
/* It only works for alpha = 1*/
double lower_bound_wpst(permutation_t *p);
/* It only works for alpha = 1*/
double lower_bound_wpsrt(permutation_t *p);

double lower_bound_wspr(permutation_t *p);
double lower_bound_wsprt(permutation_t *p);
/* It only works for alpha = 1*/
double lower_bound_wspsr(permutation_t *p);
/* It only works for alpha = 1*/
double lower_bound_wspsrt(permutation_t *p);

double lower_bound(char prob_name[], permutation_t *p);


int prt_prohibited(permutation_t *p, char type, int i, int j);
int sprt_prohibited(permutation_t *p, char type, int i, int j);


#endif /* __REARRANGEMENTS_H */

