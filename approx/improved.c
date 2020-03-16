/***********************************************************
 * Created: Sáb 14 Dez 2013 17:58:55 BRST
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************/
 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../algorithms.h"
#include "../permutations.h"
#include "../rearrangements.h"
#include "../breakpoints.h"
#include "../util.h"
#include "pr.h"
#include "psr.h"



static void fill(int *vet[5], int i, int v0, int v1, int v2, int v3, int v4) {/*{{{*/
    vet[0][i] = v0;
    vet[1][i] = v1;
    vet[2][i] = v2;
    vet[3][i] = v3;
    vet[4][i] = v4;
}/*}}}*/


static void pr_1bp_1pr(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i;

    i = p->inv_pi[p->pi[1] + 1] - 1;
    if (2 <= i && i <= p->size && is_UPR_breakpoint(p, i, i+1)) {
        (*size)++;
        fill(score, *size, i, 0, 0, 'p','r');
    }

    i = p->inv_pi[p->pi[1] - 1] - 1;
    if (2 <= i && i <= p->size && is_UPR_breakpoint(p, i, i+1)) {
        (*size)++;
        fill(score, *size, i, 0, 0, 'p','r');
    }
}/*}}}*/

static void pr_1bp_2pr(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    for (i = 1; i <= p->size; i++) {
        if (!is_UPR_breakpoint(p, i, i+1))
            continue;

        j = p->inv_pi[p->pi[i] + 1];
        if (i < j && j <= p->size && is_UPR_breakpoint(p, j, j+1)) {
            (*size)++;
            fill(score, *size, j, j-i, 0, 'p', 'r');
        }
        if (1 < i && i < j && is_UPR_breakpoint(p, j-1, j)) {
            (*size)++;
            fill(score, *size, i, j-1, 0, 'p', 'r');
        }

        j = p->inv_pi[p->pi[i] - 1];
        if (i < j && j <= p->size && is_UPR_breakpoint(p, j, j+1)) {
            (*size)++;
            fill(score, *size, j, j-i, 0, 'p', 'r');
        }
        if (1 < i && i < j && is_UPR_breakpoint(p, j-1, j)) {
            (*size)++;
            fill(score, *size, i, j-1, 0, 'p', 'r');
        }
    }
}/*}}}*/


static void psr_1bp_1pr(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i;

    i = p->inv_pi[p->pi[1] + 1] - 1;
    if (2 <= i && i < p->size && is_UPSR_breakpoint(p, i, i+1)) {
        (*size)++;
        fill(score, *size, i, 0, 0, 'p','r');
    }

    i = p->inv_pi[p->pi[1] - 1] - 1;
    if (2 <= i && i < p->size && is_UPSR_breakpoint(p, i, i+1)) {
        (*size)++;
        fill(score, *size, i, 0, 0, 'p','r');
    }
}/*}}}*/

static void psr_1bp_2pr(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    for (i = 1; i <= p->size-1; i++) {
        if (!is_UPSR_breakpoint(p, i, i+1))
            continue;

        j = p->inv_pi[p->pi[i] + 1];
        if (i < j && j <= p->size && is_UPSR_breakpoint(p, j, j+1)) {
            (*size)++;
            fill(score, *size, j, j-i, 0, 'p', 'r');
        }
        if (1 < i && i < j && j <= p->size && is_UPSR_breakpoint(p, j-1, j)) {
            (*size)++;
            fill(score, *size, i, j-1, 0, 'p', 'r');
        }

        j = p->inv_pi[p->pi[i] - 1];
        if (i < j && j <= p->size && is_UPSR_breakpoint(p, j, j+1)) {
            (*size)++;
            fill(score, *size, j, j-i, 0, 'p', 'r');
        }
        if (1 < i && i < j && j <= p->size && is_UPSR_breakpoint(p, j-1, j)) {
            (*size)++;
            fill(score, *size, i, j-1, 0, 'p', 'r');
        }
    }
}/*}}}*/

static void psr_1bp_1sr(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i;

    i = p->inv_pi[p->pi[p->size] + 1] + 1;
    if (2 <= i && i <= p->size-1 && is_UPSR_breakpoint(p, i-1, i)) {
        (*size)++;
        fill(score, *size, i, 0, 0, 's','r');
    }

    i = p->inv_pi[p->pi[p->size] - 1] + 1;
    if (2 <= i && i <= p->size-1 && is_UPSR_breakpoint(p, i-1, i)) {
        (*size)++;
        fill(score, *size, i, 0, 0, 's','r');
    }
}/*}}}*/

static void psr_1bp_2sr(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    for (j = p->size; j >= 2; j--) {
        if (!is_UPSR_breakpoint(p, j-1, j))
            continue;

        i = p->inv_pi[p->pi[j] + 1];
        if (j > i && i > 0 && is_UPSR_breakpoint(p, i-1, i)) {
            (*size)++;
            fill(score, *size, i, p->size+1-(j-i), 0, 's','r');
        }
        if (j > i && j < p->size && i > 0 && is_UPSR_breakpoint(p, i, i+1)) {
            (*size)++;
            fill(score, *size, j, i+1, 0, 's','r');
        }

        i = p->inv_pi[p->pi[j] - 1];
        if (j > i && i > 0 && is_UPSR_breakpoint(p, i-1, i)) {
            (*size)++;
            fill(score, *size, i, p->size+1-(j-i), 0, 's','r');
        }
        if (j > i && j < p->size && i > 0 && is_UPSR_breakpoint(p, i, i+1)) {
            (*size)++;
            fill(score, *size, j, i+1, 0, 's','r');
        }
    }
}/*}}}*/


static void pt_2bp_1pt(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    j = p->inv_pi[p->pi[1] - 1] + 1;
    i = p->inv_pi[p->pi[j] - 1] + 1;
    if (2 <= i && i < j) {
        (*size)++;
        fill(score, *size, i, j, 0, 'p', 't');
    }
}/*}}}*/

static void pt_1bp_1pt(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    for (i = 2; i <= p->size; i++) {
        if (!is_P_breakpoint(p, i-1, i))
            continue;

        j = p->inv_pi[p->pi[i-1] + 1];
        if (i < j) {
            (*size)++;
            fill(score, *size, i, j, 0, 'p', 't');
        }
    }

    j = p->inv_pi[p->pi[1] - 1] + 1;
    for (i = 2; i < j; i++) {
        if (is_P_breakpoint(p, i-1, i)) {
            (*size)++;
            fill(score, *size, i, j, 0, 'p', 't');
        }
    }
}/*}}}*/


static void pst_2bp_1pt(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    j = p->inv_pi[p->pi[1] - 1] + 1;
    i = p->inv_pi[p->pi[j] - 1] + 1;
    if (2 <= i && i < j && j <= p->size) {
        (*size)++;
        fill(score, *size, i, j, 0, 'p', 't');
    }
}/*}}}*/

static void pst_1bp_1pt(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    for (i = 2; i <= p->size-1; i++) {
        if (!is_PS_breakpoint(p, i-1, i))
            continue;

        j = p->inv_pi[p->pi[i-1] + 1];
        if (i < j && j <= p->size) {
            (*size)++;
            fill(score, *size, i, j, 0, 'p', 't');
        }
    }

    j = p->inv_pi[p->pi[1] - 1] + 1;
    for (i = 2; i < j; i++) {
        if (is_PS_breakpoint(p, i-1, i)) {
            (*size)++;
            fill(score, *size, i, j, 0, 'p', 't');
        }
    }
}/*}}}*/

static void pst_2bp_1st(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    i = p->inv_pi[p->pi[p->size] + 1];
    j = p->inv_pi[p->pi[i - 1] + 1];
    if (2 <= i && i < j && j <= p->size) {
        (*size)++;
        fill(score, *size, i, j, 0, 's', 't');
    }
}/*}}}*/

static void pst_1bp_1st(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    for (i = 2; i <= p->size-1; i++) {
        if (!is_PS_breakpoint(p, i-1, i))
            continue;

        j = p->inv_pi[p->pi[i-1] + 1];
        if (i < j && j <= p->size) {
            (*size)++;
            fill(score, *size, i, j, 0, 's', 't');
        }
    }

    i = p->inv_pi[p->pi[p->size] + 1];
    for (j = i+1; j <= p->size; j++) {
        if (is_PS_breakpoint(p, j-1, j)) {
            (*size)++;
            fill(score, *size, i, j, 0, 's', 't');
        }
    }
}/*}}}*/


static void prt_1bp_1pr(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i;

    i = p->inv_pi[p->pi[1] + 1] - 1;
    if (2 <= i && i <= p->size && is_UPR_breakpoint(p, i, i+1) &&
                    !prt_prohibited(p, 'r', i, 0)) {
        (*size)++;
        fill(score, *size, i, 0, 0, 'p','r');
    }

    i = p->inv_pi[p->pi[1] - 1] - 1;
    if (2 <= i && i <= p->size && is_UPR_breakpoint(p, i, i+1) &&
                    !prt_prohibited(p, 'r', i, 0)) {
        (*size)++;
        fill(score, *size, i, 0, 0, 'p','r');
    }
}/*}}}*/

static void prt_1bp_1pt(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    for (i = 2; i <= p->size; i++) {
        if (!is_UPR_breakpoint(p, i-1, i))
            continue;

        j = p->inv_pi[p->pi[i-1] + 1];
        if (j > i && j <= p->size+1 && is_UPR_breakpoint(p, j-1, j) &&
                !prt_prohibited(p, 't', i, j)) {
            (*size)++;
            fill(score, *size, i, j, 0, 'p', 't');
        }

        j = p->inv_pi[p->pi[i-1] - 1];
        if (j > i && j <= p->size+1 && is_UPR_breakpoint(p, j-1, j) &&
                !prt_prohibited(p, 't', i, j)) {
            (*size)++;
            fill(score, *size, i, j, 0, 'p', 't');
        }
    }

    j = p->inv_pi[p->pi[1] + 1] + 1;
    if (j > 2 && j <= p->size+1 && is_UPR_breakpoint(p, j-1, j)) {
        for (i = 2; i < j; i++) {
            if (is_UPR_breakpoint(p, i-1, i) && !prt_prohibited(p, 't', i, j)) {
                (*size)++;
                fill(score, *size, i, j, 0, 'p', 't');
            }
        }
    }

    j = p->inv_pi[p->pi[1] - 1] + 1;
    if (j > 2 && j <= p->size+1 && is_UPR_breakpoint(p, j-1, j)) {
        for (i = 2; i < j; i++) {
            if (is_UPR_breakpoint(p, i-1, i) && !prt_prohibited(p, 't', i, j)) {
                (*size)++;
                fill(score, *size, i, j, 0, 'p', 't');
            }
        }
    }
}/*}}}*/

static void prt_2bp_1pt(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    j = p->inv_pi[p->pi[1] + 1] + 1;
    if (j >= 0 && j < p->size+1) {
        i = p->inv_pi[p->pi[j] + 1] + 1;
        if (i > 1 && i < j && is_UPR_breakpoint(p, i-1, i) 
                && is_UPR_breakpoint(p, j-1, j) && p->pi[i] != 1) {
            (*size)++;
            fill(score, *size, i, j, 0, 'p', 't');
        }
    }
    if (j > 0 && j <= p->size+1) {
        i = p->inv_pi[p->pi[j] - 1] + 1;
        if (i > 1 && i < j && is_UPR_breakpoint(p, i-1, i) 
                && is_UPR_breakpoint(p, j-1, j) && p->pi[i] != 1) {
            (*size)++;
            fill(score, *size, i, j, 0, 'p', 't');
        }
    }

    j = p->inv_pi[p->pi[1] - 1] + 1;
    if (j >= 0 && j < p->size+1) {
        i = p->inv_pi[p->pi[j] + 1] + 1;
        if (i > 1 && i < j && is_UPR_breakpoint(p, i-1, i) 
                && is_UPR_breakpoint(p, j-1, j) && p->pi[i] != 1) {
            (*size)++;
            fill(score, *size, i, j, 0, 'p', 't');
        }
    }
    if (j > 0 && j <= p->size+1) {
        i = p->inv_pi[p->pi[j] - 1] + 1;
        if (i > 1 && i < j && is_UPR_breakpoint(p, i-1, i) 
                && is_UPR_breakpoint(p, j-1, j) && p->pi[i] != 1) {
            (*size)++;
            fill(score, *size, i, j, 0, 'p', 't');
        }
    }
}/*}}}*/


static void psrt_1bp_1pr(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i;

    i = p->inv_pi[inc_mod(p->pi[1], p->size)] - 1;
    if (2 <= i && i < p->size && is_UPSRT_breakpoint(p, i, i+1)) {
        (*size)++;
        fill(score, *size, i, 0, 0, 'p','r');
    }

    i = p->inv_pi[dec_mod(p->pi[1], p->size)] - 1;
    if (2 <= i && i < p->size && is_UPSRT_breakpoint(p, i, i+1)) {
        (*size)++;
        fill(score, *size, i, 0, 0, 'p','r');
    }
}/*}}}*/

static void psrt_2bp_1pt(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    j = p->inv_pi[dec_mod(p->pi[1], p->size)] + 1;
    i = p->inv_pi[dec_mod(p->pi[j], p->size)] + 1;
    if (2 <= i && i < j && j <= p->size &&
            is_UPSRT_breakpoint(p, i-1, i) && 
            is_UPSRT_breakpoint(p, j-1, j)) {
        (*size)++;
        fill(score, *size, i, j, 0, 'p', 't');
    }

    j = p->inv_pi[dec_mod(p->pi[1], p->size)] + 1;
    i = p->inv_pi[inc_mod(p->pi[j], p->size)] + 1;
    if (2 <= i && i < j && j <= p->size &&
            is_UPSRT_breakpoint(p, i-1, i) && 
            is_UPSRT_breakpoint(p, j-1, j)) {
        (*size)++;
        fill(score, *size, i, j, 0, 'p', 't');
    }

    j = p->inv_pi[inc_mod(p->pi[1], p->size)] + 1;
    i = p->inv_pi[dec_mod(p->pi[j], p->size)] + 1;
    if (2 <= i && i < j && j <= p->size &&
            is_UPSRT_breakpoint(p, i-1, i) &&
            is_UPSRT_breakpoint(p, j-1, j)) {
        (*size)++;
        fill(score, *size, i, j, 0, 'p', 't');
    }

    j = p->inv_pi[inc_mod(p->pi[1], p->size)] + 1;
    i = p->inv_pi[inc_mod(p->pi[j], p->size)] + 1;
    if (2 <= i && i < j && j <= p->size &&
            is_UPSRT_breakpoint(p, i-1, i) &&
            is_UPSRT_breakpoint(p, j-1, j)) {
        (*size)++;
        fill(score, *size, i, j, 0, 'p', 't');
    }
}/*}}}*/

static void psrt_1bp_1pt(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    for (i = 2; i <= p->size-1; i++) {
        if (!is_UPSRT_breakpoint(p, i-1, i))
            continue;

        j = p->inv_pi[inc_mod(p->pi[i-1], p->size)];
        if (i < j && j <= p->size && is_UPSRT_breakpoint(p, j-1, j)) {
            (*size)++;
            fill(score, *size, i, j, 0, 'p', 't');
        }

        j = p->inv_pi[dec_mod(p->pi[i-1], p->size)];
        if (i < j && j <= p->size && is_UPSRT_breakpoint(p, j-1, j)) {
            (*size)++;
            fill(score, *size, i, j, 0, 'p', 't');
        }
    }

    j = p->inv_pi[inc_mod(p->pi[1], p->size)] + 1;
    if (j == p->size+1 || (j > 2 && j < p->size+1 && is_UPSRT_breakpoint(p, j-1, j))) {
        for (i = 2; i < j; i++) {
            if (is_UPSRT_breakpoint(p, i-1, i)) {
                (*size)++;
                fill(score, *size, i, j, 0, 'p', 't');
            }
        }
    }

    j = p->inv_pi[dec_mod(p->pi[1], p->size)] + 1;
    if (j == p->size+1 || (j > 2 && j < p->size+1 && is_UPSRT_breakpoint(p, j-1, j))) {
        for (i = 2; i < j; i++) {
            if (is_UPSRT_breakpoint(p, i-1, i)) {
                (*size)++;
                fill(score, *size, i, j, 0, 'p', 't');
            }
        }
    }
}/*}}}*/

static void psrt_1bp_1sr(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i;

    i = p->inv_pi[inc_mod(p->pi[p->size], p->size)] + 1;
    if (2 <= i && i < p->size && is_UPSRT_breakpoint(p, i-1, i)) {
        (*size)++;
        fill(score, *size, i, 0, 0, 's','r');
    }

    i = p->inv_pi[dec_mod(p->pi[p->size], p->size)] + 1;
    if (2 <= i && i < p->size && is_UPSRT_breakpoint(p, i-1, i)) {
        (*size)++;
        fill(score, *size, i, 0, 0, 's','r');
    }
}/*}}}*/

static void psrt_2bp_1st(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    i = p->inv_pi[dec_mod(p->pi[p->size], p->size)];
    j = p->inv_pi[dec_mod(p->pi[i-1], p->size)];
    if (2 <= i && i < j && j <= p->size &&
           is_UPSRT_breakpoint(p, i-1, i) && 
           is_UPSRT_breakpoint(p, j-1, j)) {
        (*size)++;
        fill(score, *size, i, j, 0, 's', 't');
    }

    i = p->inv_pi[dec_mod(p->pi[p->size], p->size)];
    j = p->inv_pi[inc_mod(p->pi[i-1], p->size)];
    if (2 <= i && i < j && j <= p->size &&
           is_UPSRT_breakpoint(p, i-1, i) && 
           is_UPSRT_breakpoint(p, j-1, j)) {
        (*size)++;
        fill(score, *size, i, j, 0, 's', 't');
    }

    i = p->inv_pi[inc_mod(p->pi[p->size], p->size)];
    j = p->inv_pi[dec_mod(p->pi[i-1], p->size)];
    if (2 <= i && i < j && j <= p->size &&
           is_UPSRT_breakpoint(p, i-1, i) && 
           is_UPSRT_breakpoint(p, j-1, j)) {
        (*size)++;
        fill(score, *size, i, j, 0, 's', 't');
    }

    i = p->inv_pi[inc_mod(p->pi[p->size], p->size)];
    j = p->inv_pi[inc_mod(p->pi[i-1], p->size)];
    if (2 <= i && i < j && j <= p->size &&
           is_UPSRT_breakpoint(p, i-1, i) && 
           is_UPSRT_breakpoint(p, j-1, j)) {
        (*size)++;
        fill(score, *size, i, j, 0, 's', 't');
    }
}/*}}}*/

static void psrt_1bp_1st(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    for (i = 2; i <= p->size-1; i++) {
        if (!is_UPSRT_breakpoint(p, i-1, i))
            continue;

        j = p->inv_pi[inc_mod(p->pi[i-1], p->size)];
        if (i < j && j <= p->size && is_UPSRT_breakpoint(p, j-1, j)) {
            (*size)++;
            fill(score, *size, i, j, 0, 's', 't');
        }

        j = p->inv_pi[dec_mod(p->pi[i-1], p->size)];
        if (i < j && j <= p->size && is_UPSRT_breakpoint(p, j-1, j)) {
            (*size)++;
            fill(score, *size, i, j, 0, 's', 't');
        }
    }

    i = p->inv_pi[inc_mod(p->pi[p->size], p->size)];
    if (i == 1 || (i > 1 && i <= p->size-1 && is_UPSRT_breakpoint(p, i-1, i))) {
        for (j = i+1; j <= p->size; j++) {
            if (is_UPSRT_breakpoint(p, j-1, j)) {
                (*size)++;
                fill(score, *size, i, j, 0, 's', 't');
            }
        }
    }

    i = p->inv_pi[p->pi[p->size] - 1];
    if (i == 1 || (i > 0 && i <= p->size-1 && is_UPSRT_breakpoint(p, i-1, i))) {
        for (j = i+1; j <= p->size; j++) {
            if (is_UPSRT_breakpoint(p, j-1, j)) {
                (*size)++;
                fill(score, *size, i, j, 0, 's', 't');
            }
        }
    }
}/*}}}*/


static void spr_1bp_1pr(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i;

    i = prm_exists(p, -p->pi[1]+1);
    if (i >= 2) {
        (*size)++;
        fill(score, *size, i-1, 0, 0, 'p','r');
    }
}/*}}}*/

static void spr_1bp_2pr(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    for (i = 1; i <= p->size; i++) {
        if (!is_P_breakpoint(p, i, i+1))
            continue;

        if (i < p->size) {
            j = prm_exists(p, -p->pi[i]-1);
            if (i < j && j <= p->size) {
                (*size)++;
                fill(score, *size, j, j-i, 0, 'p', 'r');
            }
        }

        j = prm_exists(p, p->pi[i]+1);
        if (i+1 < j) {
            (*size)++;
            fill(score, *size, i, j-1, 0, 'p', 'r');
        }
    }
}/*}}}*/


static void spsr_1bp_1pr(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i;

    i = prm_exists(p, -p->pi[1]+1);
    if (i >= 2 && i <= p->size) {
        (*size)++;
        fill(score, *size, i-1, 0, 0, 'p','r');
    }
}/*}}}*/

static void spsr_1bp_2pr(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    for (i = 1; i <= p->size; i++) {
        if (!is_PS_breakpoint(p, i, i+1))
            continue;

        j = prm_exists(p, -p->pi[i]-1);
        if (j > 0 && i < j && j <= p->size) {
            (*size)++;
            fill(score, *size, j, j-i, 0, 'p', 'r');
        }

        j = prm_exists(p, p->pi[i]+1);
        if (j > 0 && i+1 < j && j <= p->size) {
            (*size)++;
            fill(score, *size, i, j-1, 0, 'p', 'r');
        }
    }
}/*}}}*/

static void spsr_1bp_1sr(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i;

    i = prm_exists(p, -p->pi[p->size]-1);
    if (i >= 1) {
        (*size)++;
        fill(score, *size, i+1, 0, 0, 's','r');
    }
}/*}}}*/

static void spsr_1bp_2sr(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    for (j = p->size; j >= 3; j--) {
        if (!is_PS_breakpoint(p, j-1, j))
            continue;

        i = prm_exists(p, -p->pi[j]+1);
        if (i >= 1 && i < j) {
            (*size)++;
            fill(score, *size, i, p->size+1-(j-i), 0, 's','r');
        }

        i = prm_exists(p, p->pi[j]-1);
        if (i >= 1 && i+1 < j) {
            (*size)++;
            fill(score, *size, j, i+1, 0, 's','r');
        }
    }
}/*}}}*/


static void sprt_1bp_1pr(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i;

    i = prm_exists(p, -p->pi[1]+1);
    if (i > 0 && is_P_breakpoint(p, i-1, i) && !sprt_prohibited(p, 'r', i-1, 0)) {
        (*size)++;
        fill(score, *size, i-1, 0, 0, 'p','r');
    }
}/*}}}*/

static void sprt_2bp_1pt(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    j = prm_exists(p, p->pi[1]-1);
    if (j > 0) {
        j++;
        i = prm_exists(p, p->pi[j]-1);
        if (2 <= i+1 && i+1 < j && is_P_breakpoint(p, i, i+1) &&
                            is_P_breakpoint(p, j-1, j) && p->pi[i+1] != 1) {
            (*size)++;
            fill(score, *size, i+1, j, 0, 'p', 't');
        }
    }
}/*}}}*/

static void sprt_1bp_1pt(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    for (i = 2; i <= p->size; i++) {
        if (!is_P_breakpoint(p, i-1, i))
            continue;

        j = prm_exists(p, p->pi[i-1]+1);
        if (i < j && is_P_breakpoint(p, j-1, j) && !sprt_prohibited(p, 't', i, j)) {
            (*size)++;
            fill(score, *size, i, j, 0, 'p', 't');
        }
    }

    j = prm_exists(p, p->pi[1]-1);
    if (j > 0 && is_P_breakpoint(p, j, j+1)) {
        j++;
        for (i = 2; i < j; i++) {
            if (is_P_breakpoint(p, i-1, i) && !sprt_prohibited(p, 't', i, j)) {
                (*size)++;
                fill(score, *size, i, j, 0, 'p', 't');
            }
        }
    }
}/*}}}*/


static void spsrt_1bp_1pr(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i;

    i = prm_exists(p, -dec_mod(p->pi[1], p->size));
    if (i >= 2 && i <= p->size && is_PSRT_breakpoint(p, i-1, i)) {
        (*size)++;
        fill(score, *size, i-1, 0, 0, 'p','r');
    }
}/*}}}*/

static void spsrt_2bp_1pt(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    j = prm_exists(p, dec_mod(p->pi[1], p->size));
    if (j > 0) {
        j++;
        i = prm_exists(p, dec_mod(p->pi[j], p->size));
        if (i >= 1 && i+1 < j && j <= p->size &&
                    is_PSRT_breakpoint(p, i, i+1) &&
                    is_PSRT_breakpoint(p, j-1, j)) {
            (*size)++;
            fill(score, *size, i+1, j, 0, 'p', 't');
        }
    }
}/*}}}*/

static void spsrt_1bp_1pt(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    for (i = 2; i <= p->size-1; i++) {
        if (!is_PSRT_breakpoint(p, i-1, i))
            continue;

        j = prm_exists(p, inc_mod(p->pi[i-1], p->size));
        if (i < j && j <= p->size && is_PSRT_breakpoint(p, j-1, j)) {
            (*size)++;
            fill(score, *size, i, j, 0, 'p', 't');
        }
    }

    j = prm_exists(p, dec_mod(p->pi[1], p->size));
    if (j == p->size || (j > 1 && j < p->size && is_PSRT_breakpoint(p, j, j+1))) {
        j++;
        for (i = 2; i < j; i++) {
            if (is_PSRT_breakpoint(p, i-1, i)) {
                (*size)++;
                fill(score, *size, i, j, 0, 'p', 't');
            }
        }
    }
}/*}}}*/

static void spsrt_1bp_1sr(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i;

    i = prm_exists(p, -inc_mod(p->pi[p->size], p->size));
    if (i >= 1 && i <= p->size-1 && is_PSRT_breakpoint(p, i, i+1)) {
        (*size)++;
        fill(score, *size, i+1, 0, 0, 's','r');
    }
}/*}}}*/

static void spsrt_2bp_1st(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    i = prm_exists(p, inc_mod(p->pi[p->size], p->size));
    if (i > 0) {
        j = prm_exists(p, inc_mod(p->pi[i-1], p->size));
        if (j > 0 && i < j && j <= p->size &&
                    is_PSRT_breakpoint(p, i-1, i) &&
                    is_PSRT_breakpoint(p, j-1, j)) {
            (*size)++;
            fill(score, *size, i, j, 0, 's', 't');
        }
    }
}/*}}}*/

static void spsrt_1bp_1st(permutation_t *p, int bps, int *score[5], int *size) {/*{{{*/
    int i, j;

    for (i = 2; i <= p->size-1; i++) {
        if (!is_PSRT_breakpoint(p, i-1, i))
            continue;

        j = prm_exists(p, inc_mod(p->pi[i-1], p->size));
        if (j > 2 && i < j && j <= p->size && is_PSRT_breakpoint(p, j-1, j)) {
            (*size)++;
            fill(score, *size, i, j, 0, 's', 't');
        }
    }

    i = prm_exists(p, inc_mod(p->pi[p->size], p->size));
    if (i == 1 || (i > 0 && i < p->size && is_PSRT_breakpoint(p, i-1, i))) {
        for (j = i+1; j <= p->size; j++) {
            if (is_PSRT_breakpoint(p, j-1, j)) {
                (*size)++;
                fill(score, *size, i, j, 0, 's', 't');
            }
        }
    }
}/*}}}*/





static int score_1pr(char prob[], permutation_t *p, int i, int bps) {/*{{{*/
    permutation_t pcpy;
    info_t info;
    init_info(&info);

    create_permutation(&pcpy, p->size, p->sig);
    prmcpy(&pcpy, p);
    op_prefix_reversal(&pcpy, i);
    alg_approx(prob, &pcpy, 0, &info);
    destroy_permutation(&pcpy);
    return info.nmoves;
}/*}}}*/

static int score_2pr(char prob[], permutation_t *p, int i, int j, int bps) {/*{{{*/
    permutation_t pcpy;
    info_t info;
    init_info(&info);

    create_permutation(&pcpy, p->size, p->sig);
    prmcpy(&pcpy, p);
    op_prefix_reversal(&pcpy, i);
    op_prefix_reversal(&pcpy, j);
    alg_approx(prob, &pcpy, 0, &info);
    destroy_permutation(&pcpy);
    return info.nmoves;
}/*}}}*/

static int score_pt(char prob[], permutation_t *p, int i, int j, int bps) {/*{{{*/
    permutation_t pcpy;
    info_t info;
    init_info(&info);

    create_permutation(&pcpy, p->size, p->sig);
    prmcpy(&pcpy, p);
    op_prefix_transposition(&pcpy, i, j);
    alg_approx(prob, &pcpy, 0, &info);
    destroy_permutation(&pcpy);
    return info.nmoves;
}/*}}}*/

static int score_1sr(char prob[], permutation_t *p, int i, int bps) {/*{{{*/
    permutation_t pcpy;
    info_t info;
    init_info(&info);

    create_permutation(&pcpy, p->size, p->sig);
    prmcpy(&pcpy, p);
    op_suffix_reversal(&pcpy, i);
    alg_approx(prob, &pcpy, 0, &info);
    destroy_permutation(&pcpy);
    return info.nmoves;
}/*}}}*/

static int score_2sr(char prob[], permutation_t *p, int i, int j, int bps) {/*{{{*/
    permutation_t pcpy;
    info_t info;
    init_info(&info);

    create_permutation(&pcpy, p->size, p->sig);
    prmcpy(&pcpy, p);
    op_suffix_reversal(&pcpy, i);
    op_suffix_reversal(&pcpy, j);
    alg_approx(prob, &pcpy, 0, &info);
    destroy_permutation(&pcpy);
    return info.nmoves;
}/*}}}*/

static int score_st(char prob[], permutation_t *p, int i, int j, int bps) {/*{{{*/
    permutation_t pcpy;
    info_t info;
    init_info(&info);

    create_permutation(&pcpy, p->size, p->sig);
    prmcpy(&pcpy, p);
    op_suffix_transposition(&pcpy, i, j);
    alg_approx(prob, &pcpy, 0, &info);
    destroy_permutation(&pcpy);
    return info.nmoves;
}/*}}}*/

static int calc_score(char prob[], permutation_t *p, int type, int rearr,
                int i, int j, int bps) {/*{{{*/
    if (type == 'p') {
        if (rearr == 'r' && j == 0)
            return score_1pr(prob, p, i, bps);
        else if (rearr == 'r' && j != 0)
            return score_2pr(prob, p, i, j, bps);
        else
            return score_pt(prob, p, i, j, bps);
    } else {
        if (rearr == 'r' && j == 0)
            return score_1sr(prob, p, i, bps);
        else if (rearr == 'r' && j != 0)
            return score_2sr(prob, p, i, j, bps);
        else
            return score_st(prob, p, i, j, bps);
    }
}/*}}}*/


int choose(char prob[], permutation_t *p, int *score[], int size, int bps) {/*{{{*/
    int min, i;

    score[2][0] = calc_score(prob, p, score[3][0], score[4][0], 
                    score[0][0], score[1][0], bps);

    min = 0;
    for (i = 1; i <= size; i++) {
        score[2][i] = calc_score(prob, p, score[3][i], score[4][i], 
                        score[0][i], score[1][i], bps);
        if (score[2][min] > score[2][i])
            min = i;
    }

    return min;
}/*}}}*/

static void alloc_score(int *score[5], char prob[], int xbps, int yop, int n) {/*{{{*/
    int i, size;

    switch (strlen(prob)) {
        case 2: /* pr, pt */
            if (prob[1] == 'r') { /* pr */
                size = 2;
                if (yop == 2)
                    size = 4*n;
            } else { /* pt */
                size = 1;
                if (xbps == 1)
                    size = 2*n;
            }
            break;

        case 3: /* prt, psr, pst */
            if (prob[1] == 'r') { /* prt */
                size = 4;
                if (xbps == 1)
                    size = 4*n + 4;
            } else if (prob[2] == 't') { /* pst */
                size = 2;
                if (xbps == 1)
                    size = 4*n;
            } else { /* psr */
                size = 4;
                if (yop == 2)
                    size = 8*n;
            }
            break;

        case 4: /* psrt */
            size = 8;
            if (xbps == 1)
                size = 8*n + 8;
            break;

        case 5: /* spsrt */
            size = 2;
            if (xbps == 1)
                size = 4*n+1;
            break;

        default:
            printf("ERRO! alloc_score(%s)\n", prob);
            exit(0);
    }

    for (i = 0; i < 5; i++)
        score[i] = Malloc(sizeof(int) * size);
}/*}}}*/

static void fill_score(char prob[], int xbps, int yop, permutation_t *p, int bps,
                int *score[], int *size) {/*{{{*/

    switch (strlen(prob)) {
        case 2: /* pr, pt */
            if (prob[1] == 'r') { /* pr */
                if (yop == 1)
                    pr_1bp_1pr(p, bps, score, size);
                else
                    pr_1bp_2pr(p, bps, score, size);
            } else { /* pt */
                if (xbps == 2)
                    pt_2bp_1pt(p, bps, score, size);
                else
                    pt_1bp_1pt(p, bps, score, size);
            }
            break;

        case 3: /* prt, psr, pst, spr */
            if (prob[0] == 's') { /* spr */
                if (yop == 1)
                    spr_1bp_1pr(p, bps, score, size);
                else
                    spr_1bp_2pr(p, bps, score, size);
            } else if (prob[1] == 'r') { /* prt */
                if (xbps == 2) {
                    prt_2bp_1pt(p, bps, score, size);
                } else {
                    prt_1bp_1pt(p, bps, score, size);
                    prt_1bp_1pr(p, bps, score, size);
                }
            } else if (prob[2] == 'r') { /* psr */
                if (yop == 1) {
                    psr_1bp_1pr(p, bps, score, size);
                    psr_1bp_1sr(p, bps, score, size);
                } else {
                    psr_1bp_2pr(p, bps, score, size);
                    psr_1bp_2sr(p, bps, score, size);
                }
            } else { /* pst */
                if (xbps == 2) {
                    pst_2bp_1pt(p, bps, score, size);
                    pst_2bp_1st(p, bps, score, size);
                } else {
                    pst_1bp_1pt(p, bps, score, size);
                    pst_1bp_1st(p, bps, score, size);
                }
            }
            break;

        case 4: /* psrt, spsr, sprt */
            if (prob[0] == 'p') { /* psrt */
                if (xbps == 2) {
                    psrt_2bp_1pt(p, bps, score, size);
                    psrt_2bp_1st(p, bps, score, size);
                } else {
                    psrt_1bp_1pt(p, bps, score, size);
                    psrt_1bp_1st(p, bps, score, size);
                    psrt_1bp_1pr(p, bps, score, size);
                    psrt_1bp_1sr(p, bps, score, size);
                }
            } else if (prob[2] == 's') { /* spsr */
                if (yop == 1) {
                    spsr_1bp_1pr(p, bps, score, size);
                    spsr_1bp_1sr(p, bps, score, size);
                } else {
                    spsr_1bp_2pr(p, bps, score, size);
                    spsr_1bp_2sr(p, bps, score, size);
                }
            } else { /* sprt */
                if (xbps == 2) {
                    sprt_2bp_1pt(p, bps, score, size);
                } else {
                    sprt_1bp_1pt(p, bps, score, size);
                    sprt_1bp_1pr(p, bps, score, size);
                }
            }
            break;

        case 5: /* spsrt */
            if (xbps == 2) {
                spsrt_2bp_1pt(p, bps, score, size);
                spsrt_2bp_1st(p, bps, score, size);
            } else {
                spsrt_1bp_1pt(p, bps, score, size);
                spsrt_1bp_1st(p, bps, score, size);
                spsrt_1bp_1pr(p, bps, score, size);
                spsrt_1bp_1sr(p, bps, score, size);
            }
            break;

        default:
            printf("ERRO! fill_score(%s)\n", prob);
            exit(0);
    }
}/*}}}*/

static int remove_Xbps_Yop(char prob[], int xbps, int yop, permutation_t *p, int bps,
                char *rear, char *type, int *iret, int *jret) {/*{{{*/
    int i, size, *score[5];
    
    alloc_score(score, prob, xbps, yop, p->size);
    
    size = -1;
    fill_score(prob, xbps, yop, p, bps, score, &size);

    if (size == 0) {
        (*iret) = score[0][0];
        (*jret) = score[1][0];
        (*type) = score[3][0];
        (*rear) = score[4][0];
    } else if (size != -1) {
        i = choose(prob, p, score, size, bps);
        (*iret) = score[0][i];
        (*jret) = score[1][i];
        (*type) = score[3][i];
        (*rear) = score[4][i];
    }

    for (i = 0; i < 5; i++)
        free(score[i]);

    if (size == -1)
        return 0;
    return 1;
}/*}}}*/




void alg_2PRx(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, bps;
    char rear, type;

    bps = nof_UPR_breakpoints(p);
    while (!is_identity(p)) {
        if (remove_Xbps_Yop("pr", 1, 1, p, bps, &rear, &type, &i, &j)) {
            bps = nofbps_after_pr(p, i, bps);
            op_prefix_reversal(p, i);
            print_verbose(p, verbose, 'r', 1, i, 0);
            info->nmoves++;
        } else if (remove_Xbps_Yop("pr", 1, 2, p, bps, &rear, &type, &i, &j)) {
            bps = nofbps_after_pr(p, i, bps);
            op_prefix_reversal(p, i);
            print_verbose(p, verbose, 'r', 1, i, 0);
            bps = nofbps_after_pr(p, j, bps);
            op_prefix_reversal(p, j);
            print_verbose(p, verbose, 'r', 1, j, 0);
            info->nmoves += 2;
        } else {
            pr_star_perm(p, verbose, info);
        }
    }

    info->nof_pr = info->nmoves;
}/*}}}*/

void alg_2PSRx(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, bps;
    char rear, type;

    bps = nof_UPSR_breakpoints(p);
    while (!is_identity(p)) {
        if (is_reverse(p)) {
            op_prefix_reversal(p, p->size);
            print_verbose(p, verbose, 'r', 1, p->size, 0);
            info->nmoves++;
            info->nof_pr++;
        } else if (remove_Xbps_Yop("psr", 1, 1, p, bps, &rear, &type, &i, &j)) {
            bps = nofbps_after_psr(p, type, i, bps);
            if (type == 'p') {
                op_prefix_reversal(p, i);
                print_verbose(p, verbose, 'r', 1, i, 0);
                info->nof_pr++;
            } else {
                op_suffix_reversal(p, i);
                print_verbose(p, verbose, 'r', i, p->size, 0);
                info->nof_sr++;
            }
            info->nmoves++;
        } else if (remove_Xbps_Yop("psr", 1, 2, p, bps, &rear, &type, &i, &j)) {
            if (type == 'p') {
                bps = nofbps_after_psr(p, type, i, bps);
                op_prefix_reversal(p, i);
                print_verbose(p, verbose, 'r', 1, i, 0);
                bps = nofbps_after_psr(p, type, j, bps);
                op_prefix_reversal(p, j);
                print_verbose(p, verbose, 'r', 1, j, 0);
                info->nof_pr += 2;
            } else {
                bps = nofbps_after_psr(p, type, i, bps);
                op_suffix_reversal(p, i);
                print_verbose(p, verbose, 'r', i, p->size, 0);
                bps = nofbps_after_psr(p, type, j, bps);
                op_suffix_reversal(p, j);
                print_verbose(p, verbose, 'r', j, p->size, 0);
                info->nof_sr += 2;
            }
            info->nmoves += 2;
        } else {
            if ((bps % 2 == 1 && p->inv_pi[1] > p->inv_pi[p->size]) ||
                        (bps % 2 == 0 && p->inv_pi[1] < p->inv_pi[p->size])) {
                op_prefix_reversal(p, p->size);
                print_verbose(p, verbose, 'r', 1, p->size, 0);
                info->nmoves++;
                info->nof_pr++;
            }
            psr_star_perm(p, verbose, info, bps);
        }
    }
}/*}}}*/

void alg_2PTx(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, bps;
    char rear, type;

    bps = nof_P_breakpoints(p);
    while (!is_identity(p)) {
        if (remove_Xbps_Yop("pt", 2, 1, p, bps, &rear, &type, &i, &j)) {
            bps = nofbps_after_pt(p, i, j, bps);
            op_prefix_transposition(p, i, j);
            print_verbose(p, verbose, 't', 1, i, j);
        } else if (remove_Xbps_Yop("pt", 1, 1, p, bps, &rear, &type, &i, &j)) {
            bps = nofbps_after_pt(p, i, j, bps);
            op_prefix_transposition(p, i, j);
            print_verbose(p, verbose, 't', 1, i, j);
        }
        info->nmoves++;
    }

    info->nof_pt = info->nmoves;
}/*}}}*/

void alg_2PSTx(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, bps;
    char rear, type;

    bps = nof_PS_breakpoints(p);
    while (!is_identity(p)) {
        if (remove_Xbps_Yop("pst", 2, 1, p, bps, &rear, &type, &i, &j)) {
            bps = nofbps_after_pst(p, type, i, j, bps);
            if (type == 'p') {
                op_prefix_transposition(p, i, j);
                print_verbose(p, verbose, 't', 1, i, j);
                info->nof_pt++;
            } else {
                op_suffix_transposition(p, i, j);
                print_verbose(p, verbose, 't', i, j, p->size+1);
                info->nof_st++;
            }
        } else if (remove_Xbps_Yop("pst", 1, 1, p, bps, &rear, &type, &i, &j)) {
            bps = nofbps_after_pst(p, type, i, j, bps);
            if (type == 'p') {
                op_prefix_transposition(p, i, j);
                print_verbose(p, verbose, 't', 1, i, j);
                info->nof_pt++;
            } else {
                op_suffix_transposition(p, i, j);
                print_verbose(p, verbose, 't', i, j, p->size+1);
                info->nof_st++;
            }
        }
        info->nmoves++;
    }
}/*}}}*/

void alg_2PRTx(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, bps;
    char rear, type;

    bps = nof_UPR_breakpoints(p);
    while (!is_identity(p)) {
        if (p->pi[1] == 1) {
            i = 1;
            while (i <= p->size && !is_UPR_breakpoint(p, i, i+1))
                i++;
            op_prefix_transposition(p, i+1, p->size+1);
            print_verbose(p, verbose, 't', 1, i+1, p->size+1);
            info->nof_pt++;
        } else if (remove_Xbps_Yop("prt", 2, 1, p, bps, &rear, &type, &i, &j)) {
            bps = nofbps_after_prt(p, rear, i, j, bps);
            op_prefix_transposition(p, i, j);
            print_verbose(p, verbose, 't', 1, i, j);
            info->nof_pt++;
        } else if (remove_Xbps_Yop("prt", 1, 1, p, bps, &rear, &type, &i, &j)) {
            bps = nofbps_after_prt(p, rear, i, j, bps);
            if (rear == 't') {
                op_prefix_transposition(p, i, j);
                print_verbose(p, verbose, 't', 1, i, j);
                info->nof_pt++;
            } else {
                op_prefix_reversal(p, i);
                print_verbose(p, verbose, 'r', 1, i, 0);
                info->nof_pr++;
            }
        }
        info->nmoves++;
    }
}/*}}}*/

void alg_2PSRTx(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, bps;
    char rear, type;

    bps = nof_UPSRT_breakpoints(p);
    while (!is_identity(p)) {
        if (is_reverse(p)) {
            op_prefix_reversal(p, p->size);
            print_verbose(p, verbose, 'r', 1, p->size, 0);
            info->nof_pr++;
        }  else if (remove_Xbps_Yop("psrt", 2, 1, p, bps, &rear, &type, &i, &j)) {
            bps = nofbps_after_psrt(p, rear, type, i, j, bps);
            if (type == 'p') {
                op_prefix_transposition(p, i, j);
                print_verbose(p, verbose, 't', 1, i, j);
                info->nof_pt++;
            } else {
                op_suffix_transposition(p, i, j);
                print_verbose(p, verbose, 't', i, j, p->size+1);
                info->nof_st++;
            }
        } else if (remove_Xbps_Yop("psrt", 1, 1, p, bps, &rear, &type, &i, &j)) {
            bps = nofbps_after_psrt(p, rear, type, i, j, bps);
            if (type == 'p' && rear == 't') {
                op_prefix_transposition(p, i, j);
                print_verbose(p, verbose, 't', 1, i, j);
                info->nof_pt++;
            } else if (type == 'p' && rear == 'r') {
                op_prefix_reversal(p, i);
                print_verbose(p, verbose, 'r', 1, i, 0);
                info->nof_pr++;
            } else if (type == 's' && rear == 't') {
                op_suffix_transposition(p, i, j);
                print_verbose(p, verbose, 't', i, j, p->size+1);
                info->nof_st++;
            } else {
                op_suffix_reversal(p, i);
                print_verbose(p, verbose, 'r', i, p->size, 0);
                info->nof_sr++;
            }
        } else {
            /* permutation has two strips */
            i = p->inv_pi[1];
            if (p->size == p->pi[i+1]) {
                op_prefix_transposition(p, i+1, p->size+1);
                print_verbose(p, verbose, 't', 1, i+1, p->size+1);
            } else {
                op_prefix_transposition(p, i, p->size+1);
                print_verbose(p, verbose, 't', 1, i, p->size+1);
            }
            info->nof_pt++;
        }
        info->nmoves++;
    }
}/*}}}*/




void alg_2SPRx(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, bps;
    char rear, type;

    bps = nof_P_breakpoints(p);
    while (!is_identity(p)) {
        if (remove_Xbps_Yop("spr", 1, 1, p, bps, &rear, &type, &i, &j)) {
            bps = nofbps_after_spr(p, i, bps);
            op_prefix_reversal(p, i);
            print_verbose(p, verbose, 'r', 1, i, 0);
            info->nmoves++;
        } else if (remove_Xbps_Yop("spr", 1, 2, p, bps, &rear, &type, &i, &j)) {
            bps = nofbps_after_spr(p, i, bps);
            op_prefix_reversal(p, i);
            print_verbose(p, verbose, 'r', 1, i, 0);
            bps = nofbps_after_spr(p, j, bps);
            op_prefix_reversal(p, j);
            print_verbose(p, verbose, 'r', 1, j, 0);
            info->nmoves += 2;
        } else {
            spr_star_perm(p, verbose, info);
        }
    }

    info->nof_pr = info->nmoves;
}/*}}}*/

void alg_2SPSRx(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, bps;
    char rear, type;

    bps = nof_PS_breakpoints(p);
    while (!is_identity(p)) {
        if (is_signed_reverse(p)) {
            op_prefix_reversal(p, p->size);
            print_verbose(p, verbose, 'r', 1, p->size, 0);
            info->nmoves++;
        } else if (remove_Xbps_Yop("spsr", 1, 1, p, bps, &rear, &type, &i, &j)) {
            info->nof_pr++;
            bps = nofbps_after_spsr(p, type, i, bps);
            if (type == 'p') {
                op_prefix_reversal(p, i);
                print_verbose(p, verbose, 'r', 1, i, 0);
                info->nof_pr++;
            } else {
                op_suffix_reversal(p, i);
                print_verbose(p, verbose, 'r', i, p->size, 0);
                info->nof_sr++;
            }
            info->nmoves++;
        } else if (remove_Xbps_Yop("spsr", 1, 2, p, bps, &rear, &type, &i, &j)) {
            if (type == 'p') {
                bps = nofbps_after_spsr(p, type, i, bps);
                op_prefix_reversal(p, i);
                print_verbose(p, verbose, 'r', 1, i, 0);
                bps = nofbps_after_spsr(p, type, j, bps);
                op_prefix_reversal(p, j);
                print_verbose(p, verbose, 'r', 1, j, 0);
                info->nof_pr += 2;
            } else {
                bps = nofbps_after_spsr(p, type, i, bps);
                op_suffix_reversal(p, i);
                print_verbose(p, verbose, 'r', i, p->size, 0);
                bps = nofbps_after_spsr(p, type, j, bps);
                op_suffix_reversal(p, j);
                print_verbose(p, verbose, 'r', j, p->size, 0);
                info->nof_sr += 2;
            }
            info->nmoves += 2;
        } else {
            if ((bps % 2 == 1 && prm_exists(p, 1) &&
                                prm_exists(p, 1) > prm_exists(p, p->size)) || 
                        (bps % 2 == 0 && prm_exists(p, -1) &&
                         prm_exists(p, -1) < prm_exists(p, -p->size))) {
                op_prefix_reversal(p, p->size);
                print_verbose(p, verbose, 'r', 1, p->size, 0);
                info->nmoves++;
                info->nof_pr++;
            }
            spsr_star_perm(p, verbose, info, bps);
        }
    }
}/*}}}*/

void alg_2SPRTx(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, bps;
    char rear, type;

    bps = nof_P_breakpoints(p);
    while (!is_identity(p)) {
        if (p->pi[1] == 1) {
            i = 1;
            while (i <= p->size && p->pi[i+1] - p->pi[i] == 1) i++;
            op_prefix_transposition(p, i+1, p->size+1);
            print_verbose(p, verbose, 't', 1, i+1, p->size+1);
            info->nof_pt++;
        } else if (remove_Xbps_Yop("sprt", 2, 1, p, bps, &rear, &type, &i, &j)) {
            bps = nofbps_after_sprt(p, 't', i, j, bps);
            op_prefix_transposition(p, i, j);
            print_verbose(p, verbose, 't', 1, i, j);
            info->nof_pt++;
        } else if (remove_Xbps_Yop("sprt", 1, 1, p, bps, &rear, &type, &i, &j)) {
            bps = nofbps_after_sprt(p, rear, i, j, bps);
            if (rear == 't') {
                op_prefix_transposition(p, i, j);
                print_verbose(p, verbose, 't', 1, i, j);
                info->nof_pt++;
            } else {
                op_prefix_reversal(p, i);
                print_verbose(p, verbose, 'r', 1, i, 0);
                info->nof_pr++;
            }
        }
        info->nmoves++;
    }
}/*}}}*/

void alg_2SPSRTx(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, bps;
    char rear, type;

    bps = nof_PSRT_breakpoints(p);
    while (!is_identity(p)) {
        if (is_signed_reverse(p)) {
            op_prefix_reversal(p, p->size);
            print_verbose(p, verbose, 'r', 1, p->size, 0);
            info->nof_pr++;
        } else if (remove_Xbps_Yop("spsrt", 2, 1, p, bps, &rear, &type, &i, &j)) {
            bps = nofbps_after_spsrt(p, rear, type, i, j, bps);
            if (type == 'p') {
                op_prefix_transposition(p, i, j);
                print_verbose(p, verbose, 't', 1, i, j);
                info->nof_pt++;
            } else {
                op_suffix_transposition(p, i, j);
                print_verbose(p, verbose, 't', i, j, p->size+1);
                info->nof_st++;
            }
        } else if (remove_Xbps_Yop("spsrt", 1, 1, p, bps, &rear, &type, &i, &j)) {
            bps = nofbps_after_spsrt(p, rear, type, i, j, bps);
            if (type == 'p' && rear == 't') {
                op_prefix_transposition(p, i, j);
                print_verbose(p, verbose, 't', 1, i, j);
                info->nof_pt++;
            } else if (type == 'p' && rear == 'r') {
                op_prefix_reversal(p, i);
                print_verbose(p, verbose, 'r', 1, i, 0);
                info->nof_pr++;
            } else if (type == 's' && rear == 't') {
                op_suffix_transposition(p, i, j);
                print_verbose(p, verbose, 't', i, j, p->size+1);
                info->nof_st++;
            } else {
                op_suffix_reversal(p, i);
                print_verbose(p, verbose, 'r', i, p->size, 0);
                info->nof_sr++;
            }
        } else {
            i = prm_exists(p, 1);
            if (i > 0) {
                op_prefix_transposition(p, i, p->size+1);
                print_verbose(p, verbose, 't', 1, i, p->size+1);
            } else {
                i = prm_exists(p, -p->size);
                op_prefix_transposition(p, i, p->size+1);
                print_verbose(p, verbose, 't', 1, i, p->size+1);
            }
            info->nof_pt++;
        }
        info->nmoves++;
    }
}/*}}}*/



