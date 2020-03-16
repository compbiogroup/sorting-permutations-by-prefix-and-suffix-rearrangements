/***********************************************************
 * Created: Mon 04 Feb 2013 03:35:22 PM BRST
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************/


#include <stdio.h>
#include <stdlib.h>
#include "../permutations.h"
#include "../rearrangements.h"
#include "../breakpoints.h"
#include "../util.h"



void alg_2PRT(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, kab, ja, jb, kcd, jc, jd, x;

    while (!is_identity(p)) {
        i = 1;
        while (i <= p->size && !is_UPR_breakpoint(p, i, i+1))
            i++;

        if (p->pi[1] == 1) {
            op_prefix_transposition(p, i+1, p->size+1);
            print_verbose(p, verbose, 't', 1, i+1, p->size+1);
            info->nof_pt++;
        } else {
            /* Tries to remove two breakpoints without bringing 1 to the
             * beginning */
            /* x indicates if this operation is made or not */
            x = 0;
            ja = jb = jc = jd = -1;

            kab = p->inv_pi[p->pi[1]-1] + 1;
            if (kab > 0 && kab <= p->size+1) {
                ja = p->inv_pi[p->pi[kab]-1] + 1;
            }
            if (kab >= 0 && kab < p->size+1) {
                jb = p->inv_pi[p->pi[kab]+1] + 1;
            }
            kcd = p->inv_pi[p->pi[1]+1] + 1;
            if (kcd > 0 && kcd <= p->size+1) {
                jc = p->inv_pi[p->pi[kcd]-1] + 1;
            }
            if (kcd >= 0 && kcd < p->size+1) {
                jd = p->inv_pi[p->pi[kcd]+1] + 1;
            }

            if (kab <= p->size+1 && is_UPR_breakpoint(p, kab-1, kab)) {
                if (ja < kab && ja > 1 && p->pi[ja] != 1 && is_UPR_breakpoint(p, ja-1, ja)) {
                    op_prefix_transposition(p, ja, kab);
                    print_verbose(p, verbose, 't', 1, ja, kab);
                    info->nof_pt++;
                    x = 1;
                } else if (jb < kab && jb > 1 && p->pi[jb] != 1 && is_UPR_breakpoint(p, jb-1, jb)) {
                    op_prefix_transposition(p, jb, kab);
                    print_verbose(p, verbose, 't', 1, jb, kab);
                    info->nof_pt++;
                    x = 1;
                }
            }
            
            if (!x && kcd <= p->size+1 && is_UPR_breakpoint(p, kcd-1, kcd)) {
                if (jc < kcd && jc > 1 && p->pi[jc] != 1 && is_UPR_breakpoint(p, jc-1, jc)) {
                    op_prefix_transposition(p, jc, kcd);
                    print_verbose(p, verbose, 't', 1, jc, kcd);
                    info->nof_pt++;
                    x = 1;
                } else if (jd < kcd && jd > 1 && p->pi[jd] != 1 && is_UPR_breakpoint(p, jd-1, jd)) {
                    op_prefix_transposition(p, jd, kcd);
                    print_verbose(p, verbose, 't', 1, jd, kcd);
                    info->nof_pt++;
                    x = 1;
                }
            }
            
            if (x == 0) {
                /* Remove one breakpoint */
                if (p->pi[1] <= p->pi[i]) {
                    x = p->inv_pi[p->pi[1]-1];
                    if (p->pi[x] == p->pi[x+1]+1) {
                        op_prefix_reversal(p, x-1);
                        print_verbose(p, verbose, 'r', 1, x-1, 0);
                        info->nof_pr++;
                    } else {
                        op_prefix_transposition(p, i+1, x+1);
                        print_verbose(p, verbose, 't', 1, i+1, x+1);
                        info->nof_pt++;
                    }
                } else {
                    x = p->inv_pi[p->pi[1]+1];
                    if (p->pi[x] == p->pi[x-1]-1) {
                        op_prefix_transposition(p, i+1, x+1);
                        print_verbose(p, verbose, 't', 1, i+1, x+1);
                        info->nof_pt++;
                    } else {
                        op_prefix_reversal(p, x-1);
                        print_verbose(p, verbose, 'r', 1, x-1, 0);
                        info->nof_pr++;
                    }
                }
            }
        }
        info->nmoves++;
    }
}/*}}}*/


void alg_2SPRT(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j;

    while (!is_identity(p)) {

        if (p->pi[1] == 1) {
            i = 1;
            while (i <= p->size && !is_P_breakpoint(p, i, i+1)) i++;
            op_prefix_transposition(p, i+1, p->size+1);
            print_verbose(p, verbose, 't', 1, i+1, p->size+1);
            info->nof_pt++;
            info->nmoves++;
            continue;
        }

        /* tries to remove 2 breakpoints */
        j = prm_exists(p, p->pi[1]-1);
        if (j > 0) {
            j++;
            i = prm_exists(p, p->pi[j]-1);
            if (2 <= i+1 && i+1 < j && is_P_breakpoint(p, i, i+1) &&
                            is_P_breakpoint(p, j-1, j) && p->pi[i+1] != 1) {
                op_prefix_transposition(p, i+1, j);
                print_verbose(p, verbose, 't', 1, i+1, j);
                info->nmoves++;
                info->nof_pt++;
                continue;
            }
        }

        /* tries to remove 1 breakpoint */
        j = prm_exists(p, p->pi[1]-1);
        if (j > 0 && is_P_breakpoint(p, j, j+1)) {
            i = 2;
            while (!is_P_breakpoint(p, i-1, i)) i++;
            if (2 <= i && i < j+1) {
                op_prefix_transposition(p, i, j+1);
                print_verbose(p, verbose, 't', 1, i, j+1);
                info->nmoves++;
                info->nof_pt++;
                continue;
            }
        }

        i = prm_exists(p, -p->pi[1]+1);
        if (i > 0 && is_P_breakpoint(p, i-1, i)) {
            op_prefix_reversal(p, i-1);
            print_verbose(p, verbose, 'r', 1, i-1, 0);
            info->nmoves++;
            info->nof_pr++;
            continue;
        }

    }
}/*}}}*/


