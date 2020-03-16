/***********************************************************
 * Created: Mon 04 Feb 2013 03:35:22 PM BRST
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include "../util.h"
#include "../permutations.h"
#include "../rearrangements.h"
#include "../breakpoints.h"


void pr_star_perm(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, l, *lengths, nofbps, new_n;

    nofbps = nof_UPR_breakpoints(p);
    lengths = Malloc(sizeof(int) * (nofbps + 1));
    lengths[0] = 0;
    new_n = 0;
    l = 2;

    /* collecting lengths of the strips and new value of n (the end of the
     * permutation might be ordered correctly) */
    for (i = 1; i <= nofbps; i++) {
        for (; p->pi[l] == p->pi[l-1]-1; l++) ;
        lengths[i] = l-1 - new_n;
        new_n += lengths[i];
        l++;
    }

    for (i = 1; i <= nofbps; i++) {
        op_prefix_reversal(p, new_n);
        print_verbose(p, verbose, 'r', 1, new_n, 0);
        op_prefix_reversal(p, new_n - lengths[i]);
        print_verbose(p, verbose, 'r', 1, new_n - lengths[i], 0);
        info->nmoves += 2;
    }

    free(lengths);
}/*}}}*/

int pr_edge_type1(int *edge, permutation_t *p) {/*{{{*/
    /* Edge (\pi_i, \pi_j) of type 1 must have i == 1 */
    int j1, j2;

    j1 = p->inv_pi[p->pi[1] - 1];
    j2 = p->inv_pi[p->pi[1] + 1];

    if (j1 > 2 && is_UPR_breakpoint(p, j1-1, j1)) {
        if (j2 > 2 && is_UPR_breakpoint(p, j2-1, j2)) {
            if (j1 > j2) {
                edge[0] = 1;
                edge[1] = j1;
                return 1;
            } else {
                edge[0] = 1;
                edge[1] = j2;
                return 1;
            }
        } else {
            edge[0] = 1;
            edge[1] = j1;
            return 1;
        }
    } else if (j2 > 2 && is_UPR_breakpoint(p, j2-1, j2)) {
        edge[0] = 1;
        edge[1] = j2;
        return 1;
    }

    return 0;
}/*}}}*/

int pr_edge_type2(int *edge, permutation_t *p) {/*{{{*/
    /* Edge (\pi_i, \pi_j) of type 2 must have i != 0 */
    int i, j;

    for (j = p->size; j >= 2; j--) {
        if (!is_UPR_breakpoint(p, j, j+1))
            continue;

        i = p->inv_pi[p->pi[j] - 1];
        if (i != 0 && j > i+1 && is_UPR_breakpoint(p, i, i+1)) {
            edge[0] = i;
            edge[1] = j;
            return 1;
        }

        i = p->inv_pi[p->pi[j] + 1];
        if (i != 0 && j > i+1 && is_UPR_breakpoint(p, i, i+1)) {
            edge[0] = i;
            edge[1] = j;
            return 1;
        }
    }

    return 0;
}/*}}}*/

int pr_edge_type3A(int *edge, permutation_t *p) {/*{{{*/
    /* Edge (\pi_i, \pi_j) of type 3A must have i == 1 */
    int j1, j2;

    if (!is_UPR_breakpoint(p, 1, 2))
        return 0;

    j1 = p->inv_pi[p->pi[1] - 1];
    j2 = p->inv_pi[p->pi[1] + 1];

    if (j1 > 2 && is_UPR_breakpoint(p, j1-1, j1)) {
        if (j2 > 2 && is_UPR_breakpoint(p, j2-1, j2)) {
            if (j1 > j2) {
                edge[0] = 1;
                edge[1] = j1;
                return 1;
            } else {
                edge[0] = 1;
                edge[1] = j2;
                return 1;
            }
        } else {
            edge[0] = 1;
            edge[1] = j1;
            return 1;
        }
    } else if (j2 > 2 && is_UPR_breakpoint(p, j2-1, j2)) {
        edge[0] = 1;
        edge[1] = j2;
        return 1;
    }

    return 0;
}/*}}}*/

int pr_edge_type3B(int *edge, permutation_t *p) {/*{{{*/
    /* Edge (\pi_i, \pi_j) of type 3B must have i > 1 */
    int i, j;

    for (j = p->size+1; j >= 3; j--) {
        if (!is_UPR_breakpoint(p, j-1, j))
            continue;

        i = p->inv_pi[p->pi[j] - 1];
        if (i > 1 && j > i+1 && is_UPR_breakpoint(p, i, i+1)) {
            edge[0] = i;
            edge[1] = j;
            return 1;
        }

        if (p->pi[j] <= p->size) {
            i = p->inv_pi[p->pi[j] + 1];
            if (i > 1 && j > i+1 && is_UPR_breakpoint(p, i, i+1)) {
                edge[0] = i;
                edge[1] = j;
                return 1;
            }
        }
    }
    return 0;
}/*}}}*/

void alg_2PR(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int edge[2], i;

    while (!is_identity(p)) {
        if (pr_edge_type1(edge, p)) {
            i = edge[1] - 1;
            op_prefix_reversal(p, i);
            print_verbose(p, verbose, 'r', 1, i, 0);
            info->nmoves++;
        } else if (pr_edge_type3A(edge, p)) {
            i = edge[1] - 1;
            op_prefix_reversal(p, i);
            print_verbose(p, verbose, 'r', 1, i, 0);
            info->nmoves++;
        } else if (pr_edge_type2(edge, p)) {
            i = edge[1];
            op_prefix_reversal(p, i);
            print_verbose(p, verbose, 'r', 1, i, 0);
            i = edge[1] - edge[0];
            op_prefix_reversal(p, i);
            print_verbose(p, verbose, 'r', 1, i, 0);
            info->nmoves += 2;
        } else if (pr_edge_type3B(edge, p)) { 
            i = edge[0];
            op_prefix_reversal(p, i);
            print_verbose(p, verbose, 'r', 1, i, 0);
            i = edge[1] - 1;
            op_prefix_reversal(p, i);
            print_verbose(p, verbose, 'r', 1, i, 0);
            info->nmoves += 2;
        } else {
            pr_star_perm(p, verbose, info);
        }
    }

    info->nof_pr = info->nmoves;
}/*}}}*/



void spr_star_perm(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, l, *lengths, nofbps, new_n;

    nofbps = nof_P_breakpoints(p);
    lengths = Malloc(sizeof(int) * (nofbps + 1));
    lengths[0] = 0;
    new_n = 0;
    l = 2;

    /* collecting lengths of the strips and new value of n (the end of the
     * permutation might be ordered correctly) */
    for (i = 1; i <= nofbps; i++) {
        for (; p->pi[l] == p->pi[l-1]+1; l++);
        lengths[i] = l-1 - new_n;
        new_n += lengths[i];
        l++;
    }

    l = 0;
    for (i = 1; i <= nofbps; i++) {
        op_prefix_reversal(p, new_n);
        print_verbose(p, verbose, 'r', 1, new_n, 0);
        op_prefix_reversal(p, new_n - lengths[i]);
        print_verbose(p, verbose, 'r', 1, new_n - lengths[i], 0);
        info->nmoves += 2;
    }

    free(lengths);
}/*}}}*/

void alg_2SPR(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, k, aux;

    while (!is_identity(p)) {

        /* tries to remove 1 breakpoint with 1 reversal */
        j = prm_exists(p, -1*p->pi[1]+1);
        if (j > 0) {
            op_prefix_reversal(p, j-1);
            print_verbose(p, verbose, 'r', 1, j-1, 0);
            info->nmoves++;
            continue;
        }

        /* otherwise, tries to remove 1 breakpoint with 2 reversals */
        aux = p->size;
        while (p->pi[aux] == aux)
            aux--;
        /* aux has the value of the greatest element (in module) out of order; 
         * now we must search an element that is smaller than aux and is the
         * greatest positive element in permutation */
        k = 0;
        for (i = 1; i <= aux; i++) {
            if (p->pi[i] > k)
                k = p->pi[i];
        }

        if (k > 0) {
            i = p->inv_pi[k];
            /* pi_i = k is the greatest positive element out of order */
            j = prm_exists(p, -1*(k+1));
            if (j > 0) {
                if (i < j) {
                    op_prefix_reversal(p, j);
                    print_verbose(p, verbose, 'r', 1, j, 0);
                    op_prefix_reversal(p, j-i);
                    print_verbose(p, verbose, 'r', 1, j-i, 0);
                } else {
                    op_prefix_reversal(p, i);
                    print_verbose(p, verbose, 'r', 1, i, 0);
                    op_prefix_reversal(p, i-j);
                    print_verbose(p, verbose, 'r', 1, i-j, 0);
                }
            } else {
                op_prefix_reversal(p, i);
                print_verbose(p, verbose, 'r', 1, i, 0);
                op_prefix_reversal(p, k);
                print_verbose(p, verbose, 'r', 1, k, 0);
            }
            info->nmoves += 2;
        } else {
            /* In this case, there is not a positive element out of order in \pi */
            i = 1;
            j = prm_exists(p, p->pi[i]+1);
            while (i < p->size && p->pi[i] < 0 && i+1 >= j) {
                i++;
                j = prm_exists(p, p->pi[i]+1);
            }

            if (i+1 < j) {
                op_prefix_reversal(p, i);
                print_verbose(p, verbose, 'r', 1, i, 0);
                op_prefix_reversal(p, j-1);
                print_verbose(p, verbose, 'r', 1, j-1, 0);
                info->nmoves += 2;
            } else {
                /* Then the permutation is equivalent to -I_n */
                spr_star_perm(p, verbose, info);
            }
        }
    }

    info->nof_pr = info->nmoves;
}/*}}}*/


