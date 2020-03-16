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


void psr_star_perm(permutation_t *p, int verbose, info_t *info, int nofbps) {/*{{{*/
    int i, l, *lengths, new_n;

    lengths = malloc(sizeof(int) * (nofbps + 3));
    lengths[0] = 0;
    lengths[nofbps + 2] = 0;
    new_n = 0;
    l = 2;

    if (nofbps % 2) {
        for (i = 1; i <= nofbps+1; i++) {
            for (; p->pi[l] == p->pi[l-1]-1; l++) ;
            lengths[i] = l-1 - new_n;
            new_n += lengths[i];
            l++;
        }

        for (i = 1; i <= nofbps+1; i += 2) {
            op_suffix_reversal(p, lengths[i]+1);
            print_verbose(p, verbose, 'r', lengths[i]+1, p->size, 0);
            op_prefix_reversal(p, p->size - lengths[i+1]);
            print_verbose(p, verbose, 'r', 1, p->size - lengths[i+1], 0);
            info->nmoves += 2;
            info->nof_sr++;
            info->nof_pr++;
        }
    } else {
        for (i = 1; i <= nofbps+1; i++) {
            for (; p->pi[l] == p->pi[l-1]+1; l++) ;
            lengths[i] = l-1 - new_n;
            new_n += lengths[i];
            l++;
        }

        for (i = nofbps+1; i >= 2; i -= 2) {
            op_prefix_reversal(p, p->size - lengths[i]);
            print_verbose(p, verbose, 'r', 1, p->size-lengths[i], 0);
            op_suffix_reversal(p, lengths[i-1]+1);
            print_verbose(p, verbose, 'r', lengths[i-1]+1, p->size, 0);
            info->nmoves += 2;
            info->nof_pr++;
            info->nof_sr++;
        }
        op_prefix_reversal(p, p->size - lengths[1]);
        print_verbose(p, verbose, 'r', 1, p->size-lengths[1], 0);
        info->nmoves++;
        info->nof_pr++;
    }

    free(lengths);
}/*}}}*/

int psr_edge_type1p(int *edge, permutation_t *p) {/*{{{*/
    /* Edge (\pi_i, \pi_j) of type 1 must have i == 1 and j < n+1 */
    int j1, j2;

    j1 = p->inv_pi[p->pi[1] - 1];
    j2 = p->inv_pi[p->pi[1] + 1];

    if (j1 > 2 && j1 <= p->size && is_UPSR_breakpoint(p, j1-1, j1)) {
        if (j2 > 2 && j2 <= p->size && is_UPSR_breakpoint(p, j2-1, j2)) {
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
    } else if (j2 > 2 && j2 <= p->size && is_UPSR_breakpoint(p, j2-1, j2)) {
        edge[0] = 1;
        edge[1] = j2;
        return 1;
    }

    return 0;
}/*}}}*/

int psr_edge_type2p(int *edge, permutation_t *p) {/*{{{*/
    /* Edge (\pi_i, \pi_j) of type 2 must have i != 0 and j <= n*/
    int i, j;
    for (j = p->size; j >= 2; j--) {
        if (!is_UPSR_breakpoint(p, j, j+1))
            continue;

        i = p->inv_pi[p->pi[j] - 1];
        if (i != 0 && j > i+1 && is_UPSR_breakpoint(p, i, i+1)) {
            edge[0] = i;
            edge[1] = j;
            return 1;
        }

        i = p->inv_pi[p->pi[j] + 1];
        if (i != 0 && j > i+1 && is_UPSR_breakpoint(p, i, i+1)) {
            edge[0] = i;
            edge[1] = j;
            return 1;
        }
    }
    return 0;
}/*}}}*/

int psr_edge_type3Ap(int *edge, permutation_t *p) {/*{{{*/
    /* Edge (\pi_i, \pi_j) of type 3A must have i == 1 and j < n+1 */
    int j1, j2;

    j1 = p->inv_pi[p->pi[1] - 1];
    j2 = p->inv_pi[p->pi[1] + 1];

    if (j1 > 2 && j1 <= p->size && is_UPSR_breakpoint(p, j1-1, j1)) {
        if (j2 > 2 && j2 <= p->size && is_UPSR_breakpoint(p, j2-1, j2)) {
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
    } else if (j2 > 2 && j2 <= p->size && is_UPSR_breakpoint(p, j2-1, j2)) {
        edge[0] = 1;
        edge[1] = j2;
        return 1;
    }

    return 0;
}/*}}}*/

int psr_edge_type3Bp(int *edge, permutation_t *p) {/*{{{*/
    /* Edge (\pi_i, \pi_j) of type 3B must have i > 1 and j < n+1 */
    int i, j;
    for (j = p->size; j >= 3; j--) {
        if (!is_UPSR_breakpoint(p, j-1, j))
            continue;

        i = p->inv_pi[p->pi[j] - 1];
        if (i > 1 && j > i+1 && is_UPSR_breakpoint(p, i, i+1)) {
            edge[0] = i;
            edge[1] = j;
            return 1;
        }

        i = p->inv_pi[p->pi[j] + 1];
        if (i > 1 && j > i+1 && is_UPSR_breakpoint(p, i, i+1)) {
            edge[0] = i;
            edge[1] = j;
            return 1;
        }
    }
    return 0;
}/*}}}*/

int psr_edge_type1s(int *edge, permutation_t *p) {/*{{{*/
    /* Edge (\pi_i, \pi_j) of type 1 must have j != n+1 and i > 0 */
    int i, j;
    for (i = 1; i <= p->size-1; i++) {
        if (!is_UPSR_breakpoint(p, i-1, i))
            continue;

        j = p->inv_pi[p->pi[i] + 1];
        if (j != p->size+1 && i+1 < j && is_UPSR_breakpoint(p, j-1, j)) {
            edge[0] = i;
            edge[1] = j;
            return 1;
        }

        j = p->inv_pi[p->pi[i] - 1];
        if (j != p->size+1 && i+1 < j && is_UPSR_breakpoint(p, j-1, j)) {
            edge[0] = i;
            edge[1] = j;
            return 1;
        }
    }
    return 0;
}/*}}}*/

int psr_edge_type2s(int *edge, permutation_t *p) {/*{{{*/
    /* Edge (\pi_i, \pi_j) of type 2 must have j == n and i > 0 */
    int i1, i2;
    
    i1 = p->inv_pi[p->pi[p->size] + 1];
    i2 = p->inv_pi[p->pi[p->size] - 1];

    if (i1 > 0 && i1 < p->size-1 && is_UPSR_breakpoint(p, i1, i1+1)) {
        if (i2 > 0 && i2 < p->size-1 && is_UPSR_breakpoint(p, i2, i2+1)) {
            if (i1 < i2) {
                edge[0] = i1;
                edge[1] = p->size;
                return 1;
            } else {
                edge[0] = i2;
                edge[1] = p->size;
                return 1;
            }
        } else {
            edge[0] = i1;
            edge[1] = p->size;
            return 1;
        }
    } else if (i2 > 0 && i2 < p->size-1 && is_UPSR_breakpoint(p, i2, i2+1)) {
        edge[0] = i2;
        edge[1] = p->size;
        return 1;
    }

    return 0;
}/*}}}*/

int psr_edge_type3As(int *edge, permutation_t *p) {/*{{{*/
    /* Edge (\pi_i, \pi_j) of type 3A must have j == n and i > 0 */
    int i1, i2;

    if (!is_UPSR_breakpoint(p, p->size-1, p->size))
        return 0;
    
    i1 = p->inv_pi[p->pi[p->size] + 1];
    i2 = p->inv_pi[p->pi[p->size] - 1];

    if (i1 > 0 && i1 < p->size-1 && is_UPSR_breakpoint(p, i1, i1+1)) {
        if (i2 > 0 && i2 < p->size-1 && is_UPSR_breakpoint(p, i2, i2+1)) {
            if (i1 < i2) {
                edge[0] = i1;
                edge[1] = p->size;
                return 1;
            } else {
                edge[0] = i2;
                edge[1] = p->size;
                return 1;
            }
        } else {
            edge[0] = i1;
            edge[1] = p->size;
            return 1;
        }
    } else if (i2 > 0 && i2 < p->size-1 && is_UPSR_breakpoint(p, i2, i2+1)) {
        edge[0] = i2;
        edge[1] = p->size;
        return 1;
    }

    return 0;
}/*}}}*/

int psr_edge_type3Bs(int *edge, permutation_t *p) {/*{{{*/
    /* Edge (\pi_i, \pi_j) of type 3B must have j < n and i > 0 */
    int i, j;
    for (i = 1; i <= p->size-2; i++) {
        if (!is_UPSR_breakpoint(p, i, i+1))
            continue;

        j = p->inv_pi[p->pi[i] + 1];
        if (j < p->size && j > i+1 && is_UPSR_breakpoint(p, j-1, j)) {
            edge[0] = i;
            edge[1] = j;
            return 1;
        }

        j = p->inv_pi[p->pi[i] - 1];
        if (j < p->size && j > i+1 && is_UPSR_breakpoint(p, j-1, j)) {
            edge[0] = i;
            edge[1] = j;
            return 1;
        }
    }
    return 0;
}/*}}}*/

void alg_2PSR(permutation_t *p, int verbose, info_t* info) {/*{{{*/
    int edge[2], i, nofbps;

    while (!is_identity(p)) {
        if (is_reverse(p)) {
            op_prefix_reversal(p, p->size);
            print_verbose(p, verbose, 'r', 1, p->size, 0);
            info->nmoves++;
            info->nof_pr++;
        } else if (psr_edge_type1p(edge, p)) {
            i = edge[1] - 1;
            op_prefix_reversal(p, i);
            print_verbose(p, verbose, 'r', 1, i, 0);
            info->nmoves++;
            info->nof_pr++;
        } else if (psr_edge_type2s(edge, p)) {
            i = edge[0] + 1;
	        op_suffix_reversal(p, i);
            print_verbose(p, verbose, 'r', i, p->size, 0);
            info->nmoves++;
            info->nof_sr++;
        } else if (psr_edge_type3Ap(edge, p)) {
            i = edge[1] - 1;
            op_prefix_reversal(p, i);
            print_verbose(p, verbose, 'r', 1, i, 0);
            info->nmoves++;
            info->nof_pr++;
        } else if (psr_edge_type3As(edge, p)) {
            i = edge[0] + 1;
	        op_suffix_reversal(p, i);
            print_verbose(p, verbose, 'r', i, p->size, 0);
            info->nmoves++;
            info->nof_sr++;
        } else if (psr_edge_type2p(edge, p)) {
            i = edge[1];
            op_prefix_reversal(p, i);
            print_verbose(p, verbose, 'r', 1, i, 0);
            i = edge[1] - edge[0];
            op_prefix_reversal(p, i);
            print_verbose(p, verbose, 'r', 1, i, 0);
            info->nmoves += 2;
            info->nof_pr += 2;
        } else if (psr_edge_type1s(edge, p)) {
            i = edge[0];
	        op_suffix_reversal(p, i);
            print_verbose(p, verbose, 'r', i, p->size, 0);
            i = p->size + 1 - (edge[1] - edge[0]);
	        op_suffix_reversal(p, i);
            print_verbose(p, verbose, 'r', i, p->size, 0);
		    info->nmoves += 2;
            info->nof_sr += 2;
        } else if (psr_edge_type3Bp(edge, p)) {
            i = edge[0];
            op_prefix_reversal(p, i);
            print_verbose(p, verbose, 'r', 1, i, 0);
            i = edge[1] - 1;
            op_prefix_reversal(p, i);
            print_verbose(p, verbose, 'r', 1, i, 0);
            info->nmoves += 2;
            info->nof_pr += 2;
        } else if (psr_edge_type3Bs(edge, p)) {
            i = edge[1];
	        op_suffix_reversal(p, i);
            print_verbose(p, verbose, 'r', i, p->size, 0);
            i = edge[0] + 1;
	        op_suffix_reversal(p, i);
            print_verbose(p, verbose, 'r', i, p->size, 0);
		    info->nmoves += 2;
            info->nof_sr += 2;
        } else {
            nofbps = nof_UPSR_breakpoints(p);
            if ((nofbps % 2 == 1 && p->inv_pi[1] > p->inv_pi[p->size]) ||
                        (nofbps % 2 == 0 && p->inv_pi[1] < p->inv_pi[p->size])) {
                op_prefix_reversal(p, p->size);
                print_verbose(p, verbose, 'r', 1, p->size, 0);
                info->nmoves++;
                info->nof_pr++;
            }
            psr_star_perm(p, verbose, info, nofbps);
        }
    }
}/*}}}*/




void spsr_star_perm(permutation_t *p, int verbose, info_t *info, int nofbps) {/*{{{*/
    int i, l, *lengths, new_n;

    lengths = malloc(sizeof(int) * (nofbps + 3));
    lengths[0] = 0;
    lengths[nofbps + 2] = 0;
    new_n = 0;
    l = 2;

    if (nofbps % 2) {
        for (i = 1; i <= nofbps+1; i++) {
            for (; p->pi[l] == p->pi[l-1]+1; l++) ;
            lengths[i] = l-1 - new_n;
            new_n += lengths[i];
            l++;
        }

        for (i = 1; i <= nofbps+1; i += 2) {
            op_suffix_reversal(p, lengths[i]+1);
            print_verbose(p, verbose, 'r', lengths[i]+1, p->size, 0);
            op_prefix_reversal(p, p->size - lengths[i+1]);
            print_verbose(p, verbose, 'r', 1, p->size - lengths[i+1], 0);
            info->nmoves += 2;
            info->nof_sr++;
            info->nof_pr++;
        }
    } else {
        for (i = 1; i <= nofbps+1; i++) {
            for (; p->pi[l] == p->pi[l-1]+1; l++) ;
            lengths[i] = l-1 - new_n;
            new_n += lengths[i];
            l++;
        }

        for (i = nofbps+1; i >= 2; i -= 2) {
            op_prefix_reversal(p, p->size - lengths[i]);
            print_verbose(p, verbose, 'r', 1, p->size-lengths[i], 0);
            op_suffix_reversal(p, lengths[i-1]+1);
            print_verbose(p, verbose, 'r', lengths[i-1]+1, p->size, 0);
            info->nmoves += 2;
            info->nof_pr++;
            info->nof_sr++;
        }
        op_prefix_reversal(p, p->size - lengths[1]);
        print_verbose(p, verbose, 'r', 1, p->size-lengths[1], 0);
        info->nmoves++;
        info->nof_pr++;
    }

    free(lengths);
}/*}}}*/

static int remove_1bp_2pr(permutation_t *p, int *iret, int *jret) {/*{{{*/
    int i, j;
    for (i = 1; i <= p->size-1; i++) {
        if (p->pi[i+1] - p->pi[i] == 1)
            continue;

        j = prm_exists(p, -p->pi[i]-1);
        if (j > 0 && i < j && j <= p->size) {
            (*iret) = j;
            (*jret) = j - i;
            return 1;
        }

        j = prm_exists(p, p->pi[i]+1);
        if (j > 0 && i+1 < j && j <= p->size) {
            (*iret) = i;
            (*jret) = j - 1;
            return 1;
        }
    }
    return 0;
}/*}}}*/

static int remove_1bp_2sr(permutation_t *p, int *iret, int *jret) {/*{{{*/
    int i, j;
    for (j = p->size; j >= 3; j--) {
        if (p->pi[j] - p->pi[j-1] == 1)
            continue;

        i = prm_exists(p, -p->pi[j]+1);
        if (i >= 1 && i < j) {
            (*iret) = i;
            (*jret) = p->size + 1 - (j - i);
            return 1;
        }

        i = prm_exists(p, p->pi[j]-1);
        if (i >= 1 && i+1 < j) {
            (*iret) = j;
            (*jret) = i + 1;
            return 1;
        }
    }
    return 0;
}/*}}}*/

void alg_2SPSR(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, nofbps;

    while (!is_identity(p)) {

        /* tries to remove 1 breakpoint with 1 prefix reversal */
        j = prm_exists(p, -p->pi[1]+1);
        if (j > 0 && j <= p->size) {
            op_prefix_reversal(p, j-1);
            print_verbose(p, verbose, 'r', 1, j-1, 0);
            info->nmoves++;
            info->nof_pr++;
            continue;
        }

        /* tries to remove 1 breakpoint with 1 suffix reversal */
        i = prm_exists(p, -p->pi[p->size]-1);
        if (i >= 1) {
            op_suffix_reversal(p, i+1);
            print_verbose(p, verbose, 'r', i+1, p->size, 0);
            info->nmoves++;
            info->nof_sr++;
            continue;
        }

        /* otherwise, tries to remove 1 breakpoint with 2 reversals */
        if (remove_1bp_2pr(p, &i, &j)) {
            op_prefix_reversal(p, i);
            print_verbose(p, verbose, 'r', 1, i, 0);
            op_prefix_reversal(p, j);
            print_verbose(p, verbose, 'r', 1, j, 0);
            info->nmoves += 2;
            info->nof_pr += 2;
            continue;
        }

        if (remove_1bp_2sr(p, &i, &j)) {
            op_suffix_reversal(p, i);
            print_verbose(p, verbose, 'r', i, p->size, 0);
            op_suffix_reversal(p, j);
            print_verbose(p, verbose, 'r', j, p->size, 0);
            info->nmoves += 2;
            info->nof_sr += 2;
            continue;
        }

        /* In this case, the permutation is equivalent to -I_n */
        if (is_signed_reverse(p)) {
            op_prefix_reversal(p, p->size);
            print_verbose(p, verbose, 'r', 1, p->size, 0);
            info->nmoves++;
            info->nof_pr++;
        } else {
            nofbps = nof_PS_breakpoints(p);
            if (nofbps % 2 == 1 && prm_exists(p, 1) && 
                            prm_exists(p, 1) > prm_exists(p, p->size)) {
                op_prefix_reversal(p, p->size);
                print_verbose(p, verbose, 'r', 1, p->size, 0);
                info->nmoves++;
                info->nof_pr++;
            } else if (nofbps % 2 == 0 && prm_exists(p, -1) &&
                            prm_exists(p, -1) < prm_exists(p, -p->size)) {
                op_prefix_reversal(p, p->size);
                print_verbose(p, verbose, 'r', 1, p->size, 0);
                info->nmoves++;
                info->nof_pr++;
            }
            spsr_star_perm(p, verbose, info, nofbps);
        }
    }
}/*}}}*/

