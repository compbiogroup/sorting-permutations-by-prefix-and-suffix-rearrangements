/***********************************************************
 * Created: Seg 17 Mar 2014 01:40:14 BRT
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
#include "wprt.h"
#include "wpr.h"
#include "wpt.h"



void alg_WPRTg(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, k, type, p1, p2;
    long double weight = 0;

    while (!is_identity(p)) {
        type = 0;
        weight = pot(p->size * p->size, info->alpha);

        k = 1;
        while (k <= p->size && !is_UPR_breakpoint(p, k, k+1))
            k++;

        if (p->pi[1] == 1) {
            op_prefix_transposition(p, k+1, p->size+1);
            info->weight += pot(p->size, info->alpha);
            print_verbose(p, verbose, 't', 1, k+1, p->size+1);
            info->nof_pt++;
            continue;
        }

        j = p->inv_pi[p->pi[1] + 1] + 1;
        if (j >= 0 && j < p->size+1) {
            i = p->inv_pi[p->pi[j] + 1] + 1;
            if (i > 1 && i < j && p->pi[i] != 1 && is_UPR_breakpoint(p, i-1, i) 
                    && is_UPR_breakpoint(p, j-1, j)) {
                weight = pot(j-1, info->alpha);
                p1 = i;
                p2 = j;
                type = 1;
            }
        }
        if (j > 0 && j <= p->size+1) {
            i = p->inv_pi[p->pi[j] - 1] + 1;
            if (i > 1 && i < j && p->pi[i] != 1 && is_UPR_breakpoint(p, i-1, i) 
                    && is_UPR_breakpoint(p, j-1, j) && pot(j-1, info->alpha) < weight) {
                weight = pot(j-1, info->alpha);
                p1 = i;
                p2 = j;
                type = 1;
            }
        }

        j = p->inv_pi[p->pi[1] - 1] + 1;
        if (j >= 0 && j < p->size+1) {
            i = p->inv_pi[p->pi[j] + 1] + 1;
            if (i > 1 && i < j && p->pi[i] != 1 && is_UPR_breakpoint(p, i-1, i) 
                    && is_UPR_breakpoint(p, j-1, j) && pot(j-1, info->alpha) < weight) {
                weight = pot(j-1, info->alpha);
                p1 = i;
                p2 = j;
                type = 1;
            }
        }
        if (j > 0 && j <= p->size+1) {
            i = p->inv_pi[p->pi[j] - 1] + 1;
            if (i > 1 && i < j && p->pi[i] != 1 && is_UPR_breakpoint(p, i-1, i) 
                    && is_UPR_breakpoint(p, j-1, j) && pot(j-1, info->alpha) < weight) {
                weight = pot(j-1, info->alpha);
                p1 = i;
                p2 = j;
                type = 1;
            }
        }

        /* remove one breakpoint with one prefix transposition */
        for (i = 2; i <= p->size; i++) {
            if (!is_UPR_breakpoint(p, i-1, i))
                continue;

            j = p->inv_pi[p->pi[i-1] + 1];
            if (j > i && j <= p->size+1 && is_UPR_breakpoint(p, j-1, j) &&
                        pot(j-1, info->alpha) < weight && !prt_prohibited(p, 't', i, j)) {
                weight = pot(j-1, info->alpha);
                p1 = i;
                p2 = j;
                type = 1;
            }

            j = p->inv_pi[p->pi[i-1] - 1];
            if (j > i && j <= p->size+1 && is_UPR_breakpoint(p, j-1, j) &&
                        pot(j-1, info->alpha) < weight && !prt_prohibited(p, 't', i, j)) {
                weight = pot(j-1, info->alpha);
                p1 = i;
                p2 = j;
                type = 1;
            }
        }

        j = p->inv_pi[p->pi[1] + 1] + 1;
        if (j > 2 && j <= p->size+1 && is_UPR_breakpoint(p, j-1, j) &&
                        pot(j-1, info->alpha) < weight) {
            for (i = 2; i < j; i++) {
                if (is_UPR_breakpoint(p, i-1, i) && !prt_prohibited(p, 't', i, j)) {
                    weight = pot(j-1, info->alpha);
                    p1 = i;
                    p2 = j;
                    type = 1;
                }
            }
        }

        j = p->inv_pi[p->pi[1] - 1] + 1;
        if (j > 2 && j <= p->size+1 && is_UPR_breakpoint(p, j-1, j) &&
                        pot(j-1, info->alpha) < weight) {
            for (i = 2; i < j; i++) {
                if (is_UPR_breakpoint(p, i-1, i) && !prt_prohibited(p, 't', i, j)) {
                    weight = pot(j-1, info->alpha);
                    p1 = i;
                    p2 = j;
                    type = 1;
                }
            }
        }

        /* remove one breakpoint with prefix reversal *//*{{{*/
        i = p->inv_pi[p->pi[1] + 1] - 1;
        if (2 <= i && i <= p->size && is_UPR_breakpoint(p, i, i+1) &&
                    pot(i, info->alpha) < weight && !prt_prohibited(p, 'r', i, 0)) {
            weight = pot(i, info->alpha);
            p1 = i;
            type = 2;
        }

        i = p->inv_pi[p->pi[1] - 1] - 1;
        if (2 <= i && i <= p->size && is_UPR_breakpoint(p, i, i+1) &&
                    pot(i, info->alpha) < weight && !prt_prohibited(p, 'r', i, 0)) {
            weight = pot(i, info->alpha);
            p1 = i;
            type = 2;
        }/*}}}*/

        if (type == 1) {
            op_prefix_transposition(p, p1, p2);
            print_verbose(p, verbose, 't', 1, p1, p2);
            info->weight += weight;
            info->nmoves++;
        } else if (type == 2) {
            op_prefix_reversal(p, p1);
            print_verbose(p, verbose, 'r', 1, p1, 0);
            info->weight += weight;
            info->nmoves++;
        }
    }
}/*}}}*/



long double WPRT(permutation_t *p, int end, int type, int verbose, info_t *info) {/*{{{*/
    int part, median;
    long double weight = 0;

    if (is_sorted_interval(p, 1, end, type))
        return 0;

    if (end <= 1)
        return 0;

    if (end == 2) {
        if ((p->pi[1] > p->pi[2] && type == INC) ||
                        (p->pi[1] < p->pi[2] && type == DEC)) {
            op_prefix_reversal(p, 2);
            weight += pot(2, info->alpha);
            print_verbose(p, verbose, 'r', 1, 2, 0);
            info->nmoves++;
        }
        return weight;
    }

    median = prm_median(p, 1, end);

    weight += partitionWPRT(p, end, 1-type, median, verbose, info);

    part = 1;
    if (type == INC) {
        while (p->pi[part] > median)
            part++;
    } else {
        while (p->pi[part] <= median)
            part++;
    }
    part--;

    weight += WPRT(p, part, 1-type, verbose, info);
    op_prefix_reversal(p, end);
    print_verbose(p, verbose, 'r', 1, end, 0);
    weight += pot(end, info->alpha);
    info->nmoves++;
    weight += WPRT(p, end - part, type, verbose, info);

    return weight;
}/*}}}*/

void alg_WPRT(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int k = p->size;

    while (k == p->pi[k] && k >= 1)
        k--;

    if (k > 0) {
        info->weight = WPRT(p, k, INC, verbose, info);
    }

    if (!is_identity(p))
        printf("ERROR! alg_WPRT()\n");
}/*}}}*/


long double partitionWPRT(permutation_t *p, int end, int type, int divisor, int verbose, info_t *info) {/*{{{*/
    int i, j, k, x, y, z, exists, half;
    long double weight = 0;

    if (end <= 1)
        return 0;

    if (end == 2) {
        if ((type == INC && p->pi[1] > divisor && p->pi[2] <= divisor) ||
                (type == DEC && p->pi[1] <= divisor && p->pi[2] > divisor)) {
            op_prefix_reversal(p, 2);
            print_verbose(p, verbose, 'r', 1, 2, 0);
            info->nmoves++;
            return pot(2, info->alpha);
        }
        return 0;
    }

    x = 0;
    y = end + 1;
    if (type == INC) {
        while (x < end && p->pi[x+1] <= divisor) x++;
        while (y > 1 && p->pi[y-1] > divisor) y--;
    } else {
        while (x < end && p->pi[x+1] > divisor) x++;
        while (y > 1 && p->pi[y-1] <= divisor) y--;
    }

    if (x == y-1)
        return 0;

    exists = 0;
    if (type == INC) {
        i = x+1;
        while (i < y && p->pi[i] > divisor) i++;
        z = i-1;
        while (i < y && p->pi[i] <= divisor) i++;
        if (i == y)
            exists = 1;
    } else {
        i = x+1;
        while (i < y && p->pi[i] <= divisor) i++;
        z = i-1;
        while (i < y && p->pi[i] > divisor) i++;
        if (i == y)
            exists = 1;
    }

    if (exists) {
        op_prefix_transposition(p, z+1, y);
        print_verbose(p, verbose, 't', 1, z+1, y);
        weight += pot(y-1, info->alpha);
        info->nmoves++;
        return weight;
    }

    end = y-1;

    half = end / 2.0;
    if (end % 2)
        half++;

    weight += partitionWPRT(p, half, type, divisor, verbose, info);
    op_prefix_transposition(p, half+1, end+1);
    print_verbose(p, verbose, 't', 1, half+1, end+1);
    weight += pot(end, info->alpha);
    info->nmoves++;
    weight += partitionWPRT(p, end - half, 1-type, divisor, verbose, info);

    i = 1;
    if (type == INC) {
        while (i <= end && p->pi[i] > divisor) i++;
        j = i;
        while (i <= end && p->pi[i] <= divisor) i++;
        k = i-1;
    } else {
        while (i <= end && p->pi[i] <= divisor) i++;
        j = i;
        while (i <= end && p->pi[i] > divisor) i++;
        k = i-1;
    }

    if (j <= k && j > 1) {
        op_prefix_reversal(p, k);
        print_verbose(p, verbose, 'r', 1, k, 0);
        weight += pot(k, info->alpha);
        info->nmoves++;
    }

    return weight;
}/*}}}*/

#if 0
long double partitionWPRT(permutation_t *p, int end, int type, int divisor, int verbose, info_t *info) {/*{{{*/
    int i, j, k, x, y, z, exists, half;
    long double weight = 0;

    if (end <= 1)
        return 0;

    if (end == 2) {
        if ((type == INC && p->pi[1] > divisor && p->pi[2] <= divisor) ||
                (type == DEC && p->pi[1] <= divisor && p->pi[2] > divisor)) {
            op_prefix_reversal(p, 2);
            print_verbose(p, verbose, 'r', 1, 2, 0);
            info->nmoves++;
            return pot(2, info->alpha);
        }
        return 0;
    }

    x = 0;
    y = end + 1;
    if (type == INC) {
        while (x < end && p->pi[x+1] <= divisor) x++;
        while (y > 1 && p->pi[y-1] > divisor) y--;
    } else {
        while (x < end && p->pi[x+1] > divisor) x++;
        while (y > 1 && p->pi[y-1] <= divisor) y--;
    }

    if (x == y-1)
        return 0;

    exists = 0;
    if (type == INC) {
        i = x+1;
        while (i < y && p->pi[i] > divisor) i++;
        z = i-1;
        while (i < y && p->pi[i] <= divisor) i++;
        if (i == y)
            exists = 1;
    } else {
        i = x+1;
        while (i < y && p->pi[i] <= divisor) i++;
        z = i-1;
        while (i < y && p->pi[i] > divisor) i++;
        if (i == y)
            exists = 1;
    }

    if (exists) {
        op_prefix_transposition(p, z+1, y);
        print_verbose(p, verbose, 't', 1, z+1, y);
        weight += pot(y-1, info->alpha);
        info->nmoves++;
        return weight;
    }

    end = y-1;

    half = end / 2.0;
    if (end % 2)
        half++;

    weight += partitionWPRT(p, half, type, divisor, verbose, info);
    op_prefix_transposition(p, half+1, end+1);
    print_verbose(p, verbose, 't', 1, half+1, end+1);
    weight += pot(end, info->alpha);
    info->nmoves++;
    weight += partitionWPRT(p, end - half, type, divisor, verbose, info);

    i = 1;
    if (type == INC) {
        while (i <= end && p->pi[i] <= divisor) i++;
        while (i <= end && p->pi[i] > divisor) i++;
        j = i-1;
        while (i <= end && p->pi[i] <= divisor) i++;
        k = i-1;
    } else {
        while (i <= end && p->pi[i] > divisor) i++;
        while (i <= end && p->pi[i] <= divisor) i++;
        j = i-1;
        while (i <= end && p->pi[i] > divisor) i++;
        k = i-1;
    }

    if (j < k && j > 1) {
        op_prefix_transposition(p, j+1, k+1);
        print_verbose(p, verbose, 't', 1, j+1, k+1);
        weight += pot(k, info->alpha);
        info->nmoves++;
    }

    return weight;
}/*}}}*/

long double partitionWPRT(permutation_t *p, int end, int type, int divisor, int verbose, info_t *info) {/*{{{*/
    int i, j, k, x, y, z, exists, half;
    long double weight = 0;

    if (end <= 1)
        return 0;

    if (end == 2) {
        if ((type == INC && p->pi[1] > divisor && p->pi[2] <= divisor) ||
                (type == DEC && p->pi[1] <= divisor && p->pi[2] > divisor)) {
            op_prefix_reversal(p, 2);
            print_verbose(p, verbose, 'r', 1, 2, 0);
            info->nmoves++;
            return pot(2, info->alpha);
        }
        return 0;
    }

    x = 0;
    y = end + 1;
    if (type == INC) {
        while (x < end && p->pi[x+1] <= divisor) x++;
        while (y > 1 && p->pi[y-1] > divisor) y--;
    } else {
        while (x < end && p->pi[x+1] > divisor) x++;
        while (y > 1 && p->pi[y-1] <= divisor) y--;
    }

    if (x == y-1)
        return 0;

    exists = 0;
    if (type == INC) {
        i = x+1;
        while (i < y && p->pi[i] > divisor) i++;
        z = i-1;
        while (i < y && p->pi[i] <= divisor) i++;
        if (i == y)
            exists = 1;
    } else {
        i = x+1;
        while (i < y && p->pi[i] <= divisor) i++;
        z = i-1;
        while (i < y && p->pi[i] > divisor) i++;
        if (i == y)
            exists = 1;
    }

    if (exists) {
        op_prefix_transposition(p, z+1, y);
        print_verbose(p, verbose, 't', 1, z+1, y);
        weight += pot(y-1, info->alpha);
        info->nmoves++;
        return weight;
    }

    end = y-1;

    half = end / 2.0;
    if (end % 2)
        half++;

    weight += partitionWPRT(p, half, 1-type, divisor, verbose, info);
    op_prefix_reversal(p, end);
    print_verbose(p, verbose, 'r', 1, end, 0);
    weight += pot(end, info->alpha);
    info->nmoves++;
    weight += partitionWPRT(p, end - half, 1-type, divisor, verbose, info);

    i = 1;
    if (type == INC) {
        while (i <= end && p->pi[i] > divisor) i++;
        j = i;
        while (i <= end && p->pi[i] <= divisor) i++;
        k = i-1;
    } else {
        while (i <= end && p->pi[i] <= divisor) i++;
        j = i;
        while (i <= end && p->pi[i] > divisor) i++;
        k = i-1;
    }

    if (j <= k && j > 1) {
        op_prefix_reversal(p, k);
        print_verbose(p, verbose, 'r', 1, k, 0);
        weight += pot(k, info->alpha);
        info->nmoves++;
    }

    return weight;
}/*}}}*/

long double partitionWPRT(permutation_t *p, int end, int type, int divisor, int verbose, info_t *info) {/*{{{*/
    int i, j, k, x, y, z, exists, half;
    long double weight = 0;

    if (end <= 1)
        return 0;

    if (end == 2) {
        if ((type == INC && p->pi[1] > divisor && p->pi[2] <= divisor) ||
                (type == DEC && p->pi[1] <= divisor && p->pi[2] > divisor)) {
            op_prefix_reversal(p, 2);
            print_verbose(p, verbose, 'r', 1, 2, 0);
            info->nmoves++;
            return pot(2, info->alpha);
        }
        return 0;
    }

    x = 0;
    y = end + 1;
    if (type == INC) {
        while (x < end && p->pi[x+1] <= divisor) x++;
        while (y > 1 && p->pi[y-1] > divisor) y--;
    } else {
        while (x < end && p->pi[x+1] > divisor) x++;
        while (y > 1 && p->pi[y-1] <= divisor) y--;
    }

    if (x == y-1)
        return 0;

    exists = 0;
    if (type == INC) {
        i = x+1;
        while (i < y && p->pi[i] > divisor) i++;
        z = i-1;
        while (i < y && p->pi[i] <= divisor) i++;
        if (i == y)
            exists = 1;
    } else {
        i = x+1;
        while (i < y && p->pi[i] <= divisor) i++;
        z = i-1;
        while (i < y && p->pi[i] > divisor) i++;
        if (i == y)
            exists = 1;
    }

    if (exists) {
        op_prefix_transposition(p, z+1, y);
        print_verbose(p, verbose, 't', 1, z+1, y);
        weight += pot(y-1, info->alpha);
        info->nmoves++;
        return weight;
    }

    end = y-1;

    half = end / 2.0;
    if (end % 2)
        half++;

    weight += partitionWPRT(p, half, 1-type, divisor, verbose, info);
    op_prefix_reversal(p, end);
    print_verbose(p, verbose, 'r', 1, end, 0);
    weight += pot(end, info->alpha);
    info->nmoves++;
    weight += partitionWPRT(p, end - half, type, divisor, verbose, info);

    i = 1;
    if (type == INC) {
        while (i <= end && p->pi[i] <= divisor) i++;
        while (i <= end && p->pi[i] > divisor) i++;
        j = i-1;
        while (i <= end && p->pi[i] <= divisor) i++;
        k = i-1;
    } else {
        while (i <= end && p->pi[i] > divisor) i++;
        while (i <= end && p->pi[i] <= divisor) i++;
        j = i-1;
        while (i <= end && p->pi[i] > divisor) i++;
        k = i-1;
    }

    if (j < k && j > 1) {
        op_prefix_transposition(p, j+1, k+1);
        print_verbose(p, verbose, 't', 1, j+1, k+1);
        weight += pot(k, info->alpha);
        info->nmoves++;
    }

    return weight;
}/*}}}*/
#endif



void alg_WSPRTg(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, type, p1, p2;
    long double weight;

    while (!is_identity(p)) {
        type = 0;
        weight = pot(p->size * p->size, info->alpha);

        if (p->pi[1] == 1) {
            i = 1;
            while (p->pi[i+1] - p->pi[i] == 1)
                i++;
            op_prefix_transposition(p, i+1, p->size+1);
            print_verbose(p, verbose, 't', 1, i+1, p->size+1);
            info->weight += pot(p->size, info->alpha);
            continue;
        }

        j = prm_exists(p, p->pi[1] - 1);
        if (j > 0) {
            j++;
            i = prm_exists(p, p->pi[j] - 1);
            if (i > 0) {
                i++;
                if (i < j && p->pi[i] != 1) {
                    weight = pot(j-1, info->alpha);
                    type = 1;
                    p1 = i;
                    p2 = j;
                }
            }
        }

        i = prm_exists(p, -p->pi[1]+1);
        if (i > 0 && is_P_breakpoint(p, i-1, i) &&
                     pot(i-1, info->alpha) < weight && !sprt_prohibited(p, 'r', i-1, 0)) {
            weight = pot(i-1, info->alpha);
            type = 2;
            p1 = i-1;
        }

        for (i = 2; i <= p->size; i++) {
            if (!is_P_breakpoint(p, i-1, i))
                continue;

            j = prm_exists(p, p->pi[i-1]+1);
            if (i < j && is_P_breakpoint(p, j-1, j) &&
                    pot(j-1, info->alpha) < weight && !sprt_prohibited(p, 't', i, j)) {
                weight = pot(j-1, info->alpha);
                type = 1;
                p1 = i;
                p2 = j;
            }
        }

        j = prm_exists(p, p->pi[1]-1);
        if (j > 0 && is_P_breakpoint(p, j, j+1) && pot(j-1, info->alpha) < weight) {
            j++;
            for (i = 2; i < j; i++) {
                if (is_P_breakpoint(p, i-1, i) && !sprt_prohibited(p, 't', i, j)) {
                    weight = pot(j-1, info->alpha);
                    type = 1;
                    p1 = i;
                    p2 = j;
                }
            }
        }

        if (type == 1) {
            op_prefix_transposition(p, p1, p2);
            print_verbose(p, verbose, 't', 1, p1, p2);
            info->weight += weight;
            info->nmoves++;
        } else if (type == 2) {
            op_prefix_reversal(p, p1);
            print_verbose(p, verbose, 'r', 1, p1, 0);
            info->weight += weight;
            info->nmoves++;
        }
    }
}/*}}}*/



long double partitionWSPRT(permutation_t *p, int end, int type, int divisor, int verbose, info_t *info) {/*{{{*/
    int i, j, k, x, y, z, exists, half;
    long double weight = 0;

    if (end <= 1)
        return 0;

    if (end == 2) {
        if ((type == INC && p->pi[1] > divisor && p->pi[2] <= divisor) ||
                (type == DEC && p->pi[1] <= divisor && p->pi[2] > divisor)) {
            op_prefix_reversal(p, 2);
            print_verbose(p, verbose, 'r', 1, 2, 0);
            info->nmoves++;
            return pot(2, info->alpha);
        }
        return 0;
    }

    x = 0;
    y = end + 1;
    if (type == INC) {
        while (x < end && p->pi[x+1] <= divisor) x++;
        while (y > 1 && p->pi[y-1] > divisor) y--;
    } else {
        while (x < end && p->pi[x+1] > divisor) x++;
        while (y > 1 && p->pi[y-1] <= divisor) y--;
    }

    if (x == y-1)
        return 0;

    exists = 0;
    if (type == INC) {
        i = x+1;
        while (i < y && p->pi[i] > divisor) i++;
        z = i-1;
        while (i < y && p->pi[i] <= divisor) i++;
        if (i == y)
            exists = 1;
    } else {
        i = x+1;
        while (i < y && p->pi[i] <= divisor) i++;
        z = i-1;
        while (i < y && p->pi[i] > divisor) i++;
        if (i == y)
            exists = 1;
    }

    if (exists) {
        op_prefix_transposition(p, z+1, y);
        print_verbose(p, verbose, 't', 1, z+1, y);
        weight += pot(y-1, info->alpha);
        info->nmoves++;
        return weight;
    }

    end = y-1;

    half = end / 2.0;
    if (end % 2)
        half++;
    if (half % 2)
        half++;

    weight += partitionWSPRT(p, half, type, divisor, verbose, info);
    op_prefix_transposition(p, half+1, end+1);
    print_verbose(p, verbose, 't', 1, half+1, end+1);
    weight += pot(end, info->alpha);
    info->nmoves++;
    weight += partitionWSPRT(p, end - half, 1-type, divisor, verbose, info);

    i = 1;
    if (type == INC) {
        while (i <= end && p->pi[i] > divisor) i++;
        j = i;
        while (i <= end && p->pi[i] <= divisor) i++;
        k = i-1;
    } else {
        while (i <= end && p->pi[i] <= divisor) i++;
        j = i;
        while (i <= end && p->pi[i] > divisor) i++;
        k = i-1;
    }

    if (j <= k && j > 1) {
        op_prefix_reversal(p, k);
        print_verbose(p, verbose, 'r', 1, k, 0);
        weight += pot(k, info->alpha);
        info->nmoves++;
    }

    return weight;
}/*}}}*/

long double WSPRT(permutation_t *p, int end, int type, int verbose, info_t *info) {/*{{{*/
    int part, median;
    long double weight = 0;

    if (is_sorted_interval(p, 1, end, type))
        return 0;

    if (end <= 1)
        return 0;

    if (end == 2) {
        if ((p->pi[1] > p->pi[2] && type == INC) ||
                        (p->pi[1] < p->pi[2] && type == DEC)) {
            op_prefix_reversal(p, 2);
            print_verbose(p, verbose, 'r', 1, 2, 0);
            weight += pot(2, info->alpha);
            info->nmoves++;
        }
        return weight;
    }

    median = prm_median(p, 1, end);
    if (median % 2)
        median++;

    weight += partitionWSPRT(p, end, 1-type, median, verbose, info);

    part = 1;
    if (type == INC) {
        while (p->pi[part] > median)
            part++;
    } else {
        while (p->pi[part] <= median)
            part++;
    }
    part--;

    weight += WSPRT(p, part, 1-type, verbose, info);
    op_prefix_reversal(p, end);
    print_verbose(p, verbose, 'r', 1, end, 0);
    weight += pot(end, info->alpha);
    info->nmoves++;
    weight += WSPRT(p, end - part, type, verbose, info);

    return weight;
}/*}}}*/

void alg_WSPRT(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int k;
    permutation_t up;

    create_permutation(&up, 2*p->size, UNSIGNED);
    prm_image(p, &up);

    k = up.size;
    while (k == up.pi[k] && k >= 1)
        k--;

    if (k > 0) {
        info->weight = WSPRT(&up, k, INC, verbose, info) / 2;
        prm_back_from_image(p, &up);
    }

    destroy_permutation(&up);

    if (!is_identity(p))
        printf("ERROR! alg_WSPRT()\n");
}/*}}}*/



