/***********************************************************
 * Created: Seg 17 Mar 2014 01:40:14 BRT
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "../util.h"
#include "../permutations.h"
#include "../rearrangements.h"
#include "../breakpoints.h"

int calls = 0;
int total_calls = 0;


void wpr_star_perm(permutation_t *p, int verbose, info_t *info) {/*{{{*/
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
        info->weight += pot(new_n, info->alpha) + pot(new_n - lengths[i], info->alpha);
    }

    free(lengths);
}/*}}}*/

void alg_WPRg(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, type, p1, p2;
    long double weight;
    int k = p->size;

    while (k == p->pi[k] && k >= 1)
        k--;

    while (!is_identity(p)) {
        type = 0;
        weight = DBL_MAX;

        /* Calcula o menor weight de remover um breakpoint com uma operação */
        i = p->inv_pi[p->pi[1]+1] - 1;
        if (i > 1 && is_UPR_breakpoint(p, i, i+1)) {
            weight = pot(i, info->alpha);
            type = 1;
            p1 = i;
        }

        i = p->inv_pi[p->pi[1]-1] - 1;
        if (i > 1 && is_UPR_breakpoint(p, i, i+1) && pot(i, info->alpha) < weight) {
            weight = pot(i, info->alpha);
            type = 1;
            p1 = i;
        }

        /* Calcula o menor weight de remover um breakpoint com duas operações */
        for (i = 1; i <= p->size-1; i++) {
            if (!is_UPR_breakpoint(p, i, i+1))
                continue;

            j = p->inv_pi[p->pi[i] + 1];
            if (j > i+1 && j <= p->size && is_UPR_breakpoint(p, j, j+1) &&
                            (pot(j, info->alpha) + pot(j-i, info->alpha)) < weight) {
                weight = pot(j, info->alpha) + pot(j-i, info->alpha);
                type = 2;
                p1 = j;
                p2 = j - i;
            }

            j = p->inv_pi[p->pi[i] - 1];
            if (j > i+1 && j <= p->size && is_UPR_breakpoint(p, j, j+1) && 
                            (pot(j, info->alpha) + pot(j-i, info->alpha)) < weight) {
                weight = pot(j, info->alpha) + pot(j-i, info->alpha);
                type = 2;
                p1 = j;
                p2 = j - i;
            }
        }

        for (i = 2; i <= p->size; i++) {
            if (!is_UPR_breakpoint(p, i, i+1))
                continue;

            j = p->inv_pi[p->pi[i] + 1];
            if (j > i+1 && is_UPR_breakpoint(p, j-1, j) && 
                            (pot(i, info->alpha) + pot(j-1, info->alpha)) < weight) {
                weight = pot(i, info->alpha) + pot(j-1, info->alpha);
                type = 2;
                p1 = i;
                p2 = j - 1;
            }

            j = p->inv_pi[p->pi[i] - 1];
            if (j > i+1 && is_UPR_breakpoint(p, j-1, j) && 
                            (pot(i, info->alpha) + pot(j-1, info->alpha)) < weight) {
                weight = pot(i, info->alpha) + pot(j-1, info->alpha);
                type = 2;
                p1 = i;
                p2 = j - 1;
            }
        }

        if (type == 1) {
            op_prefix_reversal(p, p1);
            print_verbose(p, verbose, 'r', 1, p1, 0);
            info->nmoves++;
            info->weight += weight;
        } else if (type == 2) {
            op_prefix_reversal(p, p1);
            print_verbose(p, verbose, 'r', 1, p1, 0);
            op_prefix_reversal(p, p2);
            print_verbose(p, verbose, 'r', 1, p2, 0);
            info->nmoves += 2;
            info->weight += weight;
        } else {
            wpr_star_perm(p, verbose, info);
        }
    }

    info->nof_pr = info->nmoves;
}/*}}}*/


static long double WPRm(permutation_t *p, int elem, int pos, int verbose, info_t *info) {/*{{{*/
    int elem2;
    long double weight = 0;
    int begin, end;

    if (elem < 1 || pos < 1 || elem > p->size || pos > p->size) {
        /* printf("Consec calls: %d\n", calls); */
        calls = 0;
        return weight;
    }

    find_extremes_of_strip(p, elem, &begin, &end);

    if (end > pos) {
        /* printf("Consec calls: %d\n", calls); */
        calls = 0;
        return weight;
    }

    calls++;
    total_calls++;

    elem2 = p->pi[end];
    if (p->pi[begin] < p->pi[end])
        elem2 = p->pi[begin];

    weight += WPRm(p, elem2-1, begin-1, verbose, info);

    /* bring the strip with elem to the beginning if necessary */
    if (begin != 1) {
        op_prefix_reversal(p, end);
        weight += pot(end, info->alpha);
        print_verbose(p, verbose, 'r', 1, end, 0);
        info->nmoves++;
        end = end - begin + 1;
        begin = 1;
    }

    /* reverts first strip if necessary */
    if (p->pi[end] > p->pi[begin]) {
        op_prefix_reversal(p, end);
        weight += pot(end, info->alpha);
        print_verbose(p, verbose, 'r', 1, end, 0);
        info->nmoves++;
    }

    /* move first strip to pos */
    if (pos > 1) {
        op_prefix_reversal(p, pos);
        weight += pot(pos, info->alpha);
        print_verbose(p, verbose, 'r', 1, pos, 0);
        info->nmoves++;
    }

    return weight;
}/*}}}*/

void alg_WPRm(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int out_of_order = p->size;
    int k = p->size;

    while (k == p->pi[k] && k >= 1)
        k--;

    calls = 0;
    total_calls = 0;

    while (!is_identity(p)) {
        while (p->pi[out_of_order] == out_of_order)
            out_of_order--;
        info->weight += WPRm(p, out_of_order, out_of_order, verbose, info);
    }

    if (verbose == 3 || verbose == 4) {
        print_permutation(p);
    }

    /* printf("Total Calls: %d\n", total_calls); */

    info->nof_pr = info->nmoves;
}/*}}}*/


long double partitionWPR(permutation_t *p, int end, int type, int divisor, int verbose, info_t *info) {/*{{{*/
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
        if (x > 0) {
            op_prefix_reversal(p, z);
            print_verbose(p, verbose, 'r', 1, z, 0);
            weight += pot(z, info->alpha);
            info->nmoves++;
        }
        op_prefix_reversal(p, y-1);
        print_verbose(p, verbose, 'r', 1, y-1, 0);
        weight += pot(y-1, info->alpha);
        info->nmoves++;
        return weight;
    }

    end = y-1;

    half = end / 2.0;
    if (end % 2)
        half++;

    weight += partitionWPR(p, half, 1-type, divisor, verbose, info);
    op_prefix_reversal(p, end);
    print_verbose(p, verbose, 'r', 1, end, 0);
    weight += pot(end, info->alpha);
    info->nmoves++;
    weight += partitionWPR(p, end - half, 1-type, divisor, verbose, info);

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

long double WPR(permutation_t *p, int end, int type, int verbose, info_t *info) {/*{{{*/
    int part, median;
    long double weight = 0;

    if (is_sorted_interval(p, 1, end, type))
        return 0;

    if (end <= 1)
        return 0;

    if (end == 2) {
        if ((p->pi[1] > p->pi[2] && type == INC) || (p->pi[1] < p->pi[2] && type == DEC)) {
            op_prefix_reversal(p, 2);
            weight += pot(2, info->alpha);
            print_verbose(p, verbose, 'r', 1, 2, 0);
            info->nmoves++;
        }
        return weight;
    }

    median = prm_median(p, 1, end);

    weight += partitionWPR(p, end, 1-type, median, verbose, info);

    part = 1;
    if (type == INC) {
        while (p->pi[part] > median)
            part++;
    } else {
        while (p->pi[part] <= median)
            part++;
    }
    part--;

    weight += WPR(p, part, 1-type, verbose, info);
    op_prefix_reversal(p, end);
    print_verbose(p, verbose, 'r', 1, end, 0);
    weight += pot(end, info->alpha);
    info->nmoves++;
    weight += WPR(p, end - part, type, verbose, info);

    return weight;
}/*}}}*/

void alg_WPR(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int k = p->size;

    while (k == p->pi[k] && k >= 1)
        k--;

    if (k > 0) {
        info->weight = WPR(p, k, INC, verbose, info);
        info->nof_pr = info->nmoves;
    }

    if (!is_identity(p)) 
        printf("ERROR! alg_WPR()\n");
}/*}}}*/





void wspr_star_perm(permutation_t *p, int verbose, info_t *info) {/*{{{*/
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
        info->weight += new_n + new_n - lengths[i];
    }

    free(lengths);
}/*}}}*/

void alg_WSPRg(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, type, p1, p2;
    long double weight;

    while (!is_identity(p)) {
        type = 0;
        weight = p->size * p->size;

        /* Calcula o menor weight de remover um breakpoint com uma operação */
        i = prm_exists(p, -p->pi[1]+1);
        if (i >= 2) {
            weight = i-1;
            p1 = i-1;
            type = 1;
        }

        /* Calcula o menor weight de remover um breakpoint com duas operações */
        for (i = 1; i <= p->size; i++) {
            if (!is_P_breakpoint(p, i, i+1))
                continue;

            if (i < p->size) {
                j = prm_exists(p, -p->pi[i]-1);
                if (i < j && j <= p->size && j+j-i < weight) {
                    weight = j+j-i;
                    type = 2;
                    p1 = j;
                    p2 = j-i;
                }
            }

            j = prm_exists(p, p->pi[i]+1);
            if (i+1 < j && i+j-1 < weight) {
                weight = i+j-1;
                type = 2;
                p1 = i;
                p2 = j-1;
            }
        }

        if (type == 1) {
            op_prefix_reversal(p, weight);
            print_verbose(p, verbose, 'r', 1, weight, 0);
            info->nmoves++;
            info->weight += weight;
        } else if (type == 2) {
            op_prefix_reversal(p, p1);
            print_verbose(p, verbose, 'r', 1, p1, 0);
            op_prefix_reversal(p, p2);
            print_verbose(p, verbose, 'r', 1, p2, 0);
            info->nmoves += 2;
            info->weight += weight;
        } else {
            wspr_star_perm(p, verbose, info);
        }
    }

    info->nof_pr = info->nmoves;
}/*}}}*/



long double partitionWSPR(permutation_t *p, permutation_t *sp, int end, int type, int divisor, int verbose, info_t *info) {/*{{{*/
    int i, j, k, x, y, z, exists, half;
    long double weight = 0;
        
    if (end <= 1)
        return 0;

    if (end == 2) {
        if ((type == INC && p->pi[1] > divisor && p->pi[2] <= divisor) ||
                (type == DEC && p->pi[1] <= divisor && p->pi[2] > divisor)) {
            op_prefix_reversal(p, 2);
            op_prefix_reversal(sp, 1);
            info->weight += pot(1, info->alpha);
            print_verbose(sp, verbose, 'r', 1, 1, 0);
            info->nmoves++;
            return weight;
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
        if (x > 0) {
            op_prefix_reversal(p, z);
            op_prefix_reversal(sp, z/2);
            info->weight += pot(z/2, info->alpha);
            print_verbose(sp, verbose, 'r', 1, z/2, 0);
            info->nmoves++;
        }
        op_prefix_reversal(p, y-1);
        op_prefix_reversal(sp, (y-1)/2);
        info->weight += pot((y-1)/2, info->alpha);
        print_verbose(sp, verbose, 'r', 1, (y-1)/2, 0);
        info->nmoves++;
        return weight;
    }

    end = y-1;

    half = end / 2.0;
    if (end % 2)
        half++;
    if (half % 2)
        half++;

    info->weight += partitionWSPR(p, sp, half, 1-type, divisor, verbose, info);
    op_prefix_reversal(p, end);
    op_prefix_reversal(sp, end/2);
    info->weight += pot(end/2, info->alpha);
    print_verbose(sp, verbose, 'r', 1, end/2, 0);
    info->nmoves++;
    info->weight += partitionWSPR(p, sp, end - half, 1-type, divisor, verbose, info);

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
        op_prefix_reversal(sp, k/2);
        info->weight += pot(k/2, info->alpha);
        print_verbose(sp, verbose, 'r', 1, k/2, verbose);
        info->nmoves++;
    }

    return weight;
}/*}}}*/

long double WSPR(permutation_t *p, permutation_t *sp, int end, int type, int verbose, info_t *info) {/*{{{*/
    int part, median;
    long double weight = 0;

    if (is_sorted_interval(p, 1, end, type))
        return 0;

    if (end <= 1)
        return 0;

    if (end == 2) {
        if ((p->pi[1] > p->pi[2] && type == INC) || (p->pi[1] < p->pi[2] && type == DEC)) {
            op_prefix_reversal(p, 2);
            op_prefix_reversal(sp, 1);
            info->weight += pot(1, info->alpha);
            print_verbose(sp, verbose, 'r', 1, 1, 0);
            info->nmoves++;
        }
        return weight;
    }

    median = prm_median(p, 1, end);
    if (median % 2)
        median++;

    info->weight += partitionWSPR(p, sp, end, 1-type, median, verbose, info);

    part = 1;
    if (type == INC) {
        while (p->pi[part] > median)
            part++;
    } else {
        while (p->pi[part] <= median)
            part++;
    }
    part--;

    info->weight += WSPR(p, sp, part, 1-type, verbose, info);
    op_prefix_reversal(p, end);
    op_prefix_reversal(sp, end/2);
    info->weight += pot(end/2, info->alpha);
    print_verbose(sp, verbose, 'r', 1, end/2, 0);
    info->nmoves++;
    info->weight += WSPR(p, sp, end - part, type, verbose, info);

    return weight;
}/*}}}*/

void alg_WSPR(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int k;
    permutation_t up;
    long double a;

    create_permutation(&up, 2*p->size, UNSIGNED);
    prm_image(p, &up);

    k = up.size;
    while (k == up.pi[k] && k >= 1)
        k--;

    if (k > 0) {
        a = WSPR(&up, p, k, INC, verbose, info);
        info->nof_pr = info->nmoves;
    }

    if (!is_identity(p))
        printf("ERROR! alg_WSPR()\n");

    destroy_permutation(&up);
}/*}}}*/

