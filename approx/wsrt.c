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


long double partitionWSRT(permutation_t *p, int ini, int type, int divisor, int verbose, info_t *info) {/*{{{*/
    int i, x, y, z, exists, half, n = p->size, j, k;
    long double weight = 0;

    if (ini >= n)
        return 0;

    if (ini == n-1) {
        if ((type == INC && p->pi[n-1] > divisor && p->pi[n] <= divisor) ||
                (type == DEC && p->pi[n-1] <= divisor && p->pi[n] > divisor)) {
            op_suffix_reversal(p, n-1);
            print_verbose(p, verbose, 'r', n-1, p->size+1, 0);
            info->nmoves++;
            return pot(2, info->alpha);
        }
        return 0;
    }

    x = ini-1;
    y = n+1;
    if (type == INC) {
        while (x < n && p->pi[x+1] <= divisor) x++;
        while (y > ini && p->pi[y-1] > divisor) y--;
    } else {
        while (x < n && p->pi[x+1] > divisor) x++;
        while (y > ini && p->pi[y-1] <= divisor) y--;
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
        op_suffix_transposition(p, x+1, z+1);
        weight += pot(n-x, info->alpha);
        info->nmoves++;
        print_verbose(p, verbose, 't', x+1, z+1, p->size+1);
        return weight;
    }

    ini = x+1;

    half = (n-ini+1) / 2.0;
    if ((n-ini+1) % 2)
        half++;
    half += ini-1;

    weight += partitionWSRT(p, half+1, type, divisor, verbose, info);
    op_suffix_transposition(p, ini, half+1);
    print_verbose(p, verbose, 't', ini, half+1, p->size+1);
    weight += pot(n-ini+1, info->alpha);
    info->nmoves++;
    weight += partitionWSRT(p, n-half+ini, 1-type, divisor, verbose, info);

    i = ini;
    if (type == INC) {
        while (i <= n && p->pi[i] <= divisor) i++;
        j = i;
        while (i <= n && p->pi[i] > divisor) i++;
        k = i-1;
    } else {
        while (i <= n && p->pi[i] > divisor) i++;
        j = i;
        while (i <= n && p->pi[i] <= divisor) i++;
        k = i-1;
    }

    if (j <= k && k < n) {
        op_suffix_reversal(p, j);
        weight += pot(n-j+1, info->alpha);
        info->nmoves++;
        print_verbose(p, verbose, 'r', j, p->size, 0);
    }

    return weight;
}/*}}}*/

long double WSRT(permutation_t *p, int ini, int type, int verbose, info_t *info) {/*{{{*/
    int part, median, n = p->size;
    long double weight = 0;

    if (is_sorted_interval(p, ini, n, type))
        return 0;

    if (ini >= n)
        return 0;

    if (ini == n-1) {
        if ((p->pi[n-1] > p->pi[n] && type == INC) ||
                        (p->pi[n-1] < p->pi[n] && type == DEC)) {
            op_suffix_reversal(p, n-1);
            weight += pot(2, info->alpha);
            print_verbose(p, verbose, 'r', n-1, p->size, 0);
            info->nmoves++;
        }
        return weight;
    }

    median = prm_median(p, ini, n);

    weight += partitionWSRT(p, ini, 1-type, median, verbose, info);

    part = ini;
    if (type == INC) {
        while (p->pi[part] > median)
            part++;
    } else {
        while (p->pi[part] <= median)
            part++;
    }
    part--;

    weight += WSRT(p, part+1, 1-type, verbose, info);
    op_suffix_reversal(p, ini);
    print_verbose(p, verbose, 'r', ini, p->size, 0);
    weight += pot(n-ini+1, info->alpha);
    info->nmoves++;
    weight += WSRT(p, n - part + ini, type, verbose, info);

    return weight;
}/*}}}*/

void alg_WSRT(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int k = 1;

    while (k == p->pi[k] && k <= p->size)
        k++;

    if (k <= p->size) {
        info->weight = WSRT(p, k, INC, verbose, info);
    }

    if (!is_identity(p))
        printf("ERROR! alg_SRT()\n");
}/*}}}*/




long double partitionWSSRT(permutation_t *p, int ini, int type, int divisor, int verbose, info_t *info) {/*{{{*/
    int i, j, k, x, y, z, exists, half, n = p->size;
    long double weight = 0;

    if (ini >= n)
        return 0;

    if (ini == n-1) {
        if ((type == INC && p->pi[n-1] > divisor && p->pi[n] <= divisor) ||
                (type == DEC && p->pi[n-1] <= divisor && p->pi[n] > divisor)) {
            op_suffix_reversal(p, n-1);
            info->nmoves++;
            return pot(2, info->alpha);
        }
        return 0;
    }

    x = ini-1;
    y = n+1;
    if (type == INC) {
        while (x < n && p->pi[x+1] <= divisor) x++;
        while (y > ini && p->pi[y-1] > divisor) y--;
    } else {
        while (x < n && p->pi[x+1] > divisor) x++;
        while (y > ini && p->pi[y-1] <= divisor) y--;
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
        op_suffix_transposition(p, x+1, z+1);
        weight += pot(n-x, info->alpha);
        info->nmoves++;
        return weight;
    }

    ini = x+1;

    half = (n-ini+1) / 2.0;
    if ((n-ini+1) % 2)
        half++;
    half += ini-1;
    if (half % 2)
        half++;

    weight += partitionWSSRT(p, half+1, type, divisor, verbose, info);
    op_suffix_transposition(p, ini, half+1);
    weight += pot(n-ini+1, info->alpha);
    info->nmoves++;
    weight += partitionWSSRT(p, n-half+ini, 1-type, divisor, verbose, info);

    i = ini;
    if (type == INC) {
        while (i <= n && p->pi[i] <= divisor) i++;
        j = i;
        while (i <= n && p->pi[i] > divisor) i++;
        k = i-1;
    } else {
        while (i <= n && p->pi[i] > divisor) i++;
        j = i;
        while (i <= n && p->pi[i] <= divisor) i++;
        k = i-1;
    }

    if (j <= k && k < n) {
        op_suffix_reversal(p, j);
        weight += pot(n-j+1, info->alpha);
        info->nmoves++;
    }

    return weight;
}/*}}}*/

long double WSSRT(permutation_t *p, int ini, int type, int verbose, info_t *info) {/*{{{*/
    int part, median, n = p->size;
    long double weight = 0;

    n = p->size;

    if (is_sorted_interval(p, ini, n, type))
        return 0;

    if (ini >= n)
        return 0;

    if (ini == n-1) {
        if ((p->pi[n-1] > p->pi[n] && type == INC) ||
                        (p->pi[n-1] < p->pi[n] && type == DEC)) {
            op_suffix_reversal(p, n-1);
            weight += pot(2, info->alpha);
            info->nmoves++;
        }
        return weight;
    }

    median = prm_median(p, ini, n);
    if (median % 2)
        median++;

    weight += partitionWSSRT(p, ini, 1-type, median, verbose, info);

    part = ini;
    if (type == INC) {
        while (p->pi[part] > median)
            part++;
    } else {
        while (p->pi[part] <= median)
            part++;
    }
    part--;

    weight += WSSRT(p, part+1, 1-type, verbose, info);
    op_suffix_reversal(p, ini);
    weight += pot(n-ini+1, info->alpha);
    info->nmoves++;
    weight += WSSRT(p, n - part + ini, type, verbose, info);

    return weight;
}/*}}}*/

void alg_WSSRT(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int k;
    permutation_t up;

    create_permutation(&up, 2*p->size, UNSIGNED);
    prm_image(p, &up);

    k = 1;
    while (k == up.pi[k] && k <= up.size)
        k++;

    if (k <= up.size) {
        info->weight = WSSRT(&up, k, INC, verbose, info) / 2;
        prm_back_from_image(p, &up);
    }

    destroy_permutation(&up);

    if (!is_identity(p))
       printf("ERROR! alg_WSSRT()\n");
}/*}}}*/



