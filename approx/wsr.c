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

long double partitionWSR(permutation_t *p, int ini, int type, int divisor, int verbose, info_t *info) {/*{{{*/
    int i, j, k, x, y, z, exists, half, n = p->size;
    long double weight = 0;
        
    if (ini >= n)
        return 0;

    if (ini == n-1) {
        if ((type == INC && p->pi[n-1] > divisor && p->pi[n] <= divisor) ||
                (type == DEC && p->pi[n-1] <= divisor && p->pi[n] > divisor)) {
            op_suffix_reversal(p, n-1);
            print_verbose(p, verbose, 'r', n-1, p->size, 0);
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
        if (y < n+1) {
            op_suffix_reversal(p, z+1);
            print_verbose(p, verbose, 'r', z+1, p->size, 0);
            weight += pot(n-z, info->alpha);
            info->nmoves++;
        }
        op_suffix_reversal(p, x+1);
        print_verbose(p, verbose, 'r', x+1, p->size, 0);
        weight += pot(n-x, info->alpha);
        info->nmoves++;
        return weight;
    }

    ini = x+1;

    half = (n-ini+1) / 2.0;
    if ((n-ini+1) % 2)
        half++;
    half += ini-1;

    weight += partitionWSR(p, half+1, 1-type, divisor, verbose, info);
    op_suffix_reversal(p, ini);
    print_verbose(p, verbose, 'r', ini, p->size, 0);
    weight += pot(n-ini+1, info->alpha);
    info->nmoves++;
    weight += partitionWSR(p, n - half + ini, 1-type, divisor, verbose, info);

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
        print_verbose(p, verbose, 'r', j, p->size, 0);
        weight += pot(n-j+1, info->alpha);
        info->nmoves++;
    }

    return weight;
}/*}}}*/

long double WSR(permutation_t *p, int ini, int type, int verbose, info_t *info) {/*{{{*/
    int part, median, n = p->size;
    long double weight = 0;

    if (is_sorted_interval(p, ini, n, type))
        return 0;

    if (ini >= n)
        return 0;

    if (ini == n-1) {
        if ((p->pi[n-1] > p->pi[n] && type == INC) || (p->pi[n-1] < p->pi[n] && type == DEC)) {
            op_suffix_reversal(p, n-1);
            weight += pot(2, info->alpha);
            print_verbose(p, verbose, 'r', n-1, p->size, 0);
            info->nmoves++;
        }
        return weight;
    }

    median = prm_median(p, ini, n);

    weight += partitionWSR(p, ini, 1-type, median, verbose, info);

    part = ini;
    if (type == INC) {
        while (p->pi[part] > median)
            part++;
    } else {
        while (p->pi[part] <= median)
            part++;
    }
    part--;

    weight += WSR(p, part+1, 1-type, verbose, info);
    op_suffix_reversal(p, ini);
    print_verbose(p, verbose, 'r', ini, p->size, 0);
    weight += pot(n-ini+1, info->alpha);
    info->nmoves++;
    weight += WSR(p, n - part + ini, type, verbose, info);

    return weight;
}/*}}}*/

void alg_WSR(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int k = 1;

    while (k == p->pi[k] && k <= p->size)
        k++;

    if (k < p->size+1) {
        info->weight = WSR(p, k, INC, verbose, info);
        info->nof_sr = info->nmoves;
    }

    if (!is_identity(p))
        printf("ERROR! alg_WSR()\n");
}/*}}}*/



long double partitionWSSR(permutation_t *p, permutation_t *sp, int ini, int type, int divisor, int verbose, info_t *info) {/*{{{*/
    int i, j, k, x, y, z, exists, half, n = p->size;
    long double weight = 0;
        
    if (ini >= n)
        return 0;

    if (ini == n-1) {
        if ((type == INC && p->pi[n-1] > divisor && p->pi[n] <= divisor) ||
                (type == DEC && p->pi[n-1] <= divisor && p->pi[n] > divisor)) {
            op_suffix_reversal(p, n-1);
            op_suffix_reversal(sp, sp->size);
            info->weight += 1;
            print_verbose(sp, verbose, 'r', sp->size, sp->size, 0);
            info->nmoves++;
            return weight;
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
        if (y < n+1) {
            op_suffix_reversal(p, z+1);
            op_suffix_reversal(sp, (z+1)/2 + 1);
            info->weight += pot((n-z)/2, info->alpha);
            print_verbose(sp, verbose, 'r', (z+1)/2 + 1, sp->size, 0);
            info->nmoves++;
        }
        op_suffix_reversal(p, x+1);
        op_suffix_reversal(sp, (x+1)/2 + 1);
        info->weight += pot((n-x)/2, info->alpha);
        print_verbose(sp, verbose, 'r', (x+1)/2 + 1, sp->size, 0);
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

    info->weight += partitionWSSR(p, sp, half+1, 1-type, divisor, verbose, info);
    op_suffix_reversal(p, ini);
    op_suffix_reversal(sp, ini/2 + 1);
    info->weight += pot((n-ini+1)/2, info->alpha);
    print_verbose(sp, verbose, 'r', ini/2 + 1, sp->size, 0);
    info->nmoves++;
    info->weight += partitionWSSR(p, sp, n - half + ini, 1-type, divisor, verbose, info);

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
        op_suffix_reversal(sp, j/2 + 1);
        info->weight += pot((n-j+1)/2, info->alpha);
        print_verbose(sp, verbose, 'r', j/2 + 1, sp->size, 0);
        info->nmoves++;
    }

    return weight;
}/*}}}*/

long double WSSR(permutation_t *p, permutation_t *sp, int ini, int type, int verbose, info_t *info) {/*{{{*/
    int part, median, n = p->size;
    long double weight = 0;

    if (is_sorted_interval(p, ini, n, type))
        return 0;

    if (ini >= n)
        return 0;

    if (ini == n-1) {
        if ((p->pi[n-1] > p->pi[n] && type == INC) || (p->pi[n-1] < p->pi[n] && type == DEC)) {
            op_suffix_reversal(p, n-1);
            op_suffix_reversal(sp, sp->size);
            info->weight += 1;
            print_verbose(sp, verbose, 'r', sp->size, sp->size, 0);
            info->nmoves++;
        }
        return weight;
    }

    median = prm_median(p, ini, n);
    if (median % 2)
        median++;

    info->weight += partitionWSSR(p, sp, ini, 1-type, median, verbose, info);

    part = ini;
    if (type == INC) {
        while (p->pi[part] > median)
            part++;
    } else {
        while (p->pi[part] <= median)
            part++;
    }
    part--;

    info->weight += WSSR(p, sp, part+1, 1-type, verbose, info);
    op_suffix_reversal(p, ini);
    op_suffix_reversal(sp, ini/2 + 1);
    info->weight += pot((n-ini+1)/2, info->alpha);
    print_verbose(sp, verbose, 'r', ini/2 + 1, sp->size, 0);
    info->nmoves++;
    info->weight += WSSR(p, sp, n - part + ini, type, verbose, info);

    return weight;
}/*}}}*/

void alg_WSSR(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int k;
    permutation_t up;
    long double a;

    create_permutation(&up, 2*p->size, UNSIGNED);
    prm_image(p, &up);

    k = 1;
    while (k == up.pi[k] && k <= up.size)
        k++;

    if (k < up.size+1) {
        a = WSSR(&up, p, k, INC, verbose, info);
        info->nof_sr = info->nmoves;
    }

    if (!is_identity(p)) 
        printf("ERROR! alg_WSSR()\n");

    destroy_permutation(&up);
}/*}}}*/

