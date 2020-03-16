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


void alg_2PT(permutation_t *p, int verbose, info_t* info) {/*{{{*/
    int i, j;

    while (!is_identity(p)) {
        j = p->inv_pi[p->pi[1] - 1] + 1;
        i = p->inv_pi[p->pi[j] - 1] + 1;
        if (1 < i && i < j) {
            op_prefix_transposition(p, i, j);
            print_verbose(p, verbose, 't', 1, i, j);
        } else {
            i = 1;
            while (p->pi[i] == p->pi[i+1] - 1)
                i++;

            j = p->inv_pi[p->pi[i]+1];
            op_prefix_transposition(p, i+1, j);
            print_verbose(p, verbose, 't', 1, i+1, j);
        }
        info->nmoves++;
    }
}/*}}}*/

