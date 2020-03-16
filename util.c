/***********************************************************
 * Created: Mon 04 Feb 2013 03:35:22 PM BRST
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "permutations.h"
#include "util.h"

void* Malloc(int size) {
    void *p = malloc(size);
    if (!p) {
        printf("ERROR: Malloc()\n");
        exit(0);
    }
    return p;
}

void init_info(info_t *info) {
    info->nmoves = 0;
    info->weight = 0;
    info->alpha = 0;
    info->nof_pr = 0;
    info->nof_sr = 0;
    info->nof_pt = 0;
    info->nof_st = 0;
    info->nof_r = 0;
    info->nof_t = 0;
}

void print_verbose(permutation_t *p, int verbose, char r, int i, int j, int k) {
    if (verbose == 5) {
        if (r == 'r') {
            printf("%d %d\n", i, j);
        } else if (r == 't') {
            printf("%d %d %d\n", i, j, k);
        } else {
            printf("print_verbose: Wrong rearrangement\n");
            exit(0);
        }
    }

    if (verbose == 3 || verbose == 4) {
        print_permutation(p);
    }
}

long double pot(long double x, long double y) {
    if (y == 0) return 1;
    if (y == 1) return x;
    return powl(x, y);
}

int inc_mod(int i, int n) {
    if (i == -1) return -n;
    if (i == n) return 1;
    return i+1;
}

int dec_mod(int i, int n) {
    if (i == -n) return -1;
    if (i == 1) return n;
    return i-1;
}

void exchange(int *a, int *b) {
    int aux;
    aux = *a;
    *a = *b;
    *b = aux;
}
