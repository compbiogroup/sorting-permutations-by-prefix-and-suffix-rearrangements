/***********************************************************
 * Created: Qua 26 Ago 2015 13:57:23 BRT
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************
 *
 * For a given algorithm and size, generate the approximation factors
 *
 ***********************************************************/
 

#define _LARGEFILE64_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "../../permutations.h"
#include "../../rearrangements.h"

long double calculate_lower_bound(char *perm, int size, char *problem, int sig);

int main(int argc, char *argv[]) {
    char *problem, perm[5000];
    int n, null, nmoves, sig;
    unsigned int nperms;
    FILE *fresults, *fout;
    long double lb, weight, avgapprox, approx, minapprox, maxapprox;

    if (argc != 6) {
        printf("usage: %s file_results size ", argv[0]);
        printf("problem signed? fileout\n");
        return 0;
    }

    fresults = fopen64(argv[1], "r");
    n = atoi(argv[2]);
    problem = argv[3];
    sig = atoi(argv[4]);
    fout = fopen64(argv[5], "w");

    nperms = 0;
    maxapprox = 0;
    minapprox = LDBL_MAX;
    avgapprox = 0;
    while (fscanf(fresults, "%s %d %Lf %d %d %d %d %d %d", perm, &nmoves,
                            &weight, &null, &null, &null, &null, &null, &null) != EOF) {
        lb = calculate_lower_bound(perm, n, problem, sig);
        nperms++;

        if (lb == 0) {
            approx = 1;
        } else {
            approx = (long double) nmoves / lb;
        }

        avgapprox += approx;
        if (approx > maxapprox)
            maxapprox = approx;
        if (approx < minapprox)
            minapprox = approx;

        fprintf(fout, "%Lf\n", approx);
    }

    printf("%d %Lf %Lf %Lf\n", n, avgapprox/nperms, maxapprox, minapprox);

    fclose(fresults);
    fclose(fout);
    return 0;
}

long double calculate_lower_bound(char *perm, int size, char *problem, int sig) {
    long double lb;
    permutation_t p;
    create_permutation(&p, size, sig);
    prmfill_fromstr(perm, &p);
    lb = lower_bound(problem, &p);
    destroy_permutation(&p);
    return lb;
}

