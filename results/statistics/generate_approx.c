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


long double calculate_lower_bound(char *perm, int size, char *problem, double alpha, int sig);
void fill_dists(FILE *file, long double *dists, int bin, int sig, int n);
unsigned int fat(int n);
unsigned int expo(int n);

int main(int argc, char *argv[]) {
    char *lbtype, *problem, perm[5000];
    int n, bin, null, nmoves, sig;
    unsigned int nperms;
    FILE *fresults, *fdists, *fout;
    long double alpha, weight, avgapprox, approx, minapprox, maxapprox, lb, *dists;

    if (argc != 10) {
        printf("usage: %s file_results size type_of_lower_bound(T|D)", argv[0]);
        printf(" problem alpha file_dists binary? signed? fileout\n");
        return 0;
    }

    fresults = fopen64(argv[1], "r");
    n = atoi(argv[2]);
    lbtype = argv[3];
    problem = argv[4];
    alpha = atof(argv[5]);
    if (lbtype[0] == 'D')
        fdists = fopen64(argv[6], "r");
    bin = atoi(argv[7]);
    sig = atoi(argv[8]);

    if (lbtype[0] == 'D') {
        if (sig)
            dists = malloc(sizeof(long double) * fat(n) * expo(n));
        else
            dists = malloc(sizeof(long double) * fat(n));
        fill_dists(fdists, dists, bin, sig, n);
    }

    fout = fopen64(argv[9], "w");

    nperms = 0;
    maxapprox = 0;
    minapprox = LDBL_MAX;
    while (fscanf(fresults, "%s %d %Lf %d %d %d %d %d %d", perm, &nmoves,
                            &weight, &null, &null, &null, &null, &null, &null) != EOF) {
        if (lbtype[0] == 'D')
            lb = dists[nperms];
        else
            lb = calculate_lower_bound(perm, n, problem, alpha, sig);
        
        nperms++;

        if (lb == 0) {
            approx = 1;
        } else {
            if (alpha != 0)
                approx = (long double) weight / lb;
            else
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
    if (lbtype[0] == 'D') {
        free(dists);
        fclose(fdists);
    }
    return 0;
}

long double calculate_lower_bound(char *perm, int size, char *problem, double alpha, int sig) {
    long double lb;
    permutation_t p;
    create_permutation(&p, size, sig);
    prmfill_fromstr(perm, &p);
    p.alpha = alpha;
    lb = lower_bound(problem, &p);
    destroy_permutation(&p);
    return lb;
}

void fill_dists(FILE *file, long double *dists, int bin, int sig, int n) {
    unsigned int i = 0;
    long double d;
    char c;
    permutation_t p;
    char str[100], lixo[100];

    if (bin) {
        while (fscanf(file, "%c", &c) != EOF) {
            dists[i++] = c;
            fscanf(file, "%c", &c);
            fscanf(file, "%c", &c);
            fscanf(file, "%c", &c);
            fscanf(file, "%c", &c);
        }
    } else {
        create_permutation(&p, n, sig);
        while (fscanf(file, "%s %Lf %s", str, &d, lixo) != EOF) {
            prmfill_fromstr(str, &p);
            i = rank(&p);
            dists[i] = d;
        }
        destroy_permutation(&p);
    }
}

unsigned int fat(int n) {
    if (n == 1)
        return 1;
    return n * fat(n-1);
}

unsigned int expo(int n) {
    unsigned int e = 2, i;
    for (i = 1; i < n; i++)
        e *= 2;
    return e;
}

