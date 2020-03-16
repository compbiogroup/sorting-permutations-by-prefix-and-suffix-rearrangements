/***********************************************************
 * Created: Sex 28 Ago 2015 10:29:24 BRT
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************
 *
 * Compare the actual distance when using only prefix and when using prefix +
 * suffix
 *
 ***********************************************************/
 

#define _LARGEFILE64_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include "../../permutations.h"

void fill_dists(FILE *file, long double *dists, int bin);
unsigned int fat(int n);
unsigned int expo(int n);

int main(int argc, char *argv[]) {
    FILE *fpref, *fsuff, *fout, *fperms;
    int n, bin, sig, maxn;
    char fname[1000], prmstr[100];
    long double *dists_pref, *dists_suff;
    unsigned int total, greater, i, maxdiff;
    permutation_t p;

    if (argc != 7) {
        printf("Usage: %s path_dists_pref path_dists_prefsuff bin? sig? fileout fileperms\n", argv[0]);
        return 0;
    }

    bin = atoi(argv[3]);
    sig = atoi(argv[4]);
    fout = fopen64(argv[5], "w");

    if (sig) {
        maxn = 9;
        dists_pref = malloc(sizeof(long double) * fat(9) * expo(9));
        dists_suff = malloc(sizeof(long double) * fat(9) * expo(9));
        create_permutation(&p, 9, SIGNED);
    } else {
        maxn = 10;
        dists_pref = malloc(sizeof(long double) * fat(10));
        dists_suff = malloc(sizeof(long double) * fat(10));
        create_permutation(&p, 10, UNSIGNED);
    }


    fprintf(fout, "#n qtd_greater %%greater maxdiff\n");

    for (n = 2; n <= maxn; n++) {
        p.size = n;
        sprintf(fname, "%s/%d.txt", argv[1], n);
        fpref = fopen64(fname, "r");
        sprintf(fname, "%s/%d.txt", argv[2], n);
        fsuff = fopen64(fname, "r");

        fill_dists(fpref, dists_pref, bin, maxn, sig);
        fill_dists(fsuff, dists_suff, bin, maxn, sig);

        if (sig) {
            total = fat(n) * expo(n);
        } else {
            total = fat(n);
        }

        maxdiff = greater = 0;
        for (i = 0; i < total; i++) {
            if (dists_pref[i] > dists_suff[i]) {
                greater++;
                if (maxdiff < dists_pref[i] - dists_suff[i])
                    maxdiff = dists_pref[i] - dists_suff[i];
            }
        }

        sprintf(fname, "%s.%d", argv[6], n);
        fperms = fopen64(fname, "w");
        for (i = 0; i < total; i++) {
            if (dists_pref[i] - dists_suff[i] == maxdiff) {
                unrank(i, &p);
                prm2str(&p, prmstr);
                fprintf(fperms, "%s\n", prmstr);
            }
        }

        fprintf(fout, "%d %d %Lf %d\n", n, greater, (long double) greater/total*100, maxdiff);

        fclose(fpref);
        fclose(fsuff);
        fclose(fperms);
    }

    destroy_permutation(&p);
    fclose(fout);

    return 0;
}


void fill_dists(FILE *file, long double *dists, int bin, int maxn, int sig) {
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
        create_permutation(&p, maxn, sig);
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

