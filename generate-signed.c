 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define _GNU_SOURCE
#include <getopt.h>
#include <math.h>


int permcompare(short int *permutation, int size, short int *allperms, int n) {
    int i;
    for (i = 0; i < n; i++)
        if (memcmp(permutation, allperms + i * size, size * sizeof(short int)) == 0)
            return 0;
    return 1;
}


int main(int argc, char *argv[]) {
    extern char *optarg;
    char op;
    int size, nperms, n, i, t, bp;
    short int *permutation, *allperms;

    struct option longopts[] = {
        {"size", 1, NULL, 's'},
        {"nperms", 1, NULL, 'n'},
        {"help", 0, NULL, 'h'},
    };

    while ((op = getopt_long(argc, argv, "s:n:h", longopts, NULL)) != -1) {/*{{{*/
        switch (op) {
            case 's':
		        size = atoi(optarg);
                break;
            case 'n':
		        nperms = atoi(optarg);
                break;
            case 'h':
                printf("usage: %s [options] \n", argv[0]);
                printf("options: \n");
                printf("\t--size or -s: permutation size\n");
                printf("\t--nperms or -n: number of permutations to be created\n");
                exit(0);
        }
    }/*}}}*/

    permutation = malloc(size * sizeof (short int));
    allperms = malloc(size * nperms * sizeof(short int));
    
    if (allperms == NULL) {
        printf("Out of memory!\n");
        return 1;
    }
      
    memset(allperms, 0, size  * nperms * sizeof(short int));
    
    for (n = 0; n < nperms; n++) {
        memset(permutation, 0, size * sizeof(short int));
        
        do {
            for (i = 0; i < size; i++) {
                t = rand() % (i + 1);
                permutation[i] = permutation[t];
                permutation[t] = i;
            }
            t = rand();
            for (i = 0; i < size; i++) {
                permutation[i]++;
                if (rand() % 2 == 0)
                    permutation[i] *= -1;
            }
            bp = 0;
            if (permutation[size-1] == size)
                bp++;
            if (permutation[0] == 1)
                bp++;
            for (i = 0; i < size-1 && bp == 0; i++) {
                if (permutation[i+1] - permutation[i] == 1)
                    bp++;
            }
        } while (!permcompare(permutation, size, allperms, n) || bp != 0);

        memcpy(allperms + n*size, permutation, size * sizeof(short int));

        for (i = 0; i < size-1; i++) {
            printf("%d,", permutation[i]);
        }
        printf("%d\n", permutation[size-1]);
    }
    
    return 0;
}
