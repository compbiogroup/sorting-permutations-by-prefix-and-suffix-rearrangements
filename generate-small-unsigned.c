 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define _GNU_SOURCE
#include <getopt.h>
#include <math.h>
#include "permutations.h"


unsigned long fat(int size) {
    if (size == 2 || size == 1)
        return size;
    return size * fat(size - 1);
}

unsigned long total(int size) {
    return fat(size);
}

int main(int argc, char *argv[]) {
    int size;
    permutation_t p;
    unsigned long i, max;
    char str[100];

    if (argc != 2) {
        printf("ERRO! Especifique o tamanho da permutação\n");
    }

    size = atoi(argv[1]);

    if (size > 32000) {
        printf("Only sizes until 32000\n");
        exit(0);
    }

    create_permutation(&p, size, 0);

    max = total(size);
    for (i = 0; i < max; i++) {
        unrank(i, &p);
        prm2str(&p, str);
        printf("%s\n", str);
    }

    return 0;
}
