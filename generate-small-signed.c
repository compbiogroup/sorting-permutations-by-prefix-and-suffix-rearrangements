/***********************************************************
 * Created: Dom 19 Jan 2014 00:15:09 BRST
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************/
 

#include <stdio.h>
#include <stdlib.h>
#include "permutations.h"


unsigned long fat(int size) {
    if (size == 2 || size == 1)
        return size;
    return size * fat(size - 1);
}

unsigned long expo(int x, int y) {
    int i;
    unsigned long ret = x;

    for (i = 1; i < y; i++) {
        ret = ret * x;
    }
    return ret;
}
                
unsigned long total(int size) {
    return expo(2, size) * fat(size);
}

int main(int argc, char *argv[]) {
    int size;
    permutation_t p;
    unsigned long i, max;
    char str[500];

    if (argc != 2) {
        printf("ERRO! Especifique o tamanho da permutação\n");
    }

    size = atoi(argv[1]);

    if (size > 32000) {
        printf("Only sizes until 32000\n");
        exit(0);
    }

    create_permutation(&p, size, 1);

    max = total(size);
    for (i = 0; i < max; i++) {
        unrank(i, &p);
        prm2str(&p, str);
        printf("%s\n", str);
    }

    return 0;
}
