#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>
#include "aes.h"
#include "const.c"
#include <omp.h>

int main(){

    long size;
    char *buffer;

    char filepath[] = "C:/Users/nethm/Downloads/sample1.txt";

    FILE *fp;

    fp = fopen(filepath, "r");
    if (fp == NULL) {
        printf("Failed to open file\n");
        exit(1);
    }

    fseek(fp, 0, SEEK_END);
    size = ftell(fp);
    fseek(fp, 0, SEEK_SET);

    buffer = (char*) malloc(size + 1);
    if (buffer == NULL) {
        printf("Failed to allocate memory\n");
        exit(1);
    }

    fread(buffer, size, 1, fp);
    buffer[size] = '\0';

    fclose(fp);

    //printf("%s\n", buffer);

    int charCount1 = 0;
    int charCount = 0;

    //prints the character count in the original input
    for (int i = 0; i < strlen(buffer); i++) {
        charCount1 = charCount1 + 1;
    }

    printf("\nThe character count in the original input: %d\n\n",strlen(buffer));


       return 0;
}
