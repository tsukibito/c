#ifndef MMUTIL
#define MMUTIL

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../../lib/matrix_vector.h"

matrix mm_matrix(char fileName[], int *sizeX, int *sizeY, int start) {
	FILE *fp;
	long index;
	char buf[256];
	int total;
	int i, j, x, y;
	matrix mat;

	fp = fopen(fileName, "r");

	/*do {
		index = ftell(fp) - 1L;
		fgets(buf, 256, fp);
	} while(buf[0] == '%');
	*/

	fseek(fp, index, SEEK_SET);

	fscanf(fp, "%d", sizeX);
	fscanf(fp, "%d", sizeY);
	fscanf(fp, "%d", &total);

	mat = new_matrix(*sizeY, *sizeX);

	for (i = 0; i < *sizeY; i++) {
		for (j = 0; j < *sizeX; j++) {
			mat[i][j] = 0.0;
		}
	}

	for (i = 0; i < total; i++) {
		fscanf(fp, "%d", &x);
		fscanf(fp, "%d", &y);
		fscanf(fp, "%lf", &mat[y-!start][x-!start]);
	}

	return mat;
}

void mm_text(char inputName[], char outputName[]) {
	FILE *fp;
	int n, m, i, j;
	matrix mat = mm_matrix(inputName, &n, &m, 0);

	fp = fopen(outputName, "w");
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			fprintf(fp, "%f\t", mat[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

#endif