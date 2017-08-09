#include "mm.h"

matrix mm_matrix(char *fileName, int *n, int *m, int start) {
	FILE *fp;
	long index;
	char buf[MM_MAX_LINE_BUF];
	int total;
	int i, j, x, y;
	matrix mat;

#ifdef _MSC_VER
	if ((fopen_s(&fp, fileName, "r")) != 0) {
#else
	if ((fp = fopen(fileName, "r")) == NULL) {
#endif
		printf("mm [Error]: %s は存在しないか開けません\n", fileName);
		return NULL;
	}

	index = 0L;

	while(fgetc(fp) == '%') {
		fseek(fp, index, SEEK_SET);
		fgets(buf, 256, fp);
		index = ftell(fp) - 1L;
	}

	fseek(fp, index, SEEK_SET);

#ifdef _MSC_VER
	fscanf_s(fp, "%d %d %d", n, m, &total);
#else
	fscanf(fp, "%d %d %d", n, m, &total);
#endif

	if (*m <= 0 || *n <= 0) {
		fclose(fp);

		printf("mm [Error]: %s は対応していないフォーマットです\n", fileName);
		return NULL;
	}

	mat = new_matrix(*m, *n);

	for (i = 0; i <= *m; i++) {
		for (j = 0; j <= *n; j++) {
			mat[i][j] = 0.0;
		}
	}

	for (i = 0; i < total; i++) {
#ifdef _MSC_VER
		fscanf_s(fp, "%d %d", &x, &y);
		fscanf_s(fp, "%lf", &mat[y-!start][x-!start]);
#else
		fscanf(fp, "%d %d", &x, &y);
		fscanf(fp, "%lf", &mat[y-!start][x-!start]);
#endif
	}

	return mat;
}

void mm_text(char *inputName, char *outputName) {
	FILE *fp;
	int n, m, i, j;
	matrix mat = mm_matrix(inputName, &n, &m, 0);

#ifdef _MSC_VER
	if ((fopen_s(&fp, outputName, "w")) != 0) {
#else
	if ((fp = fopen(outputName, "w")) == NULL) {
#endif
		printf("mm [Error]: %s に書き込めません\n", outputName);
		return;
	}

	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			fprintf(fp, "%.13e\t", mat[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}
