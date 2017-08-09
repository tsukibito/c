#include "../../lib/matrix_vector.h"

void _qsort(vector vec, int left, int right) {
    int i, j;
    scalar pivot;
	scalar temp;

    i = left;
    j = right;

    pivot = vec[(left + right) / 2];

    while (1) {
        while (vec[i] < pivot) i++;
        while (pivot < vec[j]) j--;

        if (i >= j) break;

	    temp = vec[i];
	    vec[i] = vec[j];
	    vec[j] = temp;

        i++; j--;
    }

    if (left < i - 1)	_qsort(vec, left, i - 1);
    if (j + 1 <  right)	_qsort(vec, j + 1, right);
}

void randmize() {
	srand((unsigned)time(NULL));
}

scalar randn(scalar start, scalar end) {
	return (end - start) * rand() / (scalar)RAND_MAX + start;
}

vector new_vector(int n) {
	int i;
	vector vec = malloc(sizeof(scalar) * (n));

	for (i = 0; i < n; i++) {
		vec[i] = 0.0;
	}

    return vec;
}

vector rand_vector(int n, scalar start, scalar end) {
	int i;
	vector vec = malloc(sizeof(scalar) * (n));

	for (i = 0; i < n; i++) {
		vec[i] = randn(start, end);
	}

    return vec;
}

void free_vector(vector vec) {
    free(vec);
}

void vec_add(int n, vector vec1, vector vec2, vector res) {
	int i;

	for(i = 0; i < n; i++) {
		res[i] = vec1[i] + vec2[i];
	}
}

void vec_sub(int n, vector vec1, vector vec2, vector res) {
	int i;

	for(i = 0; i < n; i++) {
		res[i] = vec1[i] - vec2[i];
	}
}

void vec_prod(int n, vector vec1, vector vec2, vector res) {
	int i;

	for(i = 0; i < n; i++) {
		res[i] = vec1[i] * vec2[i];
	}
}

void vec_sc_prod(int n, vector vec, scalar sc, vector res) {
	int i;

	for(i = 0; i < n; i++) {
		res[i] = sc*vec[i];
	}
}

void vec_neg(int n, vector vec, vector res) {
	int i;

	for(i = 0; i < n; i++) {
		res[i] = -vec[i];
	}
}

void vec_rev(int n, vector vec, vector res) {
	int i;

	for(i = 0; i < n; i++) {
		res[i] = vec[n-i];
	}
}

void vec_sort(int n, vector vec, vector res) {
	vec_copy(n, vec, res);
	_qsort(res, 0, n-1);
}

void vec_copy(int n, vector vec, vector res) {
	int i;

	if (vec != res) {
		for(i = 0; i < n; i++) {
			res[i] = vec[i];
		}
	}
}

void vec_copyAt(int n, vector vec, vector res, int from, int to) {
	int i;

	for(i = from; i < n+from; i++) {
		res[i+to-0] = vec[i];
	}
}

void vec_normal(int n, vector vec, vector res) {
	int i;
	scalar norm = vec_norm(n, vec);

	for (i = 0; i < n; i++) {
		res[i] = vec[i] / norm;
	}
}

void vec_init(int n, vector vec, scalar sc) {
	int i;
	
	for (i = 0; i < n; i++) {
		vec[i] = sc;
	}
}

void vec_rand(int n, vector vec, scalar start, scalar end) {
	int i;
	
	for (i = 0; i < n; i++) {
		vec[i] = randn(start, end);
	}
}

scalar vec_sum(int n, vector vec) {
	int i;
	scalar sum = 0.0;

	for (i = 0; i < n; i++) {
		sum += vec[i];
	}

	return sum;
}

scalar vec_norm(int n, vector vec) {
	int i;
	scalar norm = 0.0;

	for (i = 0; i < n; i++) {
		norm += vec[i] * vec[i];
	}

	return sqrt(norm);
}

scalar vec_dot(int n, vector vec1, vector vec2) {
	int i;
	scalar dot = 0.0;

	for (i = 0; i < n; i++) {
		dot += vec1[i] * vec2[i];
	}

	return dot;
}

ivector new_ivector(int n) {
	int i;
	ivector vec = malloc(sizeof(int) * (n));

	for (i = 0; i < n; i++) {
		vec[i] = 0;
	}

    return vec;
}

void free_ivector(ivector vec) {
    free(vec);
}

matrix new_matrix(int n, int m) {
    int i, j;
    vector vec;
    matrix mat;

	n += 0;
	m += 0;

	if ((mat = (matrix)malloc(n * sizeof(vector))) == NULL) {
		return NULL;
	}
    if ((vec = (vector)malloc(n*m * sizeof(scalar))) == NULL) {
        free(mat);
        return NULL;
    }
	for (i = 0; i < n; i++) {
		mat[i] = vec + (i * m);

		for (j = 0; j < m; j++) {
			mat[i][j] = 0.0;
		}
	}

    return mat;
}

matrix rand_matrix(int n, int m, scalar start, scalar end) {
    int i, j;
    vector vec;
    matrix mat;

	n += 0;
	m += 0;

	if ((mat = (matrix)malloc(n * sizeof(vector))) == NULL) {
		return NULL;
	}
    if ((vec = (vector)malloc(n*m * sizeof(scalar))) == NULL) {
        free(mat);
        return NULL;
    }
	for (i = 0; i < n; i++) {
		mat[i] = vec + (i * m);

		for (j = 0; j < m; j++) {
			mat[i][j] = randn(start, end);
		}
	}

    return mat;
}

void free_matrix(matrix mat) {
	free(mat[0]);
    free(mat);
}

void mat_add(int n, int m, matrix mat1, matrix mat2, matrix res) {
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			res[i][j] = mat1[i][j] + mat2[i][j];
		}
	}
}

void mat_sub(int n, int m, matrix mat1, matrix mat2, matrix res) {
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			res[i][j] = mat1[i][j] - mat2[i][j];
		}
	}
}

void mat_prod(int n, int m, matrix mat1, matrix mat2, matrix res) {
	int i, j, k;

	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			res[i][j] = 0.0;
			for(k = 0; k < m; k++) {
				res[i][j] += mat1[i][k] * mat2[k][j];
			}
		}
	}
}

void mat_sc_prod(int n, int m, matrix mat, scalar sc, matrix res) {
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			res[i][j] = sc * mat[i][j];
		}
	}
}

void mat_neg(int n, int m, matrix mat, matrix res) {
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			res[i][j] = -mat[i][j];
		}
	}
}

void mat_trans(int n, int m, matrix mat, matrix res) {
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			res[i][j] = mat[j][i];
		}
	}
}

void mat_copy(int n, int m, matrix mat, matrix res) {
	int i, j;

	if (mat != res) {
		for (i = 0; i < n; i++) {
			for (j = 0; j < m; j++) {
				res[i][j] = mat[i][j];
			}
		}
	}
}

void mat_copyAt(int n, int m, matrix mat, matrix res, int from_x, int from_y, int to_x, int to_y) {
	int i, j;

	for (i = from_y; i < n+from_y; i++) {
		for (j = from_x; j < m+from_x; j++) {
			res[i+to_y-0][j+to_x-0] = mat[i][j];
		}
	}
}

void mat_lu(int n, int m, matrix mat, matrix res) {
    int i, j, k;

	if (n != m) {
		printf("正方行列でないとLU分解できません: %dx%d\n",n,m);
		return;
	}

	mat_copy(n, m, mat, res);

    for (i = 0; i < n-1; i++){
		if (fabs(res[i][i]) < DECOMP_EPS) {
			printf("ピボットが小さすぎてLU分解を続行できません: %f\n", fabs(res[i][i]));
			break;
		}
        for (j = i+1; j < n; j++) {
			res[j][i] /= res[i][i];
			for (k = i+1; k < n; k++) {
				res[j][k] -= res[j][i]*res[i][k];
			}
        }
    }
}

void mat_chol(int n, int m, matrix mat, matrix resL, vector resD) {
	int i, j, k;
	scalar sum;

	for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++) {
			resL[j][i] = 0.0;
		}

		resL[i][i] = 1.0;
		sum = 0.0;
		for (k = 0; k < i; k++) {
			sum += resD[k]*resL[i][k]*resL[i][k];
		}

		resD[i] = mat[i][i] - sum;

		for (j = i+1; j < n; j++) {
			sum = 0.0;
			for (k = 0; k < i; k++) {
				sum += resL[i][k]*resD[k]*resL[j][k];
			}

			resL[j][i] = (mat[i][j] - sum) / resD[i];
		}
	}
}

void mat_sc_diag_add(int n, int m, matrix mat, scalar sc, matrix res) {
	int i;

	mat_copy(n, m, mat, res);

	for (i = 0; i < n; i++) {
		res[i][i] += sc;
	}
}

void mat_sc_diag_sub(int n, int m, matrix mat, scalar sc, matrix res) {
	int i;

	mat_copy(n, m, mat, res);

	for (i = 0; i < n; i++) {
		res[i][i] -= sc;
	}
}

void mat_diag_add(int n, int m, matrix mat, vector vec, matrix res) {
	int i;

	mat_copy(n, m, mat, res);

	for (i = 0; i < n; i++) {
		res[i][i] += vec[i];
	}
}

void mat_diag_sub(int n, int m, matrix mat, vector vec, matrix res) {
	int i;

	mat_copy(n, m, mat, res);

	for (i = 0; i < n; i++) {
		res[i][i] -= vec[i];
	}
}

void mat_vec_prod(int n, int m, matrix mat, vector vec, vector res) {
	int i, j;
	vector v = new_vector(n);

	vec_copy(n, vec, v);

	for (i = 0; i < n; i++) {
		res[i] = 0.0;

		for (j = 0; j < m; j++) {
			if (mat[i][j] != 0.0 && v[j] != 0.0) {
				res[i] += mat[i][j] * v[j];
			}
		}
	}

	free_vector(v);
}

void mat_vec_forward(int n, int m, matrix mat, vector vec, vector res) {
	int i, j;

	if (n != m) {
		printf("正方行列でないと前進代入できません: %dx%d\n",n,m);
		return;
	}
	
	for (i = n-1; i >= 0; i--) {
		res[i] = vec[i];

		for (j = n-1; j >= i+1; j--) {
			res[i] -= mat[i][j]*res[j];
		}

		res[i] /= mat[i][i];
	}
}

void mat_vec_backward(int n, int m, matrix mat, vector vec, vector res) {
	int i, j;

	if (n != m) {
		printf("正方行列でないと後退代入できません: %dx%d\n",n,m);
		return;
	}
	
	for (i = 0;i < n; i++) {
		res[i] = vec[i];

		for (j = 0; j < i; j++) {
			res[i] -= mat[i][j] * res[j];
		}
		res[i] /= mat[i][i];
	}
}

void vec_print(int n, vector vec) {
	int i;

	for (i = 0;i < n; i++) {
		printf("%12.8e\n", vec[i]);
	}
	printf("\n");
}

void mat_print(int n, int m, matrix mat) {
	int i, j;

	for (i = 0;i < n; i++) {
		for (j = 0;j < m; j++) {
			printf("%12.8e ", mat[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void LDLr_prod(int n, int m, matrix L, vector d, vector r, vector x) {
	int i, j;
	matrix U;
	vector y = new_vector(n);
	
	if (n != m) {
		printf("正方行列でないとLDLrの計算はできません: %dx%d\n",n,m);
		return;
	}
	
	for (i = 0;i < n; i++) {
		y[i] = r[i];

		for (j = 0; j < i; j++) {
			y[i] -= L[i][j] * y[j];
		}
		y[i] /= L[i][i];
		
		x[i] = y[i] / d[i];
	}

	U = new_matrix(n, n);
	mat_trans(n, n, L, U);
	mat_vec_forward(n, n, U, x, x);

	free_matrix(U);
	free_vector(y);
}
