#ifndef __MATRIX_VECTOR_UTIL__
#define __MATRIX_VECTOR_UTIL__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifndef true
#define true	1
#endif
#ifndef false
#define false	0
#endif

#ifndef MAT_VEC_IDX
#define MAT_VEC_IDX	0
#endif

#ifndef DECOMP_EPS
#define DECOMP_EPS	1e-8
#endif

typedef double scalar;
typedef scalar *vector;
typedef vector *matrix;

typedef int *ivector;

void randmize();
scalar randn(scalar start, scalar end);

vector new_vector(int n);
vector rand_vector(int n, scalar start, scalar end);
void free_vector(vector vec);

void vec_add(int n, vector vec1, vector vec2, vector res);
void vec_sub(int n, vector vec1, vector vec2, vector res);
void vec_prod(int n, vector vec1, vector vec2, vector res);
void vec_sc_prod(int n, vector vec, scalar sc, vector res);
void vec_neg(int n, vector vec, vector res);
void vec_rev(int n, vector vec, vector res);
void vec_sort(int n, vector vec, vector res);
void vec_copy(int n, vector vec, vector res);
void vec_copyAt(int n, vector vec, vector res, int from, int to);
void vec_normal(int n, vector vec, vector res);
void vec_init(int n, vector vec, scalar sc);
void vec_rand(int n, vector vec, scalar start, scalar end);
scalar vec_sum(int n, vector vec);
scalar vec_norm(int n, vector vec);
scalar vec_dot(int n, vector vec1, vector vec2);

ivector new_ivector(int n);
void free_ivector(ivector vec);

matrix new_matrix(int n, int m);
matrix rand_matrix(int n, int m, scalar start, scalar end);
void free_matrix(matrix mat);

void mat_add(int n, int m, matrix mat1, matrix mat2, matrix res);
void mat_sub(int n, int m, matrix mat1, matrix mat2, matrix res);
void mat_prod(int n, int m, matrix mat1, matrix mat2, matrix res);
void mat_sc_prod(int n, int m, matrix mat, scalar sc, matrix res);
void mat_neg(int n, int m, matrix mat, matrix res);
void mat_trans(int n, int m, matrix mat, matrix res);
void mat_copy(int n, int m, matrix mat, matrix res);
void mat_copyAt(int n, int m, matrix mat, matrix res, int from_x, int from_y, int to_x, int to_y);

//void mat_lu(int n, int m, matrix mat, matrix res);
//void mat_chol(int n, int m, matrix mat, matrix resL, vector resD);
void mat_cholinc(int n, int m, matrix mat, matrix resL, vector resD);

void mat_sc_diag_add(int n, int m, matrix mat, scalar sc, matrix res);
void mat_sc_diag_sub(int n, int m, matrix mat, scalar sc, matrix res);
void mat_diag_add(int n, int m, matrix mat, vector vec, matrix res);
void mat_diag_sub(int n, int m, matrix mat, vector vec, matrix res);

void mat_vec_prod(int n, int m, matrix mat, vector vec, vector res);
void mat_vec_forward(int n, int m, matrix mat, vector vec, vector res);
void mat_vec_backward(int n, int m, matrix mat, vector vec, vector res);

void vec_print(int n, vector vec);
void mat_print(int n, int m, matrix mat);

void LDLr_prod(int n, int m, matrix L, vector d, vector r, vector x);

#endif
