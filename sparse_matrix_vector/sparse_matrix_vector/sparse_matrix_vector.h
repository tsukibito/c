#ifndef __SPARSE_MATRIX_VECTOR_UTIL__
#define __SPARSE_MATRIX_VECTOR_UTIL__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void randmize();
double drand(double min, double max);

typedef struct _dvector{
	int n;
	double *elements;
};

typedef struct _dvector *dvector;

dvector new_dvector(int n);
void clear_dvector(dvector vec, double sc);
void rand_dvector(dvector vec, double min, double max);
void copy_elements_dvector(dvector source, dvector target);
dvector copy_dvector(dvector source);
void free_dvector(dvector vec);
void add_dvector(dvector vec1, dvector vec2, dvector res);
void sub_dvector(dvector vec1, dvector vec2, dvector res);
void prod_dvector(dvector vec1, dvector vec2, dvector res);
void prod_scalar_dvector(dvector vec, double sc, dvector res);
void neg_dvector(dvector vec, dvector res);
void rev_dvector(dvector vec, dvector res);
void normal_dvector(dvector vec, dvector res);
double sum_dvector(dvector vec);
double norm_dvector(dvector vec);
double dot_dvector(dvector vec1, dvector vec2);
void sort_dvector(dvector vec, dvector res);

typedef struct _dnode {
	int index;
	double value;
	struct _dnode *next;
} dnode;

typedef dnode *dnode_ptr;

typedef struct {
	dnode_ptr first;
	dnode_ptr diag;
} first_dnode;

typedef first_dnode *first_dnode_ptr;

typedef struct _dmatrix_sp {
	int n;
	first_dnode *elements;
};

typedef struct _dmatrix_sp *dmatrix_sp;

typedef struct _LDU_dmatrix_sp {
	dmatrix_sp matL;
	dmatrix_sp matU;
	dvector vecD;
};

typedef struct _LDU_dmatrix_sp *LDU_dmatrix_sp;

dmatrix_sp new_dmatrix_sp(int nt);
void clear_dmatrix_sp(dmatrix_sp mat, double sc);
void free_dmatrix_sp(dmatrix_sp mat);
dnode_ptr* set_dmatrix_sp_element(dmatrix_sp mat, int row_index, int col_index, double value, dnode_ptr *current);
void trans_dmatrix_sp(dmatrix_sp mat, dmatrix_sp res);
void add_dmatrix_sp(dmatrix_sp mat0, dmatrix_sp mat1, dmatrix_sp mat_res);
void sub_dmatrix_sp(dmatrix_sp mat0, dmatrix_sp mat1, dmatrix_sp mat_res);
void prod_dmatrix_sp(dmatrix_sp mat, dvector vec, dvector res);

void prod_dmatrix_sp_diag(dmatrix_sp mat, dvector vec, double sc, dvector res);

void forward_dmatrix_sp(dmatrix_sp mat, dvector vec, dvector res);
void backward_dmatrix_sp(dmatrix_sp mat, dvector vec, dvector res);
void LDU_prod_dmatrix_sp(LDU_dmatrix_sp matLDU, dvector vec, dvector res);
void cholinc_dmatrix_sp(dmatrix_sp mat, LDU_dmatrix_sp matLDU, double droptol, double shift);

LDU_dmatrix_sp new_LDU_dmatrix_sp(int n);

typedef struct _dvector_multi {
	int n;
	dvector *elements;
};

typedef struct _dvector_multi *dvector_multi;

dvector_multi new_dvector_multi(int n, int m);
void free_dvector_multi(dvector_multi vecs);

void print_dvector(dvector vec);
void print_dmatrix_sp(dmatrix_sp mat);

#endif
