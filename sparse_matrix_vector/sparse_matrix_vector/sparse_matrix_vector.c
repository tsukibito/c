#include "sparse_matrix_vector.h"

void _qsort(dvector vec, int left, int right) {
    int i, j;
    double pivot;
	double temp;

    i = left;
    j = right;

	pivot = vec->elements[(left + right) / 2];

    while (1) {
        while (vec->elements[i] < pivot) i++;
        while (pivot < vec->elements[j]) j--;

        if (i >= j) break;

	    temp = vec->elements[i];
	    vec->elements[i] = vec->elements[j];
	    vec->elements[j] = temp;

        i++; j--;
    }

    if (left < i - 1)	_qsort(vec, left, i - 1);
    if (j + 1 <  right)	_qsort(vec, j + 1, right);
}

void randmize() {
	srand((unsigned)time(NULL));
}

double drand(double min, double max) {
	return (max - min) * rand() / (double)RAND_MAX + min;
}

dvector new_dvector(int n) {
	dvector vec;
	int i;

	vec = (dvector)malloc(sizeof(struct _dvector));
	
	vec->n = n;
	vec->elements = (double *)malloc(sizeof(double) * n);

	for (i = 0; i < n; i++) {
		vec->elements[i] = 0.0;
	}

	return vec;
}

void clear_dvector(dvector vec, double sc) {
	int i;
	
	for (i = 0; i < vec->n; i++) {
		vec->elements[i] = sc;
	}
}

void rand_dvector(dvector vec, double min, double max) {
	int i;
	
	for (i = 0; i < vec->n; i++) {
		vec->elements[i] = drand(min, max);
	}
}

void copy_elements_dvector(dvector source, dvector target) {
	int i;

	for (i = 0; i < source->n; i++) {
		target->elements[i] = source->elements[i];
	}
}

dvector copy_dvector(dvector source) {
	dvector vec;
	int i;

	vec = (dvector)malloc(sizeof(struct _dvector));

	vec->n = source->n;
	vec->elements = (double *)malloc(sizeof(double) * source->n);

	for (i = 0; i < source->n; i++) {
		vec->elements[i] = source->elements[i];
	}

	return vec;
}

void free_dvector(dvector vec) {
    free(vec->elements);
}

void add_dvector(dvector vec1, dvector vec2, dvector res) {
	int i;

	for(i = 0; i < res->n; i++) {
		res->elements[i] = vec1->elements[i] + vec2->elements[i];
	}
}

void sub_dvector(dvector vec1, dvector vec2, dvector res) {
	int i;

	for(i = 0; i < res->n; i++) {
		res->elements[i] = vec1->elements[i] - vec2->elements[i];
	}
}

void prod_dvector(dvector vec1, dvector vec2, dvector res) {
	int i;

	for(i = 0; i < res->n; i++) {
		res->elements[i] = vec1->elements[i] * vec2->elements[i];
	}
}

void prod_scalar_dvector(dvector vec, double sc, dvector res) {
	int i;

	for(i = 0; i < res->n; i++) {
		res->elements[i] = sc * vec->elements[i];
	}
}

void neg_dvector(dvector vec, dvector res) {
	int i;

	for(i = 0; i < res->n; i++) {
		res->elements[i] = -vec->elements[i];
	}
}

void rev_dvector(dvector vec, dvector res) {
	int i;

	for(i = 0; i < res->n; i++) {
		res->elements[i] = vec->elements[res->n-i-1];
	}
}

void normal_dvector(dvector vec, dvector res) {
	int i;
	double norm = norm_dvector(vec);

	for(i = 0; i < res->n; i++) {
		res->elements[i] = vec->elements[i] / norm;
	}
}

double sum_dvector(dvector vec) {
	int i;
	double sum = 0.0;

	for (i = 0; i < vec->n; i++) {
		sum += vec->elements[i];
	}

	return sum;
}

double norm_dvector(dvector vec) {
	int i;
	double norm = 0.0;

	for (i = 0; i < vec->n; i++) {
		norm += vec->elements[i] * vec->elements[i];
	}

	return sqrt(norm);
}

double dot_dvector(dvector vec1, dvector vec2) {
	int i;
	double dot = 0.0;

	for (i = 0; i < vec1->n; i++) {
		dot += vec1->elements[i] * vec2->elements[i];
	}

	return dot;
}

void sort_dvector(dvector vec, dvector res) {
	copy_elements_dvector(vec, res);
	_qsort(res, 0, vec->n-1);
}

dmatrix_sp new_dmatrix_sp(int n) {
	dmatrix_sp mat;
	int i;

	mat = (dmatrix_sp)malloc(sizeof(struct _dmatrix_sp));

	mat->n = n;
	mat->elements = (first_dnode_ptr)malloc(sizeof(first_dnode) * n);

	for (i = 0; i < n; i++) {
		mat->elements[i].diag = (dnode_ptr)malloc(sizeof(dnode));
		mat->elements[i].diag->index = i;
		mat->elements[i].diag->value = 0.0;
		mat->elements[i].diag->next = NULL;

		mat->elements[i].first = mat->elements[i].diag;
	}

	return mat;
}

void clear_dmatrix_sp(dmatrix_sp mat, double sc) {
	int i;
	dnode_ptr current, next;

	for (i = 0; i < mat->n; i++) {
		for (current = mat->elements[i].first; current != NULL; current = next) {
			next = current->next;
			free(current);
		}
		
		mat->elements[i].diag = (dnode_ptr)malloc(sizeof(dnode));
		mat->elements[i].diag->index = i;
		mat->elements[i].diag->value = sc;
		mat->elements[i].diag->next = NULL;

		mat->elements[i].first = mat->elements[i].diag;
	}
}

void free_dmatrix_sp(dmatrix_sp mat) {
	int i;
	dnode_ptr current, next;

	for (i = 0; i < mat->n; i++) {
		for (current = mat->elements[i].first; current != NULL; current = next) {
			next = current->next;
			free(current);
		}
	}

	free(mat->elements);
}

dnode_ptr* set_dmatrix_sp_element(dmatrix_sp mat, int row_index, int col_index, double value, dnode_ptr *current) {
	int index;
	dnode_ptr new_node, free_ptr;

	if (current == NULL) {
		current = (col_index >= row_index)? &(mat->elements[row_index].diag): &(mat->elements[row_index].first);
	}

	while(*current != NULL) {
		index = (*current)->index;
		if (index < col_index) {
			current = &(*current)->next;
			continue;
		}
		
		else if (index == col_index) {
			(*current)->value = value;
		}
		else {
			new_node = (dnode_ptr)malloc(sizeof(dnode));
			new_node->index = col_index;
			new_node->value = value;
			new_node->next = *current;

			*current = new_node;
		}
		break;
	}
	if (*current == NULL) {
		*current = (dnode_ptr)malloc(sizeof(dnode));
		(*current)->index = col_index;
		(*current)->value = value;
		(*current)->next = NULL;
	}

	if (col_index != row_index && (*current)->value == 0.0) {
		free_ptr = *current;
		*current = (*current)->next;
		free(free_ptr);
	}

	return current;
}

void trans_dmatrix_sp(dmatrix_sp mat, dmatrix_sp res) {
	int i, index;
	dnode_ptr mat_next, **res_next;

	clear_dmatrix_sp(res, 0.0);
	res_next = (dnode_ptr **)malloc(sizeof(dnode_ptr *) * mat->n);

	for (i = 0; i < mat->n; i++) {
		res_next[i] = NULL;
	}

	for (i = 0; i < mat->n; i++) {
		for (mat_next = mat->elements[i].first; mat_next != NULL; mat_next = mat_next->next) {
			index = mat_next->index;
			res_next[index] = set_dmatrix_sp_element(res, index, i, mat_next->value, res_next[index]);
		}
	}

	free(res_next);
}

void add_dmatrix_sp(dmatrix_sp mat0, dmatrix_sp mat1, dmatrix_sp res) {
	int i, mat0_index, mat1_index, index;
	double value;
	dnode_ptr *next, mat0_next, mat1_next;

	if (mat0 != res && mat1 != res) {
		clear_dmatrix_sp(res, 0.0);
	}

	for (i = 0; i < res->n; i++) {
		next = NULL;
		mat0_next = mat0->elements[i].first;
		mat1_next = mat1->elements[i].first;

		while(mat0_next != NULL || mat1_next != NULL) {
			mat0_index = (mat0_next == NULL)? res->n: mat0_next->index;
			mat1_index = (mat1_next == NULL)? res->n: mat1_next->index;

			if (mat0_index < mat1_index) {
				index = mat0_index;
				value = mat0_next->value;
				
				mat0_next = (mat0_next == NULL)? NULL: mat0_next->next;
			}
			else if (mat0_index > mat1_index) {
				index = mat1_index;
				value = mat1_next->value;

				mat1_next = (mat1_next == NULL)? NULL: mat1_next->next;
			}
			else {
				index = mat0_index;
				value = mat0_next->value + mat1_next->value;

				mat0_next = (mat0_next == NULL)? NULL: mat0_next->next;
				mat1_next = (mat1_next == NULL)? NULL: mat1_next->next;
			}

			next = set_dmatrix_sp_element(res, i, index, value, next);
		}
	}
}

void sub_dmatrix_sp(dmatrix_sp mat0, dmatrix_sp mat1, dmatrix_sp res) {
	int i, mat0_index, mat1_index, index;
	double value;
	dnode_ptr *next, mat0_next, mat1_next;

	if (mat0 != res && mat1 != res) {
		clear_dmatrix_sp(res, 0.0);
	}

	for (i = 0; i < res->n; i++) {
		next = NULL;
		mat0_next = mat0->elements[i].first;
		mat1_next = mat1->elements[i].first;

		while(mat0_next != NULL || mat1_next != NULL) {
			mat0_index = (mat0_next == NULL)? res->n: mat0_next->index;
			mat1_index = (mat1_next == NULL)? res->n: mat1_next->index;

			if (mat0_index < mat1_index) {
				index = mat0_index;
				value = mat0_next->value;
				
				mat0_next = (mat0_next == NULL)? NULL: mat0_next->next;
			}
			else if (mat0_index > mat1_index) {
				index = mat1_index;
				value = -mat1_next->value;

				mat1_next = (mat1_next == NULL)? NULL: mat1_next->next;
			}
			else {
				index = mat0_index;
				value = mat0_next->value - mat1_next->value;

				mat0_next = (mat0_next == NULL)? NULL: mat0_next->next;
				mat1_next = (mat1_next == NULL)? NULL: mat1_next->next;
			}

			next = set_dmatrix_sp_element(res, i, index, value, next);
		}
	}
}

void prod_dmatrix_sp(dmatrix_sp mat, dvector vec, dvector res) {
	int i;
	dvector temp_vec;
	dnode_ptr next;
	
	temp_vec = copy_dvector(vec);
	
	for (i = 0; i < mat->n; i++) {
		res->elements[i] = 0.0;
		
		for (next = mat->elements[i].first; next != NULL; next = next->next) {
			res->elements[i] += next->value * temp_vec->elements[next->index];
		}
	}

	free_dvector(temp_vec);
}

void prod_dmatrix_sp_diag(dmatrix_sp mat, dvector vec, double sc, dvector res) {
	int i;
	dvector temp_vec;
	dnode_ptr next;
	
	temp_vec = copy_dvector(vec);
	
	for (i = 0; i < mat->n; i++) {
		res->elements[i] = 0.0;
		
		for (next = mat->elements[i].first; next != NULL; next = next->next) {
			if (next == mat->elements[i].diag) {
				res->elements[i] += (next->value - sc) * temp_vec->elements[next->index];
			}
			else {
				res->elements[i] += next->value * temp_vec->elements[next->index];
			}
		}
	}

	free_dvector(temp_vec);
}

void forward_dmatrix_sp(dmatrix_sp mat, dvector vec, dvector res) {
	int i;
	double temp;
	dnode_ptr next;
	
	for (i = 0; i < mat->n; i++) {
		temp = vec->elements[i];
		
		for (next = mat->elements[i].first; next->index < i; next = next->next) {
			temp -= next->value * res->elements[next->index];
		}

		res->elements[i] = temp / mat->elements[i].diag->value;
	}
}

void backward_dmatrix_sp(dmatrix_sp mat, dvector vec, dvector res) {
	int i;
	double temp;
	dnode_ptr next;
	
	for (i = mat->n-1; i >= 0; i--) {
		temp = vec->elements[i];
		
		for (next = mat->elements[i].diag->next; next != NULL; next = next->next) {
			temp -= next->value * res->elements[next->index];
		}

		res->elements[i] = temp / mat->elements[i].diag->value;
	}
}

void LDU_prod_dmatrix_sp(LDU_dmatrix_sp matLDU, dvector vec, dvector res) {
	int i;

	forward_dmatrix_sp(matLDU->matL, vec, res);
	
	for (i = 0; i < vec->n; i++) {
		res->elements[i] /= matLDU->vecD->elements[i];
	}

	backward_dmatrix_sp(matLDU->matU, res, res);
}

void cholinc_dmatrix_sp(dmatrix_sp mat, LDU_dmatrix_sp matLDU, double droptol, double shift) {
	int i;
	double temp;
	dnode_ptr next, i_next, j_next, **matL_next, *matU_next;

	matL_next = (dnode_ptr **)malloc(sizeof(dnode_ptr *) * mat->n);

	for (i = 0; i < mat->n; i++) {
		matL_next[i] = NULL;
	}

	clear_dmatrix_sp(matLDU->matL, 0.0);
	clear_dmatrix_sp(matLDU->matU, 0.0);
	clear_dvector(matLDU->vecD, 0.0);

	for (i = 0; i < mat->n; i++) {
		matU_next = NULL;

		for (next = mat->elements[i].diag; next != NULL; next = next->next) {
			i_next = matLDU->matL->elements[i].first;
			j_next = matLDU->matL->elements[next->index].first;

			temp = next->value - shift;
			while (i_next != NULL && j_next != NULL) {
				while (i_next->index < j_next->index) {
					i_next = i_next->next;
					if (i_next == NULL) {
						break;
					}
				}
				if (i_next == NULL) {
					break;
				}

				while (j_next->index < i_next->index) {
					j_next = j_next->next;
					if (j_next == NULL) {
						break;
					}
				}
				if (j_next == NULL) {
					break;
				}

				if (i_next->index >= i || j_next->index >= i) {
					break;
				}

				if (i_next->index == j_next->index) {
					temp -= i_next->value * matLDU->vecD->elements[i_next->index] * j_next->value;

					i_next = i_next->next;
					j_next = j_next->next;
				}
			}
			matL_next[next->index] = set_dmatrix_sp_element(matLDU->matL, next->index, i, temp, matL_next[next->index]);
			matU_next = set_dmatrix_sp_element(matLDU->matU, i, next->index, temp, matU_next);
		}

		if (fabs(matLDU->matL->elements[i].diag->value) < droptol) {
			matLDU->matL->elements[i].diag->value = droptol;
			matLDU->matU->elements[i].diag->value = droptol;
		}
		matLDU->vecD->elements[i] = 1.0 / matLDU->matL->elements[i].diag->value;
	}

	free(matL_next);
}

LDU_dmatrix_sp new_LDU_dmatrix_sp(int n) {
	LDU_dmatrix_sp matLDU;

	matLDU = (LDU_dmatrix_sp)malloc(sizeof(struct _LDU_dmatrix_sp));

	matLDU->matL = new_dmatrix_sp(n);
	matLDU->matU = new_dmatrix_sp(n);
	matLDU->vecD = new_dvector(n);

	return matLDU;
}

dvector_multi new_dvector_multi(int n, int m) {
	dvector_multi vecs;
    int i;

	vecs = (dvector_multi)malloc(sizeof(struct _dvector_multi));

	vecs->n = n;
	vecs->elements = (dvector *)malloc(sizeof(dvector) * n);

	for (i = 0; i < n; i++) {
		vecs->elements[i] = new_dvector(m);
	}

	return vecs;
}

void free_dvector_multi(dvector_multi vecs) {
	int i;

	for (i = 0; i < vecs->n; i++) {
		free_dvector(vecs->elements[i]);
	}
	free(vecs->elements);
}

void print_dmatrix_sp(dmatrix_sp mat) {
	int i;
	dnode_ptr next;

	printf("Matrix Size = %d x %d\n", mat->n, mat->n);
	for (i = 0; i < mat->n; i++) {
		for (next = mat->elements[i].first; next != NULL; next = next->next) {
			printf("(%d, %d) = %20.12e\n", i, next->index, next->value);
		}
	}
	printf("\n");
}

void print_dvector(dvector vec) {
	int i;

	printf("Vector Size = %d\n", vec->n);
	for (i = 0; i < vec->n; i++) {
		printf("(%d) = %20.12e\n", i, vec->elements[i]);
	}
	printf("\n");
}

