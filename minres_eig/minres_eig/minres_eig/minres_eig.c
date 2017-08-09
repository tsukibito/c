#pragma comment(lib, "mkl_core.lib")
#pragma comment(lib, "mkl_intel_c.lib")
#pragma comment(lib, "mkl_intel_thread.lib")
#pragma comment(lib, "libguide40.lib")

#define HELM_LEN		50
#define HELM_N			HELM_LEN*HELM_LEN

#define POSITIVE		1
#define NON_POSITIVE	0

#define EPS				1e-16
#define BI_ITR_MAX		500

#ifndef TRUE
#define TRUE			1
#endif

#ifndef FALSE
#define FALSE			0
#endif

#define PI				3.1415926535897932384626433832795028841971 

#define	MINRES
#define PREC

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mkl_types.h>
#include <mkl_spblas.h>
#include <mkl_cblas.h>
#include <mkl_lapack.h>

#include "../../lib/matrix_vector.h"

typedef struct {
	int n;
	vector value;
	ivector row;
	ivector column;
}smatrix;

void free_smatrix(smatrix S) {
	free_vector(S.value);
	free_ivector(S.row);
	free_ivector(S.column);
}

typedef struct {
	int n;
	smatrix A;
	vector b;
	smatrix L;
	vector d;
	double shift;
	double max_eig;
	double min_eig;
}input_datas;

typedef struct {
	int n;
	vector x;
	int itr;
	int sign;
}output_datas;

typedef struct {
	matrix E;
	vector d;
	int cnt;
}eig_datas;

double helm(int w, int h, int x, int y);
void helmholtz(int x, int y, vector A, double *min_eig, double *max_eig, vector trueEig);
void cholinc(int n, int m, vector L, vector d);

void bisection(input_datas *input, output_datas *output, eig_datas *eig, int max_cnt);
void minres(input_datas *input, output_datas *output, eig_datas *eig, int itr_max, double tol, int isBreak);
void pcg(input_datas *input, output_datas *output, eig_datas *eig, int itr_max, double tol, int isBreak);
void subspace(input_datas *input, output_datas *output, eig_datas *eig, int m, int itr_max, double tol);

double reyleigh(int n, smatrix A, vector x);

void mkl_LDL_prod(input_datas *input, vector x, vector y);

void dense2sparse_sym(int n, int m, matrix A, smatrix *S);
void dense2sparse_tri(int n, int m, vector A, smatrix *S);
double sp_get(smatrix S, int x, int y);
void smat_copy(smatrix A, smatrix *B);

static char uplo = 'l';

FILE *fp_lambda;
FILE *fp_bi;

int main(void) {
	input_datas input;
	output_datas output;
	eig_datas eig;
	vector A, L, T;
	int max_cnt, i, j, len, n;
	char fileName[50];
	time_t time;

	printf("nx, ny = ");
	scanf("%d", &len);
	n = len*len;

#ifdef MINRES
	sprintf(fileName, "result\\minres_%d_lambda.csv", len);
	fp_lambda = fopen(fileName, "w");
	sprintf(fileName, "result\\minres_%d_bi.csv", len);
	fp_bi = fopen(fileName, "w");
#else
	sprintf(fileName, "result\\pcg_%d_lambda.csv", len);
	fp_lambda = fopen(fileName, "w");
	sprintf(fileName, "result\\pcg_%d_bi.csv", len);
	fp_bi = fopen(fileName, "w");
#endif

	T = new_vector(100);

	input.n = n;
	//input.A = new_matrix(HELM_N, HELM_N);
	input.b = new_vector(n);
	//input.L = new_matrix(HELM_N, HELM_N);
	input.d = new_vector(n);
	input.shift = 0.0;

	A = new_vector((1+n)*n/2+1);
	helmholtz(len, len, A, &input.min_eig, &input.max_eig, T);
	dense2sparse_tri(input.n, input.n, A, &input.A);
	free_vector(A);

	L = new_vector((1+n)*n/2+1);
	cholinc(len, len, L, input.d);
	dense2sparse_tri(input.n, input.n, L, &input.L);
	free_vector(L);

	/*
	A = new_vector(n+n*(n-1)/2);
	for (i = 1; i <= n; i++) {
		for (j = 1; j <= i; j++) {
			A[i+(2*n-j)*(j-1)/2-1] = helm(len, len, j-1, i-1);
		}
	}
	*/

	vec_init(input.n, input.b, 1.0);

	output.n = input.n;
	output.x = new_vector(input.n);

	printf("求める固有値数: ");
	scanf("%d", &max_cnt);

	eig.E = new_matrix(max_cnt,input.n);
	eig.d = new_vector(max_cnt);
	eig.cnt = 0;

	time = clock();
	bisection(&input, &output, &eig, max_cnt);
	time = clock() - time;

	fprintf(fp_bi, "%f\n", time / (double)CLOCKS_PER_SEC);
	printf("time = %f\n", time / (double)CLOCKS_PER_SEC);

	for (i = 0; i < max_cnt; i++) {
		fprintf(fp_lambda, "%.16f,%.16f,%.16f\n", eig.d[i], T[i], fabs(eig.d[i]-T[i]));
	}

	free_smatrix(input.A); free_smatrix(input.L); free_vector(input.b); free_vector(input.d);
	free_vector(output.x);
	free_matrix(eig.E); free_vector(eig.d);

	fclose(fp_lambda);
	fclose(fp_bi);

	return 0;
}

double helm(int w, int h, int x, int y) {
	double dx = 1.0 / (double)(w+1);
	double dy = 1.0 / (double)(h+1);

	if (x == y) {
		return 2.0/(dx*dx) + 2.0/(dy*dy);
	}
	else if (abs(x-y) == 1) {
		if ((x+1)%w == 0) {
			return 0.0;
		}
		else {
			return -1.0/(dy*dy);
		}
	}
	else if (abs(x-y) == h) {
		return -1.0/(dx*dx);
	}
	else {
		return 0.0;
	}
}

void helmholtz(int x, int y, vector A, double *min_eig, double *max_eig, vector trueEig) {
	int i, j, l, ll, ii;
	double sum;
	double dx = 1.0 / (double)(x+1);
	double dy = 1.0 / (double)(y+1);

	for (i = 1; i <= x; i++) {
		for (j = 1; j <= y; j++) {
			l = (i-1)*y + j - 1;
			ll = (1+l)*l/2;
			A[ll+l+1] = 2.0/(dx*dx) + 2.0/(dy*dy);

			if (i < x) {
				ll = (l+y+1)*(l+y)/2;
				A[ll+l+1] = -1.0/(dx*dx);
			}
			if (j < y) {
				ll = (l+2)*(l+1)/2;
				A[ll+l+1] = -1.0/(dy*dy);
			}
		}
	}

	for (i = 0; i < x*x; i++) {
		ii = (1+i)*i/2;
		sum = 0;
		for (j = 0; j < i; j++) {
			if (i != j) sum += 2*fabs(A[ii+j+1]);
		}

		if (i == 0) {
			*min_eig = A[ii+i+1] - sum;
			*max_eig = A[ii+i+1] + sum;
		}
		else {
			if (*min_eig > A[ii+i+1] - sum)	*min_eig = A[ii+i+1] - sum;
			if (*max_eig < A[ii+i+1] + sum)	*max_eig = A[ii+i+1] + sum;
		}
	}

	for (i = 1, l = 0; i <= 10; i++) {
		for (j = 1; j <= 10; j++) {
			trueEig[l] = 2.0/(dx*dx)*(1-cos(i*PI/(x+1))) + 2.0/(dy*dy)*(1-cos(j*PI/(y+1)));
			l++;
		}
	}

	vec_sort(100, trueEig, trueEig);

	/*int i, j, l, m;
	double dx = 1.0 / (double)(x+1);
	double dy = 1.0 / (double)(y+1);

	for (i = 1; i <= x; i++) {
		for (j = 1; j <= y; j++) {
			l = (i-1)*y + j - 1;
			A[l][l] = 2.0/(dx*dx) + 2.0/(dy*dy);
			if (i < x) {
				A[l][l+y] = -1.0/(dx*dx);
				A[l+y][l] = A[l][l+y];
			}
			if (j < y) {
				A[l][l+1] = -1.0/(dy*dy);
				A[l+1][l] = A[l][l+1];
			}
		}
	}*/
}

void cholinc(int n, int m, vector L, vector d) {
	int i, j, k, ii, jj;
	scalar sum;

	for (i = 0; i < n*n; i++) {
		ii = (1+i)*i/2;
		L[ii+i+1] = 1;

		for (j = 0, sum = 0; j < i; j++) {
			sum -= d[j]*L[ii+j+1]*L[ii+j+1];
		}
		d[i] = helm(n, m, i, i) + sum;

		for (j = i+1; j < m*m; j++) {
			jj = (1+j)*j/2;

			if (helm(n, m, j, i) != 0.0) {
				for (k = 0, sum = 0; k < i; k++) {
					sum -= L[ii+k+1]*d[k]*L[jj+k+1];
				}

				L[jj+i+1] = (helm(n, m, j, i) + sum) / d[i];
			}
		}
	}

	/*for (i = 0; i < n; i++) {
		L[i][i] = 1.0;

		for (j = 0, sum = 0; j < i; j++) {
			sum -= d[j]*L[i][j]*L[i][j];
		}
		d[i] = mat[i][i] + sum;

		for (j = i+1; j < n; j++) {
			if (mat[j][i] != 0.0) {
				for (k = 0, sum = 0; k < i; k++) {
					sum -= L[i][k]*d[k]*L[j][k];
				}

				L[j][i] = (mat[j][i] + sum) / d[i];
			}
		}
	}*/
}

void bisection(input_datas *input, output_datas *output, eig_datas *eig, int max_cnt) {
	int cnt, i, j, isBreak, info, col = 1, n = input->n+input->n*(input->n-1)/2, itr;
	double min_eig, max_eig;
	double lambda, pre_lambda, dot, rho, beta, fd, pre_fd;
	double left, right, center, pre_center;
	vector u, Ax, r, x, Ap;
	ivector ipiv;
	input_datas input_minres;

	u = new_vector(input->n);
	x = new_vector(input->n);
	r = new_vector(input->n);
	Ax = new_vector(input->n);

	//Ap = new_vector(n);
	//ipiv = new_ivector(input->n);

	input_minres.n = input->n;
	input_minres.A = input->A;
	input_minres.L = input->L;
	input_minres.d = input->d;
	input_minres.b = u;

	min_eig = input->min_eig; max_eig = input->max_eig;

	left = 0;
	right = max_eig;

	itr = 0;
	for (cnt = 0; cnt < max_cnt; cnt++) {
		if (cnt > 0) {
			left = lambda;
			right = max_eig;
		}

		vec_rand(input_minres.n, u, -1.0, 1.0);
		vec_normal(input_minres.n, u, u);

		lambda = reyleigh(input_minres.n, input_minres.A, u);

		center = (left + right) / 2.0;

		while(1) {
			isBreak = TRUE;

			for (i = 1; i <= BI_ITR_MAX; i++) {
				input_minres.shift = center;
#ifdef MINRES
				minres(&input_minres, output, eig, BI_ITR_MAX, 1e-4, TRUE);
#else
				pcg(&input_minres, output, eig, BI_ITR_MAX, 1e-4, TRUE);
#endif

				pre_center = center;
				if (output->sign) {
					vec_normal(output->n, output->x, u);
					left = center;
				}
				else {
					if (fabs(vec_sum(output->n, output->x)) > EPS) {
						vec_normal(output->n, output->x, u);
					}
					right = center;
				}

				center = (right + left) / (2.0 - 0.0005);
				//lambda = center;

				fprintf(fp_bi, "%d,%d,%d,%.16f,%.16f,%.16f\n", i, output->sign, output->itr, left, center, right);
				//printf("[%d] %d,%d | %f < %f < %f\n", i, output->sign, output->itr, left, center, right);
				itr += output->itr;

				if (fabs((center-pre_center) / center) < 1e-6) break;
			}

			if (isBreak) break;
		}

		fprintf(fp_bi,"\n");

		input_minres.shift = center;
		pre_lambda = center;
		
		fd = 1;

		for (i = 1; i <= 10; i++) {
			vec_copy(input_minres.n, u, x);
			/*vec_copy(n, A, Ap);

			for (j = 1; j <= input->n; j++) {
				Ap[j+(2*input->n-j)*(j-1)/2-1] -= input_minres.shift;
			}

			dsptrf(&uplo, &input->n, Ap, ipiv, &info);
			dsptrs(&uplo, &input->n, &col, Ap, ipiv, output->x, &input->n, &info);*/
			
//#ifdef MINRES
			minres(&input_minres, output, eig, 150, 1e-8, FALSE);
			itr += output->itr;
//#else
//			pcg(&input_minres, output, eig, 120, 1e-10, FALSE);
//#endif
			for (j = 0; j < eig->cnt; j++) {
				dot = vec_dot(output->n, eig->E[j], output->x);
				vec_sc_prod(input->n, eig->E[j], dot, r);
				vec_sub(input_minres.n, output->x, r, output->x);
			}
			
			vec_normal(output->n, output->x, u);
			lambda = reyleigh(input->n, input->A, u);

			mkl_cspblas_dcoosymv(&uplo, &input->n, input->A.value, input->A.row, input->A.column, &input->A.n, output->x, Ax);
			vec_sc_prod(input->n, output->x, lambda, r);
			vec_sub(input->n, Ax, r, r);

			if (vec_norm(input->n, r) / vec_norm(input->n, Ax) < 1e-10) break;
			if (fabs(lambda-pre_lambda) / lambda < 1e-12) break;

			rho = fabs(lambda - pre_lambda);
			dot = vec_dot(input->n, r, r);
			beta = rho*dot/(rho*rho-dot);

			pre_fd = fd;
			fd = -1.0/vec_dot(input->n, x, u);

			if (fd*pre_fd < 0) {
				input_minres.shift = lambda + fabs(beta); 
			}
			else {
				input_minres.shift = lambda - fabs(beta); 
			}

			pre_lambda = lambda;
		}

		if (lambda + 1e-8 < left || lambda > right + 1e-8) {
			printf("Error![%e <> %e // %f -- %f]\n", left-lambda, lambda-right, left, right);
			//system("pause");
		}

		printf("[%d] lambda = %.12f\n\n", cnt+1, lambda);
		vec_copy(output->n, u, eig->E[cnt]);
		eig->d[cnt] = lambda;
		eig->cnt = cnt+1;
	}

	printf("itr = %d\n", itr);

	free_vector(u); free_vector(Ax); free_vector(r);
	//free_vector(ipiv);
}

void minres(input_datas *input, output_datas *output, eig_datas *eig, int itr_max, double tol, int isBreak) {
	//matrix A;
	vector v0, v1, v2, w0, w1, w2, y, z, tmp;
	double alpha, beta, beta0, delta, eta, gamma, rho1, rho2, rho3, sigma0, sigma1;
	double r, rho_bar, dbar, dot;
	int i, j;

	//A = new_matrix(input->n, input->n);
	//mat_sc_diag_sub(input->n, input->n, input->A, input->shift, A);
	vec_init(output->n, output->x, 0.0);
	
	v0 = new_vector(input->n);
	v1 = new_vector(input->n);
	v2 = new_vector(input->n);
	w0 = new_vector(input->n);
	w1 = new_vector(input->n);
	w2 = new_vector(input->n);
	y = new_vector(input->n);
	z = new_vector(input->n);

	tmp = new_vector(input->n);

	vec_copy(input->n, input->b, v1);
	//LDLr_prod(input->n, input->n, input->L, input->d, v1, z);
#ifdef PREC
		//if (isBreak) {
			vec_copy(input->n, v1, z);
		//}
		//else {
		//	mkl_LDL_prod(input, v1, z);
			//for (i = 0; i < input->n; i++) {
			//	z[i] = v1[i] / helm(50, 50, 1, 1);
			//}
		//}
#else
		vec_copy(input->n, v1, z);
#endif

	beta = sqrt(vec_dot(input->n, v1, z)); beta0 = 1;
	eta = beta;

	gamma = 1; dbar = beta;
	sigma0 = 0; sigma1 = 0;

	r = beta;
	
	output->sign = POSITIVE;
	output->itr = itr_max;

	for (i = 1; i <= itr_max; i++) {
		/* 固有値の抜き取り */
		for (j = 0; j < eig->cnt; j++) {
			dot = vec_dot(input->n, eig->E[j], z);
			vec_sc_prod(input->n, eig->E[j], dot, tmp);
			vec_sub(input->n, z, tmp, z);
		}

		vec_sc_prod(input->n, z, 1.0/beta, y);

		//mat_vec_prod(input->n, input->n, A, y, v2);
		mkl_cspblas_dcoosymv(&uplo, &input->n, input->A.value, input->A.row, input->A.column, &input->A.n, y, v2);

		vec_sc_prod(input->n, y, input->shift, tmp);
		vec_sub(input->n, v2, tmp, v2);

		vec_sc_prod(input->n, v0, beta/beta0, tmp);
		vec_sub(input->n, v2, tmp, v2);

		alpha = vec_dot(input->n, y, v2);

		delta = gamma*alpha - sigma1*dbar;

		if (isBreak) {
			//printf("[%02d] delta = %f\n", i, delta);
		}

		if (delta < 0) {
			output->sign = NON_POSITIVE;
			output->itr = i;

			if (isBreak) break;
		}

		vec_sc_prod(input->n, v1, alpha/beta, tmp);
		vec_sub(input->n, v2, tmp, v2);

		vec_copy(input->n, v1, v0);
		vec_copy(input->n, v2, v1);

		//LDLr_prod(input->n, input->n, input->L, input->d, v1, z);
		//vec_copy(input->n, v1, z);
#ifdef PREC
		if (i <= 5) {
			vec_copy(input->n, v1, z);
		}
		else {
			mkl_LDL_prod(input, v1, z);
			//for (j = 0; j < input->n; j++) {
			//	z[j] = v1[j] / helm(50, 50, 1, 1);
			//}
		}
#else
		vec_copy(input->n, v1, z);
#endif

		beta0 = beta; beta = sqrt(vec_dot(input->n, v1, z));

		rho1 = sqrt(delta*delta + beta*beta);
		rho2 = sigma1*alpha + gamma*dbar;
		rho3 = sigma0*beta0;

		dbar = gamma*beta;

		if (rho1 < EPS) rho1 = EPS;

		rho_bar = 1.0 / rho1;
		gamma = delta * rho_bar;
		sigma0 = sigma1; sigma1 = beta * rho_bar;

		vec_copy(input->n, w1, w0);
		vec_copy(input->n, w2, w1);

		vec_copy(input->n, y, w2);
		
		vec_sc_prod(input->n, w1, rho2, tmp);
		vec_sub(input->n, w2, tmp, w2);
		
		vec_sc_prod(input->n, w0, rho3, tmp);
		vec_sub(input->n, w2, tmp, w2);

		vec_sc_prod(input->n, w2, rho_bar, w2);

		vec_sc_prod(input->n, w2, gamma*eta, tmp);
		vec_add(input->n, output->x, tmp, output->x);

		eta *= -sigma1;
		r *= fabs(sigma1);

		if (!isBreak && i%20 == 0) {
			//printf("[%02d] err = %e\n", i, r);
		}

		if (r < tol) {
			output->itr = i;

			break;
		}
	}
	
	if (!isBreak) {
		printf("[%02d] err = %e\n", i, r);
	}
	
	//free_matrix(A);
	free_vector(v0); free_vector(v1); free_vector(v2);
	free_vector(w0); free_vector(w1); free_vector(w2);
	free_vector(y); free_vector(z); free_vector(tmp);
}

void minres_bak(input_datas *input, output_datas *output, eig_datas *eig, int itr_max, double tol, int isBreak) {
	vector p0, p1, w0, w1, w2, y, z, tmp, Ay;
	double alpha, beta, beta0, delta, eta, gamma0, gamma1, rho1, rho2, rho3, sigma0, sigma1;
	double gbeta, rho_bar, dot;
	int i, j;

	vec_init(output->n, output->x, 0.0);
	
	p0 = new_vector(input->n);
	p1 = new_vector(input->n);
	w0 = new_vector(input->n);
	w1 = new_vector(input->n);
	w2 = new_vector(input->n);
	y = new_vector(input->n);
	z = new_vector(input->n);
	
	Ay = new_vector(input->n);
	tmp = new_vector(input->n);

	mkl_cspblas_dcoosymv(&uplo, &input->n, input->A.value, input->A.row, input->A.column, &input->A.n, output->x, p1);
	vec_sc_prod(input->n, output->x, input->shift, tmp);
	vec_sub(input->n, p1, tmp, p1);
	vec_sub(input->n, input->b, p1, p1);

	mkl_LDL_prod(input, p1, z);

	beta = sqrt(vec_dot(input->n, p1, z)); beta0 = 1;
	eta = beta;

	gamma0 = 1; gamma1 = 1;
	sigma0 = 0; sigma1 = 0;
	
	output->sign = POSITIVE;
	output->itr = -1;

	for (i = 1; i <= itr_max; i++) {
		/* 固有値の抜き取り */
		for (j = 0; j < eig->cnt; j++) {
			dot = vec_dot(input->n, eig->E[j], z);
			vec_sc_prod(input->n, eig->E[j], dot, tmp);
			vec_sub(input->n, z, tmp, z);
		}

		vec_sc_prod(input->n, z, 1.0/beta, y);

		mkl_cspblas_dcoosymv(&uplo, &input->n, input->A.value, input->A.row, input->A.column, &input->A.n, y, Ay);
		vec_sc_prod(input->n, y, input->shift, tmp);
		vec_sub(input->n, Ay, tmp, Ay);

		alpha = vec_dot(input->n, y, Ay);

		vec_copy(input->n, p0, tmp);
		vec_copy(input->n, p1, p0);

		vec_sc_prod(input->n, tmp, beta/beta0, tmp);
		vec_sc_prod(input->n, p1, alpha/beta, p1);
		vec_sub(input->n, Ay, p1, p1);
		vec_sub(input->n, p1, tmp, p1);

		mkl_LDL_prod(input, p1, z);

		gbeta = gamma0*beta;
		delta = gamma1*alpha - sigma1*gbeta;

		if (isBreak) {
			//printf("[%02d] delta = %f\n", i, delta);
		}

		if (delta < 0) {
			output->sign = NON_POSITIVE;
			output->itr = i;

			if (isBreak) break;
		}

		rho2 = sigma1*alpha + gamma1*gbeta;
		rho3 = sigma0*beta;

		beta0 = beta; beta = sqrt(vec_dot(input->n, p1, z));
		rho1 = sqrt(delta*delta + beta*beta);

		if (rho1 < EPS) rho1 = EPS;

		rho_bar = 1.0 / rho1;
		gamma0 = gamma1; gamma1 = delta * rho_bar;
		sigma0 = sigma1; sigma1 = beta * rho_bar;

		vec_copy(input->n, w1, w0);
		vec_copy(input->n, w2, w1);

		vec_copy(input->n, y, w2);		
		vec_sc_prod(input->n, w1, rho2, tmp);
		vec_sub(input->n, w2, tmp, w2);		
		vec_sc_prod(input->n, w0, rho3, tmp);
		vec_sub(input->n, w2, tmp, w2);
		vec_sc_prod(input->n, w2, rho_bar, w2);

		vec_sc_prod(input->n, w2, gamma1*eta, tmp);
		vec_add(input->n, output->x, tmp, output->x);

		eta *= -sigma1;

		if (!isBreak && i%20 == 0) {
			//printf("[%02d] err = %e\n", i, r);
		}

		if (fabs(eta) < tol) {
			output->itr = i;
			break;
		}
	}
	
	free_vector(p0); free_vector(p1); 
	free_vector(w0); free_vector(w1); free_vector(w2);
	free_vector(y); free_vector(z);
	free_vector(tmp); free_vector(Ay);
}

void pcg(input_datas *input, output_datas *output, eig_datas *eig, int itr_max, double tol, int isBreak) {
	int i;
	vector r, p, tmp_vec, tmp_vec2;
	scalar alpha, tmp_sc, pre_rnorm = 1;
	scalar b_norm = vec_norm(input->n, input->b);

	vec_init(input->n, output->x, 0.0);
	r = new_vector(input->n);
	p = new_vector(input->n);

	tmp_vec = new_vector(input->n);
	tmp_vec2 = new_vector(input->n);

	vec_copy(input->n, input->b, r);

	output->sign = POSITIVE;

	for (output->itr = 1; output->itr <= itr_max; output->itr++) {
		/***** Cr/(r,Cr) *****/
		//LDLr_prod(n, n, L, d, r, tmp_vec);
		mkl_LDL_prod(input, r, tmp_vec);

		tmp_sc = vec_dot(input->n, r, tmp_vec);
		vec_sc_prod(input->n, tmp_vec, 1.0/tmp_sc, tmp_vec);
		vec_add(input->n, p, tmp_vec, p);

		/***** 固有値抜き取り *****/
		for (i = 0; i < eig->cnt; i++) {
			tmp_sc = vec_dot(input->n, eig->E[i], p);
			vec_sc_prod(input->n, eig->E[i], tmp_sc, tmp_vec);

			vec_sub(input->n, p, tmp_vec, p);
		}

		/***** alpha *****/
		//mat_vec_prod(n, n, A, p, tmp_vec);
		mkl_cspblas_dcoosymv(&uplo, &input->n, input->A.value, input->A.row, input->A.column, &input->A.n, p, tmp_vec);
		
		vec_sc_prod(input->n, p, input->shift, tmp_vec2);
		vec_sub(input->n, tmp_vec, tmp_vec2, tmp_vec);

		alpha = 1.0 / vec_dot(input->n, p, tmp_vec);

		if (alpha < 0) {
			output->sign = NON_POSITIVE;

			if (isBreak) break;
		}

		/***** r *****/
		vec_sc_prod(input->n, tmp_vec, alpha, tmp_vec);
		vec_sub(input->n, r, tmp_vec, r);

		/***** x *****/
		vec_sc_prod(input->n, p, alpha, tmp_vec);
		vec_add(input->n, output->x, tmp_vec, output->x);
		
		tmp_sc = vec_norm(input->n, r) / b_norm;
		if (tmp_sc < tol) {
			break;
		}

		if (output->itr%20 == 0 && !isBreak) {
			printf("[%02d] err = %e\n", output->itr, tmp_sc);
		}

		pre_rnorm = tmp_sc;
	}

	free_vector(r); free_vector(p);
	free_vector(tmp_vec); free_vector(tmp_vec2);
}

/*
void subspace(input_datas *input, output_datas *output, eig_datas *eig, int m, int itr_max, double tol) {
	int i, j, k;
	matrix X = rand_matrix(m, input->n, -1, 1);
	vector Y = new_vector(m * input->n);
	vector Y1 = new_vector(m * input->n);
	vector y = new_vector(input->n);
	matrix K = new_matrix(m, m);
	matrix M = new_matrix(m, m);

	char transa, matdescra[6] = {'s'};
	double alpha = 1, beta = 0;

	for (i = 0; i < m; i++) {
		vec_normal(input->n, X[i], X[i]);
	}

	for (i = 1; i <= itr_max; i++) {
		transa = 'n';
		for (j = 0; j < m; j++) {
			mkl_dcoosv(&transa, &input->n, &alpha, matdescra, input->A.value, input->A.row, input->A.column, &input->A.n, X[j], y);
			for (k = 0; k < input->n; k++) {
				Y[j*m + k] = y[k];
			}
		}

		mkl_dcoomm(&transa, &input->n, &m, &input->n, &alpha, matdescra, input->A.value, input->A.row, input->A.column, &input->A.n, Y, &m, &beta, Y1, &input->n);

		cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, input->n, m, alpha, Y, input->n, Y1, input->n, beta, K, m);
		cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, input->n, m, alpha, Y, input->n, Y, input->n, beta, M, m);
	}
}
*/

double reyleigh(int n, smatrix A, vector x) {
	double lambda;
	vector tmp = new_vector(n);
	
	//mat_vec_prod(n, n, A, x, tmp);
	mkl_cspblas_dcoosymv(&uplo, &n, A.value, A.row, A.column, &A.n, x, tmp);

	lambda = vec_dot(n, x, tmp);

	free_vector(tmp);

	return lambda;
}

void mkl_LDL_prod(input_datas *input, vector x, vector y) {
	int i;
	char transa = 'n', diag = 'u';

	mkl_cspblas_dcootrsv(&uplo, &transa, &diag, &input->n, input->L.value, input->L.row, input->L.column, &input->L.n, x, y);

	for (i = 0; i < input->n; i++) {
		y[i] /= input->d[i];
	}

	transa = 't';
	mkl_cspblas_dcootrsv(&uplo, &transa, &diag, &input->n, input->L.value, input->L.row, input->L.column, &input->L.n, y, y);
}

void dense2sparse_sym(int n, int m, matrix A, smatrix *S) {
	int i, j, cnt;

	for (i = 0, cnt = 0; i < n; i++) {
		for (j = 0; j <= i; j++) {
			if (A[i][j] != 0.0) cnt++;
		}
	}

	S->n = cnt;
	S->value = new_vector(cnt);
	S->row = new_ivector(cnt);
	S->column = new_ivector(cnt);

	for (i = 0, cnt = 0; i < n; i++) {
		for (j = 0; j <= i; j++) {
			if (A[i][j] != 0.0) {
				S->value[cnt] = A[i][j];
				S->row[cnt] = i;
				S->column[cnt] = j;
				cnt++;
			}
		}
	}
}

void dense2sparse_tri(int n, int m, vector A, smatrix *S){
	int i, j, cnt, ii;

	for (i = 0, cnt = 0; i < n; i++) {
		ii = (1+i)*i/2;
		for (j = 0; j <= i; j++) {
			if (A[ii+j+1] != 0.0) cnt++;
		}
	}

	S->n = cnt;
	S->value = new_vector(cnt);
	S->row = new_ivector(cnt);
	S->column = new_ivector(cnt);

	cnt = 0;

	for (i = 0; i < n; i++) {
		ii = (1+i)*i/2;
		S->value[cnt] = A[ii+i+1];
		S->row[cnt] = i;
		S->column[cnt] = i;
		cnt++;
	}

	for (i = 0; i < n; i++) {
		ii = (1+i)*i/2;
		for (j = 0; j < i; j++) {
			if (A[ii+j+1] != 0.0) {
				S->value[cnt] = A[ii+j+1];
				S->row[cnt] = i;
				S->column[cnt] = j;
				cnt++;
			}
		}
	}
}

scalar sp_get(smatrix S, int x, int y) {
	int i;
	scalar value = 0.0;

	for (i = 0; i < S.n; i++) {
		if (S.row[i] == y && S.column[i] == x) {
			value = S.value[i];
			break;
		}
	}

	return value;
}

void smat_copy(smatrix A, smatrix *B) {
	int i;

	B->n = A.n;
	B->row = new_ivector(A.n);
	B->column = new_ivector(A.n);
	B->value = new_vector(A.n);

	for (i = 0; i < A.n; i++) {
		B->row[i] = A.row[i];
		B->column[i] = A.column[i];
		B->value[i] = A.value[i];
	}
}
