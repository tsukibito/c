#include <stdio.h>
#include "sparse_matrix_vector.h"

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

#define PI 3.1415926535897932386

#define PREC
#define PROP1
#define PROP2
#define MINRES

typedef struct {
	dmatrix_sp A;
	dvector b;
	LDU_dmatrix_sp M;
	double shift;
	int precon;
} input_datas;

typedef struct {
	dvector x;
	int itr;
	int sign;
} output_datas;

typedef struct {
	dvector_multi E;
	dvector d;
	int cnt;
} eig_datas;

void helmholtz(dmatrix_sp mat, double lx, double ly, int nx, int ny, dvector eigs);
void gershgorin(dmatrix_sp A, double *min, double *max);
double reyleigh(dmatrix_sp A, dvector x);

void bisection(input_datas *input, output_datas *output, eig_datas *eig, int max_cnt);
void minres(input_datas *input, output_datas *output, eig_datas *eig, int itr_max, double tol, int isBreak);
void pcg(input_datas *input, output_datas *output, eig_datas *eig, int itr_max, double tol, int isBreak);

FILE *fp_lambda;
FILE *fp_bi;

int main() {
	int n, cnt, i;
	time_t time;
	input_datas input;
	output_datas output;
	eig_datas eig;
	dvector eigs;
	char fileName[100];

	randmize();

	printf("サイズ: "); scanf("%d", &n);
	printf("求める個数: "); scanf("%d", &cnt);

	system("mkdir result");
#ifdef MINRES
	sprintf(fileName, "result\\minres_%d_lambda.csv", n);
	fp_lambda = fopen(fileName, "w");
	sprintf(fileName, "result\\minres_%d_bi.csv", n);
	fp_bi = fopen(fileName, "w");
#else
	sprintf(fileName, "result\\pcg_%d_lambda.csv", n);
	fp_lambda = fopen(fileName, "w");
	sprintf(fileName, "result\\pcg_%d_bi.csv", n);
	fp_bi = fopen(fileName, "w");
#endif

	input.A = new_dmatrix_sp(n*n);
	input.M = new_LDU_dmatrix_sp(n*n);
	input.b = new_dvector(n*n);
	input.shift = 0;

	output.x = new_dvector(n*n);

	eigs = new_dvector(cnt);
	
	eig.E = new_dvector_multi(10, n*n);
	eig.d = new_dvector(cnt);
	eig.cnt = 0;

	helmholtz(input.A, 1.0, 1.0, n, n, eigs);
	cholinc_dmatrix_sp(input.A, input.M, 1e-4, 0);

	clear_dvector(input.b, 1.0);

	time = clock();
	bisection(&input, &output, &eig, cnt);
	printf("time = %f\n", (double)(clock() - time) / 1000.0);
	fprintf(fp_bi, "%f\n", (double)(clock() - time) / 1000.0);

	for (i = 0; i < cnt; i++) {
		fprintf(fp_lambda, "%.16f,%.16f,%.16f\n", eig.d->elements[i], eigs->elements[i], fabs(eigs->elements[i] - eig.d->elements[i]));
	}
	print_dvector(eigs);

	fclose(fp_lambda);
	fclose(fp_bi);

	free_dvector_multi(eig.E);

	return 0;
}

void helmholtz(dmatrix_sp mat, double lx, double ly, int nx, int ny, dvector eigs) {
	int i, j, l;
	double dx, dy, value;
	dvector lambda;

	dx = lx / (nx + 1.0);
	dy = ly / (ny + 1.0);

	for (i = 0; i < ny; i++) {
		for (j = 0; j < nx; j++) {
			l = j*ny + i;

			value = 2.0 / (dx*dx) + 2.0 / (dy*dy);
			set_dmatrix_sp_element(mat, l, l, value, NULL);

			if (j < nx-1) {
				value = -1.0 / (dx*dx);
				set_dmatrix_sp_element(mat, l, l+ny, value, NULL);
				set_dmatrix_sp_element(mat, l+ny, l, value, NULL);
			}
			if (i < ny-1) {
				value = -1.0 / (dy*dy);
				set_dmatrix_sp_element(mat, l, l+1, value, NULL);
				set_dmatrix_sp_element(mat, l+1, l, value, NULL);
			}
		}
	}

	lambda = new_dvector(eigs->n * eigs->n);
	for (i = 0; i < eigs->n; i++) {
		for (j = 0; j < eigs->n; j++) {
			lambda->elements[i*eigs->n + j] = 2.0/(dx*dx) * (1.0 - cos((i+1)*PI/(nx+1))) + 2.0/(dy*dy) * (1.0 - cos((j+1)*PI/(ny+1)));
		}
	}

	sort_dvector(lambda, lambda);
	for (i = 0; i < eigs->n; i++) {
		eigs->elements[i] = lambda->elements[i];
	}

	free_dvector(lambda);
}

void gershgorin(dmatrix_sp A, double *min, double *max) {
	int i;
	double sum, value;
	dnode_ptr next;

	for (i = 0; i < A->n; i++) {
		sum = 0;
		for (next = A->elements[i].first; next != NULL; next = next->next) {
			if (next != A->elements[i].diag) sum += fabs(next->value);
		}

		value = A->elements[i].diag->value;
		if (i == 0) {
			*min = value - sum;
			*max = value + sum;
		}
		else {
			if (*min > value - sum) *min = value - sum;
			if (*max < value + sum) *max = value + sum;
		}
	}
}

double reyleigh(dmatrix_sp A, dvector x) {
	double dot;
	dvector tmp;
	tmp = new_dvector(A->n);

	prod_dmatrix_sp(A, x, tmp);
	dot = dot_dvector(x, tmp);

	free_dvector(tmp);

	return dot;
}

void bisection(input_datas *input, output_datas *output, eig_datas *eig, int max_cnt) {
	int cnt, i, j, isBreak, n = input->A->n, itr, itrs, itr_minres, pre_sign, pre_prop;
	double min_eig, max_eig;
	double lambda, pre_lambda, dot, rho, beta, fd, pre_fd;
	double left, right, center, pre_center, d_center;
	dvector x, Ax, r;

	r = new_dvector(n); x = new_dvector(n); Ax = new_dvector(n);

	gershgorin(input->A, &min_eig, &max_eig);

	left = 0;
	right = max_eig;

	itrs = 0;

	for (cnt = 0; cnt < max_cnt; cnt++) {
		if (cnt > 0) {
			left = lambda;
			right = max_eig;
		}

		rand_dvector(input->b, 0.0, 1.0);
		normal_dvector(input->b, input->b);

		lambda = reyleigh(input->A, input->b);

		center = (left + right) / 2.0;

		while(1) {
			isBreak = TRUE;

			itr_minres = 0;
			d_center = 1;
			output->sign = NON_POSITIVE;
			pre_prop = 0;
			for (i = 1; i <= BI_ITR_MAX; i++) {
				pre_sign = output->sign;
				output->sign = POSITIVE;
				input->shift = center;
				itr = 0;
#ifdef PREC
#ifdef PROP2
				if (pre_sign || pre_prop) {
					input->precon = 0;
#ifdef MINRES
					minres(input, output, eig, 1, 1e-5, TRUE);
#else
					pcg(input, output, eig, 10, 1e-5, TRUE);
#endif
					/*
					if (!output->sign && fabs(sum_dvector(output->x)) > EPS) {
						//normal_dvector(output->x, &input->b);
					}
					*/

					itr = output->itr;
					pre_prop = 1;
				}
				else {
					pre_prop = 0;
				}
#endif
				if (output->sign) {
					input->precon = 1;
#ifdef MINRES
					minres(input, output, eig, BI_ITR_MAX, 1e-5, TRUE);
#else
					pcg(input, output, eig, BI_ITR_MAX, 1e-5, TRUE);
#endif
					output->itr += itr;
				}
				itr_minres += output->itr;
#else
				input->precon = 0;
				minres(input, output, eig, BI_ITR_MAX, 1e-5, TRUE);
#endif			
				
				pre_center = center;
				if (output->sign) {
					normal_dvector(output->x, input->b);
					left = center;
				}
				else {
					if (fabs(sum_dvector(output->x)) > EPS) {
						normal_dvector(output->x, input->b);
					}
					right = center;
				}

#ifdef PROP1
				/*
				if (d_center < 1e-3) {
					center = (right + left) / (2.0 + 1e-4);
					if (center < left) {
						center = (right + left) / 2;
					}
				}
				else {*/
					center = (right + left) / (2.0 - 0*1e-3);
					if (center > right) {
						center = (right + left) / 2;
					}
				//}
#else
				center = (right + left) / 2.0;
#endif

				fprintf(fp_bi, "%d,%d,%d,%.16f,%.16f,%.16f\n", i, output->sign, output->itr, left, center, right);
				//printf("[%d] %d,%d | %f < %f < %f\n", i, output->sign, output->itr, left, center, right);
				itrs += output->itr;

				d_center = fabs(center-pre_center) / center;
				if (d_center < 1e-6) break;
			}

			if (isBreak) break;
		}

		printf("minres_itr = %d\n", itr_minres);
		fprintf(fp_bi,"\n");
		printf("%20.16f < %20.16f < %20.16f\n", left, center, right);

		/*
		center = (right + left) / (2.0 + 1e-4);
		input->shift = center;
		minres(input, output, eig, 500, 1e-5, TRUE);
		if (output->sign) {
			normal_dvector(output->x, &input->b);
		}
		else {
			if (fabs(sum_dvector(output->x)) > EPS) {
				normal_dvector(output->x, &input->b);
			}
		}
		*/

		input->shift = center;
		pre_lambda = center;

		//fd = 1;

		for (i = 1; i <= 10; i++) {
			//copy_dvector(input->b, &x);
#ifdef PREC
			cholinc_dmatrix_sp(input->A, input->M, 1e-8, input->shift);
			input->precon = 1;
//#ifdef MINRES
			minres(input, output, eig, BI_ITR_MAX, 1e-10, FALSE);
//#else
//			pcg(input, output, eig, BI_ITR_MAX, 1e-10, FALSE);
//#endif
#else
			input->precon = 0;
//#ifdef MINRES
			minres(input, output, eig, BI_ITR_MAX, 1e-10, FALSE);
//#else
//			pcg(input, output, eig, BI_ITR_MAX, 1e-10, FALSE);
//#endif
#endif
			printf("itr = %d\n", output->itr);
			
			for (j = 0; j < eig->cnt; j++) {
				dot = dot_dvector(eig->E->elements[j], output->x);
				prod_scalar_dvector(eig->E->elements[j], dot, r);
				sub_dvector(output->x, r, output->x);
			}
			
			normal_dvector(output->x, input->b);
			lambda = reyleigh(input->A, input->b);

			if (output->sign) {
				left = lambda;
			}
			else {
				right = lambda;
			}

			center = (right + left) / 2;
			input->shift = center;

			//prod_dmatrix_sp(input->A, input->b, &Ax);

			//prod_scalar_dvector(input->b, lambda, &r);
			//sub_dvector(Ax, r, &r);

			//if (norm_dvector(r) / norm_dvector(Ax) < 1e-8) break;
			printf("err = %e\n", fabs(lambda-pre_lambda) / lambda);
			if (fabs(lambda-pre_lambda) / lambda < 1e-10) break;
/*
			rho = fabs(lambda - pre_lambda);
			dot = dot_dvector(r, r);
			beta = rho*dot/(rho*rho-dot);

			pre_fd = fd;
			fd = -1.0/dot_dvector(x, input->b);

			if (fd*pre_fd < 0) {
				input->shift = lambda + fabs(beta); 
			}
			else {
				input->shift = lambda - fabs(beta); 
			}
*/
			pre_lambda = lambda;
		}

		if (lambda < left - 1e-8) {
			printf("Error![%e < %e]\n", lambda, left);
			//system("pause");
		}
		if (lambda > right + 1e-8) {
			printf("Error![%e > %e]\n", lambda, right);
			//system("pause");
		}

		printf("[%d] lambda = %.12f\n\n", cnt+1, lambda);
		copy_elements_dvector(input->b, eig->E->elements[cnt]);
		eig->d->elements[cnt] = lambda;
		eig->cnt = cnt+1;
	}

	fprintf(fp_bi, "%d\n", itrs);
	printf("itrs = %d\n", itrs);

	free_dvector(r); free_dvector(x); free_dvector(Ax); 
}

void minres(input_datas *input, output_datas *output, eig_datas *eig, int itr_max, double tol, int isBreak) {
	int n = input->A->n;

	dvector p0, p1, p2, w0, w1, w2, y, z, tmp;
	double alpha, beta, beta0, delta, eta, gamma, rho1, rho2, rho3, sigma0, sigma1;
	double r, rho_bar, dbar, dot;
	int i, j;

	p0 = new_dvector(n); p1 = new_dvector(n); p2 = new_dvector(n);
	w0 = new_dvector(n); w1 = new_dvector(n); w2 = new_dvector(n);
	z = new_dvector(n); y = new_dvector(n); tmp = new_dvector(n);

	clear_dvector(output->x, 0.0);

	//prod_dmatrix_sp_diag(input->A, output->x, input->shift, p1);
	//sub_dvector(input->b, p1, p1);
	copy_elements_dvector(input->b, p1);

	if (input->precon) {
		LDU_prod_dmatrix_sp(input->M, p1, z);
	}
	else {
		copy_elements_dvector(p1, z);
	}

	beta = sqrt(dot_dvector(p1, z)); beta0 = 1;
	eta = beta;

	gamma = 1; dbar = beta;
	sigma0 = 0; sigma1 = 0;

	r = beta;
	
	output->sign = POSITIVE;
	output->itr = itr_max;

	for (i = 1; i <= itr_max; i++) {
		/* 固有値の抜き取り */
		for (j = 0; j < eig->cnt; j++) {
			dot = dot_dvector(eig->E->elements[j], z);
			prod_scalar_dvector(eig->E->elements[j], dot, tmp);
			sub_dvector(z, tmp, z);
		}

		prod_scalar_dvector(z, 1.0/beta, y);

		prod_dmatrix_sp_diag(input->A, y, input->shift, p2);

		if (i > 1) {
			prod_scalar_dvector(p0, beta/beta0, tmp);
			sub_dvector(p2, tmp, p2);
		}

		alpha = dot_dvector(y, p2);

		prod_scalar_dvector(p1, alpha/beta, tmp);
		sub_dvector(p2, tmp, p2);

		copy_elements_dvector(p1, p0);
		copy_elements_dvector(p2, p1);

		if (input->precon) {
			LDU_prod_dmatrix_sp(input->M, p1, z);
		}
		else {
			copy_elements_dvector(p1, z);
		}

		beta0 = beta; beta = sqrt(dot_dvector(p1, z));

		delta = gamma*alpha - sigma1*dbar;
		rho1 = sqrt(delta*delta + beta*beta);
		rho2 = sigma1*alpha + gamma*dbar;
		rho3 = sigma0*beta0;

		dbar = gamma*beta;

		if (delta < 0) {
			output->sign = NON_POSITIVE;
			output->itr = i;

			if (isBreak) break;
		}
		else if (itr_max == 1) {
			break;
		}

		if (rho1 < EPS) rho1 = EPS;

		rho_bar = 1.0 / rho1;
		gamma = delta * rho_bar;
		sigma0 = sigma1; sigma1 = beta * rho_bar;

		copy_elements_dvector(w1, w0);
		copy_elements_dvector(w2, w1);

		copy_elements_dvector(y, w2);
		
		prod_scalar_dvector(w1, rho2, tmp);
		sub_dvector(w2, tmp, w2);
		
		prod_scalar_dvector(w0, rho3, tmp);
		sub_dvector(w2, tmp, w2);

		prod_scalar_dvector(w2, rho_bar, w2);

		prod_scalar_dvector(w2, gamma*eta, tmp);
		add_dvector(output->x, tmp, output->x);

		eta *= -sigma1;
		r *= fabs(sigma1);

		/*
		if (!isBreak && i%20 == 0) {
			printf("[%02d] err = %e\n", i, r);
		}
		*/

		if (r < tol) {
			output->itr = i;
			break;
		}
	}
	
	free_dvector(p0); free_dvector(p1); free_dvector(p2);
	free_dvector(w0); free_dvector(w1); free_dvector(w2);
	free_dvector(y); free_dvector(z); free_dvector(tmp);
}

void pcg(input_datas *input, output_datas *output, eig_datas *eig, int itr_max, double tol, int isBreak) {
	int i, n = input->A->n;
	dvector r, p, tmp_vec, tmp_vec2;
	double alpha, tmp_sc, pre_rnorm = 1;
	double b_norm = norm_dvector(input->b);

	clear_dvector(output->x, 0.0);
	r = new_dvector(n);
	p = new_dvector(n);

	tmp_vec = new_dvector(n);
	tmp_vec2 = new_dvector(n);

	copy_elements_dvector(input->b, r);

	output->sign = POSITIVE;

	for (output->itr = 1; output->itr <= itr_max; output->itr++) {
		/***** Cr/(r,Cr) *****/
		LDU_prod_dmatrix_sp(input->M, r, tmp_vec);

		tmp_sc = dot_dvector(r, tmp_vec);
		prod_scalar_dvector(tmp_vec, 1.0/tmp_sc, tmp_vec);
		add_dvector(p, tmp_vec, p);

		/***** 固有値抜き取り *****/
		for (i = 0; i < eig->cnt; i++) {
			tmp_sc = dot_dvector(eig->E->elements[i], p);
			prod_scalar_dvector(eig->E->elements[i], tmp_sc, tmp_vec);

			sub_dvector(p, tmp_vec, p);
		}

		/***** alpha *****/
		prod_dmatrix_sp(input->A, p, tmp_vec);
		
		prod_scalar_dvector(p, input->shift, tmp_vec2);
		sub_dvector(tmp_vec, tmp_vec2, tmp_vec);

		alpha = 1.0 / dot_dvector(p, tmp_vec);

		if (alpha < 0) {
			output->sign = NON_POSITIVE;

			if (isBreak) break;
		}

		/***** r *****/
		prod_scalar_dvector(tmp_vec, alpha, tmp_vec);
		sub_dvector(r, tmp_vec, r);

		/***** x *****/
		prod_scalar_dvector(p, alpha, tmp_vec);
		add_dvector(output->x, tmp_vec, output->x);
		
		tmp_sc = norm_dvector(r) / b_norm;
		if (tmp_sc < tol) {
			break;
		}

		/*
		if (output->itr%20 == 0 && !isBreak) {
			printf("[%02d] err = %e\n", output->itr, tmp_sc);
		}
		*/

		pre_rnorm = tmp_sc;
	}

	free_dvector(r); free_dvector(p);
	free_dvector(tmp_vec); free_dvector(tmp_vec2);
}
