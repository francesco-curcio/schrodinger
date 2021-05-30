#include <iostream>
#include <vector>
#include <fstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_linalg.h>

double Potential(double x)
{
	double mu = 0.5;
	double norm = 1./(sqrt(M_PI * 2.));

	return -norm * gsl_sf_exp(gsl_pow_2(x - mu));
}

void print_matrix(const gsl_matrix_complex * m, size_t n)
{
	for(int i = 0; i < n; ++i)
	{	
		for(int j = 0; j < n; ++j)
		{
			gsl_complex z = gsl_matrix_complex_get(m, i, j);
			printf("%+5.2f %+5.2f*i\t", GSL_REAL(z), GSL_IMAG(z));
		}
		
		std::cout<<std::endl;
	}
}

int main(int argc, char* argv[])
{		
	typedef unsigned long int lu_int;
	
	double lx = 1.; //x will vary from 0 to lx
	double t_tot = 1.; //t will vary from 0 to t_tot
	
	lu_int N_POINTS;
	lu_int N_t;

	if(argc > 2)
	{
		N_POINTS = atoi(argv[1]);
		N_t		 = atoi(argv[2]);
	}
	else
	{
		N_POINTS = 100;
		N_t		 = 100;
	}
	
	//lu_int n = N_POINTS+2;
	
	gsl_vector *x 			= gsl_vector_calloc(N_POINTS);
	gsl_vector *potential 	= gsl_vector_calloc(N_POINTS);
	gsl_vector_complex *psi_n = gsl_vector_complex_calloc(N_POINTS);
	gsl_vector_complex *psi_np1 = gsl_vector_complex_calloc(N_POINTS);
	
	gsl_matrix_complex *lhs_coeffs = gsl_matrix_complex_calloc(N_POINTS, N_POINTS);
	gsl_matrix_complex *rhs_coeffs = gsl_matrix_complex_calloc(N_POINTS, N_POINTS);

	double dt = t_tot/N_t;
	double dx = lx/N_POINTS;
	double a = dt/(2.*dx*dx);

	for (lu_int i = 0; i < N_POINTS; ++i)
	{
		gsl_vector_set(x, i, i*dx);
		gsl_vector_set(potential, i, Potential(i * dx));
		gsl_vector_complex_set(psi_n, i, gsl_complex_rect(1., 0));
		gsl_vector_complex_set(psi_np1, i, gsl_complex_rect(1., 0));
		
		gsl_matrix_complex_set(lhs_coeffs, i, i, gsl_complex_rect(1., 2. * a + 0.5 * dt * gsl_vector_get(potential, i)));
		gsl_matrix_complex_set(rhs_coeffs, i, i, gsl_complex_rect(0., -2. * a - 0.5 * dt * gsl_vector_get(potential, i)));
		
		if(i < N_POINTS-1)
		{
			gsl_matrix_complex_set(lhs_coeffs, i, i+1, gsl_complex_rect(0., -a));
			gsl_matrix_complex_set(lhs_coeffs, i+1, i, gsl_complex_rect(0., -a));
			gsl_matrix_complex_set(rhs_coeffs, i, i+1, gsl_complex_rect(0., a));
			gsl_matrix_complex_set(rhs_coeffs, i+1, i, gsl_complex_rect(0., a));
		}
	}

	gsl_vector_complex_scale(psi_n, gsl_complex_rect(1./gsl_blas_dznrm2(psi_n), 0.));

/*
	std::cout<<"finished setup"<<std::endl;

	for(int i = 1; i < 10; ++i)
	{	
		for(int j = 1; j < 10; ++j)
		{
			gsl_complex z = gsl_matrix_complex_get(coeff_matrix, i, j);
			printf("%+5.2f %+5.2f*i\t", GSL_REAL(z), GSL_IMAG(z));
		}
		
		std::cout<<std::endl;
	}
*/
	//lu_int N_t = t_tot/dt;
	//lu_int N_print = 10;
	
	gsl_permutation *p = gsl_permutation_alloc(N_POINTS);
	int sign = 0;

	//print_matrix(lhs_coeffs, 5);

	
	gsl_linalg_complex_LU_decomp(lhs_coeffs, p, &sign);

	gsl_complex alpha = gsl_complex_rect(1., 0.);
	gsl_complex beta  = gsl_complex_rect(0., 0.);
	gsl_vector_complex *y = gsl_vector_complex_calloc(N_POINTS);

	//print_matrix(lhs_coeffs, 5);
	
	std::ofstream file;
	file.open("res.dat");

	for (lu_int i=0; i < N_t; ++i)
	{
		file<<i*dt<<", ";
		for (lu_int j = 0; j < N_POINTS; ++j)
		{
			file<<gsl_complex_abs2(gsl_vector_complex_get(psi_n, j))<<", ";
		}
		file<<"\n";
		
		gsl_blas_zgemv(CblasNoTrans, alpha, rhs_coeffs, psi_n, beta, y);
			
		gsl_linalg_complex_LU_solve(lhs_coeffs, p, psi_n, psi_np1);

		gsl_vector_complex_swap(psi_n, psi_np1);
	}
	file.close();
	
	return 0;
}
