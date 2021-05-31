#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_linalg.h>

typedef unsigned long int lu_int;

double Potential(double x)
{
	double mu = 7.5;
	double V_min = 1.;
	double V_max = 1.;
	double x_star = mu - sqrt(2. * V_min);
	double width = 1.;
	
	if(x > x_star)
		return 0.5 * gsl_pow_2(x - mu) - V_min;
	else if(x > x_star - width)
		return V_max;
	else 
		return 0.;

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

double psi2(gsl_vector_complex *psi, double dx)
{
	double sum = 0;

	for(lu_int i = 0; i < psi->size; i++)
	{
		sum += gsl_complex_abs2(gsl_vector_complex_get(psi, i)) * dx;
	}

	return sum;
}

int main(int argc, char* argv[])
{		

	//x will vary from 0 to lx
	double lx = 10.;

	//t will vary from 0 to t_tot
	double t_tot = 100.; 
	
	lu_int N_x;
	lu_int N_t;

	if(argc > 2)
	{
		N_x = atoi(argv[1]);
		N_t = atoi(argv[2]);
	}
	else
	{
		N_x = 100;
		N_t = 100;
	}
	
	//lu_int n = N_x+2;
	
	//x grid to initialize potential
	gsl_vector *x 			= gsl_vector_calloc(N_x);

	//potential V to be calculated on the grid	
	gsl_vector *potential 	= gsl_vector_calloc(N_x);

	//wavefunction at step n and n+1
	gsl_vector_complex *psi_n = gsl_vector_complex_calloc(N_x);
	gsl_vector_complex *psi_np1 = gsl_vector_complex_alloc(N_x);
	
	//coefficients of the matrix to be inverted to find the 
	gsl_matrix_complex *lhs_coeffs = gsl_matrix_complex_calloc(N_x, N_x);
	gsl_matrix_complex *rhs_coeffs = gsl_matrix_complex_calloc(N_x, N_x);

	double dt = t_tot/N_t;
	double dx = lx/N_x;
	double a = dt/(2.*dx*dx);

	for (lu_int i = 0; i < N_x; ++i)
	{
		gsl_vector_set(x, i, i*dx);
		gsl_vector_set(potential, i, Potential(i * dx));
		gsl_vector_complex_set(psi_n, i, gsl_complex_rect(gsl_sf_exp(-gsl_pow_2(i*dx - 1.)), 0));
		
		gsl_matrix_complex_set(lhs_coeffs, i, i, gsl_complex_rect(1., 2. * a + 0.5 * dt * gsl_vector_get(potential, i)));
		gsl_matrix_complex_set(rhs_coeffs, i, i, gsl_complex_rect(0., -2. * a - 0.5 * dt * gsl_vector_get(potential, i)));
		
		if(i < N_x-1)
		{
			gsl_matrix_complex_set(lhs_coeffs, i, i+1, gsl_complex_rect(0., -a));
			gsl_matrix_complex_set(lhs_coeffs, i+1, i, gsl_complex_rect(0., -a));
			gsl_matrix_complex_set(rhs_coeffs, i, i+1, gsl_complex_rect(0., a));
			gsl_matrix_complex_set(rhs_coeffs, i+1, i, gsl_complex_rect(0., a));
		}
	}

	gsl_permutation *p = gsl_permutation_alloc(N_x);
	int sign = 0;

	//print_matrix(lhs_coeffs, 5);

	
	gsl_linalg_complex_LU_decomp(lhs_coeffs, p, &sign);

	gsl_complex alpha = gsl_complex_rect(1., 0.);
	gsl_complex beta  = gsl_complex_rect(0., 0.);
	gsl_vector_complex *y = gsl_vector_complex_calloc(N_x);

	//print_matrix(lhs_coeffs, 5);
	
	std::ofstream file;
	file.open("res.dat");

	for (lu_int i=0; i < N_t; ++i)
	{
		
		//normalize psi at every step to ensure consistent normalization
		gsl_vector_complex_scale(psi_n, gsl_complex_rect(1./sqrt(psi2(psi_n, dx)), 0.));
		
		file<<i*dt<<", ";
		for (lu_int j = 0; j < N_x; ++j)
		{
			file<<gsl_complex_abs2(gsl_vector_complex_get(psi_n, j))<<", ";
		}
		file<<"\n";
		
		gsl_blas_zgemv(CblasNoTrans, alpha, rhs_coeffs, psi_n, beta, y);
			
		gsl_linalg_complex_LU_solve(lhs_coeffs, p, psi_n, psi_np1);

		gsl_vector_complex_swap(psi_n, psi_np1);

	}
	file.close();

	file.open("pot.dat");
	for (lu_int i = 0; i < N_x; ++i)
	{
		file<<i*dx<<", "<< gsl_vector_get(potential, i)<<"\n";
	}
	file.close();
	
	return 0;
}
