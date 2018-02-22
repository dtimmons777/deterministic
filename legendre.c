# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>

# include "legendre_polynomial.h"

int main ( );
void p_exponential_product_test ( int p, double b );
void p_integral_test ( );
void p_polynomial_coefficients_test ( );
void p_polynomial_prime_test ( );
void p_polynomial_prime2_test ( );
void p_polynomial_value_test ( );
void p_polynomial_zeros_test ( );
void p_power_product_test ( int p, int e );
void p_quadrature_rule_test ( );
void pm_polynomial_value_test ( );
void pmn_polynomial_value_test ( );
void pmns_polynomial_value_test ( );
void pn_pair_product_test ( int p );
void pn_polynomial_coefficients_test ( );
void pn_polynomial_value_test ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for LEGENDRE_POLYNOMIAL_PRB.

  Discussion:

    LEGENDRE_POLYNOMIAL_PRB tests the LEGENDRE_POLYNOMIAL library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 March 2016

  Author:

    John Burkardt
*/
{
  double b;
  int e;
  int p;

/*  timestamp ( );
  fprintf (f, "\n" );
  fprintf (f, "LEGENDRE_POLYNOMIAL_PRB:\n" );
  fprintf (f, "  C version.\n" );
  fprintf (f, "  Test the LEGENDRE_POLYNOMIAL library.\n" );
*/
  p = 15;
  b = 0.0;
//  p_exponential_product_test ( p, b );

  p = 15;
  b = 1.0;
//  p_exponential_product_test ( p, b );

//  p_integral_test ( );

  p_polynomial_coefficients_test ( );
  p_polynomial_prime_test ( );
//  p_polynomial_prime2_test ( );
//  p_polynomial_value_test ( );
  p_polynomial_zeros_test ( );

  p = 15;
  e = 0;
//  p_power_product_test ( p, e );

  p = 11;
  e = 1;
//  p_power_product_test ( p, e );

//  p_quadrature_rule_test ( );

//  pm_polynomial_value_test ( );
//  pmn_polynomial_value_test ( );
//  pmns_polynomial_value_test ( );

  p = 15;
//  pn_pair_product_test ( p );
//  pn_polynomial_coefficients_test ( );
//  pn_polynomial_value_test ( );
/*
  Terminate.
*/
/*  fprintf (f, "\n" );
  fprintf (f, "LEGENDRE_POLYNOMIAL_PRB:\n" );
  fprintf (f, "  Normal end of execution.\n" );
  fprintf (f, "\n" );
  timestamp ( );*/

  return 0;
}

/******************************************************************************/

void p_polynomial_coefficients_test ( )

/******************************************************************************/

{
  double *c;
  int i;
  int j;
  int n = 15;

FILE *f = fopen("legendre_out.txt", "w");
  fprintf (f, "\n" );
  fprintf (f, "P_POLYNOMIAL_COEFFICIENTS_TEST\n" );
  fprintf (f, "  P_POLYNOMIAL_COEFFICIENTS: coefficients of Legendre polynomial P(n,x).\n" );

  c = p_polynomial_coefficients ( n );

  for ( i = 0; i <= n; i++ )
  {
    fprintf (f, "\n" );
    fprintf (f, "  P(%d,x) =\n", i );
    fprintf (f, "\n" );
    for ( j = i; 0 <= j; j-- )
    {
      if ( c[i+j*(n+1)] == 0.0 )
      {
      }
      else if ( j == 0 )
      {
        fprintf (f, "%14.6g\n", c[i+j*(n+1)] );
      }
      else if ( j == 1 )
      {
        fprintf (f, "%14.6g * x\n", c[i+j*(n+1)] );
      }
      else
      {
        fprintf (f, "%14.6g * x^%d\n", c[i+j*(n+1)], j );
      }
    }
  }
  free ( c );
fclose(f);
  return;
}

/******************************************************************************/

void p_polynomial_prime_test ( )

/******************************************************************************/

{
  double a;
  double b;
  int i;
  int j;
  int m;
  int n;
  double *vp;
  double *x;
FILE *f = fopen("legendre_out.txt", "a");
  fprintf (f, "\n" );
  fprintf (f, "P_POLYNOMIAL_PRIME:\n" );
  fprintf (f, "  P_POLYNOMIAL_PRIME evaluates the derivative of the\n" );
  fprintf (f, "  Legendre polynomial P(n,x).\n" );
  fprintf (f, "\n" );
  fprintf (f, "                        Computed\n" );
  fprintf (f, "     N        X           P'(N,X)\n" );
  fprintf (f, "\n" );

  m = 11;
  a = - 1.0;
  b = + 1.0;
  x = r8vec_linspace_new ( m, a, b );

  n = 15;
  vp = p_polynomial_prime ( m, n, x );

  for ( i = 0; i < m; i++ )
  {
    fprintf (f, "\n" );
    for ( j = 0; j <= n; j++ )
    {
      fprintf (f, "  %4d  %12g  %24.16g\n", j, x[i], vp[i+j*m] );
    }
  }

  free ( vp );
  free ( x );
fclose(f);
  return;
}

/******************************************************************************/

void p_polynomial_zeros_test ( )

/******************************************************************************/

{
  int degree;
  double *hz;
  char title[80];
  double *z;
FILE *f = fopen("legendre_out.txt", "a");
  fprintf (f, "\n" );
  fprintf (f, "P_POLYNOMIAL_ZEROS_TEST:\n" );
  fprintf (f, "  P_POLYNOMIAL_ZEROS computes the zeros of P(n,x)\n" );
  fprintf (f, "  Check by calling P_POLYNOMIAL_VALUE there.\n" );

  for ( degree = 1; degree <= 5; degree++ )
  {
    z = p_polynomial_zeros ( degree );
    sprintf (f, title, "  Computed zeros for P(%d,z):", degree );
    r8vec_print (f, degree, z, title );

    hz = p_polynomial_value ( degree, degree, z );
    fprintf (f, title, "  Evaluate P(%d,z):", degree );
    r8vec_print (f, degree, hz+degree*degree, title );

    free ( hz );
    free ( z );
  }
fclose(f);
  return;
}
