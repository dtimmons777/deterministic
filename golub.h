# ifndef golub
# define golub

long double xExponential_Integral_Ei( long double x );
double Exponential_Integral_En( double x, int n );
double *p_polynomial_coefficients ( int n );
double *p_polynomial_zeros ( int nt );
double *p_polynomial_prime ( int m, int n, double x[] );
void p_quadrature_rule ( int nt, double t[], double wts[] );
float determinant(float [][25], float);
void cofactor(float [][25], float);
float transpose(float [][25], float [][25], float);
//static long double Power_Series_En( long double x, int n)
//static long double Continued_Fraction_En( long double x, int n, long double expx );


void function();

#endif 
