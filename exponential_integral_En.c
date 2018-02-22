////////////////////////////////////////////////////////////////////////////////
// File: exponential_integral_En.c                                            //
// Routine(s):                                                                //
//    Exponential_Integral_En                                                 //
//                                                                            //
// Required Externally Defined Routines:                                      //
//    xExponential_Integral_Ei                                                //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                    // required for expl() and logl()
#include <float.h>                   // required for DBL_MAX and LDBL_EPSILON

//                    Required Externally Defined Routines                    //
long double xExponential_Integral_Ei( long double x );

//                        Internally Defined Routines                         //
static long double Power_Series_En( long double x, int n);
static long double Continued_Fraction_En( long double x, int n,
                                                            long double expx );

////////////////////////////////////////////////////////////////////////////////
// double Exponential_Integral_En( double x, int n )                          //
//                                                                            //
//  Description:                                                              //
//     The exponential integral En(x), introduced by Schloemilch, is the      //
//     integral of                                                            //
//                            exp(-xt) / t^n dt                               //
//     where the integral extends from 1 to inf.  The parameter n is called   //
//     the order of the integral.                                             //
//     At x = 0, E0 = inf, E1 = inf and En = 1 / (n - 1) for n >= 2.          //
//     For x > 0, E0 = exp(-x) / x and E1(x) = - Ei(-x), where Ei(x) is the   //
//     exponential integral.  For n >= 2,  En(x) is then given recursively by //
//                 En(x) = ( x / (n - 1) ) [E0(x) - En-1(x)].                 //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of Schloemilch's exponential integral En().    //
//     int     n  The order of Schloemilch's exponential integral.            //
//                                                                            //
//  Return Value:                                                             //
//     The value of Schloemilch's exponential integral En evaluated at x.     //
//     If x < 0, then the integral diverges and DBL_MAX is returned.          //
//     If x = 0, then if n = 0 or n = 1, the integral diverges and DBL_MAX    //
//     is returned and if n >= 2, En(0) = 1 / (n-1).                          //
//                                                                            //
//  Example:                                                                  //
//     double y, x;                                                           //
//     int    n;                                                              //
//                                                                            //
//     ( code to initialize x and n )                                         //
//                                                                            //
//     y = Exponential_Integral_En( x, n);                                    //
////////////////////////////////////////////////////////////////////////////////

double Exponential_Integral_En( double x, int n )
{
   long double xx = (long double) x;
   long double exp_x, e1;

   if (x < 0.0) return DBL_MAX;
   if (x == 0.0) return (n < 2) ? DBL_MAX : 1.0 / (double)(n - 1);

   exp_x = expl(-xx);
   if (n == 0) return (double) (exp_x / xx);

   e1 = - xExponential_Integral_Ei(-xx);
   if ( n == 1) return (double) e1;

   if ( (xx + n) >= 20.0L ) return (double) Continued_Fraction_En(xx, n, exp_x);
   if ( xx <= 1.0L ) return (double) Power_Series_En(xx, n);
   return (double) Continued_Fraction_En(xx, n, exp_x);
}


////////////////////////////////////////////////////////////////////////////////
// static long double Power_Series_En( long double x, int n)                  //
//                                                                            //
//  Description:                                                              //
//     For 0 < x <= 1.0 and x + n < 20 the power series representation for    //
//     En(x) + ln(x) * (-x)^(n-1) / (n-1)! is used.                           //
//                                                                            //
//     The power series expansion of En + ln(x) * (-x)^(n-1) / (n-1)! is      //
//        + psi(n) * (-x)^(n-1)/(n-1)! - Sum (-x)^k / [k - (n-1)] k!, where   //
//     the Sum extends from k = 0 to inf excluding k = n-1.                   //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of the exponential integral En().         //
//     int          n  The order of Schloemilch's exponential integral.       //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral En evaluated at x.               //
////////////////////////////////////////////////////////////////////////////////

static long double Power_Series_En( long double x, int n)
{ 
   long double xn = 1.0L;
   long double Sn;
   long double Sm1 = 0.0L;
   long double g = 0.5772156649015328606065121L;
   long double factorial;
   long double psi_n;
   long double term;
   int i;

   term = 1.0L / (long double)(1-n);
   psi_n = -g;
   factorial = 1.0L;
   Sn = term;
   for (i = 1; i < (n-1); i++) {
      factorial *= (long double) i;
      psi_n += 1.0L / (long double) i;
      xn *= (-x);
      term = xn / (factorial * (long double)(i - n + 1));
      Sn += term;
   }

   factorial *= (long double) (n-1);
   psi_n += 1.0L / (long double) (n-1);
   xn *= (-x);
   Sn = xn * (-logl(x) + psi_n) / factorial - Sn;

   Sm1 = Sn; 
   for (i = n;;i++) {
      factorial *= (long double) i;
      xn *= (-x);
      term = xn / (factorial * (long double)(i - n + 1));
      Sn -= term;
      if (fabsl(term) <= LDBL_EPSILON * fabsl(Sm1) ) break;
      Sm1 = Sn;
   }
   return (double) Sn;
}


////////////////////////////////////////////////////////////////////////////////
// static long double Continued_Fraction_En( long double x, int n,            //
//                                                         long double expx)  //
//                                                                            //
//  Description:                                                              //
//     For x > 1 or x + n >= 20, the continued fraction representation of En  //
//     converges fairly rapidly.                                              //
//                                                                            //
//     The continued fraction expansion of En(x) is:                          //
//        En(x) = exp(-x) { 1/(x+n-) 1*n/(x+n+2-) 2*(n+1)/(x+n+4-)  ... }.    //
//                                                                            //
//                                                                            //
//  Arguments:                                                                //
//     long double  x     The argument of the exponential integral En().      //
//     int          n     The order of Schloemilch's exponential integral.    //
//     long double  expx  exp(-x).                                            //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral En evaluated at x.               //
////////////////////////////////////////////////////////////////////////////////

static long double Continued_Fraction_En( long double x, int n,
                                                            long double expx )
{
   long double Am1 = 1.0L;
   long double A0 = 0.0L;
   long double Bm1 = 0.0L;
   long double B0 = 1.0L;
   long double a = expx;
   long double b = x + (long double) n;
   long double Ap1 = b * A0 + a * Am1;
   long double Bp1 = b * B0 + a * Bm1;
   long double eps = 10.0L*LDBL_EPSILON;
   double lhs, rhs;
   int j = 1;

   while ( fabsl(Ap1 * B0 - A0 * Bp1) > eps * fabsl(A0 * Bp1) ) {
      if ( fabsl(Bp1) > 1.0L) {
         Am1 = A0 / Bp1;
         A0 = Ap1 / Bp1;
         Bm1 = B0 / Bp1;
         B0 = 1.0L;
      } else {
         Am1 = A0;
         A0 = Ap1;
         Bm1 = B0;
         B0 = Bp1;
      }
      a = -j * (n + j - 1);
      b += 2.0L;
      Ap1 = b * A0 + a * Am1;
      Bp1 = b * B0 + a * Bm1;
      j += 1;
      lhs = (double)fabsl(Ap1 * B0 - A0 * Bp1);
      rhs = (double)fabsl(eps * A0 * Bp1);
      lhs = (double)Ap1;
      rhs = (double)Bp1;
   }
   return (double) (Ap1 / Bp1);
}
