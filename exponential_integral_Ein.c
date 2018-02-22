////////////////////////////////////////////////////////////////////////////////
// File: exponential_integral_Ein.c                                           //
// Routine(s):                                                                //
//    Exponential_Integral_Ein                                                //
//                                                                            //
// Required Externally Defined Routines:                                      //
//    xExponential_Integral_Ei                                                //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                        // required for fabsl() and logl()     

//                    Required Externally Defined Routines                    //
long double xExponential_Integral_Ei( long double x );

////////////////////////////////////////////////////////////////////////////////
// double Exponential_Integral_Ein( double x )                                //
//                                                                            //
//  Description:                                                              //
//     The exponential integral Ein(x), called the entire exponential         //
//     integral, is the integral with integrand                               //
//                         (1 - exp(-t)) / t dt                               //
//     where the integral extends from 0 to x.                                //
//     At x = 0, Ein(x) = 0, and                                              //
//     for x != 0, Ein(x) = gamma + log |x| - Ei(-x), where gamma is the Euler//
//     constant and Ei(x) is the exponential integral.                        //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of the entire exponential integral Ein().      //
//                                                                            //
//  Return Value:                                                             //
//     The value of the entire exponential integral Ein evaluated at x.       //
//                                                                            //
//  Example:                                                                  //
//     double y, x;                                                           //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Exponential_Integral_Ein( x );                                     //
////////////////////////////////////////////////////////////////////////////////

double Exponential_Integral_Ein( double x )
{
   static long double g = 0.5772156649015328606065121L;
   long double xx = (long double) x;

   if ( x == 0.0 ) return 0.0;
   return (double) (g + logl(fabsl(xx)) - xExponential_Integral_Ei(-xx));
}
