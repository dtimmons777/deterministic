# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>
# include <iostream>
# include </usr/include/boost/math/special_functions/legendre.hpp>
# include <eigen/Eigen/Dense>
# include <eigen/Eigen/Eigenvalues> 
# define max(a,b)\
({ __typeof__ (a) _a = (a);\
    __typeof__ (b) _b = (b); \
   _a > _b ? _a : _b; })

#include <boost/math/special_functions/legendre.hpp>
# include "golub.h"
using namespace Eigen;
namespace boost{ namespace math{
}}


int main ( ) {


   int N=10;
   int pts=N;
   double A[2][2], Mass[2];
   double flux_ml[N], flux_mr[N], sclr_flux_ml[N], sclr_flux_mr[N];
   double Blk[N], Brk[N];
   double q_ml[N], q_mr[N];
   int i, j, k;
   double start=0.0;
   double end=5.0;
   double sigt=1.0, c=0.9, sigs=sigt*c, siga=sigt-sigs;
   double delta_x=(end-start)/N;
   double Q=2.0;
   double flux_boundRight=4.35;
   double flux_boundLeft=4.35;
   double flux_boundR=4.35;
   double flux_boundL=4.35;
   double flux_bound[2];

///////  Initialize numbers on x mesh /////////////////////////////////////////////////////////////////////

    for (i=0; i<N; i++){
        sclr_flux_ml[i]=1.0;
        sclr_flux_mr[i]=1.0;
        Blk[i]=(end-i*delta_x)/delta_x;
        Brk[i]=(i*delta_x-start)/delta_x;
        q_mr[i]=delta_x*Q/2.0;
        q_ml[i]=delta_x*Q/2.0;
    }

//////////////  Legendre Stuff   //////////////////////////////////////////////////////////////////////////

    
    double mu[pts], wght[pts];
    p_quadrature_rule ( pts, mu ,wght); 



/////////////////////////////////////////////////////////////////////////////////////

    for(i=0;i<pts;i++){
        if (mu[i]>0){
            A[0][0]= mu[i]/2.0+1.0/3.0*sigt*delta_x;
            A[1][1]= mu[i]/2.0+1.0/3.0*sigt*delta_x;
            A[1][0]=-mu[i]/2.0+1.0/6.0*sigt*delta_x;
            A[0][1]= mu[i]/2.0+1.0/6.0*sigt*delta_x;
            flux_bound[0]=mu[i]*flux_boundL;
            flux_bound[1]=0;
            if (i==0) flux_bound[1]=mu[i]*flux_boundleft;
        }
        else{
            A[0][0]=-mu[i]/2.0+1.0/3.0*sigt*delta_x;
            A[1][1]=-mu[i]/2.0+1.0/3.0*sigt*delta_x;
            A[1][0]=-mu[i]/2.0+1.0/6.0*sigt*delta_x;
            A[0][1]= mu[i]/2.0+1.0/6.0*sigt*delta_x;
            flux_bound[0]=0;
            flux_bound[1]=mu[i]*flux_boundR;
            if (i==pts-1) flux_bound[1]=mu[i]*flux_boundRight;
        }

/////////////////////////////////  Mass side  //////////////////////////////////////////////////// 
 
        Mass[0]=sigs/2*delta_x/3*sclr_flux_ml[i]+sigs/2*delta_x/6*sclr_flux_mr[i]+q_ml[i]+flux_bound[0];
        Mass[1]=sigs/2*delta_x/3*sclr_flux_mr[i]+sigs/2*delta_x/6*sclr_flux_ml[i]+q_mr[i]-flux_bound[1];

/////////////////////////////////////////////////////////////////////////////////////







    }  // end of mu



return 0;

}
