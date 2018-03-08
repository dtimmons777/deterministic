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
   double flux_ml[N], flux_mr[N], sclr_flux_ml[N], sclr_flux_mr[N];
   double Blk[N], Brk[N],delta_xr[N],delta_xl[N], error=5;
   double q_ml[N], q_mr[N];
   int i, j, k;
   double start=0.0;
   double end=5.0;
   double sigt=1.0, c=0.9, sigs=sigt*c, siga=sigt-sigs;
   double delta_x=(end-start)/N;
   double Q=2.0;
   double flux_boundRight=delta_x*Q/2.0;
   double flux_boundLeft=delta_x*Q/2.0;
   double flux_boundR=4.35;
   double flux_boundL=4.35;
   double flux_bound[2];
   double sclr_flux_old[N];
   VectorXd Mass(2);
    VectorXd f(2);
    VectorXd sclr_flux(N);
    MatrixXd A(2,2);

///////  Initialize numbers on x mesh /////////////////////////////////////////////////////////////////////

    for (i=0; i<N; i++){
        sclr_flux_ml[i]=1.0;
        sclr_flux_mr[i]=1.0;
        flux_ml[i]=1.0;
        flux_mr[i]=1.0;
        Blk[i]=(end-i*delta_x)/delta_x;
        Brk[i]=(i*delta_x-start)/delta_x;
        q_mr[i]=delta_x*Q/2.0;
        q_ml[i]=delta_x*Q/2.0;
        delta_xr[i]=delta_x;//(end-i*delta_x);
        delta_xl[i]=delta_x;//(i*delta_x);
        sclr_flux_old[i]=100;
    }

//////////////  Legendre Stuff   //////////////////////////////////////////////////////////////////////////

    
    double mu[pts], wght[pts];
    p_quadrature_rule ( pts, mu ,wght); 



/////////////////////////////////////////////////////////////////////////////////////
while(fabs(error)>1){
  for(i=0;i<pts;i++){
     for(j=0;j<N;j++){     
        if (mu[i]>0){
            A(0,0)= mu[i]/2.0+1.0/3.0*sigt*delta_xr[j];
            A(1,1)= mu[i]/2.0+1.0/3.0*sigt*delta_xr[j];
            A(1,0)=-mu[i]/2.0+1.0/6.0*sigt*delta_xr[j];
            A(0,1)= mu[i]/2.0+1.0/6.0*sigt*delta_xr[j];
            flux_bound[0]=mu[i]*flux_boundL;
            flux_bound[1]=0;
            if (j==0) flux_bound[0]=mu[i]*flux_boundLeft;
        }
        else{
            A(0,0)=-mu[i]/2.0+1.0/3.0*sigt*delta_xr[j];
            A(1,1)=-mu[i]/2.0+1.0/3.0*sigt*delta_xr[j];
            A(1,0)=-mu[i]/2.0+1.0/6.0*sigt*delta_xr[j];
            A(0,1)= mu[i]/2.0+1.0/6.0*sigt*delta_xr[j];
            flux_bound[0]=0;
            flux_bound[1]=mu[i]*flux_boundR;
            if (j==pts-1) flux_bound[1]=mu[i]*flux_boundRight;
        }

/////////////////////////////////  Mass side  //////////////////////////////////////////////////// 
 
        Mass(0)=sigs/2*delta_xr[j]/3*sclr_flux_ml[j]+sigs/2*delta_xr[j]/6*sclr_flux_mr[j]+q_ml[j]+flux_bound[0];
        Mass(1)=sigs/2*delta_xr[j]/3*sclr_flux_mr[j]+sigs/2*delta_xr[j]/6*sclr_flux_ml[j]+q_mr[j]-flux_bound[1];

////////////////////////////////   solve flux /////////////////////////////////////////////////
        f=A.householderQr().solve(Mass);

        if (mu[i]>0) flux_ml[j+1]=flux_boundL=f(0);
        if (mu[i]<0) flux_mr[j+1]=flux_boundR=f(1);

        if (mu[i]>0) sclr_flux(j)+=(wght[i]*f(0));
        if (mu[i]<0) sclr_flux(j)+=(wght[i]*f(1));
    }  // end of mu
   

  }  // end of Right sweep loop


////////////////////////////////  Start Left Sweep  /////////////////////////////////////////////////////
  for(i=pts-1;i>=0;i--){
    for(j=N-1;j>=0;j--){
        if (mu[i]>0){
            A(0,0)= mu[i]/2.0+1.0/3.0*sigt*delta_xl[j];
            A(1,1)= mu[i]/2.0+1.0/3.0*sigt*delta_xl[j];
            A(1,0)=-mu[i]/2.0+1.0/6.0*sigt*delta_xl[j];
            A(0,1)= mu[i]/2.0+1.0/6.0*sigt*delta_xl[j];
            flux_bound[0]=mu[i]*flux_boundL;
            flux_bound[1]=0;
            if (j==0) flux_bound[1]=mu[i]*flux_boundLeft;
        }
        else{
            A(0,0)=-mu[i]/2.0+1.0/3.0*sigt*delta_xl[j];
            A(1,1)=-mu[i]/2.0+1.0/3.0*sigt*delta_xl[j];
            A(1,0)=-mu[i]/2.0+1.0/6.0*sigt*delta_xl[j];
            A(0,1)= mu[i]/2.0+1.0/6.0*sigt*delta_xl[j];
            flux_bound[0]=0;
            flux_bound[1]=mu[i]*flux_boundR;
            if (j==pts-1) flux_bound[1]=mu[i]*flux_boundRight;
        }

/////////////////////////////////  Mass side  //////////////////////////////////////////////////// 
 
        Mass(0)=sigs/2*delta_xl[j]/3*sclr_flux_ml[j]+sigs/2*delta_xl[j]/6*sclr_flux_mr[j]+q_ml[j]+flux_bound[0];
        Mass(1)=sigs/2*delta_xl[j]/3*sclr_flux_mr[j]+sigs/2*delta_xl[j]/6*sclr_flux_ml[j]+q_mr[j]-flux_bound[1];

/////////////////////////////////////////////////////////////////////////////////////
        f=A.householderQr().solve(Mass);

        if (mu[i]>0) flux_ml[j+1]=flux_boundL=f(0);
        if (mu[i]<0) flux_mr[j+1]=flux_boundR=f(1);

        if (mu[i]>0) sclr_flux(j)+=(wght[i]*f(0));
        if (mu[i]<0) sclr_flux(j)+=(wght[i]*f(1));
    }  // end of mu
   

  }  // end of left sweep sweep loop

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (i=0;i<N;i++){
        error+=(1-sclr_flux[i]/sclr_flux_old[i]);
        sclr_flux_old[i]=sclr_flux[i];
        sclr_flux_mr[i]=sclr_flux[i];
        sclr_flux_ml[i]=sclr_flux[i];
    }
std::cout << sclr_flux << std::endl;
}   //end of convergence




return 0;

}
