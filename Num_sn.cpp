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
   double sclr_flux_ml[N], sclr_flux_mr[N];
   double delta_xm[N], error=5.0;
   double q_ml[N], q_mr[N];
   int i, j, k,ii, run=0;
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


/////////////////////////////  Legendre Stuff   ////////////////////////////////////////////////////////////

    
    double mu[pts], wght[pts];
    p_quadrature_rule ( pts, mu ,wght); 
    ii=0;
    while (mu[ii]<0) ii+=1;    //the number of mu points that are negative


/////////////////////  Initialize numbers on x mesh ////////////////////////////////////////////////////////

    double flux_ml[ii][N], flux_mr[pts-ii][N];
    
    for (i=0; i<N; i++){
        sclr_flux_ml[i]=1.0;
        sclr_flux_mr[i]=1.0;
        q_mr[i]=delta_x*Q/2.0;
        q_ml[i]=delta_x*Q/2.0;
        delta_xm[i]=delta_x;//(end-i*delta_x);
        sclr_flux_old[i]=0.1;
    }

    for (i=0; i<N; i++){     
        for (j=0; j<ii; j++) flux_ml[j][i]=1.0;
        for (j=0; j<pts-ii; j++) flux_mr[j][i]=1.0;
    }

//  For mu(m)>0 keep the left point and sweep left to right
/////////////////////////////  Start Left to Right Sweep  ////////////////////////////////////////////////

while(fabs(error)>0.1){
  for(i=0;i<pts-ii;i++){
     for(j=0;j<N;j++){     
        A(0,0)=A(1,1)= mu[ii+i]/2.0+1.0/3.0*sigt*delta_xm[j];
        A(1,0)=-mu[ii+i]/2.0+1.0/6.0*sigt*delta_xm[j];
        A(0,1)= mu[ii+i]/2.0+1.0/6.0*sigt*delta_xm[j];
        flux_bound[0]=mu[ii+i]*flux_boundL;
        flux_bound[1]=0;
        //if (j==0) flux_bound[0]=mu[i]*flux_boundLeft;
 
        Mass(0)=sigs/2*delta_xm[j]/3*sclr_flux_ml[j]+sigs/2*delta_xm[j]/6*sclr_flux_mr[j]+q_ml[j]+flux_bound[0];
        Mass(1)=sigs/2*delta_xm[j]/3*sclr_flux_mr[j]+sigs/2*delta_xm[j]/6*sclr_flux_ml[j]+q_mr[j]-flux_bound[1];

//------------------solve flux------------------
        f=A.householderQr().solve(Mass);

        flux_ml[i][j]=flux_boundL=f(0);
       
        
    }  // end of left to Right sweep loop
    flux_boundL=flux_boundLeft;

  }  // end of mu
 
  

////////////////////////////////  Start Right to Left Sweep  /////////////////////////////////////////////////////
  for(i=0;i<ii;i++){
    for(j=N-1;j>=0;j--){
        
        A(0,0)= A(1,1)= -mu[i]/2.0+1.0/3.0*sigt*delta_xm[j];
        A(1,0)=-mu[i]/2.0+1.0/6.0*sigt*delta_xm[j];
        A(0,1)= mu[i]/2.0+1.0/6.0*sigt*delta_xm[j];
        flux_bound[0]=0;
        flux_bound[1]=mu[i]*flux_boundR;
        //if (j==N-1) flux_bound[1]=mu[i]*flux_boundRight;
 
        Mass(0)=sigs/2*delta_xm[j]/3*sclr_flux_ml[j]+sigs/2*delta_xm[j]/6*sclr_flux_mr[j]+q_ml[j]+flux_bound[0];
        Mass(1)=sigs/2*delta_xm[j]/3*sclr_flux_mr[j]+sigs/2*delta_xm[j]/6*sclr_flux_ml[j]+q_mr[j]-flux_bound[1];

//---------------solve flux-------------------------
        f=A.householderQr().solve(Mass);

        flux_mr[i][j]=flux_boundR=f(1);
        
    }  // end of right to left sweep  loop
   
    flux_boundR=flux_boundRight;

  }  // end of mu
  

////////////////////////////////////////////  Scalar Flux   /////////////////////////////////////////////
  for(i=0;i<N;i++){ 
      for(j=0;j<ii;j++){
          sclr_flux_mr[i]+=wght[j]*flux_mr[j][i];
      }
      for(j=0;j<pts-ii;j++){
          sclr_flux_ml[i]+=wght[j+ii]*flux_ml[j][i];
      }

  }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    error=0;
    for (i=0;i<N;i++){
        sclr_flux(i)=sclr_flux_ml[i]+sclr_flux_mr[i];
        error+=fabs(sclr_flux(i)-sclr_flux_old[i]);
        sclr_flux_old[i]=sclr_flux(i);
    }
std::cout << error << "  "<< run <<std::endl;
run++;
}   //end of convergence
std::cout << sclr_flux << "  "<< run <<std::endl;

return 0;

}
