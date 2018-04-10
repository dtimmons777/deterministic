# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>
#include <vector>
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


  void MMS(int N, int pts, double c, double end, double error_max, double flux[]  ) {

   int run_max=900;
   FILE *fp;
//   fp=fopen("Sn_num51210.txt","w");
   std::vector<double> sclr_flux_ml(N), sclr_flux_mr(N), sclr_flux_ml_old(N), sclr_flux_mr_old(N),s_kl(N),s_kr(N);
   double error=5.0;
   std::vector<double> q_ml(N,pts), q_mr(N,pts), pos(N), delta_xm(N);
   int i, j, k,ii, run=0;
   double start=0.0;
   double eps=pow(0.5, 0.0);
   double sigt=2.0/eps, siga=(sigt-sigt*c)*eps, sigs=sigt*c,a;
   double delta_x=(end-start)/(N);
   double Q[N][pts];
   double flux_boundRight=0.0;
   double flux_boundLeft=0.0;
   double flux_boundR=0.0;
   double flux_boundL=0.0;
   double flux_bound[2];
   double sclr_flux_old[N];
   double err_old=1.0,err=1.0 ; 
   VectorXd Mass(2);
   VectorXd f(2),r(4*N),d(4*N);
   VectorXd sclr_flux(N);
   MatrixXd A(2,2), P(4*N,4*N), PP(4*N,4*N);


/////////////////////////////  Legendre Stuff   ////////////////////////////////////////////////////////////

    
    double mu[pts], wght[pts];
    p_quadrature_rule ( pts, mu ,wght); 
    ii=0;
    while (mu[ii]<0) ii+=1;    //the number of mu points that are negative


/////////////////////  Initialize numbers on x mesh ////////////////////////////////////////////////////////

    double flux_ml[pts][N], flux_mr[pts][N];
    
    for (i=0; i<N; i++){
        pos[i]=i*delta_x+delta_x/2.0;
        sclr_flux_ml[i]=0.2;
        sclr_flux_mr[i]=0.2;
        
        for(j=0;j<pts;j++){
            a=(.25+1.5*mu[j]*0.75*mu[j]*mu[j])
            Q[i][j]=mu[j]*a*(end*2.0*pos[i])+sigt*a*(pos[i]*(end-pos[i]))-sigs/2.0*(pos[i]*(end-pos[i]));

            q_mr(i,j)=Q[i]*delta_x/4.0;
            q_ml(i,j)=Q[i]*delta_x/4.0;
        }
        delta_xm[i]=delta_x;
        sclr_flux_old[i]=0.1;
//printf("%lf \n", sclr_flux_mr[i]);
    }

    for (i=0; i<N; i++){     
        for (j=0; j<pts; j++) flux_ml[j][i]=1.0;
        for (j=0; j<pts; j++) flux_mr[j][i]=1.0;
    }

/////////////////////// P //////////////////////////////////
 

  for(k=0;k<N;k++){
      for (i=0;i<N;i++){
          P(k*4,i*4)=P(k*4+1,i*4+1)=P(k*4+2,i*4+2)=P(k*4+3,i*4+3)=0.0;
      }
  }
     P(0,0)=P(2,2)=0.25+siga*delta_x/3.0;  
     P(1,1)=P(3,3)=1.0+sigt*delta_x;
     P(0,2)=P(2,0)=siga*delta_x/6.0;
     P(0,3)=P(1,2)=0.5;
     P(2,1)=P(3,0)=-0.5;
     P(1,3)=P(3,1)=sigt*delta_x/2.0;
     P(2,5)=P(3,4)=0.5;
     P(2,4)=-0.25;
     P(3,5)=-1.0;

  for(k=1;k<N-1;k++){

      P(4*k,4*k-2)=P(4*k+2,4*k+4)=-0.25;
      P(4*k,4*k-1)=P(4*k+1,4*k-2)=P(4*k+3,4*k)=P(4*k+2,4*k+1)=-0.5;
      P(4*k,4*k)=P(4*k+2,4*k+2)=0.25+siga*delta_x/3.0;
      P(4*k,4*k+2)=P(4*k+2,4*k)=siga*delta_x/6.0;
      P(4*k,4*k+3)=P(4*k+1,4*k+2)=P(4*k+2,4*k+5)=P(4*k+3,4*k+4)=0.5;
      P(4*k+1,4*k-1)=P(4*k+3,4*k+5)=-1.0;
      P(4*k+1,4*k+1)=P(4*k+3,4*k+3)=1.0+sigt*delta_x;
      P(4*k+1,4*k+3)=P(4*k+3,4*k+1)=sigt*delta_x/2.0;
  } 
     P(4*N-4,4*N-4)=P(4*N-2,4*N-2)=0.25+siga*delta_x/3.0;  
     P(4*N-3,4*N-3)=P(4*N-1,4*N-1)=1.0+sigt*delta_x;
     P(4*N-4,4*N-2)=P(4*N-2,4*N-4)=siga*delta_x/6.0;
     P(4*N-4,4*N-1)=P(4*N-3,4*N-2)=0.5;
     P(4*N-2,4*N-3)=P(4*N-1,4*N-4)=-0.5;
     P(4*N-3,4*N-1)=P(4*N-1,4*N-3)=sigt*delta_x/2.0;
     P(4*N-4,4*N-5)=P(4*N-3,4*N-6)=-0.5;
     P(4*N-4,4*N-6)=-0.25;
     P(4*N-3,4*N-5)=-1.0;
     PP=P.inverse();
//std::cout<< P << std::endl;

///////////////////////////////////////////////////////////////////////////////////
//  For mu(m)>0 keep the left point and sweep left to right
//                  Start Left to Right Sweep 
///////////////////////////////////////////////////////////////////////////////////
     bool converge;      
     converge=true;

while(converge){
  for(i=0;i<pts-ii;i++){
     for(j=0;j<N;j++){     
        A(0,0)=A(1,1)= mu[ii+i]/2.0+1.0/3.0*sigt*delta_xm[j];
        A(1,0)=-mu[ii+i]/2.0+1.0/6.0*sigt*delta_xm[j];
        A(0,1)= mu[ii+i]/2.0+1.0/6.0*sigt*delta_xm[j];
        flux_bound[0]=mu[ii+i]*flux_boundL;

 
        Mass(0)=sigs/6.0*delta_x*(sclr_flux_ml[j]+sclr_flux_mr[j]/2.0)+q_ml[j]+flux_bound[0];
        Mass(1)=sigs/6.0*delta_x*(sclr_flux_mr[j]+sclr_flux_ml[j]/2.0)+q_mr[j];

//------------------solve for angular flux------------------

        f=A.householderQr().solve(Mass);

        flux_ml[ii+i][j]=f(0);  // fills the lower half of the matrix
        flux_mr[ii+i][j]=flux_boundL=f(1);
       

    }  // end of left to Right sweep loop
    flux_boundL=flux_boundLeft;

  }  // end of mu
//////////////////////////////////////////////////////////////////////////////////
//                   mu(m)<0
//           Start Right to Left Sweep  
//////////////////////////////////////////////////////////////////////////////////

  for(i=0;i<ii;i++){
    for(j=N-1;j>=0;j--){
        
        A(0,0)= A(1,1)= -mu[i]/2.0+1.0/3.0*sigt*delta_xm[j];
        A(1,0)=-mu[i]/2.0+1.0/6.0*sigt*delta_xm[j];
        A(0,1)= mu[i]/2.0+1.0/6.0*sigt*delta_xm[j];
        flux_bound[1]=mu[i]*flux_boundR;
         
        Mass(0)=sigs/6.0*delta_xm[j]*(sclr_flux_ml[j]+sclr_flux_mr[j]/2.0)+q_ml[j];
        Mass(1)=sigs/6.0*delta_xm[j]*(sclr_flux_mr[j]+sclr_flux_ml[j]/2.0)+q_mr[j]-flux_bound[1];

//---------------solve for angular flux-------------------------
 

        f=A.householderQr().solve(Mass);

        flux_mr[i][j]=f(1);
        flux_ml[i][j]=flux_boundR=f(0);
       
    }  // end of right to left sweep  loop
   
    flux_boundR=flux_boundRight;

  }  // end of mu
  
//if(run==0) std::cout << A <<std::endl;

//////////////////////  Scalar Flux   ///////////////////////////////////////////

  for(i=0;i<N;i++){
          sclr_flux_mr_old[i]=sclr_flux_mr[i];
          sclr_flux_ml_old[i]=sclr_flux_ml[i]; 
          sclr_flux_mr[i]=0.0;
          sclr_flux_ml[i]=0.0;
  }

  for(i=0;i<N;i++){ 
      for(j=0;j<pts;j++){
          sclr_flux_mr[i]+=wght[j]*flux_mr[j][i];
          sclr_flux_ml[i]+=wght[j]*flux_ml[j][i];
      }
  }

//////////////////////   Residual  ////////////////////////////////////////////////
  for(k=0;k<N;k++){
      s_kl[k]=sigs/3.0*delta_x*(sclr_flux_ml[k]-sclr_flux_ml_old[k]+(sclr_flux_mr[k]-sclr_flux_mr_old[k])/2.0);
      s_kr[k]=sigs/3.0*delta_x*(sclr_flux_mr[k]-sclr_flux_mr_old[k]+(sclr_flux_ml[k]-sclr_flux_ml_old[k])/2.0);

      r(4*k)=s_kl[k];
      r(4*k+2)=s_kr[k];
      r(4*k+1)=r(4*k+3)=0.0;
  }

/////////////////////// solve for d and New scalar flux /////////////////////////////////////

    d=PP*r;
    for(k=0;k<N;k++){
        sclr_flux_ml[k]+=d(4*k);
        sclr_flux_mr[k]+=d(4*k+2);
//std::cout<< sclr_flux_ml[k] << std::endl;

    }
//std::getchar();
///////////////////////  Scalar flux error ////////////////////////////////////////
    double rho;  
    error=0.0;

    for (i=0;i<N;i++){
        error+=std::pow(sclr_flux_ml[i]+sclr_flux_mr[i]-sclr_flux_ml_old[i]-sclr_flux_mr_old[i],2.0);
    }
    error = std::sqrt(error);
    rho=(error/err_old);
    if((1.0-rho)/rho*error_max>=error) converge=false;
    err_old=error;
//std::cout<< error << " "<<std::endl;   
    run++;
    if (run%(run_max)==0) std::cout<<"error= "<<error<<std::endl;


}   //end of convergence

//////////////////////////////////////////////////////////////////////////////

    for (i=0;i<N;i++){
        sclr_flux(i)=(sclr_flux_ml[i]+sclr_flux_mr[i])/2.0;
        sclr_flux_old[i]=sclr_flux(i);
        flux[i]=sclr_flux_old[i];
    }


/////////////////////////////  Balance  /////////////////////////////////////

double AA,SS,BB,JL,JR;
    AA=SS=JL=JR=0;

    for(i=0;i<N;i++){ 

        AA+=siga*delta_x/2.0*(sclr_flux_mr[i]+sclr_flux_ml[i]);
        SS+=delta_x*(Q[i]);
//printf("%d\n", run);
    }

    for(i=0;i<ii;i++){

        JL+=wght[i]*mu[i]*flux_ml[i][0];
        JR+=wght[ii+i]*mu[ii+i]*flux_mr[ii+i][N-1];
    }

    BB=SS-AA-(JR-JL);


//std::cout << JR<<"  "<<JL<< "  "<<  AA<< " " << SS <<" "<< BB <<std::endl;
//std::cout<<sclr_flux<<std::endl;    

    

/////////////////////////////////////////////////////////

return;

}
