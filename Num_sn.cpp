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

   int N=512;
   int pts=8;
   int run_max=90000;
   FILE *fp;
   fp=fopen("Sn_num51210.txt","w");
   double sclr_flux_ml[N], sclr_flux_mr[N];
   double delta_xm[N], error=5.0;
   double q_ml[N], q_mr[N], pos[N];
   int i, j, k,ii, run=0;
   double start=0.0;
   double end=1.0;
   double eps=pow(0.5,10.0);
   double sigt=1.0/eps, c=0.90, siga=(0.1)*eps, sigs=sigt-siga;
   double delta_x=(end-start)/(N);
   double Q[N];
   double flux_boundRight=0.0;
   double flux_boundLeft=0.0;
   double flux_boundR=0.0;
   double flux_boundL=0.0;
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

    double flux_ml[pts][N], flux_mr[pts][N];
    
    for (i=0; i<N; i++){
        pos[i]=i*delta_x+delta_x/2.0;
        sclr_flux_ml[i]=0.2;
        sclr_flux_mr[i]=0.2;
        
        Q[i]=eps;//*(2.0*(pos[i]/end)*(2.0*(end-pos[i])/end));
        q_mr[i]=Q[i]*delta_x/4.0;//-1/3*pow(delta_x,2)*(-2.0*end*(2.5*pos[i]+delta_x/2.0)+3.0*pow(pos[i],2)*2.0*x);
        q_ml[i]=Q[i]*delta_x/4.0;
        delta_xm[i]=delta_x;
        sclr_flux_old[i]=0.1;
//printf("%lf \n", sclr_flux_mr[i]);
    }

    for (i=0; i<N; i++){     
        for (j=0; j<pts; j++) flux_ml[j][i]=1.0;
        for (j=0; j<pts; j++) flux_mr[j][i]=1.0;
    }

//  For mu(m)>0 keep the left point and sweep left to right
/////////////////////////////  Start Left to Right Sweep  ////////////////////////////////////////////////

while(error > 1e-8 && run<run_max){
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
//std::cout<<f<<std::endl;
//std::getchar();         

    }  // end of left to Right sweep loop
    flux_boundL=flux_boundLeft;

  }  // end of mu

///////////////////////////////////////  mu(m)<0
////////////////////////////////  Start Right to Left Sweep  /////////////////////////////////////////////////////

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

////////////////////////////////////////////  Scalar Flux   /////////////////////////////////////////////

  for(i=0;i<N;i++){ 
          sclr_flux_mr[i]=0.0;
          sclr_flux_ml[i]=0.0;
  }

  for(i=0;i<N;i++){ 
      for(j=0;j<pts;j++){
          sclr_flux_mr[i]+=wght[j]*flux_mr[j][i];
          sclr_flux_ml[i]+=wght[j]*flux_ml[j][i];

      }

  }

///////////////////////  Scalar flux error ////////////////////////////////////////////////////////////

    error=0;

    for (i=0;i<N;i++){
        sclr_flux(i)=(sclr_flux_ml[i]+sclr_flux_mr[i])/2.0;
        error+=fabs(sclr_flux(i)-sclr_flux_old[i]);
        sclr_flux_old[i]=sclr_flux(i);
        
    }

    run++;
    if (run%run_max/10.0==0) std::cout<<"error= "<<error<<std::endl;
//std::cout<<sclr_flux<<std::endl;
//std::getchar(); 
}   //end of convergence
//printf("%d \n", run);

/////////////////////////////  Balance  /////////////////////////////////////

double AA,SS,BB,JL,JR;
    AA=SS=JL=JR=0;

    for(i=0;i<N;i++){ 

        AA+=siga*delta_x/2.0*(sclr_flux_mr[i]+sclr_flux_ml[i]);
        SS+=delta_x*(Q[i]);
//printf("%lf    %lf\n", sclr_flux_mr[i], sclr_flux_ml[i]);
    }

    for(i=0;i<ii;i++){

        JL+=wght[i]*mu[i]*flux_ml[i][0];
        JR+=wght[ii+i]*mu[ii+i]*flux_mr[ii+i][N-1];
    }

    BB=SS-AA-(JR-JL);

//std::cout << JR<<"  "<<JL<< "  "<<  AA<< " " << SS <<" "<< BB <<std::endl;
/////////////////////////////////////////////////////////////////////////////////////////////
//std::cout<<sclr_flux<<std::endl;    

    for(i=0;i<N;i++){
        fprintf(fp, "%lf %14.12f \n", pos[i], sclr_flux(i)) ;
    }
    fclose(fp);

/////////////////////////////////////////////////////////

    


////////////////////////////////////////////////////////
return 0;

}
