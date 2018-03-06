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
int main ( )

{

  int N=31, xbin=100;
  double c=0.9;
  double g[N+1], quad[N+1], wght[N+1];
  int i,j,k,n,m;
  double start=0, length=5;
  double  flux_particular[xbin+1], pos[xbin+1] ;
  double sigt[xbin+1],sigs[xbin+1],q0[xbin+1],siga[xbin+1];

////////////////////////////////////////////////////////////////////////////////////

    for (i=0;i<xbin+1;i++){
        sigt[i]=1.0;
        sigs[i]=c/1.0;
        siga[i]=sigt[i]-sigs[i];
        q0[i]=2.0;
        flux_particular[i]=q0[i]/(2.0*siga[i]);
        pos[i]=(length-start)/(xbin)*i;
//std::cout << pos[i] << std::endl;
    }


/////////////////////////////////  V and W Matrix //////////////////////////////////////////////////
    MatrixXd V(N+1,N+1), W(N+1,N+1);

    for (i=0;i<=N;i++){
//std::cout << i << std::endl;
        for (j=0;j<=N;j++){
            if (i==j){
                if (i==0){
                    V(i,j)=siga[i];
                    W(i,j)=0.0;
                }
                else{
                    V(i,j)=sigt[i]+2*i*sigt[i];
                    W(i,j)=0.0;
                }
            }
            else if (j==i+1){
                V(i,j)=0.0;
                W(i,j)=sigt[i]*(i+1);
            }  
            else if(i>0 && j==i-1){
                V(i,j)=0.0;
                W(i,j)=sigt[i]*(i);
             
            }
            else{
                V(i,j)=0.0;
                W(i,j)=0.0;
            }
        }
    }
//std::cout << W << std::endl;

/////////////////////////////////////////////////////////////////////////////////////
    MatrixXd W_inv(N+1,N+1), gnk(N+1,N+1), WV(N+1,N+1);
    VectorXd eivals(N);

    //W_inv=W.inverse();
    WV=-W.inverse()*V;



    EigenSolver<MatrixXd> es(WV);
    eivals=es.eigenvalues().real();
    gnk=es.eigenvectors().real();

//std::cout << eivals << std::endl;

//////////////////////////////////// Normalize Eigenvectors ///////////////////////////////////////////////
    
    double Nk[N+1], gprime[N+1][N+1];

    p_quadrature_rule ( N+1, quad ,wght); 

    for (k=0; k<=N;k++){
  
        for (m=0;m<=N;m++){
            Nk[k]+=c/2.0*wght[m]*gnk.row(m)(k);           
        
        }
        for (m=0;m<=N;m++){
            gprime[k][m]=gnk.row(k)(m)/Nk[k];
        } 
    }
//std::cout << gnk << std::endl;
///////////////////////////////////  Construct Ba=f B matrix  //////////////////////////////////////////////////

    double B[N+1][N+1],fun;

    for (k=0;k<=N;k++){
        for (n=0;n<=N;n++){
           fun=0.0;
           for (j=0;j<=N;j++){   
               fun+=(2.0*j+1.0)/2.0*gnk(j,n)*boost::math::legendre_p(j,quad[k]);
            }
            if (quad[k]>0){
                B[k][n]=fun*exp(eivals(n)*start*sigt[k]);
            }
            else{
                B[k][n]=fun*exp(eivals(n)*length*sigt[k]);

            }
        }

    }

///////////////////////////////////// Solve for Coefficients ///////////////////////////////////////////////////


    VectorXd a(N+1);
    VectorXd f(N+1);
    MatrixXd BB(N+1,N+1);
    for (i=0;i<=N;i++){
        for (j=0;j<=N;j++){
            BB(i,j)=B[i][j];
        }
        f(i)=-flux_particular[i];
    }
     a=BB.householderQr().solve(f);
//std::cout << BB << std::endl;
//std::cout << a << std::endl;
//////////////////////////////////// Solve for Flux at all points /////////////////////////////////////////////

    double sclr_flux[N+1];

    for (i=0;i<=xbin;i++){
        fun=0;
        for (k=0;k<=N;k++){
            fun+=a(k)*gnk(0,k)*exp(sigt[i]*pos[i]*eivals(k));

         }
         sclr_flux[i]=q0[i]/siga[i]+fun;
//printf("%lf \n", sclr_flux[i]);
    }


///////////////////////////////////////////////////   output   ///////////////////////////////////////////////////////////
FILE *fp;
fp=fopen("Pn32.txt","w");
for(i=0;i<=xbin;i++){
fprintf(fp, "%lf %lf \n", pos[i], sclr_flux[i]);
}

fclose(fp);




return 0;

}
