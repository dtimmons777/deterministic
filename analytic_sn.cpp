# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>
# include <iostream>
# include <eigen/Eigen/Dense>
#include <eigen/Eigen/Eigenvalues> 
# define max(a,b)\
({ __typeof__ (a) _a = (a);\
    __typeof__ (b) _b = (b); \
   _a > _b ? _a : _b; })

# include "golub.h"
using namespace Eigen;

int main ( )

{
    int N=64;

    double z[N], A[N][N], Nk[N],B[N][N],zz[N][N],w[N],zero[N],mu[N],zprime[N][N], hi[N],a[N],aold[N],zi[N],lambda[N],lam[N] ;
    double c=0.9;
    double length=5.0, start=0.0;
    int xbins=100;
    double xpos[xbins], fluxold[xbins], flux[xbins];
    double sigt=1.0,siga=sigt*(1-c),sigs=sigt-siga;
    int i, j, k, m ;
    double q0=2.0, flux_part; 
    double error=1.0e-8, err=10.0, fun, errz=10, nmax=0;

// mu is the weights , zero is the quadrature points[-1,1]
    p_quadrature_rule ( N, zero ,mu); 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (i=0;i<=xbins;i++){
        xpos[i]=i*(length-start)/xbins;
        //printf("%lf \n", xpos[i]);   
        
        
        fluxold[i]=15.0;
        flux[i]=10.0;
    }
flux_part=q0/(2.0*siga);
//printf("%lf",flux_part);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (i=0; i<N;i++){
        w[i]=(length-start)/2.0*mu[i];
        zi[i]=(length-start)/2.0*zero[i]+(length+start)/2.0;
        z[i]=10.0;
        a[i]=10.0;
        aold[i]=15.0;
        lam[i]=10;
        Nk[i]=0;
//printf("%lf \n", w[i]);
    }


///////////////////////////////////////////////   Amj matrix   ///////////////////////////////////////////////////////
 double lad[N][N];   
    for (m=0; m<N;m++){
        
        for (j=0;j<N;j++){
            if (j == m){
                A[m][j]=c/(2.0*zero[m])*mu[m]-1.0/zero[m];
            }
            else{
                A[m][j]=c/(2.0*zero[m])*mu[j];
            }
            zz[i][j]=10;
            lad[i][j]=5; //printf("%lf ", A[m][j]);
        }
        //printf("\n "); 
    }

////////////////////////////////////////////////////  zk   ////////////////////////////////////////////////////////////
    MatrixXd mm(N,N);
    VectorXd eivals(N);
    MatrixXd zkm(N,N);
        
    for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            mm(i,j)=A[i][j];
         }
    }
EigenSolver<MatrixXd> es(mm);
 eivals=es.eigenvalues().real();
 zkm=es.eigenvectors().real();



////////////////////////////////////////////////////  z'k   ///////////////////////////////////////////////////////////
    for (k=0; k<N;k++){
  
        for (m=0;m<N;m++){
            Nk[k]+=c/2.0*mu[m]*zkm.row(m)(k);
//printf(" %lf ", zkm(k,m));            
        
        }
        for (m=0;m<N;m++){
            zprime[k][m]=zkm.row(m)(k)/Nk[k];
//std::cout << zprime[k][m]<<" ";
        } 
//std::cout << std::endl;
    }

////////////////////////////////////////////////   Bmk and fm   ////////////////////////////////////////////////////////
    for (m=0; m<N;m++){
        
        for (k=0;k<N;k++){
            if (zero[m]>0){
                B[m][k]=zprime[k][m]*exp(sigt*start*eivals[k]);
            }
            else{
                B[m][k]=zprime[k][m]*exp(sigt*length*eivals[k]);
            }
         //std::cout << B[m][k]<<" ";
        }
//std::cout << std::endl;
        
    }


////////////////////////////////////////////   Jacobi for a coefficients  //////////////////////////////////////////////
  VectorXd aa(N);
  VectorXd f(N);
  MatrixXd BB(N,N);
  for (i=0;i<N;i++){
     for (j=0;j<N;j++){
         BB(i,j)=B[i][j];
     }
  f(i)=-flux_part;
  }
  aa=BB.householderQr().solve(f);

//std::cout<< aa <<std::endl;     
//printf("%lf \n", a[i]);
     
//////////////////////////////////////////////////   flux   //////////////////////////////////////////////////////////////
    double funn;
    for (i=0;i<=xbins;i++){
       funn=0.0;
        for (m=0;m<N;m++){
            fun=0.0;
            for (k=0;k<N;k++){
                fun+=aa(k)*zprime[k][m]*exp(sigt*xpos[i]*eivals[k]);
            }
            funn+=mu[m]*(fun+flux_part);

        }
        flux[i]=funn;

    }
///////////////////////////////////////////////////   output   ///////////////////////////////////////////////////////////
FILE *fp;
fp=fopen("sn64.txt","w");
for(i=0;i<=xbins;i++){
fprintf(fp, " %lf  %lf \n", xpos[i], flux[i]);
}

fclose(fp);




return 0;
}



