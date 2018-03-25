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
    int N=16;
    int xbins=2048;
    FILE *fp;
    fp=fopen("sn11.txt","w");
    double z[N], A[N][N], Nk[N],B[N][N],w[N],mu[N],wght[N],zprime[N][N], hi[N],a[N],aold[N],zi[N],lambda[N],lam[N] ;
    double c=0.95;
    double length=2.0, start=0.0;
    double xpos[xbins], fluxold[xbins], flux[xbins];
    double sigt=2.0,siga=sigt*(1.0-c),sigs=sigt-siga;
    int i, j, k, m ;
    double q0=1.0, flux_part,c0,c1,c2; 
    double error=1.0e-12, err=10.0, fun, errz=10, nmax=0,delta_x=(length-start)/xbins;

// wght is the weights , mu is the quadrature points[-1,1]
    p_quadrature_rule ( N, mu ,wght); 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    c0=-4.0*length*start/pow((length-start),2);
    c1=4.0*(length+start)/pow((length-start),2);
    c2=-4.0/pow((length-start),2);
 // for (j=0;j<N;j++){
    for (i=0;i<xbins;i++){
        xpos[i]=i*(length-start)/xbins;
        
        flux_part=q0/(2.0*siga);//1/siga*(c0/2.0+c1/2.0*(xpos[i]-mu[j]/sigt)+ c2*(pow(mu[j]/sigt,2)+pow(xpos[i],2)/2+sigs/(3.0*siga*pow(sigt,2))-xpos[i]*wght[j]/sigt));//
        fluxold[i]=15.0;
        flux[i]=10.0;

    }

 // }
//printf("%lf",flux_part);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for (i=0; i<N;i++){
        w[i]=(length-start)/2.0*wght[i];
        zi[i]=(length-start)/2.0*mu[i]+(length+start)/2.0;
        z[i]=10.0;
        a[i]=10.0;
        aold[i]=15.0;
        lam[i]=10.0;
        Nk[i]=0;
    }



///////////////////////////////////////////////   Amj matrix   ///////////////////////////////////////////////////////
   
    for (m=0; m<N;m++){
        for (j=0;j<N;j++){
            if (j == m){
                A[m][j]=c/(2.0*mu[m])*wght[m]-1.0/mu[m];
            }
            else{
                A[m][j]=c/(2.0*mu[m])*wght[j];
            }
 //printf("%lf ", A[m][j]);
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
            Nk[k]+=c/2.0*wght[m]*zkm.row(m)(k);
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
            if (mu[m]>0){
                if (eivals(k)<0){
                    B[m][k]=zprime[k][m];
                }
                else{
                    B[m][k]=zprime[k][m]*exp(sigt*-length*eivals[k]);
                }
            }
            else{
                if (eivals(k)<0){
                    B[m][k]=zprime[k][m]*exp(sigt*length*eivals[k]);
                }
                else{
                    B[m][k]=zprime[k][m];
                }
  
            }
            //std::cout << B[m][k]<<" ";
        }
        //std::cout << std::endl;
        
    }


////////////////////////////////////////////   QR Method for a coefficients  //////////////////////////////////////////////
    VectorXd aa(N);
    VectorXd f(N);
    MatrixXd BB(N,N);
    for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            BB(i,j)=B[i][j];
        }
        f(i)=0.0;
        for (j=0;j<xbins;j++){
            f(i)+=-flux_part;
        }
    }
    aa=BB.householderQr().solve(f);

//std::cout<< aa <<std::endl;     
//printf("%lf \n", a[i]);
     
//////////////////////////////////////////////////   flux   //////////////////////////////////////////////////////////////

    double funn;
    for (i=0;i<xbins;i++){
       funn=0.0;
        for (m=0;m<N;m++){
            fun=0.0;
            for (k=0;k<N;k++){
                if (eivals(k)>0){
                    fun+=aa(k)*zprime[k][m]*(exp(sigt*(xpos[i]+delta_x-length)*eivals(k))-exp(sigt*(xpos[i]-length)*eivals(k)))*1.0/(eivals(k)*sigt*delta_x);
                }
                else{
                    fun+=aa(k)*zprime[k][m]*(exp(sigt*(xpos[i]+delta_x)*eivals(k))-exp(sigt*(xpos[i])*eivals(k)))*1.0/(eivals(k)*sigt*delta_x);
                }
            }
            funn+=wght[m]*(fun/(xbins)+flux_part);

        }
        flux[i]=funn;

    }
///////////////////////////////////////////////////   output   ///////////////////////////////////////////////////////////

    
    for(i=0;i<xbins;i++){
        fprintf(fp, " %lf  %14.12f \n", i*xpos[1]+xpos[1]/2.0, flux[i]);
    }
    fclose(fp);


return 0;
}



