# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>
# define max(a,b)\
({ __typeof__ (a) _a = (a);\
    __typeof__ (b) _b = (b); \
   _a > _b ? _a : _b; })

# include "golub.h"

int main ( )

{
    int N=8;
    double z[N], A[N][N], Nk[N],B[N][N],zz[N][N],f[N],w[N],zero[N],mu[N],zprime[N], hi[N],a[N],aold[N],zi[N],lambda[N],lam[N] ;
    double c=0.9;
    double length=5.0, start=0.0;
    int xbins=20;
    double xpos[xbins], fluxold[xbins], flux[xbins];
    double sigt=1.0,siga=sigt*(1-c),sigs=sigt-siga;
    int i, j, k, m ;
    double q0=2.0, flux_part[xbins]; 
    double error=1.0e-8, err=10.0, fun, errz=10, nmax=0;

// mu is the weights , zero is the quadrature points[-1,1]
    p_quadrature_rule ( N, zero ,mu); 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (i=0;i<=xbins;i++){
        xpos[i]=i*(length-start)/xbins;
        if(xpos[i]==0 || xpos[i]==length){
        flux_part[i]=0;
        }
        else{
        flux_part[i]=q0/(2.0*siga);
        }
        fluxold[i]=15.0;
        flux[i]=10.0;
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (i=0; i<N;i++){
        w[i]=(length-start)/2.0*mu[i];
        zi[i]=(length-start)/2.0*zero[i]+(length+start)/2.0;
        z[i]=10.0;
        a[i]=10.0;
        aold[i]=15.0;
        lam[i]=10;
    }


///////////////////////////////////////////////   Amj matrix   ///////////////////////////////////////////////////////
 double lad[N][N];   
    for (m=0; m<N;m++){
        
        for (j=0;j<N;j++){
            if (j == m){
                A[m][j]=c/(2.0*zero[m])*w[m]-1.0/zero[m];
            }
            else{
                A[m][j]=c/(2.0*zero[m])*w[j];
            }
            zz[i][j]=10;
            lad[i][j]=5;
        }
 
    }

////////////////////////////////////////////////////  zk   ////////////////////////////////////////////////////////////
   while(error*1000<fabs(errz)){
        for (i=0;i<N;i++){  
            for (j=0;j<N;j++){
                lambda[i]+=z[j]*A[i][j];
                nmax=max(nmax,lambda[i]);
            }
            errz=max(error,fabs(1-lambda[i]/lam[i]));
	   // printf("%lf \n", fabs(1-lambda[i]/lam[i])); 
        }
        for(i=0;i<N;i++){ 
            z[i]=lambda[i]/nmax;
            lam[i]=lambda[i];
        }
    }
    
////////////////////////////////////////////////  zkm vector  /////////////////////////////////////////////////////////
    errz=10.0;
    
    for (i=0;i<N;i++){
        while(errz>error*1000){
            fun=0;
            for (j=0;j<N;j++){
                fun+=A[i][j]*zz[i][j];
            }
            for (j=0;j<N;j++){
                lad[i][j]=fun/lambda[j];
                errz=max(error,1.0-lad[i][j]/zz[i][j]); 
                zz[i][j]=lad[i][j];
printf("%lf \n", zz[i][j]);
            }
        }
    }

////////////////////////////////////////////////////  z'k   ///////////////////////////////////////////////////////////
    for (k=0; k<N;k++){
        
        for (m=0;m<N;m++){
            Nk[k]+=c/2.0*w[m]*zz[m][k];
            
        }
        zprime[k]=z[k]/Nk[k];        
//printf("%lf \n", zprime[k]);
    }

////////////////////////////////////////////////   Bmk and fm   ////////////////////////////////////////////////////////
    for (m=0; m<N;m++){
        
        for (k=0;k<N;k++){
            if (zero[m]>0){
                B[m][k]=zprime[k]*exp(sigt*start*A[m][k]);
            }
            else{
                B[m][k]=zprime[k]*exp(sigt*length*A[m][k]);
            }
        }
        f[m]=-flux_part[m];
    }


////////////////////////////////////////////   Jacobi for a coefficients  //////////////////////////////////////////////
    
    while (fabs(err)>error){
        for(i=0;i<N; i++){
            fun=0.0;
            for (j=0; j<N; j++){
                if (i !=j){ 
                    fun+=aold[j]*B[i][j];

                }
            }
        a[i]=1.0/B[i][i]*(f[m]-fun);
        err=max(error, a[i]-aold[i]);
        aold[i]=a[i];
printf("%lf \n", a[i]);
        }
    }

//////////////////////////////////////////////////   flux   //////////////////////////////////////////////////////////////
    double funn;
    for (i=0;i<=xbins;i++){
        for (m=0;m<N;m++){
            fun=0.0;
            for (k=0;k<N;k++){
                fun+=a[k]*zprime[k]*exp(sigt*xpos[i]*A[m][k])+flux_part[m];
            }
            funn+=w[m]*fun;

        }
        flux[i]=funn;

    }
///////////////////////////////////////////////////   output   ///////////////////////////////////////////////////////////
FILE *fp;
fp=fopen("sn.txt","w");
for(i=0;i<=xbins;i++){
fprintf(fp, "At %lf flux is %lf \n", xpos[i], flux[i]);
}

fclose(fp);




return 0;
}



