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
  double anginc, fun, diff, s  ;
  int i, j, n, k;
  int N=64;
  double c=0.9;
  double length=5.0;
  double qx=1.0;
  double sigt=1.0;
  double u=0;
  double v=length+u; 
  double sx[N],sv[N];
  double z[N];               //angular discretization
  double hi[N];             //spacial descretization
  i=0;
  j=0;
  anginc=2.0/(N+2);
  double x[N];
  double rx[N];
  double w[N];
  float A[N][N], A_inv[N], flux[N];
  double error=1e-8, errflux=20;

 p_quadrature_rule ( N, z, hi );

//printf("%lf \n", z);
//Set up the Angular and Spatial points for the problem
  for (i=0; i<N; i++) {
    
    x[i]= (v-u)/2.0*(z[i])+(v+u)/2.0;
    sx[i]=(pow(2.0*x[i]/5.0,2.0)*pow(2.0*(5.0-x[i])/5.0,2.0));
    w[i]= (v-u)/2.0* hi[i];
//    flux[i]=1;
printf("%lf\n", x[i]);
  }

//x=1.0;
n=1;
i=0;
  for (i=0;i<N;i++){
    
    rx[i]=2.0-Exponential_Integral_En( x[i], n+1)-Exponential_Integral_En( length-x[i], n+1);
    
    for (j=0;j<N; j++){
      if (i != j){
        A[i][j]=-c/2.0*Exponential_Integral_En( fabs(x[i]-x[j]), n)*w[j];
        
      }
      else{
        fun=0;
        
        for (k=0;k<N;k++){    // sums up exponential integrals when i != k
          if (k!=i){
            fun=fun+Exponential_Integral_En(fabs(x[i]-x[k]),n)*w[k];
            
          }
        }
        A[i][j]=1.0-c/2.0*rx[i]+c/2.0*fun;
        

      }

    } 
//printf("%lf\n", rx[i]);
  } 

  for (i=0;i<N;i++){
    s=0.0;
    for (k=0;k<N; k++){
      if (i != k){

        s+=Exponential_Integral_En(fabs(x[i]-x[k]),n)*w[k]*(sx[k]-sx[i]);

      }
    }

sv[i]=sx[i]*rx[i]+s; 

  }
////////////////////////////////////////////////////////////////////////////////////////////////
double newflux[N], count=0;
double oldflux[N];

  for (i=0;i<N;i++){
    newflux[i]=1.0;
    oldflux[i]=1.0;
  }
while(errflux>error){
  for (i=0;i<N;i++){
    fun=0;
    for (j=0;j<N;j++){
      if (i!=j){
         fun +=oldflux[j]*A[i][j];
      }
    }   
      newflux[i]=1.0/A[i][i]*(sv[i]-fun);

      errflux=max(error, newflux[i]-oldflux[i]);
      oldflux[i]=newflux[i];

  }  
//count++;
//printf("%lf \n",count);
      
}
/////////////////////////////////////////////////////////////////////////////////////////////////
FILE *fp;
fp=fopen("Peierls64q.txt","w");
for(i=0;i<N;i++){
fprintf(fp, "For ordinate %lf flux is %lf \n", x[i], newflux[i]);
}

fclose(fp);




  return 0;
}


