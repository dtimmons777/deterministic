# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>


main ( )
{

int N=3;
double A[N][N], determ, I[N][N], v[N],d[N];
int i,j,k,l;

for (i=0;i<N;i++){
   for(j=0;j<N;j++){
      if (i==j){
         I[i][j]=1;
         A[i][j]=2;
      }
      else{
         I[i][j]=0; 
         A[i][j]=1;
      }
   }
}
////////////////////////////////////////////////////////////
double fun,funn, funnn, P[N][N][N], R[N][N],D;
while(.0001<fabs(A[0][N-1])){
   for (i=0;i<N;i++){
      fun=0;
      for (j=0;j<N;j++){
         fun+=pow(A[j][i],2.0);
      }
      for (j=0;j<N;j++){//d[j]
         d[j]=A[j][i]/sqrt(fun);
      }
      j=0;
      if(i>j){
         while(i>j){//v[j] for j<i ==0
             v[j]=0;
             j=j+1;
         }
      }
      funn=0;   
      for (j=i;j<N;j++){// determinate
          funn+=pow(d[j],2.0);

      }
      D=-A[i][i]/fabs(A[i][i])*sqrt(funn);
      for (j=i;j<N;j++){
         if(i==j){
             v[j]=sqrt(0.5*(1.0-d[j]/D));
         }
         else{
             v[j]=d[j]/(2.0*-D*v[i]);
         }
      }
      for (j=0;j<N;j++){
          for (k=0;k<N;k++){
              P[i][j][k]=I[j][k]-2.0*v[j]*v[k];
          }
      }
      for (j=0;j<N;j++){
          for (k=0;k<N;k++){
              funnn=0;
              for (l=0;l<N;l++){
                  funnn+=P[i][j][l]*A[l][k];
              }
              R[j][k]=funnn;
printf(" %lf ", R[j][k]);
          }
printf("\n");
      }

      for (j=0;j<N;j++){
          for (k=0;k<N;k++){      
             A[j][k]=R[j][k];

          }
      } 
   }
}


return 0;

}
