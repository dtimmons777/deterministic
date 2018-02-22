# include <math.h>
# include <iostream>
# include <stdlib.h>
# include <string.h>
# include <time.h>
# include <eigen/Eigen/Dense>
#include <eigen/Eigen/Eigenvalues> 

using namespace Eigen;


int main ()
{
    int N=2;
    MatrixXd mm(N,N);
    VectorXd eivals(N);
    MatrixXd zkm(N,N);
mm(0,0)=14.0;
mm(1,1)=-15.0;
mm(0,1)=-4.0/7.0;
mm(1,0)=10.0;

 EigenSolver<MatrixXd> es(mm);
 eivals=es.eigenvalues().real();
 zkm=es.eigenvectors().real();

std::cout << mm << std::endl;
std::cout << zkm << std::endl;
std::cout << eivals << std::endl;

}

