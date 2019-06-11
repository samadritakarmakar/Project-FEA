#ifndef BUILDHEXAHEDRALMATRIX_HPP
#define BUILDHEXAHEDRALMATRIX_HPP
#include<armadillo>
using namespace arma;
/// Builds Matrix for Hexahedral.
/// Input order of mesh, x coordinates, y coordinates, z coordinates, NumOfNodes per element,
/// and an empty dense matrix that will be used to store the matrix
void BuildHexahedralMatrix(int order, mat x, mat y, mat z, int NumOfNodes, mat &Mat)
{
  int lengthx, pos,i, j, k;
  lengthx = std::max(x.n_rows,x.n_cols);
  Mat = arma::zeros<mat>(lengthx, NumOfNodes) ;
  //Mat.col(0) = arma::ones<vec>(lengthx) ;
  pos = 0;
  for (i=0; i<=order; i++)
  {
    for (j=0; j<=order; j++)
    {
        for(k=0; k<=order; k++)
        {
            Mat.col(pos) = arma::pow(x, i)%arma::pow(y, j)%arma::pow(z, k) ; //x^i*y^j*z^k
            pos = pos+1 ;
        }
    }
  }
}

#endif // BUILDHEXAHEDRALMATRIX_HPP
