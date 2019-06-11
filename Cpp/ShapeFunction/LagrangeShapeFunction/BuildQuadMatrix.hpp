// Automatically translated using m2cpp 2.0 on 2019-04-13 23:47:05

#ifndef BUILDQUADMATRIX_M_HPP
#define BUILDQUADMATRIX_M_HPP
#include <armadillo>
using namespace arma ;
/// Builds Matrix for Quadrilateral.
/// Input order of mesh, x coordinates, y coordinates, NumOfNodes per element,
/// and an empty dense matrix that will be used to store the matrix
void BuildQuadMatrix(int order, mat x, mat y, int NumOfNodes, mat &Mat)
{
  int lengthx, pos,i, j;;
  lengthx = std::max(x.n_rows,x.n_cols);
  Mat = arma::zeros<mat>(lengthx, NumOfNodes) ;
  Mat.col(0) = arma::ones<vec>(lengthx) ;
  pos = 1;
  for (i=1; i<=order; i++)
  {
    Mat.col(pos) = arma::pow(x, i) ; //Builds x; x^2;... etc
    Mat.col(pos+1) = arma::pow(y, i) ; //Builds y; y^2;... etc
    pos = pos+2 ;
    for (j=1; j<=order; j++)
    {
      Mat.col(pos) = arma::pow(x, i)%arma::pow(y, j) ; //Builds xy; xy, xy^2, x^2y, x^2y^2;.. etc
      pos = pos+1 ;
    }
  }
}
#endif
