#ifndef BUILDLINEMATRIX_HPP
#define BUILDLINEMATRIX_HPP

#endif // BUILDLINEMATRIX_HPP

#include <armadillo>
#include <algorithm>
using namespace arma ;
/// Builds Matrix for Line.
/// Input order of mesh, x coordinates, NumOfNodes per element,
/// and an empty dense matrix that will be used to store the matrix
void BuildLineMatrix(int order, mat x, const int& NumOfNodes, mat& Mat)
{
  int i, j, lengthx, pos ;
  lengthx = std::max(x.n_rows,x.n_cols) ;
  Mat = zeros(lengthx, NumOfNodes) ;
  Mat.col(0) = ones(lengthx,1) ;
  pos = 0 ;
  for (i=1; i<=order; i++)
  {
    pos = i;
    Mat.col(pos) = arma::pow(x, i) ;
  }
}
