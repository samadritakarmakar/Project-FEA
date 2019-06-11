#ifndef BUILDTETRAHEDRALMATRIX_HPP
#define BUILDTETRAHEDRALMATRIX_HPP

#endif // BUILDTETRAHEDRALMATRIX_HPP

#include<armadillo>
using namespace arma;
/// Builds Matrix for Tetrahedral.
/// Input order of mesh, x coordinates, y coordinates, z coordinates, NumOfNodes per element,
/// and an empty dense matrix that will be used to store the matrix
void BuildTetrahedralMatrix(int order, mat x, mat y, mat z, int NumOfNodes, mat &Mat)
{
  int lengthx, pos,i, j, k, orderCount;
  lengthx = std::max(x.n_rows,x.n_cols);
  Mat = arma::zeros<mat>(lengthx, NumOfNodes) ;
  //cout<<"Num of Columns of Mat ="<<NumOfNodes<<"\n";
  //Mat.col(0) = arma::ones<vec>(lengthx) ;
  pos = 0;
  for (orderCount=0; orderCount<=order; orderCount++)
  {
    for (i=0; i<=orderCount; i++)
    {
        for(j=0; (i+j)<=orderCount; j++)
        {
            k=orderCount-(i+j);
            Mat.col(pos) = arma::pow(x, i)%arma::pow(y, j)%arma::pow(z, k) ;
            //std::cout<<"x^"<<i<<" y^"<<j<<" z^"<<k<<"\n";
            //std::cout<<pos<<"\n";
            pos = pos+1 ;
        }
    }
  }
}
