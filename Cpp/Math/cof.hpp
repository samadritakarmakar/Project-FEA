// Automatically translated using m2cpp 2.0 on 2019-04-11 23:48:24

#ifndef COF_M_HPP
#define COF_M_HPP

#include <armadillo>
using namespace arma ;
/// Returns the cofator of a matrix. Accepts the reqd. matrix as input.
 mat cof(mat matrix)
{

  int r = matrix.n_rows;
  int c = matrix.n_cols;
  int i,j;
  mat matrx  = arma::zeros<mat>(r-1, c-1) ;
  mat cofactors = arma::zeros<mat>(r, c) ;
  // Left topmost Corner of matrix
  i=0; j=0;
  matrx =matrix(span(i+1,r-1),span(j+1,c-1));
  cofactors(i, j) = det(matrx)*pow((-1), (i+j));

  // Right topmost Corner of matrix
  i=0; j=c-1;
  matrx=matrix(span(i+1,r-1),span(0,j-1));
  cofactors(i, j) = det(matrx)*pow((-1), (i+j));

  // Bottom leftmost corner
  i=r-1; j=0;
  matrx = matrix(span(0,i-1),span(j+1,c-1));
  cofactors(i, j) = det(matrx)*pow((-1), (i+j));
  // Bottom rightmost corner
  i=r-1; j=c-1;
  matrx=matrix(span(0,i-1),span(0,j-1));
  cofactors(i, j) = det(matrx)*pow((-1), (i+j));

  // Topmost and bottom most row of other than  the corners
  for (j=1; j<c-1; j++)
  {
      i=0;
      matrx=join_horiz(matrix(span(i+1,r-1),span(0,j-1)), matrix(span(i+1,r-1),span(j+1,c-1)));
      cofactors(i, j) = det(matrx)*pow((-1), (i+j));
      i=r-1;
      matrx=join_horiz(matrix(span(0,i-1),span(0,j-1)), matrix(span(0,i-1),span(j+1,c-1)));
      cofactors(i, j) = det(matrx)*pow((-1), (i+j));
  }

  // Leftmost and Rightmost columns other than the corners
  for (i=1; i<r-1; i++)
  {
      j=0;
      matrx=join_vert(matrix(span(0,i-1),span(j+1,c-1)),matrix(span(i+1,r-1),span(j+1,c-1)));
      cofactors(i, j) = det(matrx)*pow((-1), (i+j));
      j=c-1;
      matrx=join_vert(matrix(span(0,i-1),span(0,j-1)),matrix(span(i+1,r-1),span(0,j-1)));
      cofactors(i, j) = det(matrx)*pow((-1), (i+j));
  }

  for (i=1; i<r-1; i++)
  {
    for (j=1; j<c-1; j++)
    {
      matrx =join_vert( join_horiz(matrix(span(0,i-1),span(0,j-1)), matrix(span(0,i-1),span(j+1,c-1))),
                        join_horiz(matrix(span(i+1,r-1),span(0,j-1)), matrix(span(i+1,r-1),span(j+1,c-1))));
      cofactors(i, j) = det(matrx)*pow((-1), (i+j)) ;
    }
  }
  return cofactors ;
}

 mat cofTopRow(mat matrix)
{

  int r = matrix.n_rows;
  int c = matrix.n_cols;
  int i,j;
  mat matrx  = arma::zeros<mat>(r-1, c-1) ;
  mat cofactors = arma::zeros<mat>(1, c) ;
  // Left topmost Corner of matrix
  i=0; j=0;
  matrx =matrix(span(i+1,r-1),span(j+1,c-1));
  cofactors(i, j) = det(matrx)*pow((-1), (i+j));

  // Right topmost Corner of matrix
  i=0; j=c-1;
  matrx=matrix(span(i+1,r-1),span(0,j-1));
  cofactors(i, j) = det(matrx)*pow((-1), (i+j));

  // Topmost and bottom most row of other than  the corners
  for (j=1; j<c-1; j++)
  {
      i=0;
      matrx=join_horiz(matrix(span(i+1,r-1),span(0,j-1)), matrix(span(i+1,r-1),span(j+1,c-1)));
      cofactors(i, j) = det(matrx)*pow((-1), (i+j));
  }
  return cofactors;
}
#endif
