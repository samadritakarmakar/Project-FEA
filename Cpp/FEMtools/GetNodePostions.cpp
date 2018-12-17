// Automatically translated using m2cpp 2.0 on 2018-12-13 19:39:15

#ifndef GETNODEPOSTIONS_M_HPP
#define GETNODEPOSTIONS_M_HPP

#include "mconvert.h"
#include <armadillo>
using namespace arma ;

TYPE GetNodePostions(TYPE ElementNodes, TYPE vectorLevel)
{
  TYPE NodePositons, column, i, j, row ;
  row = size(ElementNodes, 1) ;
  column = size(ElementNodes, 2) ;
  NodePositons = arma::zeros<mat>(row, vectorLevel*column) ;
  for (i=1; i<=column; i++)
  {
    for (j=0; j<=vectorLevel-1; j++)
    {
      NodePositons.cols(vectorLevel*i-j) = vectorLevel*ElementNodes.cols(i)-j*arma::ones<vec>(size(row, 2)) ;
    }
  }
  return NodePositons ;
}
#endif