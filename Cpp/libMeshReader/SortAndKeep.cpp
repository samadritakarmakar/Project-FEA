// Automatically translated using m2cpp 2.0 on 2018-10-29 16:56:45

#ifndef SORTANDKEEP_M_HPP
#define SORTANDKEEP_M_HPP
#include <iostream>
#include "libMeshReader.h"
#include <armadillo>
using namespace arma ;

mat MeshReader::SortAndKeep(mat ElementNodes)
{
  mat ContainsNodes ;
  int rows = ElementNodes.n_rows;
  int columns = ElementNodes.n_cols;

  ContainsNodes = arma::zeros(1,rows*columns) ;
  for (int row=1; row<=rows; row++)
  {
      //ContainsNodes(1,columns*row-(columns-1):columns*row)=ElementNodes(row,:);
    ContainsNodes(span(0,0),span(columns*row-(columns-1)-1,columns*row-1)) = ElementNodes.row(row-1) ;
  }
  ContainsNodes=unique(ContainsNodes);
  return ContainsNodes;
}
#endif
