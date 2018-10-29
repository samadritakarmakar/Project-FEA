// Automatically translated using m2cpp 2.0 on 2018-10-29 16:56:45

#ifndef SORTANDKEEP_M_HPP
#define SORTANDKEEP_M_HPP

#include "mconvert.h"
#include <armadillo>
using namespace arma ;

TYPE SortAndKeep(TYPE ElementNodes)
{
  TYPE ContainsNodes, columns, row, rows ;
  rows = ElementNodes.n_rows;
  columns = ElementNodes.n_cols;
  
  ContainsNodes = arma::zeros<rowvec>(rows+columns) ;
  for (row=1; row<=rows; row++)
  {
    ContainsNodes.col(0).rows(arma::fspan(columns*row-(columns-1), 1, columns*row)) = ElementNodes.rows(row) ;
  }
  arma::unique(ContainsNodes, ContainsNodes);
  ContainsNodes = ContainsNodes(ContainsNodes!=0) ;
  return ContainsNodes ;
}
#endif
