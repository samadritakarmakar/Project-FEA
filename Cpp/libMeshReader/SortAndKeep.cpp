// Automatically translated using m2cpp 0.5 on 2018-10-28 02:38:18
            
#ifndef SORTANDKEEP_M_HPP
#define SORTANDKEEP_M_HPP

#include "/usr/include/armadillo"
using namespace arma ;
TYPE SortAndKeep(TYPE ElementNodes)
{
  mat ContainsNodes, column, columns, first, flag, last, middle, row, rows ;
  rows = ElementNodes.n_rows ;
  columns = ElementNodes.n_cols ;
  for (row=1; row<=rows; row++)
  {
    for (column=1; column<=columns; column++)
    {
      if ((row==1 && column==1))
      {
        ContainsNodes(1) = ElementNodes(row, column) ;
        ContainsNodes(2) = ElementNodes(row, column+1) ;
      }
      flag = 1 ;
      first = 1 ;
      last = m2cpp::length(ContainsNodes) ;
      while (flag)
      {
        middle = int32((first+last)/2.0) ;
        if ((ContainsNodes(middle-1)<ElementNodes(row, column) && ContainsNodes(middle)<ElementNodes(row, column)))
        {
          first = middle ;
          if ((middle==last))
          {
            ContainsNodes = {arma::join_rows(ContainsNodes, ElementNodes(row, column))} ;
            flag = 0 ;
          }
        }
        else if ((ContainsNodes(middle-1)>ElementNodes(row, column) && ContainsNodes(middle)>ElementNodes(row, column)))
        {
          last = middle ;
          if ((middle-1==first))
          {
            ContainsNodes = {arma::join_rows(ElementNodes(row, column), ContainsNodes)} ;
            flag = 0 ;
          }
        }
        else if ((ContainsNodes(middle-1)<ElementNodes(row, column) && ContainsNodes(middle)>ElementNodes(row, column)))
        {
          ContainsNodes = {arma::join_rows(arma::join_rows(ContainsNodes(m2cpp::fspan(1, 1, middle-1)), ElementNodes(row, column)), ContainsNodes(m2cpp::fspan(middle, 1, ContainsNodes.n_rows)))} ;
          flag = 0 ;
        }
        else if ((ContainsNodes(middle-1)==ElementNodes(row, column) || ContainsNodes(middle)==ElementNodes(row, column)))
        {
          flag = 0 ;
        }
      }
    }
  }
  return ContainsNodes ;
}
#endif
