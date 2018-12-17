// Automatically translated using m2cpp 2.0 on 2018-12-13 19:11:36

#ifndef VECTORIZEGRADU_M_HPP
#define VECTORIZEGRADU_M_HPP

#include "mconvert.h"
#include <armadillo>
using namespace arma ;

TYPE VectorizeGradu(TYPE graduTemp, TYPE vectorLevel)
{
  TYPE NoOfGraduTempColumns, NoOfGraduTempRows, gradu, i, j, k, startingColumn ;
  NoOfGraduTempRows = size(graduTemp, 1) ;
  NoOfGraduTempColumns = size(graduTemp, 2) ;
  gradu = spalloc(vectorLevel*NoOfGraduTempRows, vectorLevel*NoOfGraduTempColumns, NoOfGraduTempRows*NoOfGraduTempColumns) ;
  j = 1 ;
  k = 1 ;
  for (i=1; i<=vectorLevel*NoOfGraduTempRows; i++)
  {
    startingColumn = k*vectorLevel-(vectorLevel-1) ;
    gradu(i, m2cpp::fspan(startingColumn, 1, startingColumn+(NoOfGraduTempColumns-1))) = graduTemp.rows(j) ;
    j = j+(mod(i, vectorLevel)==0)*1 ;
    k = k*(mod(i, vectorLevel)!=0) ;
    k = k+1 ;
  }
  return gradu ;
}
#endif