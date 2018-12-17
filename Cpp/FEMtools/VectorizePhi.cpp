// Automatically translated using m2cpp 2.0 on 2018-12-13 19:11:54

#ifndef VECTORIZEPHI_M_HPP
#define VECTORIZEPHI_M_HPP

#include "mconvert.h"
#include <armadillo>
using namespace arma ;

TYPE VectorizePhi(TYPE phiTemp, TYPE vectorLevel)
{
  TYPE NoOfPhiColumns, endColumn, i, phi, startColumn ;
  NoOfPhiColumns = size(phiTemp, 2) ;
  phi = spalloc(vectorLevel, NoOfPhiColumns*vectorLevel, NoOfPhiColumns*vectorLevel) ;
  for (i=1; i<=NoOfPhiColumns; i++)
  {
    startColumn = i*vectorLevel-(vectorLevel-1) ;
    endColumn = startColumn+vectorLevel-1 ;
    phi(m2cpp::fspan(1, 1, vectorLevel), m2cpp::fspan(startColumn, 1, endColumn)) = sparse(bsxfun(@times, speye(vectorLevel), phiTemp(i))) ;
  }
  return phi ;
}
#endif