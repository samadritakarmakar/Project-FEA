// Automatically translated using m2cpp 2.0 on 2019-04-11 23:33:14

#ifndef NUMOFTRINODES_M_HPP
#define NUMOFTRINODES_M_HPP

#include <armadillo>
using namespace arma ;

int NumOfTriNodes(int order)
{
  int NumOfNodes ;
  NumOfNodes = ((order+1)*(order+2))/2.0 ;
  return NumOfNodes ;
}
#endif
