// Automatically translated using m2cpp 2.0 on 2019-04-13 23:38:41

#ifndef NUMOFQUADNODES_M_HPP
#define NUMOFQUADNODES_M_HPP

#include <armadillo>

using namespace arma ;

int NumOfQuadNodes(int order)
{
  int NumOfNodes ;
  NumOfNodes = pow((order+1), 2) ;
  return NumOfNodes;
}
#endif
