#ifndef NUMOFLINENODES_HPP
#define NUMOFLINENODES_HPP
#include <armadillo>
using namespace arma ;

int NumOfLineNodes(int order)
{
  int NumOfNodes ;
  NumOfNodes = (order+1) ;
  return NumOfNodes ;
}

#endif // NUMOFLINENODES_HPP
