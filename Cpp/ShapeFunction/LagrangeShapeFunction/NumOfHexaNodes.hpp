#ifndef NUMOFHEXANODES_HPP
#define NUMOFHEXANODES_HPP
#include <armadillo>
using namespace arma ;
int NumOfHexaNodes(int order)
{
  int NumOfNodes ;
  NumOfNodes = pow((order+1), 3) ;
  return NumOfNodes;
}
#endif // NUMOFHEXANODES_HPP
