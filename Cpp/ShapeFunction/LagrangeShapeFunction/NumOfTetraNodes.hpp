// Automatically translated using m2cpp 2.0 on 2019-04-15 21:08:17

#ifndef NUMOFTETRANODES_M_HPP
#define NUMOFTETRANODES_M_HPP

int NumOfTetraNodes(int order)
{
  int NumOfNodes ;
   //(n(n+1)(2n+1)/6 + 3n(n+1) + 2n)/2 where n=order
  NumOfNodes = (order*(order+1)*(2*order+1))/6+3*(order*(order+1))/2+2*order;
  NumOfNodes = NumOfNodes/2;
  NumOfNodes = NumOfNodes+1;
  return NumOfNodes ;
}
#endif
