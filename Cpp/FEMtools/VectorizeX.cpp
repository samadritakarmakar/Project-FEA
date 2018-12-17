// Automatically translated using m2cpp 2.0 on 2018-12-13 19:12:00

#ifndef VECTORIZEX_M_HPP
#define VECTORIZEX_M_HPP

#include "mconvert.h"
#include <armadillo>
using namespace arma ;

struct _Otherdata
{
  TYPE Property, xCoords, yCoords, zCoords ;
} ;

TYPE VectorizeX(_Otherdata OtherData, TYPE vectorLevel)
{
  TYPE Property, VectorizedX, VectorizedX_Temp, i, xCoords, yCoords, zCoords ;
  _Otherdata OtherData ;
  Property = OtherData.Property ;
  xCoords = OtherData.xCoords ;
  yCoords = OtherData.yCoords ;
  zCoords = OtherData.zCoords ;
  for (i=1; i<=m2cpp::length(xCoords); i++)
  {
    VectorizedX_Temp = {xCoords(i), yCoords(i), zCoords(i)} ;
    VectorizedX.row(m2cpp::fspan(vectorLevel*i-(vectorLevel-1), 1, vectorLevel*i)).cols(0) = VectorizedX_Temp(m2cpp::fspan(1, 1, vectorLevel), 1) ;
  }
  return VectorizedX ;
}
#endif