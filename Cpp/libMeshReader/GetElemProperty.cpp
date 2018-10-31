// Automatically translated using m2cpp 2.0 on 2018-10-30 17:05:13

#ifndef GETELEMPROPERTY_M_HPP
#define GETELEMPROPERTY_M_HPP

#include <armadillo>
#include<string>
#include "libMeshReader.h"
using namespace arma ;
using namespace std;
void GetElemProperty(int GmshElementNum, string& Type, string& Degree, int& NumOfElementNodes, int& NumofDimensions, int& Supported)
{
  switch (GmshElementNum)
  {
  case 1:

    Type = "1D" ;
    Degree = "1" ;
    NumOfElementNodes = 2 ;
    NumofDimensions = 1 ;
    Supported = 1 ;
    break;

  case 2:

    Type = "Triangle" ;
    Degree = "1" ;
    NumOfElementNodes = 3 ;
    NumofDimensions = 2 ;
    Supported = 1 ;
    break;

  case 3:

    Type = "Quadrilateral" ;
    Degree = "1" ;
    NumOfElementNodes = 4 ;
    NumofDimensions = 2 ;
    Supported = 1 ;
    break;

  case 4:

    Type = "Tetrahedral" ;
    Degree = "1" ;
    NumOfElementNodes = 4 ;
    NumofDimensions = 3 ;
    Supported = 1 ;
    break;

  case 5:
    Type = "Hexahedral" ;
    Degree = "1" ;
    NumOfElementNodes = 8 ;
    NumofDimensions = 3 ;
    Supported = 1 ;
    break;

  case 8:

    Type = "1D" ;
    Degree = "2" ;
    NumOfElementNodes = 3 ;
    NumofDimensions = 1 ;
    Supported = 1 ;
    break;

  case 9:

    Type = "Triangle" ;
    Degree = "2" ;
    NumOfElementNodes = 6 ;
    NumofDimensions = 2 ;
    Supported = 1 ;
    break;

  case 10:

    Type = "Quadrilateral" ;
    Degree = "Biquadratic" ;
    NumOfElementNodes = 9 ;
    NumofDimensions = 2 ;
    Supported = 1 ;
    break;

  case 11:

    Type = "Tetrahedral" ;
    Degree = "2" ;
    NumOfElementNodes = 10 ;
    NumofDimensions = 3 ;
    Supported = 1 ;
    break;

  case 16:

    Type = "Quadrilateral" ;
    Degree = "Serendipity" ;
    NumOfElementNodes = 8 ;
    NumofDimensions = 2 ;
    Supported = 1 ;
    break;

  case 21:

    Type = "Triangle" ;
    Degree = "3" ;
    NumOfElementNodes = 10 ;
    NumofDimensions = 2 ;
    Supported = 1 ;
    break;

  case 26:

    Type = "1D" ;
    Degree = "3" ;
    NumOfElementNodes = 4 ;
    NumofDimensions = 1 ;
    Supported = 1 ;
    break;

  default:

    Type = "" ;
    Degree = "" ;
    NumOfElementNodes = 0 ;
    NumofDimensions = 0 ;
    Supported = 0 ;
    break;

  }
}
#endif
