// Automatically translated using m2cpp 2.0 on 2018-12-13 19:10:32

#ifndef LOADGAUSSFILE_M_HPP
#define LOADGAUSSFILE_M_HPP

#include <armadillo>
using namespace arma ;

struct _Property
{
  TYPE Type, degree ;
} ;

TYPE LoadGaussFile(_Property Property)
{
  TYPE GaussDegree, file ;
  _Property Property ;
  if ((strcmp(Property.Type, "1D")))
  {
    GaussDegree = 2*str2num(Property.degree) ;
    file = strcat("Data2/n", num2str(GaussDegree)) ;
  }
  else if ((strcmp(Property.Type, "Triangle")))
  {
    TYPE _var_TYPE = (Property.degree) ;
    if ("1" == _var_TYPE)
    {
      GaussDegree = 1 ;
    }
    else if ("2" == _var_TYPE)
    {
      GaussDegree = 3 ;
    }
    else if ("3" == _var_TYPE)
    {
      GaussDegree = 4 ;
    }
    file = strcat("DataTriangle2/n", num2str(GaussDegree)) ;
  }
  else if ((strcmp(Property.Type, "Quadrilateral")))
  {
    TYPE _var_TYPE = (Property.degree) ;
    if ("1" == _var_TYPE)
    {
      GaussDegree = 2 ;
    }
    else if ("Biquadratic" == _var_TYPE)
    {
      GaussDegree = 3 ;
    }
    else if ("Serendipity" == _var_TYPE)
    {
      GaussDegree = 3 ;
    }
    file = strcat("Data2/n", num2str(GaussDegree)) ;
  }
  else if ((strcmp(Property.Type, "Tetrahedral")))
  {
    TYPE _var_TYPE = (Property.degree) ;
    if ("1" == _var_TYPE)
    {
      GaussDegree = 1 ;
    }
    else if ("2" == _var_TYPE)
    {
      GaussDegree = 4 ;
    }
    else if ("3" == _var_TYPE)
    {
      GaussDegree = 5 ;
    }
    file = strcat("DataTetrahedral2/n", num2str(GaussDegree)) ;
  }
  else if ((strcmp(Property.Type, "Hexahedral")))
  {
    TYPE _var_TYPE = (Property.degree) ;
    if ("1" == _var_TYPE)
    {
      GaussDegree = 2 ;
    }
    else if ("2" == _var_TYPE)
    {
      GaussDegree = 3 ;
    }
    else if ("3" == _var_TYPE)
    {
      GaussDegree = 5 ;
    }
    file = strcat("Data2/n", num2str(GaussDegree)) ;
  }
  return file ;
}
#endif