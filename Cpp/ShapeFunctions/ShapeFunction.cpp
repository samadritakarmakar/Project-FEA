// Automatically translated using m2cpp 2.0 on 2018-12-13 19:10:47

#ifndef SHAPEFUNCTION_M_HPP
#define SHAPEFUNCTION_M_HPP

#include <armadillo>
#include <string>
using namespace arma ;

struct _Property
{
  std::string Type, degree ;
} ;

vec ShapeFunction(vec epsilon, _Property Property)
{
  vec   phi;
  double gamma, zeta;
  std::string degree, type;
  type = Property.Type ;
  degree = Property.degree ;
  if ("1D" == type)
  {
    if ("1" == degree)
    {
      phi = arma::zeros<colvec>(2) ;
      phi(1) = -0.5*(epsilon(1)-1) ;
      phi(2) = 0.5*(epsilon(1)+1) ;
    }
    else if ("2" == degree)
    {
      phi = arma::zeros<colvec>(3) ;
      phi(1) = 0.5*epsilon(1)*(epsilon(1)-1) ;
      phi(2) = -(epsilon(1)+1)*(epsilon(1)-1) ;
      phi(3) = 0.5*epsilon(1)*(epsilon(1)+1) ;
    }
    else if ("3" == degree)
    {
      phi = arma::zeros<colvec>(4) ;
      phi(1) = -9/16.0*(epsilon(1)+1/3.0)*(epsilon(1)-1/3.0)*(epsilon(1)-1) ;
      phi(2) = 27/16.0*(epsilon(1)+1)*(epsilon(1)-1/3.0)*(epsilon(1)-1) ;
      phi(3) = -27/16.0*(epsilon(1)+1)*(epsilon(1)+1/3.0)*(epsilon(1)-1) ;
      phi(4) = -9/16.0*(epsilon(1)+1)*(epsilon(1)+1/3.0)*(epsilon(1)-1/3.0) ;
    }
    else
    {
      printf("Shape Function Degree specified not found") ;
    }
  }
  else if ("Triangle" == type)
  {
    zeta = 1-(epsilon(1)+epsilon(2)) ;
    if ("1" == degree)
    {
      phi = arma::zeros<colvec>(3) ;
      phi(1) = epsilon(1) ;
      phi(2) = epsilon(2) ;
      phi(3) = (1-zeta) ;
    }
    else if ("2" == degree)
    {
      phi = arma::zeros<colvec>(6) ;
      phi(1) = epsilon(1)*(2*epsilon(1)-1) ;
      phi(2) = epsilon(2)*(2*epsilon(2)-1) ;
      phi(3) = zeta*(2*zeta-1) ;
      phi(4) = 4*epsilon(1)*epsilon(2) ;
      phi(5) = 4*epsilon(2)*zeta ;
      phi(6) = 4*epsilon(1)*zeta ;
    }
    else if ("3" == degree)
    {
      phi = arma::zeros<colvec>(10) ;
      phi(1) = 0.5*epsilon(1)*(3*epsilon(1)-1)*(3*epsilon(1)-2) ;
      phi(2) = 0.5*epsilon(2)*(3*epsilon(2)-1)*(3*epsilon(2)-2) ;
      phi(3) = 0.5*zeta*(3*zeta-1)*(3*zeta-2) ;
      phi(4) = 4.5*epsilon(1)*epsilon(2)*(3*epsilon(1)-1) ;
      phi(5) = 4.5*epsilon(1)*epsilon(2)*(3*epsilon(2)-1) ;
      phi(6) = 4.5*epsilon(2)*zeta*(3*epsilon(2)-1) ;
      phi(7) = 4.5*epsilon(2)*zeta*(3*zeta-1) ;
      phi(8) = 4.5*epsilon(1)*zeta*(3*zeta-1) ;
      phi(9) = 4.5*epsilon(1)*zeta*(3*epsilon(1)-1) ;
      phi(10) = 4.5*epsilon(1)*epsilon(2)*zeta ;
    }
    else
    {
      printf("Shape Function Degree specified not found") ;
    }
  }
  else if ("Quadrilateral" == type)
  {
    if ("1" == degree)
    {
      phi = arma::zeros<colvec>(4) ;
      phi(1) = 0.25*(1-epsilon(1))*(1-epsilon(2)) ;
      phi(2) = 0.25*(1+epsilon(1))*(1-epsilon(2)) ;
      phi(3) = 0.25*(1+epsilon(1))*(1+epsilon(2)) ;
      phi(4) = 0.25*(1-epsilon(1))*(1+epsilon(2)) ;
    }
    else if ("Biquadratic" == degree)
    {
      phi = arma::zeros<colvec>(9) ;
      phi(1) = 0.25*epsilon(1)*epsilon(2)*(epsilon(1)-1)*(epsilon(2)-1) ;
      phi(2) = -0.25*epsilon(1)*epsilon(2)*(epsilon(1)+1)*(epsilon(2)-1) ;
      phi(3) = 0.25*epsilon(1)*epsilon(2)*(epsilon(1)+1)*(epsilon(2)+1) ;
      phi(4) = -0.25*epsilon(1)*epsilon(2)*(epsilon(1)-1)*(epsilon(2)+1) ;
      phi(5) = -0.5*epsilon(2)*(epsilon(2)-1)*(pow(epsilon(1), 2)-1) ;
      phi(6) = -0.5*epsilon(1)*(epsilon(1)+1)*(pow(epsilon(2), 2)-1) ;
      phi(7) = -0.5*epsilon(2)*(epsilon(2)+1)*(pow(epsilon(1), 2)-1) ;
      phi(8) = -0.5*epsilon(1)*(epsilon(1)-1)*(pow(epsilon(2), 2)-1) ;
      phi(9) = (1-pow(epsilon(1), 2))*(1-pow(epsilon(2), 2)) ;
    }
    else if ("Serendipity" == degree)
    {
      phi = arma::zeros<colvec>(8) ;
      phi(1) = -0.25*(epsilon(1)-1)*(epsilon(2)-1)*(epsilon(1)+epsilon(2)+1) ;
      phi(2) = -0.25*(epsilon(1)+1)*(epsilon(2)-1)*(epsilon(1)-epsilon(2)-1) ;
      phi(3) = -0.25*(epsilon(1)+1)*(1+epsilon(2))*(1-epsilon(1)-epsilon(2)) ;
      phi(4) = -0.25*(1-epsilon(1))*(1+epsilon(2))*(1+epsilon(1)-epsilon(2)) ;
      phi(5) = 0.5*(1-epsilon(1))*(1+epsilon(1))*(1-epsilon(2)) ;
      phi(6) = 0.5*(1+epsilon(1))*(1+epsilon(2))*(1-epsilon(2)) ;
      phi(7) = 0.5*(1-epsilon(1))*(1+epsilon(1))*(1+epsilon(2)) ;
      phi(8) = 0.5*(1-epsilon(1))*(1+epsilon(2))*(1-epsilon(2)) ;
    }
    else
    {
      printf("Shape Function Degree specified not found") ;
    }
  }
  else if ("Tetrahedral" == type)
  {
    gamma = 1-(epsilon(1)+epsilon(2)+epsilon(3)) ;
    if ("1" == degree)
    {
      phi = arma::zeros<colvec>(4) ;
      phi(1) = epsilon(1) ;
      phi(2) = epsilon(2) ;
      phi(3) = epsilon(3) ;
      phi(4) = gamma ;
    }
    else if ("2" == degree)
    {
      phi = arma::zeros<colvec>(10) ;
      phi(1) = epsilon(1)*(2*epsilon(1)-1) ;
      phi(2) = epsilon(2)*(2*epsilon(2)-1) ;
      phi(3) = epsilon(3)*(2*epsilon(3)-1) ;
      phi(4) = gamma*(2*gamma-1) ;
      phi(5) = 4*epsilon(1)*epsilon(2) ;
      phi(6) = 4*epsilon(2)*epsilon(3) ;
      phi(7) = 4*epsilon(3)*epsilon(1) ;
      phi(8) = 4*epsilon(1)*gamma ;
      phi(9) = 4*epsilon(2)*gamma ;
      phi(10) = 4*epsilon(3)*gamma ;
    }
    else
    {
      printf("Shape Function Degree specified not found") ;
    }
  }
  else if ("Hexahedral" == type)
  {
    if ("1" == degree)
    {
      phi = arma::zeros<colvec>(8) ;
      phi(1) = 1/8.0*(1-epsilon(1))*(1-epsilon(2))*(1-epsilon(3)) ;
      phi(2) = 1/8.0*(1+epsilon(1))*(1-epsilon(2))*(1-epsilon(3)) ;
      phi(3) = 1/8.0*(1+epsilon(1))*(1+epsilon(2))*(1-epsilon(3)) ;
      phi(4) = 1/8.0*(1-epsilon(1))*(1+epsilon(2))*(1-epsilon(3)) ;
      phi(5) = 1/8.0*(1-epsilon(1))*(1-epsilon(2))*(1+epsilon(3)) ;
      phi(6) = 1/8.0*(1+epsilon(1))*(1-epsilon(2))*(1+epsilon(3)) ;
      phi(7) = 1/8.0*(1+epsilon(1))*(1+epsilon(2))*(1+epsilon(3)) ;
      phi(8) = 1/8.0*(1-epsilon(1))*(1+epsilon(2))*(1+epsilon(3)) ;
    }
    else
    {
      printf("Shape Function Degree specified not found") ;
    }
  }
  return phi ;
}
#endif
