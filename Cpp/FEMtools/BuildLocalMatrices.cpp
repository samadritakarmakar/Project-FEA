// Automatically translated using m2cpp 2.0 on 2018-12-13 19:09:50

#ifndef BUILDLOCALMATRICES_M_HPP
#define BUILDLOCALMATRICES_M_HPP

#include <armadillo>
using namespace arma ;

struct _Localmatrices
{
  TYPE LHSmatrixLocal, RHSmatrixLocal, RHSvectorLocal ;
} ;

struct _Localgaussmatrices
{
  TYPE LHSmatrixLocalGauss, RHSmatrixLocalGauss, RHSvectorLocalGauss ;
} ;

void BuildLocalMatrices(_Localgaussmatrices LocalGaussMatrices, _Localmatrices LocalMatrices, TYPE F, TYPE w, TYPE& LHSmatrixLocal, TYPE& RHSmatrixLocal, TYPE& RHSvectorLocal)
{
  TYPE LHSmatrixLocalGauss, RHSmatrixLocalGauss, RHSvectorLocalGauss ;
  _Localgaussmatrices LocalGaussMatrices ;
  _Localmatrices LocalMatrices ;
  LHSmatrixLocal = LocalMatrices.LHSmatrixLocal ;
  RHSmatrixLocal = LocalMatrices.RHSmatrixLocal ;
  RHSvectorLocal = LocalMatrices.RHSvectorLocal ;
  LHSmatrixLocalGauss = LocalGaussMatrices.LHSmatrixLocalGauss ;
  RHSmatrixLocalGauss = LocalGaussMatrices.RHSmatrixLocalGauss ;
  RHSvectorLocalGauss = LocalGaussMatrices.RHSvectorLocalGauss ;
  LHSmatrixLocal = w(1)*w(2)*w(3)*LHSmatrixLocalGauss*det(F)+LHSmatrixLocal ;
  RHSmatrixLocal = w(1)*w(2)*w(3)*RHSmatrixLocalGauss*det(F)+RHSmatrixLocal ;
  RHSvectorLocal = w(1)*w(2)*w(3)*RHSvectorLocalGauss*det(F)+RHSvectorLocal ;
}
#endif