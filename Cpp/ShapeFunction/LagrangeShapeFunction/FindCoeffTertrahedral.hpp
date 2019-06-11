#ifndef FINDCOEFFTERTRAHEDRAL_HPP
#define FINDCOEFFTERTRAHEDRAL_HPP

#include <armadillo>
#include "GmshNodeListGen.h"
#include "NumOfTetraNodes.hpp"
#include "BuildTetrahedralMatrix.hpp"
#include "cof.hpp"
#include "CheckSupportedOrder.hpp"

using namespace arma;

/// This code calculates the coefficients for each node equation
/// Eg. Coeff(1,:) gives the set of Coefficients of node 1.
/// Eg. Coeff(2,:) gives the set of Coefficients of node 2.
///
/// To use the resulting coefficients you may do the following.
/// Step 1: Determine the points x & y you want to find the shape functions for.
/// Step 2: Determine the order of the shape function.
/// Step 3: Call function "Mat=Build<Triangle, Quad, etc>Matrix(order, x, y, NumOfNodes)"
/// Step 4: Multiply Mat with Transpose(CoeffMatrix) as ,  N=Mat*CoeffMatrix.t()
///          N(1,1), N(1,2), N(1,3), ... with represent the different shape functions
/// If COLUMN vectors are passed in x and y then N(:,1), N(:,2), N(:,3), ... will
/// represent the different shape functions at different values of x and y.

mat FindCoeffTetrahedral(int order)
{
  mat CoeffMatrix, CofMatrix, Mat, NodeList, x, y, z ;
  int orderLimit=5;
  CheckSupportedOrder(order, orderLimit);
  NodeList = GmshNodeListTetrahedral(order) ;
  x = NodeList.col(0) ;
  y = NodeList.col(1) ;
  z = NodeList.col(2) ;
  int NumOfNodes = NumOfTetraNodes(order) ;
  cout<<"Num of Tertra Nodes="<<NumOfNodes<<"\n";
  BuildTetrahedralMatrix(order, x, y, z, NumOfNodes, Mat) ;
  CofMatrix = cof(Mat) ;
  CoeffMatrix = speye(NumOfNodes,NumOfNodes)*cof(Mat)/det(Mat) ;
  return CoeffMatrix ;
}
#endif // FINDCOEFFTERTRAHEDRAL_HPP
