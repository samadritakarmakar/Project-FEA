// Automatically translated using m2cpp 2.0 on 2019-04-11 23:33:14

#ifndef FINDCOEFFTRIANGLE_M_HPP
#define FINDCOEFFTRIANGLE_M_HPP
#include "GmshNodeListGen.h"
#include "cof.hpp"
#include "BuildTriangleMatrix.hpp"
#include "NumOfTriNodes.hpp"
#include "CheckSupportedOrder.hpp"
#include <armadillo>
using namespace arma ;

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

mat FindCoeffTriangle(int order)
{
  mat CoeffMatrix, CofMatrix, Mat, NodeList, x, y ;
  int NumOfNodes;
  int orderLimit=9;
  CheckSupportedOrder(order, orderLimit);
  NodeList = GmshNodeListTriangle(order) ;
  x = NodeList.col(0) ;
  y = NodeList.col(1) ;
  NumOfNodes = NumOfTriNodes(order) ;
  BuildTriangleMatrix(order, x, y, NumOfNodes, Mat) ;
  CofMatrix = cof(Mat) ;
  CoeffMatrix = speye(NumOfNodes,NumOfNodes)*cof(Mat)/det(Mat) ;
  return CoeffMatrix ;
}
#endif
