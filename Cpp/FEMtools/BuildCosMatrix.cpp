// Automatically translated using m2cpp 2.0 on 2018-12-13 19:09:24

#ifndef BUILDCOSMATRIX_M_HPP
#define BUILDCOSMATRIX_M_HPP


#include <armadillo>
#include "libGmshReader.h"
using namespace arma ;


void BuildCosMatrix(mat xCoords, mat yCoords, mat zCoords, libGmshReader::ElementData ElementData1, int i, int vectorLevel, mat& cosMatrix, mat& xVector)
{
    mat Coords, CoordsTempStore, DiffVector, ElementLength, SquareOfDiffVector, SumOfSquare, cos_xyz, xDiff, yDiff, zDiff ;
    int Row;
    xDiff = xCoords(xCoords.n_rows-1,0)-xCoords(0,0) ;
    yDiff = yCoords(yCoords.n_rows-1,0)-yCoords(0,0) ;
    zDiff = zCoords(zCoords.n_rows-1,0)-zCoords(0,0) ;
    DiffVector=join_cols( xDiff,join_rows(yDiff,zDiff));
    SquareOfDiffVector = arma::square(DiffVector) ;
    SumOfSquare = sum(SquareOfDiffVector(span(0,vectorLevel), 0)) ;
    ElementLength = sqrt(SumOfSquare) ;
    cos_xyz = trans(DiffVector(0,span(0, vectorLevel)))/ElementLength ;
    cosMatrix = zeros<mat>(ElementData1.NumOfElementNodes[i], vectorLevel*ElementData1.NumOfElementNodes[i]) ;
    xVector = zeros<vec>(vectorLevel*ElementData1.NumOfElementNodes[i]) ;
    CoordsTempStore = {join_rows(join_rows(xCoords, yCoords), zCoords)} ;
    for (Row=0; Row<ElementData1.NumOfElementNodes[i]; Row++)
    {
      Coords = arma::trans(CoordsTempStore(Row, span( 0, vectorLevel))) ;
      cosMatrix(Row, span((vectorLevel*Row-(vectorLevel-1)), (vectorLevel*Row))) = cos_xyz ;
      xVector(span((vectorLevel*Row-(vectorLevel-1)), (vectorLevel*Row)),0)= Coords ;
    }

}
#endif
