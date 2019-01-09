#ifndef FEMVARIABLES_H
#define FEMVARIABLES_H
#include <iostream>
#include <armadillo>

using namespace arma;

class FemVariables
{
public:
    mat *LHSmatrix, *RHSmatrix, *RHSvector;
    int *BuildMatrixList;
    FemVariables();
};

#endif // FEMVARIABLES_H
