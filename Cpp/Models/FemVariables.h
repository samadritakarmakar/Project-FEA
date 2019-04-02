#ifndef FEMVARIABLES_H
#define FEMVARIABLES_H
#include <iostream>
#include <armadillo>

using namespace arma;

class FemVariables
{
public:
    FemVariables();
    ~FemVariables();
    sp_mat *LHSmatrix, *RHSmatrix, *RHSvector;
    int *BuildMatrixList;
    void SetNoOfLHSmatrix(int instances);
    void SetNoOfRHSmatrix(int instances);
    void SetNoOfRHSvector(int instances);
};

#endif // FEMVARIABLES_H
