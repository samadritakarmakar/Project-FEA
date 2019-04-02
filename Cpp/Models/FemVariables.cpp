#include "FemVariables.h"

FemVariables::FemVariables()
{

}
FemVariables::~FemVariables()
{
    delete []LHSmatrix;
    delete []RHSmatrix;
    delete []RHSvector;
}

void FemVariables::SetNoOfLHSmatrix(int instances)
{
    LHSmatrix = new sp_mat [instances];
}

void FemVariables::SetNoOfRHSmatrix(int instances)
{
    RHSmatrix = new sp_mat [instances];
}

void FemVariables::SetNoOfRHSvector(int instances)
{
    RHSvector = new sp_mat [instances];
}
