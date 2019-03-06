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
    LHSmatrix = new mat [instances];
}

void FemVariables::SetNoOfRHSmatrix(int instances)
{
    RHSmatrix = new mat [instances];
}

void FemVariables::SetNoOfRHSvector(int instances)
{
    RHSvector = new mat [instances];
}
