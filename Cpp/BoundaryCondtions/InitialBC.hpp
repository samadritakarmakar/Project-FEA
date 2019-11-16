#ifndef INITIALBC_HPP
#define INITIALBC_HPP
#include "DirichletBC.hpp"
class InitialBC : protected DirichletBC
{
public:
    InitialBC(TrialFunction &u, int PhysicalGroupNumber, umat& boolInitialBCNodes):
        DirichletBC (u, PhysicalGroupNumber, boolInitialBCNodes)
    {}
    void SetInitialBCExpression(Expression& DeclaredExprssnClss)
    {
        SetDirichletBCExpression(DeclaredExprssnClss);
    }
    void ApplyInitialBC(mat& Sol_of_u)
    {
        for (int ElementType=0; ElementType<Msh->NumOfElementTypes; ElementType++)
        {
            umat uniquePositions;
            SetBC_on_b_n_GetUniquePositions(Sol_of_u, uniquePositions, ElementType);
        }

    }
};


#endif // INITIALBC_HPP
