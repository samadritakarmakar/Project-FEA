#ifndef EXPRESSION_HPP
#define EXPRESSION_HPP
#include <armadillo>
#include "TrialFunction.hpp"
using namespace arma;
class Expression
{
public:
    int vctrLvlInternal;
    Expression(int vectorLevel):
        vctrLvlInternal(vectorLevel)
    {
    }

    virtual vec Eval(vec& x)
    {
        vec Exprssn=zeros(vctrLvlInternal,1);
        return Exprssn;
    }
private:

};

#endif // EXPRESSION_HPP
