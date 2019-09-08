#ifndef VARIABLE_HPP
#define VARIABLE_HPP
#include <armadillo>
using namespace arma;
class VariableMatrix
{
public:
    unsigned int nOfRows, nOfCols;
    VariableMatrix(int NoOfInpndntTerms)
    {
        Matrix=std::vector<std::vector<sp_mat>>(NoOfInpndntTerms);
        nOfRows=NoOfInpndntTerms;
        nOfCols=NoOfInpndntTerms;
        for (int row=0; row<Matrix.size(); row++)
        {
            Matrix[row]=std::vector<sp_mat>(NoOfInpndntTerms);
        }
    }
    std::vector<std::vector<sp_mat>> Matrix;
};

class VariableVector
{
public:
    unsigned int nOfRows;
    VariableVector(int NoOfInpndntTerms)
    {
        Vector=std::vector<mat>(NoOfInpndntTerms);
        nOfRows=NoOfInpndntTerms;
    }
    std::vector<mat> Vector;
};

template<class GenericMatrix>
class VariableGeneric
{
 public:
    unsigned int nOfTerms;
    VariableGeneric(int NoOfInpndntTerms)
    {
        Term = std::vector<GenericMatrix>(NoOfInpndntTerms);
        nOfTerms=NoOfInpndntTerms;
    }
    std::vector<GenericMatrix> Term;

    sp_mat operator[](int index)
    {
        return Term[index];
    }
};

#endif // VARIABLE_HPP
