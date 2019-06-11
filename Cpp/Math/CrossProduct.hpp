#ifndef CROSSPRODUCT_HPP
#define CROSSPRODUCT_HPP
#include <armadillo>
#include "cof.hpp"

using namespace arma;
 mat CrossProduct(mat Matrix1, mat Matrix2)
{
    vec mtrx1=vectorise(Matrix1);
    vec mtrx2=vectorise(Matrix2);
    if (mtrx1.n_rows!=mtrx2.n_rows && mtrx1.n_cols!=mtrx2.n_cols)
    {
        cout<<"Matrix1 and Matrix2 should have the same size for Cross Product\n";
        throw;
    }
    vec ons =vectorise(ones(mtrx1.n_rows, mtrx1.n_cols));
    mat temp=join_vert(ons.t(),mtrx1.t());
    mat crossProduct=cofTopRow(join_vert(temp,mtrx2.t()));
    return crossProduct;
}
#endif // CROSSPRODUCT_HPP
