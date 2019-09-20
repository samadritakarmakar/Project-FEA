#ifndef VECNORM_HPP
#define VECNORM_HPP
#include <armadillo>

using namespace arma;

vec vecnorm(mat Matrix)
{
    vec Vecnorm(Matrix.n_cols);
    for (int colmn=0;colmn<Matrix.n_cols;colmn++)
    {
        Vecnorm.row(colmn)=norm(Matrix.col(colmn));
    }
    return Vecnorm;
}



#endif
