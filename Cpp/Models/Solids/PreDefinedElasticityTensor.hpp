#ifndef PREDEFINEDELASTICITYTENSOR_HPP
#define PREDEFINEDELASTICITYTENSOR_HPP
#include <armadillo>
using namespace arma;
void Set_C_3D(sp_mat& C, double E, double nu=0.3)
{
    double mu=E/(2.0*(1.0+nu));
    double lambda=(nu*E)/((1.0+nu)*(1.0-2.0*nu));
    C.set_size(6,6);
    C(0, 0) = lambda + 2.0*mu;
    C(1, 1) = lambda + 2.0*mu;
    C(2, 2) = lambda + 2.0*mu;
    C(3, 3) = mu;
    C(4, 4) = mu;
    C(5, 5) = mu;
    C(0, 1) = lambda;
    C(0, 2) = lambda;
    C(1, 0) = lambda;
    C(1, 2) = lambda;
    C(2, 0) = lambda;
    C(2, 1) = lambda;

}
#endif // PREDEFINEDELASTICITYTENSOR_HPP
