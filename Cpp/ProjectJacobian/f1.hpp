#ifndef F1_HPP
#define F1_HPP
#include <armadillo>
using namespace arma;
class otherdata
{
    //---do nothing----
};

mat f1(mat x, otherdata junk)
{
    mat fnc(3,1);
    vec x1=vectorise(x);
    fnc(0,0)=3.0*pow(x(0),2.0)+2.0*x(1)+x(2);
    fnc(1,0)=pow(x(0),3.0)+3.0*pow(x(1),3.0);
    fnc(2,0)=3.0*x(0)+2.0*(x(1))+pow(x(2),4.0);
    return fnc;
}
#endif // F1_HPP
