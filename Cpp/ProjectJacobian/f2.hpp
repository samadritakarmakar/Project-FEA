#ifndef F2_HPP
#define F2_HPP
#include <armadillo>
using namespace arma;
class Otherdata
{
    //---do nothing----
};

mat f2(mat x, Otherdata junk);
sp_mat f3(mat x, Otherdata junk);

mat f2(mat x, Otherdata junk)
{
    mat fnc(3,3);
    vec x1=vectorise(x);

    fnc(0,0)=3.0*pow(x(0),2.0)+2.0*x(1)+x(2); fnc(0,1)=4*pow(x(0),3)+16*pow(x(1),2)+3*x(2);
    fnc(0,2)=5*pow(x(0),3)+20*x(1)+pow(x(2),3);

    fnc(1,0)=pow(x(0),3.0)+3.0*pow(x(1),3.0); fnc(1,1)=2*pow(x(0),4)+13*x(1)+x(2);
    fnc(1,2)=13*pow(x(0),2)+16*pow(x(1),2)+pow(x(2),2);

    fnc(2,0)=3.0*x(0)+2.0*(x(1))+pow(x(2),4.0); fnc(2,1)=21*pow(x(0),2)+6*pow(x(1),3)+14*pow(x(2),2);
    fnc(2,2)=16*pow(x(0),3)+4*pow(x(1),2)+pow(x(2),4);
    return fnc;
}

sp_mat f3(mat x, Otherdata junk)
{
    sp_mat fnc(3,3);
    vec x1=vectorise(x);
    fnc(0,0)=3.0*pow(x(0),2.0)+2.0*x(1)+x(2);
    fnc(1,1)=2*pow(x(0),4)+13*x(1)+x(2);
    fnc(2,2)=16*pow(x(0),3)+4*pow(x(1),2)+pow(x(2),4);
    return fnc;
}
#endif // F2_HPP
