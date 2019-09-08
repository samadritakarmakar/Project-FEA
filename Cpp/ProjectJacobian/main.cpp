#include <iostream>
#include "f1.hpp"
#include "f2.hpp"
#include "Jacobian.hpp"
int main()
{
    otherdata junk;
    mat x={{1,2,3}};
    mat (*p)(mat, otherdata)=f1 ;
    std::cout<<Jacobian(p,x,junk)<<"\n";
    //mat (*p2)(mat, Otherdata)=f2;
    //Otherdata junk2;
    //std::cout<<JacobianMatrix(p2, x, junk2)<<"\n";
    //sp_mat (*p3)(mat, Otherdata)=f3;
    //std::cout<<mat(JacobianMatrix(p3, x, junk2))<<"\n";
    std::cout<<"d2f_by_dx2\n"<<d2f_by_dx2(p,x,junk)<<"\n";
    return 0;
}
