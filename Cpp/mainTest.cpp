/*#include<stdio.h>
#include <string>
#include <armadillo>
//#include "libMeshReader/libMeshReader.h"
//#include "libGmshReader/libGmshReader.h"
#include "libGmshReader.h"
#include "FemModule.h"*/

#include "ProjectFea.h"

//using namespace arma;
int main()
{
    mat b,c,d;
    d=zeros(3,1);
    b={{1,2,3},{5,6,7},{7,10,9}};
    c<<1<<endr
    <<2<<endr
    <<3<<endr;
    d=solve(b,c);
    std::cout<<d;

    std::cout<<"Enter File Name: ";
    std::string s;
    std::cin>>s;
    std::cout<<"Enter Dimension: ";
    int D;
    std::cin>>D;
    libGmshReader::MeshReader Test(s,D);
    //int success=Test.GetNodeData();
    //Test.GetElementData();
    if(Test.success)
    {
        std::cout<<"Number of Nodes= "<<Test.NumOfNodes<<"\n";
        // std::cout<<"Node Tags= \n"<<Test.NodeTag;
        //cout<<"Element Nodes = \n"<<Test.ElementNodes[0];
        //cout<<"Gmsh Element Nodes = \n"<<Test.GmshNodeTag[0];
        Models Model1("poisson2D");
        FemModule FEMmodel1(Test, Model1);
        Model1.RunFunction();
    }
    return 0;
}
