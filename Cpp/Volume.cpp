#include "ProjectFEA.hpp"
#include <math.h>

class FindSize : public LocalIntegrator<TrialFunction>
{
public:
    FindSize(FormMultiThread<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v):
        LocalIntegrator<TrialFunction> (a, u, v){}
    mat scalar_integration(FormMultiThread<TrialFunction>& a, TrialFunction& u, int thread)
    {
        mat dx;
        dx<<a[thread].dX(u)<<endr;
        return dx;
    }
};
int main(int argc, char *argv[])
{
    if(argc==1 || argc>3)
    {
        std::cout<<"Usage: ./SizeEx <.msh Filename> <Dimension>\n";
        return 0;
    }
    std::string FileName(argv[1]);
    int Dimension=*argv[2]-'0';
    int vectorLevel=1;
    libGmshReader::MeshReader Mesh(FileName, Dimension);
    TrialFunction u(Mesh, vectorLevel);
    TestFunctionGalerkin<TrialFunction> v(u);
    FormMultiThread<TrialFunction> a;
    mat h;
    SystemAssembler<FindSize,TrialFunction> SclrIntrgtn(a, u ,v);
    FindSize SizeOfElement(a,u,v);
    SclrIntrgtn.RunScalarIntegration(SizeOfElement, h);
    h=abs(h);
    cout<<"Total Volume = "<<sum(h)<<"\n";
    GmshWriter WriteStrss(u, "Volume.pos");
    WriteStrss.viewName="Volume";
    WriteStrss.SetDataType_to_ElementData();
    WriteStrss.WriteToGmsh(h);
}
