#include "models.h"
void Models::InitializePoission(std::string ModelName)
{
    coupleLevel=0;
    NoOfConstants=1;
    VectorLevel.set_size(coupleLevel+1);
    ModelFunction=&Models::poisson;
    if (ModelName.compare("poisson1D")==0)
    {
        VectorLevel(0)=1;
    }
    else if (ModelName.compare("poisson2D")==0)
    {
        VectorLevel(0)=2;
    }
    else if (ModelName.compare("poisson3D")==0)
    {
            VectorLevel(0)=3;
    }
    else
    {
        std::cout<<"Does not match any of the available Poisson models.\n";
        std::cout<<"Please check your Model Input or contact your developer.\n";
        throw;
    }
    std::cout<<"Model Poission Initialized"<<"\n";
}

void Models::poisson()
{
    std::cout<<"Possion model Working!\n";
}
