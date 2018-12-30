#include "models.h"

Models:: Models(std::string ModelName)
{
    if (ModelName.compare("poisson1D")==0)
    {
        coupleLevel=0;
        NoOfConstants=1;
        VectorLevel.set_size(coupleLevel+1);
        VectorLevel(0)=1;
        ModelFunction=&Models::poisson;
    }
    else if (ModelName.compare("poisson2D")==0)
    {
        coupleLevel=0;
        NoOfConstants=1;
        VectorLevel.set_size(coupleLevel+1);
        VectorLevel(0)=2;
        ModelFunction=&Models::poisson;
    }
    else if (ModelName.compare("poisson3D")==0)
    {
        coupleLevel=0;
        NoOfConstants=1;
        VectorLevel.set_size(coupleLevel+1);
        VectorLevel(0)=3;
        ModelFunction=&Models::poisson;
    }
    constants.set_size(NoOfConstants);
}

void Models::RunFunction()
{
    (this->*ModelFunction)();
}
