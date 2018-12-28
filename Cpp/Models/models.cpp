#include "models.h"

Models:: Models(std::string ModelName)
{
    if (ModelName.compare("poisson1D")==0)
    {
        NoOfConstants=1;
        VectorLevel=1;
        ModelFunction=&Models::poisson;
    }
    else if (ModelName.compare("poisson2D")==0)
    {
        NoOfConstants=1;
        VectorLevel=3;
        ModelFunction=&Models::poisson;
    }
    else if (ModelName.compare("poisson3D")==0)
    {
        NoOfConstants=1;
        VectorLevel=3;
        ModelFunction=&Models::poisson;
    }
    constants.set_size(NoOfConstants);
}

void Models::RunModel()
{
    (this->*ModelFunction)();
}
