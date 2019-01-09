#include "models.h"

Models:: Models(std::string ModelName)
{
    if (ModelName.find("poisson")!=std::string::npos)
    {
        InitializePoission(ModelName);
    }
    constants.set_size(NoOfConstants);
}

void Models::RunFunction()
{
    (this->*ModelFunction)();
}
