#ifndef MODELS_H
#define MODELS_H
#include <armadillo>
#include <string>
using namespace arma;

class Models
{
public:
    ///Initializes the Mathmatical Model.
    ///Pass "poisson1D" to use Poisson 1D Model.
    ///Pass "poisson2D" to use Poisson 2D Model.
    /// Pass "poisson3D" to use Poisson 3D Model.
    Models(std::string ModelName);
    void RunModel();
    std::string ModelName;
    mat constants;
    int NoOfConstants, VectorLevel;
protected:
    void (Models::*ModelFunction)();
    void poisson();
};



#endif // MODELS_H
