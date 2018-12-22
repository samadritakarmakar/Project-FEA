#ifndef LOADGAUSSFILE_H
#define LOADGAUSSFILE_H
#include <armadillo>
#include <string>
#include <cmath>
#include <iostream>
#include "libGmshReader.h"
namespace FEMtools
{
std::string LoadGaussFile(const libGmshReader::MeshReader &Mesh, const int & ElementType);
void TensorGauss(int &p, float &n, std::string &GaussFileName);
}


#endif // LOADGAUSSFILE_H
