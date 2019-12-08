#ifndef PROJECTFEA_HPP
#define PROJECTFEA_HPP
//#define ARMA_USE_SUPERLU 1
#include "armadillo"
#include "LagrangeShapeFunctionAllElementTypes.hpp"
#include "TrialFunction.hpp"
#include "TestFunctionGalerkin.hpp"
#include "Form.hpp"
#include "FormMultiThread.hpp"
#include "TrialFunctionNeumannSurface.hpp"
#include "TrialFunctionNeumannLine.hpp"
#include "LocalIntegration.hpp"
#include "SystemAssembly.hpp"
#include "Dynamic/Dynamic.hpp"
#include "InitialBC.hpp"
//#include "DirichletBC.hpp"
#include "Expression.hpp"
#include "Variable.hpp"
#include "GmshWriter.hpp"
#include "vecnorm.hpp"
#include "ProjectFEA_Math.hpp"
#include "GetLengthAlongVector.hpp"
//#include "PreDefinedElasticityTensor.hpp"
#include "Models.hpp"


#endif // PROJECTFEA_HPP
