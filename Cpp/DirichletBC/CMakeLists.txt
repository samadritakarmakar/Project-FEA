project(DirichletBC)

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_INCLUDE_CURRENT_DIR ON)
#set (InternalDev GmshNodeListGen.h cof.hpp
#    BuildTriangleMatrix.hpp FindCoeffTriangle.hpp NumOfTriNodes.hpp
#    FindCoeffQuad.hpp NumOfQuadNodes.hpp BuildQuadMatrix.hpp
#    BuildHexahedralMatrix.hpp FindCoeffHexahedral.hpp NumOfHexaNodes.hpp
#    BuildTetrahedralMatrix.hpp FindCoeffTertrahedral.hpp NumOfTetraNodes.hpp
#    BuildLineMatrix.hpp FindCoeffLine.hpp NumOfLineNodes.hpp
#    CheckSupportedOrder.hpp)

#set (PreDev libGmshReader GmshApi)

#set (PreInstalled armadillo lapack openblas superlu)

add_library(${PROJECT_NAME} DirichletBC.hpp Expression.hpp)
    
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)
