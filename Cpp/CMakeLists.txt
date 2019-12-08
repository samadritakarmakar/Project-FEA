cmake_minimum_required(VERSION 2.8)
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

project(ProjectFEA)

#set(CMAKE_EXE_LINKER_FLAGS "-s") #Strip Binary
set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
find_package(Armadillo REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(OpenMP)
#if(OPENMP_FOUND)
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

include_directories("armadillo")
add_subdirectory(armadillo/)
include_directories ("libGmshReader/GmshApi")
add_subdirectory(libGmshReader/GmshApi/)
include_directories ("libGmshReader")
add_subdirectory(libGmshReader/)
include_directories("Math")
add_subdirectory(Math/)
include_directories ("Variable")
add_subdirectory(Variable/)
include_directories ("ShapeFunction/LagrangeShapeFunction")
add_subdirectory(ShapeFunction/LagrangeShapeFunction/)
include_directories ("ShapeFunction")
add_subdirectory(ShapeFunction/)
include_directories ("FEMtools")
add_subdirectory(FEMtools/)
include_directories ("ProjectJacobian")
add_subdirectory(ProjectJacobian/)
include_directories ("TestAndTrialFunction")
add_subdirectory(TestAndTrialFunction/)
include_directories ("Integrator")
add_subdirectory(Integrator/)
include_directories ("Dynamic")
add_subdirectory(Dynamic/)
include_directories ("BoundaryCondtions")
add_subdirectory(BoundaryCondtions/)
include_directories("Models")
add_subdirectory(Models/)
include_directories("GmshWriter")
add_subdirectory(GmshWriter/)


#set (PreInstalled armadillo blas lapack)

set (InternalDev TestAndTrialFunction)

set (PreDev libGmshReader GmshApi LagrangeShapeFunction ShapeFunction FEMtools ProjectJacobian)
if(OPENMP_FOUND)
set (PreInstalled armadillo lapack openblas superlu OpenMP::OpenMP_CXX)
else()
    set (PreInstalled armadillo lapack openblas superlu)
endif()

add_library(${PROJECT_NAME} ProjectFEA.hpp)
target_link_libraries (${PROJECT_NAME} ${PreDev} ${PreInstalled} ${InternalDev} ${CMAKE_THREAD_LIBS_INIT}
   ${CMAKE_SOURCE_DIR}/libGmshReader/GmshApi/libgmsh.a)

add_executable(ElasticityEx ElasticityExample.cpp)
target_link_libraries (ElasticityEx ProjectFEA)

add_executable(Poisson Poisson.cpp)
target_link_libraries (Poisson ProjectFEA)

#add_executable(OneElmntVrfy OneElmntVrfy.cpp ProjectFEA.hpp)
#target_link_libraries (OneElmntVrfy ${PreDev} ${PreInstalled} ${InternalDev} ${CMAKE_THREAD_LIBS_INIT}
#   ${CMAKE_SOURCE_DIR}/libGmshReader/GmshApi/libgmsh.a)

add_executable(OneElmntVrfy OneElmntVrfy.cpp)
target_link_libraries (OneElmntVrfy ProjectFEA)

add_executable(DiffusionEx DiffusionEx.cpp)
target_link_libraries (DiffusionEx ProjectFEA)

add_executable(VolumeEx Volume.cpp)
target_link_libraries (VolumeEx ProjectFEA)

add_executable(TestDynamic DynamicTest.cpp)
target_link_libraries(TestDynamic ProjectFEA)

add_executable(DynamicHeatDiff DynamicHeatDiff.cpp)
target_link_libraries(DynamicHeatDiff ProjectFEA)

add_executable(ElasticVibration ElasticVibration.cpp)
target_link_libraries(ElasticVibration ProjectFEA)

#add_executable(Trial Trial.cpp)
#target_link_libraries (Trial ProjectFEA)

set (COPY_FILES BaseGeom GaussData OneElmntGeom SampleGmsh)
foreach(loop_var ${COPY_FILES})
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${loop_var} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
endforeach(loop_var)
