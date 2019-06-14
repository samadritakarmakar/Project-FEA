cmake_minimum_required(VERSION 2.8)

project(ProjectFEA)

#set(CMAKE_EXE_LINKER_FLAGS "-s") #Strip Binary
set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
find_package(Armadillo REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
if(OPENMP_FOUND)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories ("${PROJECT_SOURCE_DIR}/libGmshReader/GmshApi/")
add_subdirectory(libGmshReader/GmshApi/)
include_directories ("${PROJECT_SOURCE_DIR}/libGmshReader/")
add_subdirectory(libGmshReader/)
include_directories("${PROJECT_SOURCE_DIR}/Math/")
add_subdirectory(Math/)
include_directories ("${PROJECT_SOURCE_DIR}/ShapeFunction/LagrangeShapeFunction/")
add_subdirectory(ShapeFunction/LagrangeShapeFunction/)
include_directories ("${PROJECT_SOURCE_DIR}/ShapeFunction/")
add_subdirectory(ShapeFunction/)
include_directories ("${PROJECT_SOURCE_DIR}/FEMtools/")
add_subdirectory(FEMtools/)
include_directories ("${PROJECT_SOURCE_DIR}/ProjectJacobian/")
add_subdirectory(ProjectJacobian/)
include_directories ("${PROJECT_SOURCE_DIR}/TestAndTrialFunction/")
add_subdirectory(TestAndTrialFunction/)
include_directories ("${PROJECT_SOURCE_DIR}/Variable")
add_subdirectory(Variable/)
include_directories ("${PROJECT_SOURCE_DIR}/Integrator/")
add_subdirectory(Integrator/)
include_directories ("${PROJECT_SOURCE_DIR}/DirichletBC/")
add_subdirectory(DirichletBC/)
include_directories("${PROJECT_SOURCE_DIR}/Models")
add_subdirectory(Models/)
include_directories("${PROJECT_SOURCE_DIR}/GmshWriter")
add_subdirectory(GmshWriter/)


#set (PreInstalled armadillo blas lapack)

set (InternalDev TestAndTrialFunction)

set (PreDev libGmshReader GmshApi LagrangeShapeFunction ShapeFunction FEMtools ProjectJacobian)

set (PreInstalled armadillo lapack openblas superlu)

add_executable(${PROJECT_NAME} ProjectFEA.hpp main.cpp)
target_link_libraries (${PROJECT_NAME} ${PreDev} ${PreInstalled} ${InternalDev} ${CMAKE_THREAD_LIBS_INIT}
    ${CMAKE_SOURCE_DIR}/libGmshReader/GmshApi/libgmsh.a)

add_executable(ElasticityEx ProjectFEA.hpp ElasticityExample.cpp)
target_link_libraries (ElasticityEx ${PreDev} ${PreInstalled} ${InternalDev} ${CMAKE_THREAD_LIBS_INIT}
    ${CMAKE_SOURCE_DIR}/libGmshReader/GmshApi/libgmsh.a)

add_executable(OneElmntVrfy ProjectFEA.hpp OneElmntVrfy.cpp)
target_link_libraries (OneElmntVrfy ${PreDev} ${PreInstalled} ${InternalDev} ${CMAKE_THREAD_LIBS_INIT}
    ${CMAKE_SOURCE_DIR}/libGmshReader/GmshApi/libgmsh.a)

set (COPY_FILES BaseGeom GaussData OneElmntGeom SampleGmsh)
foreach(loop_var ${COPY_FILES})
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${loop_var} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
endforeach(loop_var)
