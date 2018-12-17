cmake_minimum_required(VERSION 2.8.11)

project(Project-FEA2)

#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_BUILD_TYPE Debug)
## Compiler flags
#Comment out if building on UNIX
#    set(CMAKE_CXX_FLAGS "-O3 -mwindows")        ## Optimize
#Uncomment if building on UNIX
    set(CMAKE_CXX_FLAGS "-O3")

   # set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

set(CMAKE_AUTOUIC ON)

#set Gmsh cmake Module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

#set(QT_USE_QTWEBKIT TRUE)
# Find the QtWidgets library
#find_package(Qt5Widgets REQUIRED)

find_package(Armadillo REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
#find_package(SLU REQUIRED)
#find_package(Gmsh REQUIRED)
find_package(SuperLU)
find_package(OpenMP)

#include_directories ("${PROJECT_SOURCE_DIR}/libMeshReader")
#add_subdirectory(libMeshReader)
include_directories ("${PROJECT_SOURCE_DIR}/libGmshReader/GmshApi/")
add_subdirectory(libGmshReader/GmshApi/)
include_directories ("${PROJECT_SOURCE_DIR}/libGmshReader/")
add_subdirectory(libGmshReader/)
#include_directories ("${PROJECT_SOURCE_DIR}/FEMtools/")
#add_subdirectory(FEMtools/)
include_directories ("${PROJECT_SOURCE_DIR}/FEM/")
add_subdirectory(FEM/)

add_executable(Prj-FEA2 mainTest.cpp)

#target_link_libraries (Prj-FEA2 libMeshReader armadillo lapack openblas gmsh superlu)
target_link_libraries (Prj-FEA2 libGmshReader armadillo lapack openblas superlu ${CMAKE_SOURCE_DIR}/libGmshReader/GmshApi/libgmsh.a GmshApi)
target_link_libraries(Prj-FEA2 PRIVATE "${OpenMP_CXX_FLAGS}")
target_compile_options(Prj-FEA2 PRIVATE "${OpenMP_CXX_FLAGS}")
# Use the Widgets module from Qt 5.
#target_link_libraries(GATE_calc Qt5::Widgets Qt5WebKitWidgets)

#file(COPY ${YOUR_SRC_IMAGE_DIR}/Calc/* DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Calc/)
#add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/Calc
#     COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/Calc ${CMAKE_BINARY_DIR}/Calc)
# add_custom_target(Calc ALL
#     DEPENDS ${CMAKE_BINARY_DIR}/Calc)
