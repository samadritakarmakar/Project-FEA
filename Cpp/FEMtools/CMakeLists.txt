cmake_minimum_required(VERSION 2.8.11)

project(FEMtools CXX)

set(CMAKE_BUILD_TYPE Release)
## Compiler flags
#Comment out if building on UNIX
#    set(CMAKE_CXX_FLAGS "-O2 -mwindows")        ## Optimize
#Uncomment if building on UNIX
    set(CMAKE_CXX_FLAGS "-O3")
    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
#set(CMAKE_AUTOMOC ON)

#set(CMAKE_AUTOUIC ON)

#set(QT_USE_QTWEBKIT TRUE)
# Find the QtWidgets library
#find_package(Qt5Widgets REQUIRED)

#Find Armadillo 
set (LoadGaussFile LoadGaussFile.h LoadGaussFile.cpp GetNodePostions.cpp)

# Tell CMake to create the helloworld executable
#add_executable(K_assmbly k_spring.cpp k_truss.cpp k_assmbly_file.cpp)
add_library(${PROJECT_NAME} FEMtools.h ${LoadGaussFile})#BuildCosMatrix.cpp VectorizeGradu.cpp VectorizeX.cpp BuildLocalMatrices.cpp  GetNodePostions.cpp TestFunctionDot.cpp  VectorizePhi.cpp)

