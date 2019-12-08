project(armadillo)

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_INCLUDE_CURRENT_DIR ON)

#set (PreDev libGmshReader GmshApi)

#set (PreInstalled armadillo lapack openblas superlu)

add_library(${PROJECT_NAME} armadillo)
    
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)
