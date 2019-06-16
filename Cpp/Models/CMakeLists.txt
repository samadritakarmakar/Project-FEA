project(Models)

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories("${PROJECT_SOURCE_DIR}/Solids/")
add_subdirectory(Solids/)

add_library(${PROJECT_NAME} Models.hpp)
    
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)
