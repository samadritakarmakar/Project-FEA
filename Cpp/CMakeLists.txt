cmake_minimum_required(VERSION 2.8.11)

project(prj)

set(CMAKE_BUILD_TYPE Release)
## Compiler flags
#Comment out if building on UNIX
#    set(CMAKE_CXX_FLAGS "-O2 -mwindows")        ## Optimize
#Uncomment if building on UNIX
    set(CMAKE_CXX_FLAGS "-O2")

    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

set(CMAKE_AUTOUIC ON)

#set(QT_USE_QTWEBKIT TRUE)
# Find the QtWidgets library
#find_package(Qt5Widgets REQUIRED)

# Tell CMake to create the helloworld executable
add_executable(GATE_calc main.cpp widget.cpp)

# Use the Widgets module from Qt 5.
target_link_libraries(GATE_calc Qt5::Widgets Qt5WebKitWidgets)

#file(COPY ${YOUR_SRC_IMAGE_DIR}/Calc/* DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Calc/)
add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/Calc
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/Calc ${CMAKE_BINARY_DIR}/Calc)
add_custom_target(Calc ALL
    DEPENDS ${CMAKE_BINARY_DIR}/Calc)
