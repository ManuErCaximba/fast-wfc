cmake_minimum_required(VERSION 3.13)

# Project name and C++ standard
project(MyProject VERSION 1.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Debug and Release compiling options
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Headers paths
include_directories(src/external)
include_directories(src/include)

# Pugixml library
add_library(pugixml STATIC src/external/pugixml/pugixml.cpp)
target_include_directories(pugixml PUBLIC src/external/pugixml)

# Add all .cpp files
file(GLOB_RECURSE LIB_SOURCES src/lib/*.cpp src/external/pugixml/pugixml.cpp)

# Make executable on root folder
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

# Executable
add_executable(wfc_demo main.cpp ${LIB_SOURCES})

# Link .exe with pugixml
target_link_libraries(wfc_demo PRIVATE pugixml)

target_compile_options(wfc_demo PRIVATE $<$<CONFIG:DEBUG>:-g>)
target_compile_options(wfc_demo PRIVATE $<$<CONFIG:RELEASE>:-O3>)