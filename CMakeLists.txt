cmake_minimum_required(VERSION 3.13)

# Nombre del proyecto y estándar de C++
project(MyProject VERSION 1.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Agrega opciones de compilación para Debug y Release
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Incluye el directorio de cabeceras
include_directories(src/external)
include_directories(src/include)

# Biblioteca para pugixml
add_library(pugixml STATIC src/external/pugixml/pugixml.cpp)
target_include_directories(pugixml PUBLIC src/external/pugixml)

# configure_file(${CMAKE_SOURCE_DIR}/samples.xml ${CMAKE_BINARY_DIR}/samples.xml COPYONLY)

# Agrega todos los archivos cpp en src/lib
file(GLOB_RECURSE LIB_SOURCES src/lib/*.cpp src/external/pugixml/pugixml.cpp)

# Configura el directorio de salida del ejecutable en el directorio raíz
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

# Ejecutable del proyecto
add_executable(wfc_demo main.cpp ${LIB_SOURCES})

# Vincular el ejecutable con pugixml
target_link_libraries(wfc_demo PRIVATE pugixml)

# Establece propiedades de depuración
target_compile_options(wfc_demo PRIVATE $<$<CONFIG:DEBUG>:-g>)
target_compile_options(wfc_demo PRIVATE $<$<CONFIG:RELEASE>:-O3>)