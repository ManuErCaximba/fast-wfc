# Nombre del ejecutable
TARGET = main.exe

# Directorios
SRC_DIR = src\lib
BUILD_DIR = build
DEBUG_DIR = $(BUILD_DIR)\Debug
RELEASE_DIR = $(BUILD_DIR)\Release

# Archivos fuente y objetos
SOURCES = main.cpp $(wildcard $(SRC_DIR)/*.cpp)
DEBUG_OBJECTS = $(patsubst %.cpp,$(DEBUG_DIR)/%.o,$(notdir $(SOURCES)))
RELEASE_OBJECTS = $(patsubst %.cpp,$(RELEASE_DIR)/%.o,$(notdir $(SOURCES)))

# Compilador y opciones de compilación
CXX = g++
CXXFLAGS = -Wall -std=c++11

# Opciones de compilación para Debug y Release
DEBUG_FLAGS = -g -O0 -DDEBUG
RELEASE_FLAGS = -O3 -DNDEBUG

# Objetivo predeterminado
all: debug

# Objetivo para compilar en modo Debug
debug: CXXFLAGS += $(DEBUG_FLAGS)
debug: $(DEBUG_DIR)/$(TARGET)

# Objetivo para compilar en modo Release
release: CXXFLAGS += $(RELEASE_FLAGS)
release: $(RELEASE_DIR)/$(TARGET)

# Regla para compilar el ejecutable en modo Debug
$(DEBUG_DIR)/$(TARGET): $(DEBUG_OBJECTS) | $(DEBUG_DIR)
	$(CXX) $(CXXFLAGS) $(DEBUG_OBJECTS) -o $@

# Regla para compilar el ejecutable en modo Release
$(RELEASE_DIR)/$(TARGET): $(RELEASE_OBJECTS) | $(RELEASE_DIR)
	$(CXX) $(CXXFLAGS) $(RELEASE_OBJECTS) -o $@

# Regla para compilar archivos de objeto en modo Debug
$(DEBUG_DIR)/%.o: %.cpp | $(DEBUG_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(DEBUG_DIR)/%.o: $(SRC_DIR)/%.cpp | $(DEBUG_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Regla para compilar archivos de objeto en modo Release
$(RELEASE_DIR)/%.o: %.cpp | $(RELEASE_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(RELEASE_DIR)/%.o: $(SRC_DIR)/%.cpp | $(RELEASE_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Crear directorios si no existen
$(DEBUG_DIR):
	mkdir $(DEBUG_DIR)

$(RELEASE_DIR):
	@if not exist $(RELEASE_DIR) mkdir $(RELEASE_DIR)

# Limpieza de archivos generados
clean:
	@if exist $(DEBUG_DIR) rmdir /S /Q $(DEBUG_DIR)
	@if exist $(RELEASE_DIR) rmdir /S /Q $(RELEASE_DIR)

# Phony targets para evitar conflictos con archivos reales
.PHONY: all debug release clean
