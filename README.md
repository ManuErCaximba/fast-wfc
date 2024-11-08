# fast-wfc-cpp11

A re-implementation of [fast-wfc](https://github.com/math-fehr/fast-wfc) to make the algorithm run on C++11 compilers and some minor changes.

# Requirements

You need a C++-11 compatible compiler, and CMake installed.

# Main changes from fast-wfc

* I changed std::optional (C++17 or newer) and std::nullopt to pointers and nullptr
* I removed color.hpp and all their references and changed them to use uint32_t instead
* I changed the project original paths and CMakeFiles.txt to make it work on Windows and VSCode specifically (It should still works in Unix system as well due to I'm using g++/gdb)
* I added a CMakePresets.json to specify Debug and Release builds
* I changed RapidXML to pugiXML library to read the XML
* I created a single header with all the algorithm code for an easy implementation in other projects

# Run the examples

## Debug

```
cmake --preset debug
cmake --build --preset build-debug
./wfc_demo
```

## Release

```
cmake --preset release
cmake --build --preset build-release
./wfc_demo
```

It will execute WFC on the examples defined in `samples.xml`, and will put the results in `results/`.

# Third-parties library

The files in `src/external/` comes from:
* pugiXML [https://github.com/zeux/pugixml](https://github.com/zeux/pugixml)
* stb Library [https://github.com/nothings/stb](https://github.com/nothings/stb)

# Image samples

The image samples come from [https://github.com/mxgmn/WaveFunctionCollapse](https://github.com/mxgmn/WaveFunctionCollapse)

# Licence 

Copyright (c) 2018-2019 Mathieu Fehr and Nathanaëlle Courant.

MIT License, see `LICENSE` for further details.

