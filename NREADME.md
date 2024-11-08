# fast-wfc-c++11

A reimplementation of [fast-wfc](https://github.com/math-fehr/fast-wfc) (which is an implementation of [Wave Function Collapse](https://github.com/mxgmn/WaveFunctionCollapse)) with some changes, the most important is the compatibility with c++-11 compilers.

## Changes from the original fast-wfc

- Remove all std::optionals, which are a class implemented for c++17 and beyond.
- Remove Color class, being remplaced by the type 'uint32_t'
- Changes on folders and files paths
- Remove the requirement of CMake, being only necesary the make command
- Replace of rapidxml to pugixml library

## Requirements

You need a C++-11 compatible compiler

## Install the library and run examples

After clone this repository run

### Release
```
make release
./build/Release/main.exe
```

### Debug
```
make debug
./build/Debug/main.exe
```

will execute WFC on the examples defined in `samples.xml`, and will put the results in `outputs/`.

## Third-parties library

The files in `example/src/include/external/` come from:
* RapidXML [https://github.com/dwd/rapidxml](https://github.com/dwd/rapidxml)
* stb Library [https://github.com/nothings/stb](https://github.com/nothings/stb)

## Image samples

The image samples come from [https://github.com/mxgmn/WaveFunctionCollapse](https://github.com/mxgmn/WaveFunctionCollapse)

## Licence 

Copyright (c) 2018-2019 Mathieu Fehr and Nathanaëlle Courant.

MIT License, see `LICENSE` for further details.
