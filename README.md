# GraphHarmonicMap
implement harmonic map from arbitrary surface to its pants decomposition graph

## Prerequisites
* C++ standard: C++11, some of C++14
* Compiler: 
    * Visual Studio 2015+
    * Clang/LLVM: 3.8+ (not recommended, clang/llvm shipped with Xcode does not support openmp by default)
    * gcc: 6+
* cmake: 3.0+
* openmp: any version (highly recommended)

## Installation
* clone source code
```
git clone --recursive https://github.com/cfwen/GraphHarmonicMap.git
```
* cmake to generate project file or makefile
```
cd GraphHarmonicMap
mkdir build
cd build
cmake ..
```
You can also generator for certain toolchains, for example on macOS:
```
cmake .. -G XCode
```
or generate with x64 toolchains (on Windows default is x86, which is not recommended by me)
```
cmake .. -G "Visual Studio 14 2015 Win64"
```
* if cmake goes well, then make 
```
make
```
or open generated solution file with IDE (Visual Studio, XCode, etc..) and build the solution

if everything ok, then
```
make install
```
or build "install" project in IDE, binary will be installed in ~/bin, where ~ is the root of GraphHarmonicMap source 

# Use
* first compute harmonic map to graph
```
./bin/GraphHarmonic harmonic data/eight.m data/eight.graph data/eight.cut data/eight.harmonic.m
```
* then decompose surface into cylinders
```
./bin/GraphHarmonic decompose data/eight.harmonic.m data/eight.graph data/eight.cut data/eight.harmonic.2.m
```
* visualize
```
./bin/G3DOGL data/eight.harmonic.m -texturemap ./bin/check_128.jpg
```
