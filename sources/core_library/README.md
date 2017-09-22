# LibTri library #

The core of the IA data structure for triangle meshes.
This repository provides only the portable library.

### Features ###

+ Topological relations extraction
    * single relation
    * batched relations extraction

### How to compile ###

The library requires only the standard template library installed in your system.

Once in the root of the repository type from the command line
```
cmake CMakeList.txt
```
and once configured
```
make
```
This command generates a portable library file located into *lib* folder.

The compilation has been test on linux and mac systems.
