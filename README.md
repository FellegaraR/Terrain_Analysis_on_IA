# Terrain Analysis on the IA data structure #

A 2D implementation of the IA data structure for the representation 
and analysis triangle meshes. We have developed a tool based on the LibTri for
terrain analysis, which includes state-of-the-art estimators for slope
and curvature, and for the extraction of critical points. 

### Reference Paper ###

Riccardo Fellegara, Federico Iuricich, and Leila De Floriani. 
*Efficient representation and analysis of triangulated terrains*.
In Proceedings of SIGSPATIAL’17, Los Angeles Area, CA, USA, November 7–10, 2017, 4 pages.
[doi](https://doi.org/10.1145/3139958.3140050)

### Features ###

+ Topological relations extraction
    * single relation
    * batched relations extraction
+ Terrain Features
    * Triangle/Edges slope computation
    * Critical Points extraction
+ Curvature computation ([reference1](http://dl.acm.org/citation.cfm?id=1463498)and [reference2](http://www.umiacs.umd.edu/~deflo/papers/2010grapp/2010grapp.pdf))
    * Concentrated curvature
    * Mean Curvature
    * Mean and Gaussian CCurvature 
+ Dual quad mesh computation

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
This command generates a portable library file, located into *lib* folder, as well as some executables in `dist` folder.

The compilation has been test on linux and mac systems.

### Supported Input Files ###

The library supports two input formats for the triangulated irregular network:
+ off
+ tri

For a detailed description of the input formats refer the [wiki](https://github.com/FellegaraR/Terrain_Trees/wiki/Supported-Input-Formats) page.

### Use the main library ###

For a complete list of the command line options refer the [wiki](https://github.com/FellegaraR/Terrain_Analysis_on_IA/wiki/Command-line-parameters) page.
