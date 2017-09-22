cmake_minimum_required(VERSION 2.8)

set(INCLUDE_DIRS
		ia
		utilities
	)
	
INCLUDE_DIRECTORIES( ${INCLUDE_DIRS} )

file(GLOB IA_SRC
        ia/*.h
   )

file(GLOB UTILITIES_SRC 
        utilities/*.cpp
        utilities/*.h
	)

set(ALL_SRC
        ${IA_SRC}
        ${UTILITIES_SRC}
	)
	
add_library(_IA2D ${ALL_SRC})


