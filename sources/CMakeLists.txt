cmake_minimum_required(VERSION 2.8)

include_directories(core_library)
add_subdirectory(core_library)

set(INCLUDE_DIRS		
		curvature
		terrain_features
	)
	
INCLUDE_DIRECTORIES( ${INCLUDE_DIRS} )

file(GLOB CURVATURE_SRC
        curvature/*.cpp
        curvature/*.h
	)
		
file(GLOB TERRAIN_FEATURES_SRC
        terrain_features/*.h
        terrain_features/*.cpp
    ) 

set(ALL_SRC
        ${CURVATURE_SRC}
        ${TERRAIN_FEATURES_SRC}
	)
	
add_library(_terrain_analysis ${ALL_SRC})
target_link_libraries(_terrain_analysis _IA2D)

