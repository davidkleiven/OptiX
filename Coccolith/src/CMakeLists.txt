set( SOURCES coccolithSim.cpp voxelMaterial.cpp fieldMonitors.cpp steadyCoccoSim.cpp
refractiveIndexMaterial.cpp stokesParameters.cpp cylindricalBirefringentMaterial.cpp
voxelMaterial2D.cpp bandStructureSimulation2D.cpp)
add_library( sourcefiles STATIC ${SOURCES} )
target_link_libraries( sourcefiles ${LIBMEEP} ${LIBHARM})
