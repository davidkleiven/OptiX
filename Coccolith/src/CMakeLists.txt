set( SOURCES coccolithSim.cpp voxelMaterial.cpp fieldMonitors.cpp steadyCoccoSim.cpp )
add_library( sourcefiles STATIC ${SOURCES} )
target_link_libraries( sourcefiles ${LIBMEEP} ${LIBHARM})
