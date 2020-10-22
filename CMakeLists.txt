cmake_minimum_required(VERSION 3.1)

option (FULL_BUILD "Build all libraries and programs." ON)
option (WITH_PC "Build the point cloud computers for pc2grid." ON)
option (WITH_ANN "Build the ANN KD-Tree implementation." OFF)
option (WITH_GUI "Build the GUIs." ON)

project (geotools)

enable_language(C)
enable_language(CXX)
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set (CMAKE_INSTALL_RPATH "/usr/local/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(${WIN32})
	message (STATUS "Building for Windows")
	set (CMAKE_CXX_FLAGS "/Zc:alignedNew /Gz /EHsc /Zc:__cplusplus /std:c++17" CACHE STRING "compile flags" FORCE)	
	set (CMAKE_CXX_FLAGS_DEBUG "/Od /Zi /GS- /guard:cf-" CACHE STRING "compile flags" FORCE)	
	add_definitions (-DBOOST_ALL_NO_LIB)
	set (BOOST_ROOT "C:\\dev\\boost_1_67_0")
	set (DEV_DIR "C:\\dev")
	set (GIT_DIR "${DEV_DIR}\\git")
	set (GIS_DIR "${DEV_DIR}\\OSGeo4W64")
	set (GEOS_LIBRARYDIR "${GIS_DIR}\\lib;${GIS_DIR}\\bin") #set (GEOS_LIBRARYDIR "${DEV_DIR}\\geos\\lib;${DEV_DIR}\\geos\\bin")
	set (GEOS_INCLUDEDIR "${GIS_DIR}\\include") #set (GEOS_INCLUDEDIR "${DEV_DIR}\\geos\\include")
	set (EIGEN3_INCLUDE_DIR "${GIS_DIR}\\include\\eigen3")
	set (LIBLAS_INCLUDE_DIR "${GIS_DIR}\\include")
	set (GDAL_INCLUDE_DIR "${GIS_DIR}\\include")
	set (GDAL_LIBRARYDIR "${GIS_DIR}\\lib;${GIS_DIR}\\bin")
	set (GIS_INCLUDE_DIR "${GIS_DIR}\\include")
	set (GIS_LIBRARYDIR "${GIS_DIR}\\lib;${GIS_DIR}\\bin")
	set (CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${GIS_DIR}\\bin;${GIS_DIR}\\lib;${DEV_DIR}\\Qt\\5.14.2\\msvc2017_64\\lib\\cmake")
	set (SQLITE_LIBRARY sqlite3_i)
	set (SPATIALITE_LIBRARY spatialite_i)
	set (GDAL_LIBRARY gdal_i)
	set (LIBLAS_LIBRARY liblas)
	set (LIBLAS_INCLUDE_DIRS "${GIS_DIR}\\include")
	set (LIBLAS_LIBDIR "${GIS_DIR}\\lib")
	find_package (Boost 1.67)
else()
	if (WITH_GUI)
		set (CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}:/usr/lib/x86_64-linux-gnu/cmake/Qt5")
	endif (WITH_GUI)
	set (CMAKE_CXX_COMPILER_LAUNCHER ccache)
	set (CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unknown-pragmas -fopenmp -frounding-math -D_GLIBCXX_PARALLEL" CACHE STRING "compile flags" FORCE)
	set (CMAKE_CXX_FLAGS_DEBUG "-g -O0" CACHE STRING "compile flags" FORCE)
	set (CMAKE_CXX_FLAGS_RELEASE "-O3" CACHE STRING "compile flags" FORCE)
	set (FS_LIB stdc++fs) # For filesystem features in c++17.
	find_package (Eigen3 REQUIRED)
	find_package (PDAL REQUIRED)
	find_package (PythonLibs 3)
	find_library (LIBLAS_LIBRARY las)
endif()


set (CMAKE_AUTOUIC ON)
set (CMAKE_AUTOMOC ON)
set (CMAKE_INCLUDE_CURRENT_DIR ON)
set (CMAKE_AUTOMOC_OPTIONS "-Iinclude -Isrc")

set (THREADS_PREFER_PTHREAD_FLAG ON)

find_package (GDAL 2 REQUIRED)
find_package (Qt5 COMPONENTS Core Widgets Gui REQUIRED)
find_package (Threads REQUIRED)
find_package (PDAL REQUIRED)

find_library (GEOS_LIBRARY geos_c)


add_subdirectory (libgeo)

# Configure directories###################################################################################

include_directories(
	./include 
	./libgeo/include 
	./libgeo/ann/include 
	${GDAL_INCLUDE_DIR} 
	${EIGEN3_INCLUDE_DIR} 
	${PYTHON_INCLUDE_DIRS}
	${LIBLAS_INCLUDE_DIR}
	${PDAL_INCLUDE_DIRS}
	${Boost_INCLUDE_DIRS}
	${Qt5_DIR}/include ${Qt5Core_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS} ${Qt5Gui_INCLUDE_DIRS}
) 

link_directories(
	bin/Release
	lib/Release
	${GDAL_LIBRARYDIR}
	${LIBLAS_LIBDIR}
	${PDAL_LIBRARY_DIRS}
	${GMP_LIBDIR}
	${Boost_LIBRARY_DIRS}
)

# Build support libraries ################################################################################

if (FULL_BUILD)
	message (STATUS "Full Build...")
	
	add_library (geotools_plot SHARED src/plotter.cpp)
	target_link_libraries(geotools_plot ${PYTHON_LIBRARIES})

	add_library (geotools_stats SHARED src/stats.cpp)
	target_link_libraries(geotools_stats)

	add_library (geotools_reader SHARED src/reader.cpp)
	target_link_libraries(geotools_reader geoann ${GDAL_LIBRARY})

	add_library (geotools_writer SHARED src/writer.cpp)
	target_link_libraries(geotools_writer geotools_stats ${GDAL_LIBRARY})

	add_executable (contrem src/contrem_app.cpp src/contrem.cpp src/ui/contrem_ui.cpp)
	target_include_directories(contrem PUBLIC contrem_autogen/include)
	target_link_libraries (contrem ${GEOS_LIBRARY} Qt5::Widgets geoutil geoann geotools_plot geotools_reader geotools_writer)

	add_executable (convolve src/convolve_app.cpp src/convolve.cpp src/ui/convolve_ui.cpp)
	target_include_directories(convolve PUBLIC convolve_autogen/include)
	target_link_libraries (convolve geoutil geotools_reader geotools_writer Threads::Threads Qt5::Widgets)

	add_executable (reflectance src/reflectance.cpp src/ui/reflectance_ui.cpp)
	target_link_libraries (reflectance Qt5::Widgets geoutil geogrid geotools_reader ${GEOS_LIBRARY})

	add_executable (edgereplace src/edgereplace.cpp)

	add_executable (sbet2csv src/sbet2csv.cpp)
	target_link_libraries (sbet2csv geoutil ${PROJ_LIBRARY} ${GEOS_LIBRARY})

	add_executable (splinesmooth src/splinesmooth.cpp)
	target_link_libraries (splinesmooth geoutil geogrid fitpack_mod ${GDAL_LIBRARIES} ${GEOS_LIBRARY})

	add_executable (idw src/idw.cpp)
	target_link_libraries (idw geoutil geogrid ${GDAL_LIBRARIES} ${GEOS_LIBRARY})

	add_executable (rastermerge src/rastermerge.cpp)
	target_link_libraries (rastermerge geoutil geogrid ${GDAL_LIBRARIES} ${GEOS_LIBRARY})

	add_executable (polygonize src/polygonize.cpp)
	target_link_libraries (polygonize geoutil geogrid ${Boost_LIBRARIES})

endif (FULL_BUILD)

add_executable (pc2grid src/pc2grid.cpp libgeo/src/geo.cpp)
target_compile_features(pc2grid PUBLIC cxx_std_17)
target_link_libraries (pc2grid geoutil geogrid geopc ${LIBLAS_LIBRARY} ${GEOS_LIBRARY} ${PDAL_LIBRARIES})

add_executable (pcnorm src/pcnorm.cpp)
target_compile_features(pcnorm PUBLIC cxx_std_17)
target_link_libraries (pcnorm geoutil ${LIBLAS_LIBRARY})

add_executable (voidfill src/voidfill.cpp libgeo/src/geo.cpp)
target_compile_features(voidfill PUBLIC cxx_std_17)
target_link_libraries (voidfill geoutil geogrid ${GDAL_LIBRARY})

add_executable (map2raw src/map2raw.cpp)
target_compile_features(map2raw PUBLIC cxx_std_17)
target_link_libraries (map2raw ${GDAL_LIBRARIES})

add_executable (srtmlakes src/srtmlakes.cpp)
target_compile_features(srtmlakes PUBLIC cxx_std_17)
target_link_libraries (srtmlakes geoutil geogrid)

if (FULL_BUILD)

	install(TARGETS voidfill pc2grid rastermerge pcnorm convolve geotools_plot geotools_reader geotools_writer geotools_stats RUNTIME DESTINATION bin LIBRARY DESTINATION lib) 

else ()

	install(TARGETS voidfill pc2grid pcnorm map2raw RUNTIME DESTINATION bin LIBRARY DESTINATION lib) 

endif ()