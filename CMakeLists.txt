cmake_minimum_required(VERSION 3.1)

project (geotools)

enable_language(C)
enable_language(CXX)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_CXX_FLAGS "-Wall -Wextra") # -D_GLIBCXX_PARALLEL")
set (CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set (CMAKE_CXX_FLAGS_RELEASE "-O3")
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set (CMAKE_INSTALL_RPATH "/usr/local/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(APPLE)
	# This does not use GCC.
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++ -Wno-unknown-pragmas -Wall -D_GLIBCXX_PARALLEL")
else()
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -fopenmp -Wno-unknown-pragmas -Wall -D_GLIBCXX_PARALLEL")
endif()

set (CMAKE_AUTOUIC ON)
set (CMAKE_AUTOMOC ON)
set (CMAKE_INCLUDE_CURRENT_DIR ON)
set (CMAKE_AUTOMOC_OPTIONS "-Iinclude -Isrc")

set (THREADS_PREFER_PTHREAD_FLAG ON)
set (CMAKE_PREFIX_PATH "/usr/local/share/cmake/Modules;${CMAKE_PREFIX_PATH}")
set (GDAL_DIR /usr/local)

find_package (GDAL 2 REQUIRED)
find_package (Qt5 COMPONENTS Core Widgets Gui REQUIRED)
find_package (Eigen3 REQUIRED)
find_package (PythonLibs 2.7)
find_package (PDAL REQUIRED)
find_package (PCL REQUIRED)
find_package (CGAL REQUIRED)

find_library (GEOS_LIBRARY geos_c)

add_subdirectory (libgeo)

# Configure directories###################################################################################

include_directories(BEFORE SYSTEM ./libgeo/include ./libgeo/ann/include ./include /usr/local/include /usr/include 
	${GDAL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${PYTHON_INCLUDE_DIRS}
	${Qt5_DIR}/include ${Qt5Core_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS} ${Qt5Gui_INCLUDE_DIRS}
	${PDAL_INCLUDE_DIRS} ${PCL_INCLUDE_DIRS} ${CGAL_INCLUDE_DIRS}) 

link_directories(BEFORE SYSTEM ${PDAL_LIBRARY_DIRS})

# Build support libraries ################################################################################

#add_library (geotools_util SHARED src/geo_util.cpp src/reader.cpp src/writer.cpp src/stats.cpp src/raster.cpp src/plotter.cpp)
#target_link_libraries(geotools_util ${GDAL_LIBRARY} ${PYTHON_LIBRARIES} geoann)

add_library (geotools_plot SHARED src/plotter.cpp)
target_link_libraries(geotools_plot ${PYTHON_LIBRARIES})

add_library (geotools_stats SHARED src/stats.cpp)
target_link_libraries(geotools_stats)

add_library (geotools_reader SHARED src/reader.cpp)
target_link_libraries(geotools_reader geoann ${GDAL_LIBRARY})

add_library (geotools_writer SHARED src/writer.cpp)
target_link_libraries(geotools_writer geotools_stats ${GDAL_LIBRARY})

#add_executable (contrem src/contrem_app.cpp src/contrem.cpp src/ui/contrem_ui.cpp)
#target_include_directories(contrem PUBLIC contrem_autogen/include)
#target_link_libraries (contrem ${GEOS_LIBRARY} Qt5::Widgets geoutil geoann geotools_plot geotools_reader geotools_writer)

add_executable (convolve src/convolve_app.cpp src/convolve.cpp src/ui/convolve_ui.cpp)
target_include_directories(convolve PUBLIC convolve_autogen/include)
target_link_libraries (convolve Qt5::Widgets geoutil geotools_reader geotools_writer)

add_executable (reflectance src/reflectance.cpp src/ui/reflectance_ui.cpp)
target_link_libraries (reflectance Qt5::Widgets geoutil geogrid geotools_reader ${GEOS_LIBRARY})

add_executable (edgereplace src/edgereplace.cpp)

add_executable(boresight src/boresight.cpp)

add_executable (pc2grid src/pc2grid.cpp)
target_link_libraries (pc2grid geoutil geogrid geopc ${GEOS_LIBRARY})

add_executable (pcnorm src/pcnorm.cpp)
target_link_libraries (pcnorm ${PDAL_LIBRARIES})

add_executable (pcmatch src/pcmatch.cpp)
target_link_libraries (pcmatch geoutil ${PDAL_LIBRARIES})

add_executable (pcmatch3 src/pcmatch3.cpp)
target_link_libraries (pcmatch3 geoutil ${PDAL_LIBRARIES})

add_executable (sbet2csv src/sbet2csv.cpp)
target_link_libraries (sbet2csv geoutil ${PROJ_LIBRARY} ${GEOS_LIBRARY})

add_executable (splinesmooth src/splinesmooth.cpp)
target_link_libraries (splinesmooth geoutil geogrid fitpack_mod ${GDAL_LIBRARIES} ${GEOS_LIBRARY})

add_executable (idw src/idw.cpp)
target_link_libraries (idw geoutil geogrid ${GDAL_LIBRARIES} ${GEOS_LIBRARY})

add_executable (rastermatch src/rastermatch.cpp)
target_link_libraries (rastermatch geoutil geogrid ${GDAL_LIBRARIES} ${GEOS_LIBRARY})

add_executable (rastermerge src/rastermerge.cpp)
target_link_libraries (rastermerge geoutil geogrid ${GDAL_LIBRARIES} ${GEOS_LIBRARY} ${PCL_LIBRARIES})

#add_executable (polygonize src/polygonize.cpp)
#target_link_libraries (polygonize geoutil geogrid ${Boost_LIBRARIES})


# Install ################################################################################################

install(TARGETS pc2grid rastermerge pcnorm pcmatch3 convolve reflectance splinesmooth geotools_plot geotools_reader geotools_writer geotools_stats RUNTIME DESTINATION bin LIBRARY DESTINATION lib) 
