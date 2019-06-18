cmake_minimum_required(VERSION 3.1)

enable_language(C)
enable_language(CXX)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)

if(APPLE)
	# This does not use GCC.
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++ -Wno-unknown-pragmas -Wall -D_GLIBCXX_PARALLEL")
else()
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -fopenmp -Wno-unknown-pragmas -Wall -D_GLIBCXX_PARALLEL")
endif()

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set (CMAKE_INSTALL_RPATH "/usr/local/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set (CMAKE_CXX_COMPILER_LAUNCHER ccache)

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

find_library (GEOS_LIBRARY geos)

project (contrem)

# Configure directories###################################################################################

include_directories(BEFORE SYSTEM ann/include ./include /usr/local/include /usr/include 
	${GDAL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${PYTHON_INCLUDE_DIRS}
	${Qt5_DIR}/include ${Qt5Core_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS} ${Qt5Gui_INCLUDE_DIRS}) 
link_directories(. /usr/local/lib /usr/lib)

# Build support libraries ################################################################################

set (geoannsrc ann/src/ANN.cpp ann/src/brute.cpp 
	ann/src/kd_tree.cpp ann/src/kd_util.cpp 
	ann/src/kd_split.cpp ann/src/kd_dump.cpp ann/src/kd_search.cpp ann/src/kd_pr_search.cpp 
	ann/src/kd_fix_rad_search.cpp ann/src/bd_tree.cpp ann/src/bd_search.cpp ann/src/bd_pr_search.cpp 
	ann/src/bd_fix_rad_search.cpp ann/src/perf.cpp)
add_library (geoannobj OBJECT ${geoannsrc})
set_target_properties (geoannobj PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_definitions (geoannobj PRIVATE DLL_EXPORTS=1)
add_library (geoann SHARED $<TARGET_OBJECTS:geoannobj>)

add_library (geotools_util src/util.cpp src/reader.cpp src/writer.cpp src/stats.cpp src/raster.cpp src/plotter.cpp)
target_link_libraries(geotools_util ${GDAL_LIBRARY} ${PYTHON_LIBRARIES} geoann)

add_executable (contrem src/contrem_app.cpp src/contrem.cpp src/ui/contrem_ui.cpp)
target_include_directories(contrem PUBLIC contrem_autogen/include)
target_link_libraries (contrem ${GEOS_LIBRARY} Qt5::Widgets geotools_util geoann)

add_executable (convolve src/convolve_app.cpp src/convolve.cpp src/ui/convolve_ui.cpp)
target_include_directories(convolve PUBLIC convolve_autogen/include)
target_link_libraries (convolve Qt5::Widgets geotools_util)

add_executable (reflectance src/reflectance.cpp src/ui/reflectance_ui.cpp)
target_link_libraries (reflectance Qt5::Widgets geotools_util)

add_executable (edgereplace src/edgereplace.cpp)

add_executable(boresight src/boresight.cpp)


# Install ################################################################################################

install(TARGETS convolve reflectance RUNTIME DESTINATION bin LIBRARY DESTINATION lib) 
