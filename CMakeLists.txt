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

find_package (GDAL 2.2 REQUIRED)
find_package (Qt5Widgets REQUIRED)
find_library (GEOS_LIBRARY geos PATHS /opt/lib /usr/local/lib /usr/lib NO_CMAKE_SYSTEM_PATH)

project (contrem)

# Configure directories###################################################################################

include_directories(BEFORE SYSTEM ./include /opt/include /usr/include/gdal) 
link_directories(. /opt/lib /usr/local/lib /usr/lib)

# Build support libraries ################################################################################

add_executable (contrem src/contrem.cpp src/ui/contrem_ui.cpp src/reader.cpp src/writer.cpp src/contrem_proc.cpp src/stats.cpp)
target_include_directories(contrem PUBLIC contrem_autogen/include  
	${Qt5_DIR}/include ${Qt5Core_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS} ${Qt5Gui_INCLUDE_DIRS})
target_link_libraries (contrem gdal ${GEOS_LIBRARY} Qt5::Widgets)

add_executable (convolve src/convolve.cpp src/convolver.cpp src/ui/convolve_ui.cpp)
target_include_directories(convolve PUBLIC convolve_autogen/include  
	${Qt5_DIR}/include ${Qt5Core_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS} ${Qt5Gui_INCLUDE_DIRS})
target_link_libraries (convolve Qt5::Widgets)

add_executable (nano_timesync src/nano_timesync.cpp src/reader.cpp src/raster.cpp src/ui/nano_timesync_ui.cpp)
target_link_libraries (nano_timesync gdal Qt5::Widgets)

add_executable (edgereplace src/edgereplace.cpp)


# Install ################################################################################################

install(TARGETS contrem convolve RUNTIME DESTINATION bin LIBRARY DESTINATION lib) 
