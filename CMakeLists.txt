cmake_minimum_required(VERSION 3.1)

enable_language(C)
enable_language(CXX)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
if(APPLE)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -Wno-unknown-pragmas -Wall -D_GLIBCXX_PARALLEL")
else()
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -fopenmp -Wno-unknown-pragmas -Wall -D_GLIBCXX_PARALLEL")
endif()

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set (CMAKE_INSTALL_RPATH "/usr/local/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set (CMAKE_CXX_COMPILER_LAUNCHER ccache)

set (THREADS_PREFER_PTHREAD_FLAG ON)

find_package (GDAL 2.2 REQUIRED)
find_library (GEOS_LIBRARY geos PATHS /opt/lib /usr/lib /usr/local/lib)

project (contrem)

# Configure directories###################################################################################

include_directories(BEFORE SYSTEM ./include /opt/include /usr/include/gdal) 
link_directories(. /opt/lib /usr/lib)

# Build support libraries ################################################################################

add_executable (contrem src/contrem.cpp src/reader.cpp src/writer.cpp src/processor.cpp src/stats.cpp)
target_link_libraries (contrem gdal ${GEOS_LIBRARY})

# Install ################################################################################################

install(TARGETS contrem RUNTIME DESTINATION bin LIBRARY DESTINATION lib) 