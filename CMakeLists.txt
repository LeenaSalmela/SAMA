# CMake compatibility issues: don't modify this, please!
cmake_minimum_required(VERSION 2.6.3)

project(MAGA)

# project version
set(${PROJECT_NAME}_MAJOR_VERSION 0)
set(${PROJECT_NAME}_MINOR_VERSION 0)
set(${PROJECT_NAME}_PATCH_LEVEL 0)

# set some definitions
add_definitions("-D${PROJECT_NAME}_MAJOR_VERSION=${${PROJECT_NAME}_MAJOR_VERSION}")
add_definitions("-D${PROJECT_NAME}_MINOR_VERSION=${${PROJECT_NAME}_MINOR_VERSION}")
add_definitions("-D${PROJECT_NAME}_PATCH_LEVEL=${${PROJECT_NAME}_PATCH_LEVEL}")

# set the default configuration to Release
if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# uncomment the portion below to disable assertions
if (CMAKE_BUILD_TYPE STREQUAL Release)
        add_definitions(-DNDEBUG)
else (CMAKE_BUILD_TYPE STREQUAL Release)
        add_definitions(-DDEBUG)
endif (CMAKE_BUILD_TYPE STREQUAL Release)

# set the module path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# check if zlib is present
find_package(ZLIB)
if (ZLIB_FOUND)
    add_definitions(-DHAVE_ZLIB)
    include_directories(${ZLIB_INCLUDE_DIR})
endif(ZLIB_FOUND)

set(CMAKE_VERBOSE_MAKEFILE ON)

# set include path for Google's sparse hash table
find_package(SparseHash)
if (SPARSEHASH_FOUND)
    include_directories(${SPARSEHASH_INCLUDE_DIR})
else (SPARSEHASH_FOUND)
    message(FATAL_ERROR "\nFATAL ERROR: The required Google SparseHash package"
            " could not be found on this system.  Please refer to the MAGA "
            "manual for the Google Sparsehash installation instructions.  If "
            "you installed Google Sparsehash in a non-standard location "
            "(e.g. somewhere in your homedir), you can point cmake to the "
            "installation location as follows: \ncmake "
            "-DSPARSEHASH_INCLUDE_DIR=<path-to-sparsehash>/include .")
endif(SPARSEHASH_FOUND)

set(CMAKE_CXX_FLAGS "-Wall -g")

add_subdirectory(src)
