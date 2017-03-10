cmake_minimum_required(VERSION 3.0)

project(GraphHarmonicMap)
set (GraphHarmonicMap_VERSION_MAJOR 1)
set (GraphHarmonicMap_VERSION_MINOR 0)

set(PROJECT_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${PROJECT_ROOT}/bin CACHE PATH "" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)


# Enable folders for projects in Visual Studio
# if (CMAKE_GENERATOR MATCHES "Visual Studio")
#  set_property(GLOBAL PROPERTY USE_FOLDERS ON)
# endif()

OPTION (USE_OpenMP "Use OpenMP" ON)
if(USE_OpenMP)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        set(CMAKE_C_COMPILER clang-omp CACHE STRING "C compiler" FORCE)
        set(CMAKE_CXX_COMPILER clang-omp++ CACHE STRING "C++ compiler" FORCE)
    endif()

    find_package(OpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

# set vc runtime library as /MT
if(${MSVC})
    foreach(flag
            CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
            CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
    if(${flag} MATCHES "/MD")
        string(REGEX REPLACE "/MD" "/MT" ${flag} "${${flag}}")
    endif(${flag} MATCHES "/MD")
    endforeach(flag)
    # pass to external projects
    set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_RELEASE} /MTd")
    set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
    set (CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} /MT")
    set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} /MT")
endif()

include_directories(third-party/MeshLib/base)
include_directories(third-party/MeshLib/core)
include_directories(third-party/eigen)
#include_directories(third-party/lemon/include)

set(LEMON_PATH "${PROJECT_ROOT}/third-party/lemon/")
set(LEMON_INSTALL_DIR "${PROJECT_ROOT}/third-party/lemon/")
include(ExternalProject)
ExternalProject_Add(liblemon
    PREFIX "${LEMON_PATH}"
    BINARY_DIR "${LEMON_PATH}/build"
    STAMP_DIR  "${LEMON_PATH}/build"
    SOURCE_DIR "${LEMON_PATH}"
    INSTALL_DIR "${LEMON_INSTALL_DIR}"
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${LEMON_INSTALL_DIR} -DBUILD_SHARED_LIBS=OFF 
        -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
        -DCMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}
        -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}
        -DCMAKE_CXX_FLAGS_MINSIZEREL=${CMAKE_CXX_FLAGS_MINSIZEREL}
        -DCMAKE_CXX_FLAGS_RELWITHDEBINFO=${CMAKE_CXX_FLAGS_RELWITHDEBINFO}
)
include_directories(third-party/lemon/include)
add_library(lemon STATIC IMPORTED)
set_property(TARGET lemon PROPERTY IMPORTED_LOCATION ${LEMON_INSTALL_DIR}/lib/lemon.*)
add_dependencies(lemon liblemon)
#ExternalProject_Get_Property(lemon install_dir)

# detect if x64
#if(CMAKE_SIZEOF_VOID_P EQUAL 8)
#endif()
#link_directories(third-party/lemon)

add_executable(GraphHarmonicMap 
    source/GraphHarmonicMap.cpp
    source/Graph.cpp
)
set_property(TARGET GraphHarmonicMap PROPERTY CXX_STANDARD 14)
set_property(TARGET GraphHarmonicMap PROPERTY CXX_STANDARD_REQUIRED ON)

target_link_libraries(GraphHarmonicMap lemon)

install (TARGETS GraphHarmonicMap DESTINATION ${PROJECT_ROOT}/bin)
