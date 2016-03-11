cmake_minimum_required(VERSION 2.8.11)
project(ExplicitBeamFEA)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

option(EXPLICIT_FEA_BUILD_UNIT_TESTS "Build unit tests" ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -DEIGEN_NO_DEBUG -fopenmp")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(EXT_EIGEN_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/ext/eigen3")
set(GTEST_SUBDIR "${CMAKE_CURRENT_SOURCE_DIR}/ext/googletest/googletest")
set(RAPIDJSON_SUBDIR "${CMAKE_CURRENT_SOURCE_DIR}/ext/rapidjson")
set(BOOST_SUBDIR "${CMAKE_CURRENT_SOURCE_DIR}/ext/boost_1_59_0")
set(TCLAP_SUBDIR "${CMAKE_CURRENT_SOURCE_DIR}/ext/tclap")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
                    ${EXT_EIGEN_ROOT}
                    ${RAPIDJSON_SUBDIR}/include
                    ${BOOST_SUBDIR}
                    ${TCLAP_SUBDIR}/include)

set(BOOST_SOURCES ${BOOST_SUBDIR}/libs/smart_ptr/src/sp_collector.cpp
                  ${BOOST_SUBDIR}/libs/smart_ptr/src/sp_debug_hooks.cpp)

add_subdirectory(src)

if (EXPLICIT_FEA_BUILD_UNIT_TESTS)
    enable_testing()
    add_subdirectory(${GTEST_SUBDIR})
endif(EXPLICIT_FEA_BUILD_UNIT_TESTS)