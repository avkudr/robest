cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(robest CXX)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_BUILD_TYPE Release)

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    set(CMAKE_CXX_FLAGS         "-Wall -Wno-unknown-pragmas -Wno-sign-compare -Woverloaded-virtual -Wwrite-strings -Wno-unused")
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()

include_directories(
    ${PROJECT_SOURCE_DIR}/src
)

add_library(robest)
target_include_directories(robest INTERFACE src/)

set(GOOGLETEST_ROOT 3rdparty/googletest/googletest CACHE STRING "Google Test source root")

include_directories(
    ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT}
    ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT}/include
)

set(GOOGLETEST_SOURCES
    ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT}/src/gtest-all.cc
    ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT}/src/gtest_main.cc
)

foreach(_source ${GOOGLETEST_SOURCES})
    set_source_files_properties(${_source} PROPERTIES GENERATED 1)
endforeach()

add_library(googletest ${GOOGLETEST_SOURCES})

find_package(OpenMP QUIET) # Find the package
if (OpenMP_FOUND)
    message("Enabling OpenMP support")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    target_link_libraries(${PROJECT_NAME} ${OpenMP_CXX_LIBRARIES}) # Link against it for C++
endif()

add_executable(
    unit_tests
    tests/LineFitting/LineFitting.hpp
    tests/LineFitting/LineFitting.cpp
    tests/CircleFitting/CircleFitting.hpp
    tests/CircleFitting/CircleFitting.cpp
    tests/SphereFitting/SphereFitting.hpp
    tests/SphereFitting/SphereFitting.cpp
    tests/PlaneFitting/PlaneFitting.hpp
    tests/PlaneFitting/PlaneFitting.cpp
    tests/test_iterEstimation.cpp
    tests/test_LineFitting.cpp
    tests/test_CircleFitting.cpp
    tests/test_SphereFitting.cpp
    tests/test_PlaneFitting.cpp
    tests/main.cpp
)

add_dependencies(unit_tests googletest)

target_link_libraries(
    unit_tests
    googletest
    robest
    )

include(CTest)
enable_testing()

add_test(unit ${PROJECT_BINARY_DIR}/unit_tests)
