project(GALCONDA C CXX)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_VERBOSE_MAKEFILE off)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_C_FLAGS "-O2")
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "-O2 -ftree-vectorize -std=c++11")
set(GALAXIONS_INCLUDES)
set(GALAXIONS_LIBRARIES)

# YMW16 code (provided)
add_subdirectory(libs/ymw16)
list(APPEND GALCONDA_LIBRARIES ymw16)
list(APPEND GALCONDA_INCLUDES libs/ymw16)

# include
include_directories("${CMAKE_SOURCE_DIR}/include")

# add header files
set(GALCONDA_HDRS
        include/cgs.h
        include/constants.h
        include/damping.h
        include/galconda.h
        include/gas.h
        include/JF12.h
        include/los.h
        include/magneticfield.h
        include/mymatrix.h
        include/Pshirkov.h
        include/solvers.h)

# add source files
set(GALCONDA_SRCS
        src/gas.cpp
        src/galconda.cpp
        src/Ferriere07.cpp
        src/JF12.cpp
        src/magneticfield.cpp
        src/mymatrix.cpp
        src/Pshirkov.cpp
        src/solvers.cpp
        )

# make library
include_directories(include ${GALCONDA_INCLUDES})
add_library(GALCONDA ${GALCONDA_SRCS} ${GALCONDA_HDRS} ${GALCONDA_INCLUDES})

# make executables
add_executable(los src/main.cpp)
target_link_libraries(los GALCONDA ${GALCONDA_LIBRARIES})

#add_executable(map src/calculateMaps.cpp)
#target_link_libraries(map GALCONDA ${GALCONDA_LIBRARIES})

# show summary
message(STATUS "CMAKE_SYSTEM .......... = ${CMAKE_SYSTEM}")
message(STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message(STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message(STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")
message(STATUS "GALCONDA_INCLUDES ..... = ${GALCONDA_INCLUDES}")
message(STATUS "GALCONDA_LIBRARIES .... = ${GALCONDA_LIBRARIES}")
