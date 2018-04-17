project(GALAXIONS CXX)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_VERBOSE_MAKEFILE off)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "-O2 -ftree-vectorize -std=c++11 -ferror-limit=3")
set(GALAXIONS_INCLUDES)
set(GALAXIONS_LIBRARIES)

# YMW16 code (provided)
#add_subdirectory(libs/ymw16)
#list(APPEND GALAXIONS_LIBRARIES ymw16) 
#list(APPEND GALAXIONS_INCLUDES libs/ymw16)

# include
include_directories( "${CMAKE_SOURCE_DIR}/include" )

# add source files
set( GALAXIONS_SRCS
     src/Farrar.cpp
     src/Pshirkov.cpp
     src/Ferriere07.cpp
     src/galaxions.cpp
     src/magnfield.cpp
     src/mymatrix.cpp
     src/solvers.cpp
)

# make library
include_directories (include ${GALAXIONS_INCLUDES})
add_library( GALAXIONS ${GALAXIONS_SRCS} ${GALAXIONS_INCLUDES} )

# make executables
add_executable(los src/main.cpp)
target_link_libraries (los GALAXIONS ${GALAXIONS_LIBRARIES})

add_executable(map src/calculateMaps.cpp )
target_link_libraries (map GALAXIONS ${GALAXIONS_LIBRARIES})

# show summary
message (STATUS "CMAKE_SYSTEM .......... = ${CMAKE_SYSTEM}")
message (STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message (STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")
message (STATUS "GALAXIONS_INCLUDES ....... = ${GALAXIONS_INCLUDES}")
message (STATUS "GALAXIONS_LIBRARIES ...... = ${GALAXIONS_LIBRARIES}")
