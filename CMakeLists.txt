cmake_minimum_required (VERSION 2.6)
project (grid CXX C)

include(CheckCXXCompilerFlag)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_subdirectory (grid)

option(DEBUG "Add debug" ON)

if(DEBUG)
        add_definitions(-g -DDEBUG)
endif(DEBUG)

add_subdirectory (lac)
