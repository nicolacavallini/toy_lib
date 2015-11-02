cmake_minimum_required (VERSION 2.6)
project(time_testing)

include(CheckCXXCompilerFlag)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories ("$ENV{TOY_LIB_PATH}/include")
set (EXTRA_LIBS ${EXTRA_LIBS} $ENV{TOY_LIB_PATH}/bin/libgrid.a)
set (EXTRA_LIBS ${EXTRA_LIBS} $ENV{TOY_LIB_PATH}/bin/liblac.a)

find_package( Boost 1.54 COMPONENTS program_options REQUIRED system filesystem)
include_directories( ${Boost_INCLUDE_DIR} )

add_executable(time_testing main.cpp)
target_link_libraries(time_testing ${EXTRA_LIBS} ${Boost_LIBRARIES})
