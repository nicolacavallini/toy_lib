cmake_minimum_required (VERSION 2.6)
project(tests)

message("Configuring unit tests")

include(CheckCXXCompilerFlag)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#--------------------------------------------------------------

include_directories ("$ENV{TOY_LIB_PATH}/include")
set (EXTRA_LIBS ${EXTRA_LIBS} $ENV{TOY_LIB_PATH}/bin/libgrid.a)
set (EXTRA_LIBS ${EXTRA_LIBS} $ENV{TOY_LIB_PATH}/bin/liblac.a)

set (test_dirs grid base lac)

foreach(dir ${test_dirs})
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir})
  file(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.cpp")
  foreach(filename ${files})
    get_filename_component(name ${filename} NAME_WE)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir}/${name})
    add_executable(${dir}/${name}.exe EXCLUDE_FROM_ALL ${dir}/${name}.cpp)
    target_link_libraries (${dir}/${name}.exe ${EXTRA_LIBS})
  endforeach()
endforeach()
#--------------------------------------------------------------



# Create the tests 
if (submit_test)
  include(CTest)
endif()
enable_testing()
foreach(dir ${test_dirs})
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir})
  file(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.cpp")
  foreach(filename ${files})
    get_filename_component(name ${filename} NAME_WE)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir}/${name}/)
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/${name}
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/${dir}/
      FILES_MATCHING PATTERN "*.xml")
    add_test(
      NAME ${dir}-${name}
      COMMAND ${CMAKE_COMMAND}
      -DEXE_TARGET=${dir}/${name}.exe
      -DTEST_PROG=${CMAKE_CURRENT_BINARY_DIR}/${dir}/${name}.exe
      -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/${dir}/${name}
      -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/${dir}/${name}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/runtest.cmake)
  endforeach()
endforeach()

#--------------------------------------------------------------
# Add a target for each .cpp file in test_dirs/*.cpp
#

set (blas_test_dirs lac)

find_package(BLAS)

if(BLAS_FOUND)
    foreach(dir ${blas_test_dirs})
      set(dir_blas ${dir}-blas)
      file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir_blas})
      file(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.cpp")
      foreach(filename ${files})
        get_filename_component(name ${filename} NAME_WE)
        file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir_blas}/${name})
        add_executable(${dir_blas}/${name}.exe EXCLUDE_FROM_ALL ${dir}/${name}.cpp)
        target_link_libraries (${dir_blas}/${name}.exe ${EXTRA_LIBS} blas)
        add_definitions( -DUSEBLAS)
      endforeach()
    endforeach()

    # Create the tests 
    if (submit_test)
      include(CTest)
    endif()
    enable_testing()
    foreach(dir ${blas_test_dirs})
      set(dir_blas ${dir}-blas)
      file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir_blas})
      file(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.cpp")
      foreach(filename ${files})
        get_filename_component(name ${filename} NAME_WE)
        file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir_blas}/${name}/)
        file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/${name}
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/${dir_blas}/
          FILES_MATCHING PATTERN "*.xml")
        add_test(
          NAME ${dir_blas}-${name}
          COMMAND ${CMAKE_COMMAND}
          -DEXE_TARGET=${dir_blas}/${name}.exe
          -DTEST_PROG=${CMAKE_CURRENT_BINARY_DIR}/${dir_blas}/${name}.exe
          -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/${dir}/${name}
          -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/${dir_blas}/${name}
          -P ${CMAKE_CURRENT_SOURCE_DIR}/runtest.cmake)
      endforeach()
    endforeach()    
endif(BLAS_FOUND)
