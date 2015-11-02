set(expected_level 3)


message ("******************************************************") 
set(level 0)
# Structure check, expected directories and files exist
if(NOT EXISTS ${WORK_DIR})
    message(FATAL_ERROR "** Work directory does not exist.")
endif()
if(NOT EXISTS ${SOURCE_DIR}/expected.txt)    
    message(FATAL_ERROR "** expected.txt does not exist.")
endif()

message ("** Compiling and linking test program")
execute_process(
  COMMAND         make ${EXE_TARGET} 
  RESULT_VARIABLE compile_error)

set(runtime_error "Not run")
if(NOT compile_error)
  set(level 2)
  message ("** Executing test program.")
  execute_process(
    COMMAND           ${TEST_PROG} 
    OUTPUT_FILE       stdout.txt
    ERROR_FILE        stdout.txt
    WORKING_DIRECTORY ${WORK_DIR} 
    RESULT_VARIABLE   runtime_error)
endif(NOT compile_error)


if(NOT runtime_error)
  set(level 3)
endif(NOT runtime_error)

set(level_ok FALSE)
if (level EQUAL expected_level)
  set(level_ok TRUE)
endif()

set(output_error "Not compared")
if (level_ok)
  message ("** Comparing test output with expected one.")

  set(output_file  '')
  if(EXISTS ${WORK_DIR}/output.txt)
    set(output_file ${WORK_DIR}/output.txt)
  elseif(EXISTS ${WORK_DIR}/stdout.txt)
    set(output_file ${WORK_DIR}/stdout.txt)
  endif()

  execute_process(
    COMMAND ${CMAKE_COMMAND} -E compare_files ${output_file} ${SOURCE_DIR}/expected.txt
    RESULT_VARIABLE output_error)
endif()

message("******************************************************") 
message("            Test results")
message("            ------------")
message("Compiling and linking: ${compile_error}")
message("Runnig program:        ${runtime_error}")
message("Expected output:       ${output_error}\n")


if(output_error)
  message ("Test status: FAILED.")
  message ("Aditional Information:")
  if(level_ok)
    message("  Output different from expected.")
    message("  Run something like the following for further details:\n")
    message("  tkdiff ${output_file} ${SOURCE_DIR}/expected.txt\n")
  endif()
  if(compile_error)
    message("  Compile or linking Error. See the test output above for details.\n")
  elseif(runtime_error)
    message("  Run Time Error.")
    message("  You can try to run the executable with:\n   cd ${WORK_DIR} & ${TEST_PROG}\n")
  endif()
  message(FATAL_ERROR " ")
endif()

message ("Test status: PASSED.")
message ("******************************************************")
