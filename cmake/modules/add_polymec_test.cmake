include_directories(${PROJECT_SOURCE_DIR}/tests;${PROJECT_BINARY_DIR}/include)
link_directories(${PROJECT_BINARY_DIR}/lib)

if (TRAVIS_CI)
  # On Travis CI, we use the number of logical cores for MPI testing.
  set(NUMBER_OF_TEST_CORES ${NUMBER_OF_CORES}) 

  # Travis CI and docker run everything as root, and our MPI environment there
  # (OpenMPI) really hates being run as root, so we have to ask it nicely to 
  # do so.
  if (MPIEXEC_PREFLAGS)
    set(MPIEXEC_PREFLAGS "${MPIEXEC_PREFLAGS} --allow-run-as-root")
  else()
    set(MPIEXEC_PREFLAGS "--allow-run-as-root")
  endif()
else()
  # When we're not running Travis CI, we use the number of physical cores
  # for MPI testing.
  set(NUMBER_OF_TEST_CORES ${NUMBER_OF_PHYSICAL_CORES}) 
endif()

# This function adds a parallel unit test executable to be built using cmocka,
# linking against the specified libraries. Other arguments are source files 
# and numbers of processes (in no particular order, but usually with process 
# counts following source files by convention). 1 test run will be generated 
# for each processor number value.
function(add_mpi_polymec_test_with_libs exe libs)
  foreach (arg ${ARGN})
    if (arg MATCHES ".c")
      list(APPEND sources ${arg})
    else()
      list(APPEND procs ${arg})
    endif()
  endforeach()
  add_executable(${exe} ${sources})
  target_link_libraries(${exe} cmocka ${libs})
  set_target_properties(${exe} PROPERTIES COMPILE_FLAGS "-DCMAKE_CURRENT_SOURCE_DIR=\\\"${CMAKE_CURRENT_SOURCE_DIR}\\\"")
  set_target_properties(${exe} PROPERTIES FOLDER Tests)
  if (HAVE_MPI)
    foreach (proc ${procs})
      if (NOT ${proc} GREATER ${NUMBER_OF_TEST_CORES})
        add_test(${exe}_${proc}_proc ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${proc} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${exe} ${MPIEXEC_POSTFLAGS})
        set_tests_properties(${exe}_${proc}_proc PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

        # If this test needs more than half of the available cores, set it
        # to run by itself.
        math(EXPR HALF_OF_CORES "${NUMBER_OF_TEST_CORES} / 2") 
        if (${proc} GREATER ${HALF_OF_CORES})
          set_tests_properties(${exe}_${proc}_proc PROPERTIES RUN_SERIAL ON)
        endif()
      endif()
    endforeach()
  else()
    # We only add a single-process test case when MPI is not present.
    add_test(${exe}_1_proc ${CMAKE_CURRENT_BINARY_DIR}/${exe})
    set_tests_properties(${exe}_1_proc PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  endif()
endfunction()

# This function adds a (serial) unit test executable to be built using cmocka,
# linking against the specified libraries.
function(add_polymec_test_with_libs exe libs)
  if (DEFINED BATCH_SYSTEM)
    add_mpi_polymec_test_with_libs(${exe} ${libs} ${ARGN} 1)
  else()
    add_executable(${exe} ${ARGN})
    target_link_libraries(${exe} cmocka ${libs})
    set_target_properties(${exe} PROPERTIES COMPILE_FLAGS "-DCMAKE_CURRENT_SOURCE_DIR=\\\"${CMAKE_CURRENT_SOURCE_DIR}\\\"")
    set_target_properties(${exe} PROPERTIES FOLDER Tests)
    add_test(${exe} ${exe})
    set_tests_properties(${exe} PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  endif()
endfunction()

# This function adds a (serial) simulation test run.
function(add_polymec_test_run exe input)
  file(GLOB_RECURSE sim_exe ${PROJECT_BINARY_DIR}/${exe})
  add_test(${exe}_run_${input} ${sim_exe} run ${input} ${ARGN})
  set_tests_properties(${exe}_run_${input} PROPERTIES FAIL_REGULAR_EXPRESSION "Fatal error:")
endfunction()

