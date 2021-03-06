include(add_polymec_executable)

# Travis CI's Docker environment is pretty spare, with 1 single physical core
# for testing. 
if (TRAVIS_CI)
  # On Travis, we use the number of logical cores for MPI testing.
  set(NUMBER_OF_TEST_CORES ${NUMBER_OF_CORES}) 
else()
  # ...otherwise we use the number of physical cores.
  set(NUMBER_OF_TEST_CORES ${NUMBER_OF_PHYSICAL_CORES}) 
endif()

function(add_polymec_driver_with_libs driver_name libs driver_source)
  add_polymec_executable_with_libs(${driver_name}_exe "${libs}" ${driver_source})
endfunction()

function(add_polymec_driver_test test_name driver_name test_script)
  if (HAVE_MPI)
    foreach (arg ${ARGN})
      if (arg MATCHES "=")
        set(options ${options} ${arg})
      else()
        list(APPEND procs ${arg})
      endif()
    endforeach()

    list(LENGTH procs procs_n)
    if (procs_n EQUAL 0)
      add_test(${test_name} ${CMAKE_CURRENT_BINARY_DIR}/${driver_name}_exe ${test_script} ${options})
      set_tests_properties(${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION "${test_script}:")
    else()
      foreach (proc ${procs})
        if (NOT ${proc} GREATER ${NUMBER_OF_TEST_CORES})
          set(run_test_name ${test_name}_${proc}_proc)
          add_test(${run_test_name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${proc} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${driver_name}_exe ${test_script} ${options} ${MPIEXEC_POSTFLAGS})
          set_tests_properties(${run_test_name} PROPERTIES FAIL_REGULAR_EXPRESSION "${test_script}:")

          # If this test needs more than half of the available cores, set it
          # to run by itself.
          math(EXPR HALF_OF_CORES "${NUMBER_OF_CORES} / 2") 
          if (${proc} GREATER ${HALF_OF_CORES})
            set_tests_properties(${run_test_name} PROPERTIES RUN_SERIAL ON)
          endif()
        endif()
      endforeach()
    endif()
  else()
    foreach (arg ${ARGN})
      if (arg MATCHES "=")
        set(options "${options} ${arg}")
      endif()
    endforeach()

    add_test(${test_name} ${CMAKE_CURRENT_BINARY_DIR}/${driver_name}_exe ${test_script} ${options})
    set_tests_properties(${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION "${test_script}:")
  endif()
endfunction()

