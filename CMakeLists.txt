# Minimum CMake version -- currently determined by HDF5.
cmake_minimum_required (VERSION 3.10.0)

# Adjust CMake's module path.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

# Set compilers. This must be done before enabling languages.
include(set_up_platform)
set_up_platform()
enable_language(C)
enable_language(CXX)

# We declare the project here.
project (polymec)

message(STATUS "Generating project files in build directory: ${PROJECT_BINARY_DIR}")
message(STATUS "C compiler is ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION})")
message(STATUS "C++ compiler is ${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION})")

# We require compilers that are pretty C11-compliant at this point.
if (CMAKE_C_COMPILER_ID STREQUAL "GNU" AND CMAKE_C_COMPILER_VERSION STRLESS "4.9")
  message(FATAL_ERROR "GCC v${CMAKE_C_COMPILER_VERSION} detected. v4.9+ is required.")
elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel" AND CMAKE_C_COMPILER_VERSION STRLESS "16.0")
  message(FATAL_ERROR "Intel v${CMAKE_C_COMPILER_VERSION} detected. v16.0+ is required.")
endif()

# Version numbers.
set (POLYMEC_MAJOR_VERSION 2)
set (POLYMEC_MINOR_VERSION 5)
set (POLYMEC_PATCH_VERSION 1)

set (POLYMEC_VERSION "${POLYMEC_MAJOR_VERSION}.${POLYMEC_MINOR_VERSION}.${POLYMEC_PATCH_VERSION}")
message(STATUS "Configuring polymec v${POLYMEC_VERSION}")

if (POLYMEC_PRECISION STREQUAL "single")
  set(HAVE_SINGLE_PRECISION 1)
  set(HAVE_DOUBLE_PRECISION 0)
  set(POLYMEC_REAL_TYPE float)
  set(POLYMEC_MPI_REAL_TYPE MPI_FLOAT)
  set(POLYMEC_MPI_COMPLEX_TYPE MPI_C_FLOAT_COMPLEX)
elseif(POLYMEC_PRECISION STREQUAL "double")
  set(HAVE_SINGLE_PRECISION 0)
  set(HAVE_DOUBLE_PRECISION 1)
  set(POLYMEC_REAL_TYPE double)
  set(POLYMEC_MPI_REAL_TYPE MPI_DOUBLE)
  set(POLYMEC_MPI_COMPLEX_TYPE MPI_C_DOUBLE_COMPLEX)
else()
  message(FATAL_ERROR "Unknown precision: ${POLYMEC_PRECISION}")
endif()
message(STATUS "Real number representation is ${POLYMEC_REAL_TYPE}")

# Report shared library status.
if (BUILD_SHARED_LIBS)
  message(STATUS "Building shared libraries.")
  set(HAVE_SHARED_LIBS 1)
else()
  set(HAVE_SHARED_LIBS 0)
  message(STATUS "Building static libraries.")
endif()

# Add a custom target (always out of date) that generates polymec_version.h.
add_custom_target(update_version_h ALL
                  ${BASH} ${PROJECT_SOURCE_DIR}/tools/update_version_h.sh ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/core/polymec_version.h
                  BYPRODUCTS ${PROJECT_BINARY_DIR}/core/polymec_version.h)
set_target_properties(update_version_h PROPERTIES FOLDER Misc)

# Figure out the system type.
if (APPLE)
  set(SYS_FLAGS "-DAPPLE=1")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Accelerate")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -framework Accelerate")

  # Find the number of cores available on the machine.
  execute_process(COMMAND sysctl -n hw.ncpu OUTPUT_VARIABLE NUMBER_OF_CORES)
  string(REPLACE "\n" "" NUMBER_OF_CORES ${NUMBER_OF_CORES})
  if (NUMBER_OF_CORES LESS 0)
    set(NUMBER_OF_CORES 1)
  endif()

  # Physical cores.
  execute_process(COMMAND sysctl -n hw.physicalcpu OUTPUT_VARIABLE NUMBER_OF_PHYSICAL_CORES)
  string(REPLACE "\n" "" NUMBER_OF_PHYSICAL_CORES ${NUMBER_OF_PHYSICAL_CORES})
  if (NUMBER_OF_PHYSICAL_CORES LESS 0)
    set(NUMBER_OF_PHYSICAL_CORES 1)
  endif()

  message(STATUS "Configuring on Mac with ${NUMBER_OF_CORES} logical (${NUMBER_OF_PHYSICAL_CORES} physical) cores.")
else ()
  if (LINUX)
    set(SYS_FLAGS "-DLINUX=1")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")

    # Find the number of cores available on the machine.
    execute_process(COMMAND nproc OUTPUT_VARIABLE NUMBER_OF_CORES)
    string(REPLACE "\n" "" NUMBER_OF_CORES ${NUMBER_OF_CORES})
    if (NUMBER_OF_CORES LESS 0)
      set(NUMBER_OF_CORES 1)
    endif()

    # Physical cores: tricky in general. We'll parse the lscpu command to
    # find out the number of threads per core and then just divide
    # NUMBER_OF_CORES by that to get NUMBER_OF_PHYSICAL_CORES.
    execute_process(COMMAND lscpu OUTPUT_VARIABLE lscpu_output)
    string(FIND ${lscpu_output} "Thread(s) per core:" pos)
    string(SUBSTRING ${lscpu_output} ${pos} -1 lscpu_output) # chop the front
    string(FIND ${lscpu_output} "\n" pos)
    string(SUBSTRING ${lscpu_output} 0 ${pos} lscpu_output) # chop the back
    string(REPLACE "Thread(s) per core: " "" lscpu_output ${lscpu_output})
    string(REPLACE " " "" THREADS_PER_CORE ${lscpu_output})
    math(EXPR NUMBER_OF_PHYSICAL_CORES "${NUMBER_OF_CORES} / ${THREADS_PER_CORE}")

    message(STATUS "Configuring on Linux with ${NUMBER_OF_CORES} logical (${NUMBER_OF_PHYSICAL_CORES} physical) cores.")
  else()
    message(FATAL_ERROR "Unknown system type! Currently only Mac and Linux are supported.")
  endif ()
endif ()

# General C compiler flags.
if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -Wall -pedantic-errors -Werror -Wextra -Werror-implicit-function-declaration")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wshadow -Wfloat-equal -Wmissing-noreturn -Wmissing-prototypes")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast -Wno-unused-result")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unknown-pragmas -Wno-format -Wno-maybe-uninitialized -Wno-missing-noreturn")
  if (HAVE_DOUBLE_PRECISION)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wfloat-conversion")
  endif()

  # Strip out impertinent warnings for older compilers.
  if (CMAKE_C_COMPILER_VERSION STRLESS "5")
    # Not available before GCC 5.
    string(REPLACE "-Wfloat-conversion" " " CMAKE_C_FLAGS ${CMAKE_C_FLAGS})
  endif()

  if (LINUX)
    # Counter some of GCC's more recent stinginess on Linux.
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_POSIX_C_SOURCE=200809L")# -D_BSD_SOURCE")
    # Pass some more needed flags to the compiler.
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
  endif()
elseif (CMAKE_C_COMPILER_ID MATCHES "Clang")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -Weverything -Werror -pedantic-errors -Werror-implicit-function-declaration -fno-builtin")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare -Wno-unused-parameter -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast -Wno-unused-function")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-reserved-id-macro -Wno-unused-macros -Wno-sign-conversion -Wno-padded -Wno-vla")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-format-nonliteral -Wno-cast-qual -Wno-switch-enum -Wno-unknown-pragmas -Wno-source-uses-openmp")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unknown-warning-option -Wno-documentation-deprecated-sync -Wno-missing-noreturn")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-documentation -Wno-documentation-unknown-command")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-complex-component-init") # For CMPLX macro
  if (NOT HAVE_DOUBLE_PRECISION)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-double-promotion -Wno-conversion")
  endif()

  # Clang proper behaves differently from Apple Clang, and requires
  # more suppressions.
  if (NOT CMAKE_C_COMPILER_ID MATCHES "AppleClang")
    # tolower() is a recursive macro, which causes issues.
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-disabled-macro-expansion")
  endif()
elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -Wall")
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SYS_FLAGS}")

# C++ compiler flags.
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare -Wno-unused-parameter -Wno-int-to-pointer-cast -Wno-unused-function")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare -Wno-unused-parameter -Wno-int-to-pointer-cast -Wno-unused-function")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SYS_FLAGS}")

# If the lld or gold linkers are available, use them (in that order).
find_program(LLD lld)
if (LLD STREQUAL "LLD-NOTFOUND")
  find_program(GOLD gold)
  if (NOT GOLD STREQUAL "GOLD-NOTFOUND")
    message(STATUS "Found gold linker. Using instead of ld.")
    set(CMAKE_LINKER ${GOLD})
  endif()
else()
  message(STATUS "Found lld linker. Using instead of ld.")
  set(CMAKE_LINKER ${LLD})
endif()

# Address Sanitizer
if (USE_ASAN)
  message(STATUS "Enabling Address Sanitizer.")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
endif()

# Considerations for shared libraries, including library suffixes.
if (BUILD_SHARED_LIBS)
  if (LINUX)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  endif()

  # Here we make sure CMake-installed binaries use the correct runpath, and
  # that the path is not stripped during installation.
  set(CMAKE_MACOSX_RPATH TRUE)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

  set(LIB_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
  set(LIB_SUFFIX .a)
endif()

# Support for code coverage.
if (COVERAGE)
  message(STATUS "Enabling code coverage diagnostics.")
  if (CMAKE_C_COMPILER_ID MATCHES "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
  elseif (CMAKE_C_COMPILER_ID MATCHES "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-instr-generate -fcoverage-mapping")
  endif()

  # Set up a coverage target.
  find_program(GCOV gcov)
  if (GCOV STREQUAL "GCOV-NOTFOUND")
    message(FATAL_ERROR "gcov not found--can't measure code coverage!")
  endif()

  find_program(LCOV lcov)
  if (LCOV STREQUAL "LCOV-NOTFOUND")
    message(FATAL_ERROR "lcov not found--can't report code coverage!")
  endif()
  find_program(GENHTML genhtml)
  if (GENHTML STREQUAL "GENHTML-NOTFOUND")
    message(FATAL_ERROR "genhtml not found--can't format code coverage report!")
  endif()

  add_custom_target(coverage ${LCOV} --gcov-tool=${GCOV} --directory ${PROJECT_BINARY_DIR} -c -q -o lcov.info
                    COMMAND ${LCOV} --remove lcov.info '*valgrind.h' '*/tests/*' -o lcov.info
                    COMMAND ${GENHTML} -o coverage lcov.info
                    COMMENT "Writing code coverage report to ${PROJECT_BINARY_DIR}/coverage.")
endif()

# Give a warning about not setting the install prefix.
message(STATUS "Installation prefix is ${CMAKE_INSTALL_PREFIX}")

# Basic libraries to be linked in.
set(POLYMEC_LIBRARIES m)
if (NEED_LAPACK)
  include(FindBLAS)
  include(FindLAPACK)
  find_package(BLAS REQUIRED QUIET)
  if (NOT BLAS_FOUND)
    message(FATAL_ERROR "Couldn't find BLAS. Please make sure it's installed.")
  endif()
  find_package(LAPACK REQUIRED QUIET)
  if (NOT LAPACK_FOUND)
    message(FATAL_ERROR "Couldn't find LAPACK. Please make sure it's installed.")
  endif()
  if (${LAPACK_LIBRARY_DIR})
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${LAPACK_LIBRARY_DIR}")
  endif()
  if (${BLAS_LIBRARY_DIR})
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${BLAS_LIBRARY_DIR}")
  endif()
  set(POLYMEC_LIBRARIES ${POLYMEC_LIBRARIES};${LAPACK_LIBRARIES};${BLAS_LIBRARIES})
endif()

# NOTE: HDF5 requires dynamic loading as of v1.8.11, so we need libdl.
set(POLYMEC_LIBRARIES ${POLYMEC_LIBRARIES};dl)

# Figure out MPI.
if (HAVE_MPI)
  message(STATUS "MPI is enabled.")
  set(HAVE_MPI 1) # Normalize.

  # NOTE: Disable C++ bindings for MPI, since they have never worked for anyone.
  set(NO_MPI_CXX_FLAGS "-DMPICH_SKIP_MPICXX -UHAVE_MPI_CPP -DLAM_WANT_MPI2CPP=0 -DLAM_BUILDING=1 -DOMPI_WANT_CXX_BINDINGS=0 -DOMPI_BUILDING=1")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${NO_MPI_CXX_FLAGS}")
else()
  message(STATUS "MPI is disabled.")
  set(HAVE_MPI 0) # Normalize.

  # Include our own serial implementation of MPI.
  add_subdirectory(mpi_serial)
  include_directories("${PROJECT_SOURCE_DIR}/mpi_serial")
endif ()

# OpenMP support?
if (USE_OPENMP)
  find_package(OpenMP QUIET)
  if (OPENMP_FOUND)
    message(STATUS "OpenMP is enabled.")
    set(HAVE_OPENMP 1)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  else()
    message(STATUS "OpenMP is not supported.")
    set(HAVE_OPENMP 0)
  endif()
else()
  message(STATUS "OpenMP is disabled.")
  set(HAVE_OPENMP 0)
endif()

# Other third-party libraries.
add_subdirectory(3rdparty)

# Record the libraries we've gathered so far as the "base" libraries that
# Polymec needs. We use this for efficiently linking unit tests.
set(POLYMEC_BASE_LIBRARIES ${POLYMEC_LIBRARIES})

# Include the binary directory in the header file search path,
# since it's where we place the third-party libraries.
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_BINARY_DIR}/include")
link_directories("${PROJECT_BINARY_DIR}/lib")
include_directories(${POLYMEC_INCDIRS})

# Support for valgrind -- Linux only.
if (LINUX)
  find_package(Valgrind QUIET)
  if (VALGRIND_FOUND)
    set(VALGRIND_FOUND 1) # regularize this value
    include_directories(${VALGRIND_INCLUDE_DIR})
    set(MEMORYCHECK_COMMAND ${VALGRIND_PROGRAM})
    set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --show-leak-kinds=all --errors-for-leak-kinds=definite,possible --track-origins=yes --error-exitcode=1 --trace-children=yes --suppressions=${PROJECT_SOURCE_DIR}/tools/valgrind/polymec.supp" CACHE STRING "Options passed to Valgrind." FORCE)

    # make memcheck target
    add_custom_target(memcheck ctest -T memcheck -j ${NUMBER_OF_CORES} USES_TERMINAL)
  else()
    set(VALGRIND_FOUND 0)
  endif()
else()
  # Valgrind doesn't work on Macs.
  set(VALGRIND_FOUND 0)
endif()

# Unit testing.
include(CTest)
enable_testing()

# Libraries.
include_directories("${PROJECT_SOURCE_DIR}")
add_subdirectory(core)
add_subdirectory(geometry)
add_subdirectory(solvers)
add_subdirectory(model)
add_subdirectory(io)

# Add a custom target that generates polymec.h and friends at the top level.
add_custom_target(generate_headers ALL
                  ${BASH} ${PROJECT_SOURCE_DIR}/tools/generate_headers.sh ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include core geometry solvers model io)
set_target_properties(generate_headers PROPERTIES FOLDER Misc)

# Now that we have gathered all our libraries, generate a polymec.cmake
# file that contains all the vital information.
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/templates/polymec.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/polymec.cmake"
  @ONLY
)

# Installation.
file(GLOB all_headers "${PROJECT_BINARY_DIR}/include/*")
foreach(header ${all_headers})
  get_filename_component(filename ${header} NAME)
  if (filename MATCHES "polymec")
    list(APPEND polymec_headers ${header})
  elseif(filename MATCHES ".h")
    list(APPEND tp_headers ${header})
  else()
    install(DIRECTORY ${header} DESTINATION include)
  endif()
endforeach()
set(polymec_headers ${polymec_headers} CACHE STRING "generated polymec headers" )

# Install all third-party library Polymec source headers.
install(FILES ${tp_headers} DESTINATION include)

# Install all third-party libraries.
file(GLOB tp_libs "${PROJECT_BINARY_DIR}/lib/*.a")
if (BUILD_SHARED_LIBS)
  file(GLOB shared_tp_libs "${PROJECT_BINARY_DIR}/lib/*${LIB_SUFFIX}*")
  list(APPEND tp_libs ${shared_tp_libs})
endif()
install(FILES ${tp_libs} DESTINATION lib)

# Install all third-party executables.
file(GLOB tp_exes "${PROJECT_BINARY_DIR}/bin/*")
install(FILES ${tp_exes} DESTINATION bin
        PERMISSIONS OWNER_READ;OWNER_WRITE;OWNER_EXECUTE;
                    GROUP_READ;GROUP_EXECUTE;WORLD_READ;WORLD_EXECUTE)

# Install all generated Polymec source headers.
install(FILES ${polymec_headers} DESTINATION include/polymec)

# Install miscellaneous build/test files.
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/polymec.cmake DESTINATION share/polymec)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/add_polymec_executable.cmake DESTINATION share/polymec)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/add_polymec_library.cmake DESTINATION share/polymec)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/add_polymec_test.cmake DESTINATION share/polymec)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/add_polymec_driver_test.cmake DESTINATION share/polymec)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/FindValgrind.cmake DESTINATION share/polymec)

# Additional targets.

find_package(Python COMPONENTS Interpreter)
if (Python_FOUND AND Python_VERSION_MAJOR EQUAL 3)
  message(STATUS "Found Python. Generating 'stats', 'prepend-license' targets.")
  # make stats
  add_custom_target(stats ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/gather_stats.py ${PROJECT_SOURCE_DIR})
  # make prepend-license
  add_custom_target(prepend-license ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/prepend_license.py ${PROJECT_SOURCE_DIR})
endif()

# Documentation.

find_package(Doxygen OPTIONAL_COMPONENTS dot QUIET)
if (DOXYGEN_FOUND)
  message(STATUS "Found Doxygen. Generating 'doc' target.")
  set(DOXYGEN_EXTRACT_LOCAL_CLASSES NO)
  set(DOXYGEN_FILE_PATTERNS "*.h")
  set(DOXYGEN_GENERATE_HTML YES)
  set(DOXYGEN_GENERATE_MAN YES)
  set(DOXYGEN_GENERATE_LATEX NO)
  set(DOXYGEN_GENERATE_DOCBOOK NO)
  set(DOXYGEN_GENERATE_XML NO)
  set(DOXYGEN_USE_MATHJAX YES)
  set(DOXYGEN_HIDE_UNDOC_CLASSES YES)
  set(DOXYGEN_INCLUDE_PATH "${PROJECT_SOURCE_DIR};${PROJECT_BINARY_DIR}/include")
  set(DOXYGEN_JAVADOC_AUTOBRIEF YES)
  set(DOXYGEN_OPTIMIZE_OUTPUT_FOR_C YES)
  set(DOXYGEN_QUIET YES)
  set(DOXYGEN_WARN_LOGFILE Doxygen_warnings.log)
  set(DOXYGEN_SHOW_FILES NO)
  set(DOXYGEN_SHOW_NAMESPACES NO)
  set(DOXYGEN_CALL_GRAPH YES)
  set(DOXYGEN_ALIASES
      [[collective="\xrefitem collective \"Collective functions\" \"Collective (MPI) functions\""]]
      [[refcounted="\xrefitem refcounted \"Reference counted resources\" \"Reference counted resources\""]]
      [[experimental="\xrefitem experimental \"Experimental code: not extensively tested!\" \"Experimental code: not extensively tested!\""]]
     )
  set(DOXYGEN_VERBATIM_VARS DOXYGEN_ALIASES)
  doxygen_add_docs(doc ${PROJECT_BINARY_DIR}/core/polymec.h;core;geometry;solvers;model;io
                   COMMENT "Generating Documentation")
endif()
