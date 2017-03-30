CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(RoberODE)

include(${PROJECT_SOURCE_DIR}/cmake/ToolchainOptions.cmake)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT BINARY_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
else()
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/${BINARY_OUTPUT_DIRECTORY})
endif()
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
message(STATUS "Binary output directory: " ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})


set(EXTERNALPACKAGES ${PROJECT_SOURCE_DIR}/externalpackages)
add_subdirectory(${EXTERNALPACKAGES})

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${CVODE_INSTALL}/include)
link_directories(${CVODE_INSTALL}/lib)
include_directories(${CODIPACK_INSTALL}/include)
include_directories(${MATPLOTLIB_INSTALL})


set(solver_src
  src/SolverCVode.cpp)

add_library(ad_ode ${solver_src})
add_dependencies(ad_ode cvode codipack matplotlib-cpp)

target_link_libraries(ad_ode ${ODE_LIBS_CVODE} ${PYTHON_LIBRARIES})

add_format_target(format-all
  "Formats ad_ode_solve source files"
)

add_tidy_target(tidy-run
  "Clang-tidy run on ad_ode_solve translation units"
  SOURCES ${solver_src}
  OTHER -header-filter='^${PROJECT_SOURCE_DIR}/include/.*h'
)

add_tidy_fix_target(tidy-core-fix
  "Clang-tidy core run with fixes on ad_ode_solve translation units"
  SOURCES ${solver_src}
  OTHER -checks=-*,modernize-*,llvm-header-guard,llvm-namespace-comment,google-explicit-constructor
        -header-filter='^${PROJECT_SOURCE_DIR}/include/.*h'
)

enable_testing()
add_test(NAME RunTests COMMAND test_runner${EXE_SUFFIX})

add_subdirectory(examples)
add_subdirectory(test)
