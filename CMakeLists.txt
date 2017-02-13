CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(RoberODE)

include(${PROJECT_SOURCE_DIR}/cmake/ToolchainOptions.cmake)

set(CMAKE_VERBOSE_MAKEFILE TRUE)

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

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${CVODE_INSTALL}/include)
link_directories(${CVODE_INSTALL}/lib)
include_directories(${CODIPACK_INSTALL}/include)

add_library(ad_ode src/SolverCVode.cpp)
add_dependencies(ad_ode cvode codipack)
target_compile_options(ad_ode PRIVATE -fno-rtti)
target_link_libraries(ad_ode ${ODE_LIBS_CVODE})

add_subdirectory(examples)
add_subdirectory(test)

enable_testing()
add_test(NAME RunTests COMMAND test_runner${EXE_SUFFIX})
