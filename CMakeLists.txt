CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
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

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${EXTERNALPACKAGES}/cvode/install/include)
include_directories(${EXTERNALPACKAGES}/codipack)

link_directories(${EXTERNALPACKAGES}/cvode/install/lib)

set(ODE_LIBS
  sundials_cvode
  sundials_nvecserial)

#add_executable(rober${EXE_SUFFIX} src/main.cpp)
#target_link_libraries(rober${EXE_SUFFIX} ${ODE_LIBS})
