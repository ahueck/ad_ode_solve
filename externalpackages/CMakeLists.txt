find_package(Git)
if(NOT GIT_FOUND)
  message(WARNING "git not found")
endif()

include (ExternalProject)

#############
# SubModules
function(make_subm MNAME)
  add_custom_target(
    initSubModules_${MNAME}
    COMMAND ${GIT_EXECUTABLE} submodule update --init -- externalpackages/${MNAME}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  )
  set(submodule_target initSubModules_${MNAME} PARENT_SCOPE)
  add_dependencies(${MNAME} ${submodule_target})
endfunction(make_subm)

#############
# CVODE
set(CVODE_HOME "${EXTERNALPACKAGES}/cvode")
ExternalProject_Add(cvode
  PREFIX ${CVODE_HOME}
  SOURCE_DIR ${CVODE_HOME}
  BINARY_DIR ${CVODE_HOME}/build
  STAMP_DIR ${CVODE_HOME}/build
  TMP_DIR ${CVODE_HOME}/build
#  INSTALL_DIR ${CVODE_HOME}/install
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CVODE_HOME}/install
)
make_subm(cvode)

include_directories(${EXTERNALPACKAGES}/cvode/install/include)
link_directories(${EXTERNALPACKAGES}/cvode/install/lib)

set(ODE_LIBS_CVODE
  sundials_cvode
  sundials_nvecserial)

#############
# CoDiPack
set(CODIPACK_HOME "${EXTERNALPACKAGES}/codipack")
include_directories(${CODIPACK_HOME})
add_library(codipack INTERFACE)
make_subm(codipack)

