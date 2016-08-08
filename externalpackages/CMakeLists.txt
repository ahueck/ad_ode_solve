find_package(Git)
if(NOT GIT_FOUND)
  message(WARNING "git executable not found")
endif()

set(BUILD_EXT_DIR ${PROJECT_SOURCE_DIR})

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
  add_dependencies(${MNAME} initSubModules_${MNAME})
endfunction(make_subm)

#############
# CVODE
set(CVODE_HOME ${EXTERNALPACKAGES}/cvode)
set(CVODE_INSTALL_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/cvode)
ExternalProject_Add(cvode
  SOURCE_DIR ${CVODE_HOME}
  CMAKE_ARGS ${CVODE_HOME} -DCMAKE_INSTALL_PREFIX:PATH=${CVODE_INSTALL_DIR}
  EP_BASE ${BUILD_EXT_DIR}
)
make_subm(cvode)
set(ODE_LIBS_CVODE
  sundials_cvode
  sundials_nvecserial 
  PARENT_SCOPE)
set(CVODE_INSTALL ${CVODE_INSTALL_DIR} PARENT_SCOPE)

#############
# CoDiPack
set(CODIPACK_INSTALL ${EXTERNALPACKAGES}/codipack PARENT_SCOPE)
add_library(codipack INTERFACE)
make_subm(codipack)

