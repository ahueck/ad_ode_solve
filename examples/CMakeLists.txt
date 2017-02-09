#############
# Rober ODE
add_executable(ex_rober${EXE_SUFFIX} rober/main.cpp)
add_dependencies(ex_rober${EXE_SUFFIX} cvode codipack)
target_link_libraries(ex_rober${EXE_SUFFIX} ${ODE_LIBS_CVODE})
