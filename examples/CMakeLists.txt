#############
# Rober ODE (stiff)
add_executable(ex_rober${EXE_SUFFIX} rober/main.cpp)
target_link_libraries(ex_rober${EXE_SUFFIX} ad_ode)


#############
# Van der Pol ODE (stiff)
add_executable(ex_vdp${EXE_SUFFIX} van_der_pol/main.cpp)
target_link_libraries(ex_vdp${EXE_SUFFIX} ad_ode)