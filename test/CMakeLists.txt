include_directories(${CATCH_INSTALL}/single_include)
include_directories(${PROJECT_SOURCE_DIR}/test/include)

set(unit_test_names 
  src/test_ad.cpp
  src/test_util.cpp  
  src/test_ode.cpp
  src/test_cvode.cpp
)

add_executable(test_runner${EXE_SUFFIX} src/main.cpp ${unit_test_names})
add_dependencies(test_runner${EXE_SUFFIX} catch)
target_link_libraries(test_runner${EXE_SUFFIX} ad_ode)
