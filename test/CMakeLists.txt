include_directories(${CATCH_INSTALL}/single_include)
include_directories(${PROJECT_SOURCE_DIR}/test/include)

set(unit_test_names 
  src/test_ad.cpp 
)

add_executable(test_runner${EXE_SUFFIX} src/main.cpp ${unit_test_names})