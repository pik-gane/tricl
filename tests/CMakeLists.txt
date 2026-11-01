# C++ unit tests (run via ctest or tests/run_tests.py)

add_executable(test_sigmoid test_sigmoid.cpp)
target_link_libraries(test_sigmoid tricl_core)
add_test(NAME sigmoid COMMAND test_sigmoid)
