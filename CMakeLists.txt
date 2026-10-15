cmake_minimum_required(VERSION 3.14)
project(stl_star_monitor CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(STL_STAR_BUILD_TESTS    "Build the test programs" ON)
option(STL_STAR_BUILD_EXAMPLES "Build the example programs" ON)

if (MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# ---------------------------------------------------------------------------
# Library
# ---------------------------------------------------------------------------
add_library(stl_star
    src/trace.cpp
    src/formula.cpp
    src/formula_parser.cpp
    src/interval_operations.cpp
    src/signal_constraint_updater.cpp
    src/monitor_common.cpp
    src/boolean_monitor_one.cpp
    src/boolean_monitor_general.cpp
    src/robustness_monitor_one.cpp
    src/robustness_monitor_general.cpp
    src/naive_monitors.cpp          # reference monitors, used to cross-check the fast ones
)
target_include_directories(stl_star PUBLIC include)

# ---------------------------------------------------------------------------
# Examples
# ---------------------------------------------------------------------------
if (STL_STAR_BUILD_EXAMPLES AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/examples)
    add_executable(examples examples/examples.cpp)
    target_link_libraries(examples PRIVATE stl_star)
endif()

# ---------------------------------------------------------------------------
# Tests (each program returns non-zero on any failure; run them with ctest)
# ---------------------------------------------------------------------------
if (STL_STAR_BUILD_TESTS AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    set(STL_STAR_TESTS
        test_parser_and_api                # parser, trace validation, one-call API
        test_boolean_monitor_one           # STL*_1: oracle vs naive vs fast, operators, scaling
        test_boolean_monitor_general       # general STL*: oracle vs naive vs fast, scaling
        test_robustness_monitor_one        # STL*_1 robustness vs exact robustness
        test_robustness_monitor_general    # general robustness vs exact robustness
    )
    enable_testing()
    foreach(t ${STL_STAR_TESTS})
        add_executable(${t} tests/${t}.cpp)
        target_link_libraries(${t} PRIVATE stl_star)
        target_include_directories(${t} PRIVATE tests)
        add_test(NAME ${t} COMMAND ${t})
    endforeach()
endif()
