cmake_minimum_required (VERSION 3.10)

project (tricl CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# default to an optimized build (which also disables assertions via NDEBUG);
# use -DCMAKE_BUILD_TYPE=Debug for a build with assertions and debug symbols:
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (Release, Debug, RelWithDebInfo)" FORCE)
endif()
message(STATUS "tricl build type: ${CMAKE_BUILD_TYPE}")

add_compile_options(-Wall)

find_package(ZLIB REQUIRED)
find_package(Boost REQUIRED)      # only header-only parts are used (container, iostreams)
find_package(yaml-cpp REQUIRED)

add_subdirectory(src)
