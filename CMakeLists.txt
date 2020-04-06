cmake_minimum_required (VERSION 2.6)

project (tricl)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
#set(CMAKE_CXX_FLAGS "-O3 -Wall -std=gnu++17")
#set(CMAKE_CXX_FLAGS "-O3 -Wall -std=c++2a")
set(CMAKE_CXX_FLAGS "-O3 -Wall -std=c++17")

add_subdirectory(src)
