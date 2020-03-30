cmake_minimum_required (VERSION 2.6)

project (tricl)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-O3 -Wall -std=gnu++17")

add_executable(tricl 
    tricl.cpp 
    constants.cpp
    probability.cpp 
    entity.cpp 
    link.cpp
    angle.cpp
    event.cpp
    config.cpp 
    io.cpp 
    gexf.cpp 
    init.cpp
    finish.cpp)

target_link_libraries(tricl yaml-cpp)
