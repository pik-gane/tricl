add_executable(tricl 
    tricl.cpp 
    3rdparty/tinyexpr.c
    debugging.cpp
    constants.cpp
    probability.cpp 
    entity.cpp 
    link.cpp
    angle.cpp
    event.cpp
    config.cpp 
    io.cpp 
    graphviz.cpp
    gexf.cpp 
    init.cpp
    simulate.cpp
    finish.cpp)

target_link_libraries(tricl yaml-cpp boost_iostreams)
