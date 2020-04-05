add_executable(tricl 
    tricl.cpp 
    debugging.cpp
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
    simulate.cpp
    finish.cpp)

target_link_libraries(tricl yaml-cpp)
