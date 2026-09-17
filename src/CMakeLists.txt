add_executable(tricl
    tricl.cpp
    3rdparty/tinyexpr.c
    3rdparty/gzip.cpp
    3rdparty/zlib.cpp
    debugging.cpp
    constants.cpp
    probability.cpp
    entity.cpp
    link.cpp
#    angle.cpp  # currently empty
    event.cpp
    config.cpp
    io.cpp
    graphviz.cpp
    gexf.cpp
    init.cpp
    simulate.cpp
    replay.cpp
    finish.cpp)

# the bundled tinyexpr triggers a false-positive array-bounds warning in newer gcc versions:
set_source_files_properties(3rdparty/tinyexpr.c PROPERTIES COMPILE_OPTIONS "-Wno-array-bounds")

# yaml-cpp exports different target names depending on its version:
if(TARGET yaml-cpp::yaml-cpp)
    set(TRICL_YAML_CPP_TARGET yaml-cpp::yaml-cpp)
else()
    set(TRICL_YAML_CPP_TARGET yaml-cpp)
endif()

target_include_directories(tricl PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(tricl ${TRICL_YAML_CPP_TARGET} ZLIB::ZLIB)
