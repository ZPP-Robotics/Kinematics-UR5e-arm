set(CMAKE_CXX_STANDARD 17)

add_library(analytical_ik analytical_ik.h analytical_ik.cpp)
set_target_properties(analytical_ik PROPERTIES LINKER_LANGUAGE CXX)
