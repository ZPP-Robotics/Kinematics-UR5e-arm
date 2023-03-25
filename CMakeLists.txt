cmake_minimum_required(VERSION 3.18)
project(Kinematics_UR5e_arm)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(src)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  option(Kinematics_UR5e_arm_BUILD_EXAMPLES "Build examples" ON)

  if(Kinematics_UR5e_arm_BUILD_EXAMPLES)
    add_subdirectory(examples)  
  endif()

endif()
