include(FetchContent)

FetchContent_Declare(
  mujoco
  GIT_REPOSITORY https://github.com/deepmind/mujoco.git
  GIT_TAG 2.3.0
)

set(MUJOCO_BUILD_EXAMPLES OFF)
set(MUJOCO_BUILD_TESTS OFF)
set(MUJOCO_TEST_PYTHON_UTIL OFF)

FetchContent_Declare(
        glfw
        GIT_REPOSITORY https://github.com/glfw/glfw
)

FetchContent_GetProperties(glfw)
if(NOT glfw_POPULATED)
    FetchContent_Populate(glfw)

    set(GLFW_BUILD_EXAMPLES OFF CACHE INTERNAL "Build the GLFW example programs")
    set(GLFW_BUILD_TESTS OFF CACHE INTERNAL "Build the GLFW test programs")
    set(GLFW_BUILD_DOCS OFF CACHE INTERNAL "Build the GLFW documentation")
    set(GLFW_INSTALL OFF CACHE INTERNAL "Generate installation target")

    add_subdirectory(${glfw_SOURCE_DIR} ${glfw_BINARY_DIR})
endif()

FetchContent_Declare(
    armadillo
    GIT_REPOSITORY https://gitlab.com/conradsnicta/armadillo-code.git
    GIT_SHALLOW TRUE
)


FetchContent_MakeAvailable(mujoco)
FetchContent_MakeAvailable(GLFW)

include_directories(${Kinematics_UR5e_arm_SOURCE_DIR}/src)

add_executable(inverse_kinematics_mujoco main.cpp)
target_link_libraries(inverse_kinematics_mujoco analytical_ik)
target_link_libraries(inverse_kinematics_mujoco glfw)
target_link_libraries(inverse_kinematics_mujoco mujoco)
target_link_libraries(inverse_kinematics_mujoco armadillo)

add_executable(check_jacob check_jacob.cpp)
target_link_libraries(check_jacob analytical_ik)
target_link_libraries(check_jacob armadillo)

add_executable(mujoco_jac mujoco_jac.cpp)
target_link_libraries(mujoco_jac mujoco)


