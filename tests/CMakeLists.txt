include(FetchContent)

include_directories(${Kinematics_UR5e_arm_SOURCE_DIR}/src)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        origin/main
)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)

add_executable(tests test.cpp)
target_link_libraries(tests analytical_ik gtest gmock gtest_main)