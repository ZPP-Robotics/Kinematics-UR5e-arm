#include "gtest/gtest.h"

#include "analytical_ik.h"

const double THRESHOLD = 1e-3;

TEST(AnalyticalIKTest, forwardKinematics) {
    double q[6] = {0, 0, 0, 0, 0, 0};
    auto [x, y, z] = forward_kinematics(q);
    EXPECT_NEAR(x, -0.8172, THRESHOLD);
    EXPECT_NEAR(y, -0.2329, THRESHOLD);
    EXPECT_NEAR(z, 0.0628, THRESHOLD);
}

TEST(AnalyticalIKTest, forwardKinematics6Back) {
    double q[6] = {0, 0, 0, 0, 0, 0};
    auto [x, y, z] = forward_kinematics_6_back(q);
    EXPECT_NEAR(x, -0.8172, THRESHOLD);
    EXPECT_NEAR(y, -0.0835, THRESHOLD);
    EXPECT_NEAR(z, 0.0628, THRESHOLD);
}

TEST(AnalyticalIKTest, forwardKinematics4) {
    double q[6] = {0, 0, 0, 0, 0, 0};
    auto [x, y, z] = forward_kinematics_4(q);
    EXPECT_NEAR(x, -0.8172, THRESHOLD);
    EXPECT_NEAR(y, -0.1333, THRESHOLD);
    EXPECT_NEAR(z, 0.1625, THRESHOLD);
}

TEST(AnalyticalIKTest, forwardKinematics3) {
    double q[6] = {0, 0, 0, 0, 0, 0};
    auto [x, y, z] = forward_kinematics_3(q);
    EXPECT_NEAR(x, -0.8172, THRESHOLD);
    EXPECT_NEAR(y, 0., THRESHOLD);
    EXPECT_NEAR(z, 0.1625, THRESHOLD);
}

TEST(AnalyticalIKTest, forwardKinematicsElbowJoint) {
    double q[6] = {0, 0, 0, 0, 0, 0};
    auto [x, y, z] = forward_kinematics_elbow_joint(q);
    EXPECT_NEAR(x, -0.425, THRESHOLD);
    EXPECT_NEAR(y, 0., THRESHOLD);
    EXPECT_NEAR(z, 0.1625, THRESHOLD);
}