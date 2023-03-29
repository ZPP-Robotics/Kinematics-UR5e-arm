#include "gtest/gtest.h"

#include "analytical_ik.h"

const double THRESHOLD_FK = 1e-3;
const double THRESHOLD_JAC = 1e-2;

TEST(AnalyticalIKTest, forwardKinematics) {
    double q[6] = {0, 0, 0, 0, 0, 0};
    auto [x, y, z] = forward_kinematics(q);
    EXPECT_NEAR(x, -0.8172, THRESHOLD_FK);
    EXPECT_NEAR(y, -0.2329, THRESHOLD_FK);
    EXPECT_NEAR(z, 0.0628, THRESHOLD_FK);
}

TEST(AnalyticalIKTest, forwardKinematics6Back) {
    double q[6] = {0, 0, 0, 0, 0, 0};
    auto [x, y, z] = forward_kinematics_6_back(q);
    EXPECT_NEAR(x, -0.8172, THRESHOLD_FK);
    EXPECT_NEAR(y, -0.0835, THRESHOLD_FK);
    EXPECT_NEAR(z, 0.0628, THRESHOLD_FK);
}

TEST(AnalyticalIKTest, forwardKinematics4) {
    {
        double q[6] = {0, 0, 0, 0, 0, 0};
        auto [x, y, z] = forward_kinematics_4(q);
        EXPECT_NEAR(x, -0.8172, THRESHOLD_FK);
        EXPECT_NEAR(y, -0.1333, THRESHOLD_FK);
        EXPECT_NEAR(z, 0.1625, THRESHOLD_FK);
    }
    {
        double q[6] = {-1.51, -1.76, 1.04, -2.14, 1.82, -1.45};
        auto [x, y, z] = forward_kinematics_4(q);
        EXPECT_NEAR(x, -0.1461, THRESHOLD_FK);
        EXPECT_NEAR(y, 0.2064, THRESHOLD_FK);
        EXPECT_NEAR(z, 0.8385, THRESHOLD_FK);
    }
    {
        double q[6] = {2.51, -1.19, -1.89, -2.14, -3.39, 3.46};
        auto [x, y, z] = forward_kinematics_4(q);
        EXPECT_NEAR(x, -0.1097, THRESHOLD_FK);
        EXPECT_NEAR(y, 0.2454, THRESHOLD_FK);
        EXPECT_NEAR(z, 0.5811, THRESHOLD_FK);
    }
    
}

TEST(AnalyticalIKTest, forwardKinematics3) {
    {
        double q[6] = {0, 0, 0, 0, 0, 0};
        auto [x, y, z] = forward_kinematics_3(q);
        EXPECT_NEAR(x, -0.8172, THRESHOLD_FK);
        EXPECT_NEAR(y, 0., THRESHOLD_FK);
        EXPECT_NEAR(z, 0.1625, THRESHOLD_FK);
    }
    {
        double q[6] = {-1.51, -1.76, 1.04, -2.14, 1.82, -1.45};
        auto [x, y, z] = forward_kinematics_3(q);
        EXPECT_NEAR(x, -0.0130, THRESHOLD_FK);
        EXPECT_NEAR(y, 0.2145, THRESHOLD_FK);
        EXPECT_NEAR(z, 0.8385, THRESHOLD_FK);
    }
    {
        double q[6] = {2.51, -1.19, -1.89, -2.14, -3.39, 3.46};
        auto [x, y, z] = forward_kinematics_3(q);
        EXPECT_NEAR(x, -0.1884, THRESHOLD_FK);
        EXPECT_NEAR(y, 0.1378, THRESHOLD_FK);
        EXPECT_NEAR(z, 0.5811, THRESHOLD_FK);
    }
}

TEST(AnalyticalIKTest, forwardKinematicsElbowJoint) {
    {
        double q[6] = {0, 0, 0, 0, 0, 0};
        auto [x, y, z] = forward_kinematics_elbow_joint(q);
        EXPECT_NEAR(x, -0.425, THRESHOLD_FK);
        EXPECT_NEAR(y, 0., THRESHOLD_FK);
        EXPECT_NEAR(z, 0.1625, THRESHOLD_FK);
    }
    {
        double q[6] = {-1.51, -1.76, 1.04, -2.14, 1.82, -1.45};
        auto [x, y, z] = forward_kinematics_elbow_joint(q);
        EXPECT_NEAR(x, 0.0048, THRESHOLD_FK);
        EXPECT_NEAR(y, -0.0797, THRESHOLD_FK);
        EXPECT_NEAR(z, 0.5799, THRESHOLD_FK);
    }
    {
        double q[6] = {2.51, -1.19, -1.89, -2.14, -3.39, 3.46};
        auto [x, y, z] = forward_kinematics_elbow_joint(q);
        EXPECT_NEAR(x, 0.1274, THRESHOLD_FK);
        EXPECT_NEAR(y, -0.0932, THRESHOLD_FK);
        EXPECT_NEAR(z, 0.5570, THRESHOLD_FK);
    }
}

TEST(AnalyticalIKTest, joint_jacobian) {
    double q[6] = {0, 0, 0, 0, 0, 0};
    double jac[18];
    joint_jacobian(jac, q);
    double jac_gt[18] = {0.235301, 0.0999999, 0.0999999, 0.0999999, -0.0999999, 0, -0.816626, 0.000159265, 0.000159265, 0.000159265, -0.000159265, 0, 0, -0.817, -0.392, 0, 0, 0};
    for (int i = 0; i < 18; i++) {
        EXPECT_NEAR(jac[i], jac_gt[i], THRESHOLD_JAC);
    }
}

TEST(AnalyticalIKTest, joint_jacobian_4) {
    {
        double q[6] = {0, 0, 0, 0, 0, 0};
        double jac[18];
        joint_jacobian_4(jac, q);
        double jac_gt[18] = {0.135301, 0, 0, 0, 0, 0, -0.816786, 0, 0, 0, 0, 0, 0, -0.817, -0.392, 0, 0, 0};
        for (int i = 0; i < 18; i++) {
            EXPECT_NEAR(jac[i], jac_gt[i], THRESHOLD_JAC);
        }
    }
    {
        double q[6] = {-1.51, -1.76, 1.04, -2.14, 1.82, -1.45};
        double jac[18];
        joint_jacobian_4(jac, q);
        double jac_gt[18] = {-0.206003, -0.042141, -0.0161158, 0, 0, 0, -0.14713, 0.674579, 0.257976, 0, 0, 0, 0, -0.214775, -0.294708, 0, 0, 0};
        for (int i = 0; i < 18; i++) {
            EXPECT_NEAR(jac[i], jac_gt[i], THRESHOLD_JAC);
        }
    }
    {
        double q[6] = {2.51, -1.19, -1.89, -2.14, -3.39, 3.46};
        double jac[18];
        joint_jacobian_4(jac, q);
        double jac_gt[18] = {-0.245724, 0.33831, 0.0194969, 0, 0, 0, -0.109568, -0.246667, -0.0142155, 0, 0, 0, 0, 0.233301, 0.391257, 0, 0, 0};
        for (int i = 0; i < 18; i++) {
            EXPECT_NEAR(jac[i], jac_gt[i], THRESHOLD_JAC);
        }
    }
    {
        double q[6] = {1.81, 0.99, 2.0, 1.57, -1.0, -1.0};
        double jac[18];
        joint_jacobian_4(jac, q);
        double jac_gt[18] = {-0.181815, -0.0988502, -0.0141171, 0, 0, 0, 0.0933343, 0.402549, 0.0574891, 2.77556e-17, 0, 0, 0, 0.154311, 0.387504, 0, 0, 0};
        for (int i = 0; i < 18; i++) {
            EXPECT_NEAR(jac[i], jac_gt[i], THRESHOLD_JAC);
        }
    }
}

TEST(AnalyticalIKTest, joint_jacobian_3_Test) {
    {
        double q[6] = {0, 0, 0, 0, 0, 0};
        double jac[18];
        joint_jacobian_3(jac, q);
        double jac_gt[18] = {0.00830119, 0, 0, 0, 0, 0, -0.816988, 0, 0, 0, 0, 0, 0, -0.817, -0.392, 0, 0, 0};
        for (int i = 0; i < 18; i++) {
            EXPECT_NEAR(jac[i], jac_gt[i], THRESHOLD_JAC);
        }
    }
    {
        double q[6] = {-1.51, -1.76, 1.04, -2.14, 1.82, -1.45};
        double jac[18];
        joint_jacobian_3(jac, q);
        double jac_gt[18] = {-0.213921, -0.042141, -0.0161158, 0, 0, 0, -0.0203773, 0.674579, 0.257976, 0, 0, 0, 0, -0.214775, -0.294708, 0, 0, 0};
        for (int i = 0; i < 18; i++) {
            EXPECT_NEAR(jac[i], jac_gt[i], THRESHOLD_JAC);
        }
    }
    {
        double q[6] = {2.51, -1.19, -1.89, -2.14, -3.39, 3.46};
        double jac[18];
        joint_jacobian_3(jac, q);
        double jac_gt[18] = {-0.143104, 0.33831, 0.0194969, 0, 0, 0, -0.18439, -0.246667, -0.0142155, 0, 0, 0, 0, 0.233301, 0.391257, 0, 0, 0};
        for (int i = 0; i < 18; i++) {
            EXPECT_NEAR(jac[i], jac_gt[i], THRESHOLD_JAC);
        }
    }
    {
        double q[6] = {1.81, 0.99, 2.0, 1.57, -1.0, -1.0};
        double jac[18];
        joint_jacobian_3(jac, q);
        double jac_gt[18] = {-0.151528, -0.0988502, -0.0141171, 0, 0, 0, -0.0300015, 0.402549, 0.0574891, 0, 0, 0, 0, 0.154311, 0.387504, 0, 0, 0};
        for (int i = 0; i < 18; i++) {
            EXPECT_NEAR(jac[i], jac_gt[i], THRESHOLD_JAC);
        }
    }
}