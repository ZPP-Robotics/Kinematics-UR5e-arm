#include <iostream>

#include "analytical_ik.h"

int main() {
  // double q[6] = {0, 0, 0, 0, 0, 0};
  double q[6] = {2.26, -0.754, 0.723, 2.2, -1.45, 2.2};
  auto forwards = {forward_kinematics, forward_kinematics_gripper,
                   forward_kinematics_6_back, forward_kinematics_4,
                   forward_kinematics_3, forward_kinematics_elbow_joint};
  for (auto f : forwards) {
    auto [x, y, z] = f(q);
    std::cout << x << " " << y << " " << z << std::endl;
  }
}