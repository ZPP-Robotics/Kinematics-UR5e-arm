#include "analytical_ik.h"
#include <armadillo>
#include <iostream>

auto eps = 0.02;

int main() {
  double q_start[6] = {0.44, 1.0, 0.2, 2.5, 1.4, 1.7};
  double q_diff[6] = {eps, 0.0, 0.0, 0.0, 0.0, 0.0};
  double q_end[6];
  for (int i = 0; i < 6; ++i) {
    q_end[i] = q_start[i] + q_diff[i];
  }

  double jac[18];
  joint_jacobian(jac, q_start);

  arma::mat jac_diff = arma::mat{jac, 3, 6} * arma::mat{q_diff, 6, 1};

  auto [x_gt, y_gt, z_gt] = forward_kinematics(q_end);
  auto [x, y, z] = forward_kinematics(q_start);

  std::cout << "gt xyz: " << x_gt << " " << y_gt << " " << z_gt << "\n";
  std::cout << "xyz: " << x + jac_diff[0] << " " << y+ jac_diff[0] << " " << z+ jac_diff[0] << "\n";
  std::cout << "start xyz: " << x << " " << y << " " << z << "\n";
  std::cout << "jac_diff: " << jac_diff[0] << " " << jac_diff[1] << " " << jac_diff[2] << "\n";

}