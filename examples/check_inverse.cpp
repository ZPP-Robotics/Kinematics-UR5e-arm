#include <iostream>

#include "analytical_ik.h"

int main()
{
  // double qrel[6]{1.81, 0.99, 2.0, 1.57, -1.0, -1.0};
  // double qrel[6]{-1.0509, -2.10413, 1.73413, -2.77159, 1.0509, -3.14159};
  // double qrel[6]{0.5341, -0.6780, -0.3220, -0.0283, 0.0245, 0};
  // double qrel[6]{-3.141562,0,0,0,0,0};
  // auto [x, y, z] = forward_kinematics(qrel);

  double q[6 * 8]{0};
  double x, y, z;
  x = 0.4;
  y = 0.3;
  z = 0.3;

  int num_sols = inverse_kinematics(q, x, y, z);
  std::cout << "num_sols: " << num_sols << "\n";
  for (int i = 0; i < num_sols; ++i)
  {
    std::cout << "q" << i << ": " << q[6 * i] << " " << q[6 * i + 1] << " " << q[6 * i + 2] << " " << q[6 * i + 3] << " " << q[6 * i + 4] << " " << q[6 * i + 5] << "\n";
  }
//    2.6     
//   -2.2     
//    1.7     
//   -2.2     
//    2.5     
//   -2.4   
}