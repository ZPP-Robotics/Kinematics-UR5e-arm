#ifndef ANALYTICAL_IK_H
#define ANALYTICAL_IK_H

#include <tuple>

std::tuple<double, double, double> forward_kinematics(double *q);

std::tuple<double, double, double> forward_kinematics_elbow_joint(double *q);

int inverse_kinematics_2PI(double *q_sols, double x, double y, double z);

int inverse_kinematics(double *q_sols, double x, double y, double z);

void joint_jacobian(double *jacobian, double *q);

void jacobian_elbow_joint(double *jacobian, double *q);

#endif  // ANALYTICAL_IK_H