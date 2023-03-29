#ifndef ANALYTICAL_IK_H
#define ANALYTICAL_IK_H

#include <tuple>

std::tuple<double, double, double> forward_kinematics(double *q);

std::tuple<double, double, double> forward_kinematics_6_back(double *q);

std::tuple<double, double, double> forward_kinematics_4(double *q);

std::tuple<double, double, double> forward_kinematics_3(double *q);

std::tuple<double, double, double> forward_kinematics_elbow_joint(double *q);

int inverse_kinematics_2PI(double *q_sols, double x, double y, double z);

int inverse_kinematics(double *q_sols, double x, double y, double z);

void joint_jacobian(double *jacobian, double *q);

void joint_jacobian_6_back(double *jacobian, double *q);

void joint_jacobian_4(double *jacobian, double *q);

void joint_jacobian_3(double *jacobian, double *q);

void jacobian_elbow_joint(double *jacobian, double *q);

#endif  // ANALYTICAL_IK_H