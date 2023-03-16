# Kinematics for a UR5e arm

## Run simulation
 - <code>
    mkdir budild </br>
    cd budild </br>
    cmake .. </br>
  </code>

 - execute 
    <code>
    ./examples/inverse_kinematics_mujoco [ur5e_model_file_path] [target coordinates in R^3]
    </code>


## Overview
This project contains a forward kinematics and an inverse kinematics solution for a UR5e.

The library  `analytical_ik.h` contains three functions:
 
 - std::tuple<double, double, double> forward_kinematics(double *q);

 - int inverse_kinematics_2PI(double *q_sols, double x, double y, double z);

 - int inverse_kinematics(double *q_sols, double x, double y, double z);

The `forward_kinematics()` function returns the coordinates of a point in the R^3 space of the position of the end-effector or the UR5e arm, for a given configuration of it's joint positions.

The other functions count the inverse kinematics for the UR5e arm - given a point in the R^3 space, the functions will return (in the parameter q_sols) the joint positions which will position the arm, so that it's end-effector is located in that point. 

Both functions return an integer indicating the number of inverse kinematics solutions.


The `inverse_kinematics()` function counts the joint positions in the interval `[-PI ; PI]` and the `inverse_kinematics_2PI()` function counts the joint positions in the interval `[0 ; 2 * PI]`


We base our analytical solution of `forward` and `inverse` functions on Kelsey P. Hawkins' paper "Analytic Inverse Kinematics for the Universal Robots UR-5/UR-10 Arms" and his implementation of it https://github.com/ros-industrial/universal_robot/blob/kinetic-devel/ur_kinematics/src/ur_kin.cpp.
