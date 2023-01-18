// Copyright 2021 DeepMind Technologies Limited
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <cstdio>
#include <cstring>
#include <armadillo>
#include "analytical_ik.h"

#include <GLFW/glfw3.h>
#include <mujoco/mujoco.h>

// MuJoCo data structures
mjModel* m = NULL;                  // MuJoCo model
mjData* d = NULL;                   // MuJoCo data
mjvCamera cam;                      // abstract camera
mjvOption opt;                      // visualization options
mjvScene scn;                       // abstract scene
mjrContext con;                     // custom GPU context

// mouse interaction
bool button_left = false;
bool button_middle = false;
bool button_right =  false;
double lastx = 0;
double lasty = 0;

// keyboard callback
void keyboard(GLFWwindow* window, int key, int scancode, int act, int mods) {
    // backspace: reset simulation
    if (act==GLFW_PRESS && key==GLFW_KEY_BACKSPACE) {
        mj_resetData(m, d);
        mj_forward(m, d);
    }
}

// mouse button callback
void mouse_button(GLFWwindow* window, int button, int act, int mods) {
    // update button state
    button_left = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT)==GLFW_PRESS);
    button_middle = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE)==GLFW_PRESS);
    button_right = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT)==GLFW_PRESS);

    // update mouse position
    glfwGetCursorPos(window, &lastx, &lasty);
}

// mouse move callback
void mouse_move(GLFWwindow* window, double xpos, double ypos) {
    // no buttons down: nothing to do
    if (!button_left && !button_middle && !button_right) {
        return;
    }

    // compute mouse displacement, save
    double dx = xpos - lastx;
    double dy = ypos - lasty;
    lastx = xpos;
    lasty = ypos;

    // get current window size
    int width, height;
    glfwGetWindowSize(window, &width, &height);

    // get shift key state
    bool mod_shift = (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT)==GLFW_PRESS ||
                      glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT)==GLFW_PRESS);

    // determine action based on mouse button
    mjtMouse action;
    if (button_right) {
        action = mod_shift ? mjMOUSE_MOVE_H : mjMOUSE_MOVE_V;
    } else if (button_left) {
        action = mod_shift ? mjMOUSE_ROTATE_H : mjMOUSE_ROTATE_V;
    } else {
        action = mjMOUSE_ZOOM;
    }

    // move camera
    mjv_moveCamera(m, action, dx/height, dy/height, &scn, &cam);
}

// scroll callback
void scroll(GLFWwindow* window, double xoffset, double yoffset) {
    // emulate vertical mouse motion = 5% of window height
    mjv_moveCamera(m, mjMOUSE_ZOOM, 0, -0.05*yoffset, &scn, &cam);
}

// arm control to target
// bool arm_controller(const mjModel_* mm, mjData_* dd, std::tuple<double, double, double> target) {
//     // set target
//     mm->geom_pos[0] = std::get<0>(target);
//     mm->geom_pos[1] = std::get<1>(target);
//     mm->geom_pos[2] = std::get<2>(target);

//     double q_sols[8 * 6];
//     int num_sols = inverse_kinematics(q_sols, std::get<0>(target), std::get<1>(target), std::get<2>(target));

//     if (num_sols == 0)
//         return false;

//     printf("Inverse solutions:\n");
//     for (int i=0; i < num_sols;i++) 
//         printf("%1.6f %1.6f %1.6f %1.6f %1.6f %1.6f\n", 
//              q_sols[i*6+0], q_sols[i*6+1], q_sols[i*6+2], q_sols[i*6+3], q_sols[i*6+4], q_sols[i*6+5]);

//     for (int i = 0; i < 6; ++i) {
//         dd->ctrl[i] = q_sols[i];
//     }

//     return true;
// }

void arm_controller(const mjModel_* mm, mjData_* dd) {
    double jactmp[18] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double q[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    for (int i = 0; i < 6; ++i) {
        q[i] = dd->ctrl[i];
    }

    // double jactmp[18] = {};

    auto res = joint_jacobian(jactmp, q);
    // mj_jacSite(mm, dd, jactmp, nullptr, 0);

    // auto test = arma::mat{jactmp, 6, 3};
    // for (int i = 0; i < ; i++) {
    //     cout << 
    // }

    arma::mat pjac = pinv(arma::mat{jactmp, 6, 3}.t());
    arma::vec target{0.1, 0.5, 0.6};
    arma::vec pos = {dd->site_xpos[0], dd->site_xpos[1], dd->site_xpos[2]};

    arma::vec err = (target - pos);
    if (norm(err) < 0.05) return;

    err = err / norm(err) * 0.005;
    arma::vec diff = pjac * err;

    for (int i = 0; i < 6; ++i) {
        dd->ctrl[i] += diff[i];
    }
}


// main function
int main(int argc, const char** argv) {
    
    double q_sols[8 * 6];
    int num_sols = inverse_kinematics(q_sols, 0.1, 0.5, 0.6);

    printf("Inverse solutions:\n");
    for (int i=0; i < num_sols;i++) 
        printf("%1.6f %1.6f %1.6f %1.6f %1.6f %1.6f\n", 
             q_sols[i*6+0], q_sols[i*6+1], q_sols[i*6+2], q_sols[i*6+3], q_sols[i*6+4], q_sols[i*6+5]);


    double q[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < 6; ++i) {
        q[i] = q_sols[i];
        std::cout << q[i] << " ";
    }
    std::cout << "\n";

    auto [x, y, z] = forward_kinematics(q);

    std::cout << x << " " << y << " " << z << "\n";


    // default targe position
    double target_x = 0.0;
    double target_y = 0.0;
    double target_z = 0.0;

    // check command-line arguments
    if (argc == 5) {    
        target_x = atof(argv[2]);
        target_y = atof(argv[3]);
        target_z = atof(argv[4]);
    }
    else if (argc != 2) {
        std::printf(" USAGE: inverse_kinematics_mujoco [modelfile] [target postion <x, y, z>]\n");
        return 0;
    }

    // load and compile model
    char error[1000] = "Could nif (strtod())ot load binary model";
    if (std::strlen(argv[1])>4 && !std::strcmp(argv[1]+std::strlen(argv[1])-4, ".mjb")) {
        m = mj_loadModel(argv[1], 0);
    } else {
        m = mj_loadXML(argv[1], 0, error, 1000);
    }
    if (!m) {
        mju_error_s("Load model error: %s", error);
    }

    // make data
    d = mj_makeData(m);

    // init GLFW
    if (!glfwInit()) {
        mju_error("Could not initialize GLFW");
    }

    // create window, make OpenGL context current, request v-sync
    GLFWwindow* window = glfwCreateWindow(900, 900, "Inverse Kinematics", NULL, NULL);
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // initialize visualization data structures
    mjv_defaultCamera(&cam);
    mjv_defaultOption(&opt);
    mjv_defaultScene(&scn);
    mjr_defaultContext(&con);

    // create scene and context
    mjv_makeScene(m, &scn, 2000);
    mjr_makeContext(m, &con, mjFONTSCALE_150);

    // install GLFW mouse and keyboard callbacks
    glfwSetKeyCallback(window, keyboard);
    glfwSetCursorPosCallback(window, mouse_move);
    glfwSetMouseButtonCallback(window, mouse_button);
    glfwSetScrollCallback(window, scroll);

    // m->opt.collision = mjCOL_PAIR;
    // // count inverse kinematics for our arm
    // if (!arm_controller(m, d, {target_x, target_y, target_z})) {
    //     std::cerr << "No inverse kinematics solution!\n";
    //     return 0;
    // }
    // set target
    m->geom_pos[0] = target_x;
    m->geom_pos[1] = target_y;
    m->geom_pos[2] = target_z;

    mjcb_control = arm_controller;
    m->opt.collision = mjCOL_PAIR; // DYNAMIC - with collisions

    // run main loop, target real-time simulation and 60 fps rendering
    while (!glfwWindowShouldClose(window)) {
        //  advance interactive simulation for 1/60 sec
        //  Assuming MuJoCo can simulate faster than real-time, which it usually can,
        //  this loop will finish on time for the next frame to be rendered at 60 fps.
        //  Otherwise add a cpu timer and exit this loop when it is time to render.
        mjtNum simstart = d->time;
        while (d->time - simstart < 1.0/60.0) {
            mj_step(m, d);
        }

        // get framebuffer viewport
        mjrRect viewport = {0, 0, 0, 0};
        glfwGetFramebufferSize(window, &viewport.width, &viewport.height);

        // update scene and render
        mjv_updateScene(m, d, &opt, NULL, &cam, mjCAT_ALL, &scn);
        mjr_render(viewport, &scn, &con);

        // swap OpenGL buffers (blocking call due to v-sync)
        glfwSwapBuffers(window);

        // process pending GUI events, call GLFW callbacks
        glfwPollEvents();
    }

    //free visualization storage
    mjv_freeScene(&scn);
    mjr_freeContext(&con);

    // free MuJoCo model and data
    mj_deleteData(d);
    mj_deleteModel(m);

    // terminate GLFW (crashes with Linux NVidia drivers)
#if defined(__APPLE__) || defined(_WIN32)
    glfwTerminate();
#endif

    return 1;
}
