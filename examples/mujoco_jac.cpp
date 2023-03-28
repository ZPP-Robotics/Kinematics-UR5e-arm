#include <cstdio>
#include <cstring>
#include <iostream>

#include <mujoco/mujoco.h>

// MuJoCo data structures
mjModel* m = NULL;                  // MuJoCo model
mjData* d = NULL;                   // MuJoCo data

int main(int argc, const char** argv) {

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

    d = mj_makeData(m);

    for(int i = 0; i < 6; i++) {
      d->ctrl[i] = 0.0;
    }

    if(argc == 8) {
      for (int i = 0; i < 6; ++i) {
          d->qpos[i] = atof(argv[i+2]);
      }
    } else {
      for (int i = 0; i < 6; ++i) {
          d->qpos[i] = 0.0;
      }
    }
    
    mj_step(m, d);

    double jac[18] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    mj_jacBody(m, d, jac, NULL, 5);
    // mj_jacSite(m, d, jac, NULL, 0);
    // mj_jacBodyCom(m, d, jac, NULL, 4);
    for(int i = 0; i < 18; i++) {
      std::cout << jac[i] << ", ";
    }
}
