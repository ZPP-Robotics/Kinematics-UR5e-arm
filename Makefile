COMMON=-O2 -I ../../../mujoco/mujoco-2.3.0/include -I ../armadillo-11.4.2/include/ -L ../../../mujoco/mujoco-2.3.0/lib -L ../armadillo-11.4.2  -std=c++17 -pthread -Wl,-no-as-needed -Wl,-rpath,'$$ORIGIN'/../../../../mujoco/mujoco-2.3.0/lib

.PHONY: all
all:
	$(CXX) $(COMMON) -std=c++17 -c src/analytical_ik.cpp 
	$(CXX) $(COMMON) -std=c++17 src/main.cpp analytical_ik.o -lmujoco -lglfw -larmadillo -lopenblas -o bin/inverse_kinematics_mujoco
	rm analytical_ik.o