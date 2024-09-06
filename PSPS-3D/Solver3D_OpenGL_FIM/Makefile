CXX = g++ -std=c++11
CXXFLAGS = -lGLEW -lglfw -lGL -lX11 -lXi -lXrandr -lXxf86vm -lXinerama -lXcursor -lrt -lm -pthread -ldl -lnetcdf
OBJ  = PSPS_OpenGL_FIM_3D_ComShader.o PSPS_OpenGL_FIM_3D_Tools.o PSPS_OpenGL_FIM_3D_Solver.o PSPS_OpenGL_FIM_3D_Program.o PSPS_OpenGL_FIM_3D_Main.o

%.o: %.c $(DEPS)
	$(CXX) -o $@ $< $(CXXFLAGS)

Solver3D_OpenGL_FIM.exe: $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS)

.PHONY : clean
clean :
	rm -f *.exe $(OBJ)
