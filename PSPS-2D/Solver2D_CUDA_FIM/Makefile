CXX = nvcc -std=c++11
CXXFLAGS = -lnetcdf
OBJ  = PSPS_CUDA_FIM_2D_ComShader.o PSPS_CUDA_FIM_2D_Tools.o PSPS_CUDA_FIM_2D_Solver.o PSPS_CUDA_FIM_2D_Program.o PSPS_CUDA_FIM_2D_Main.o

%.o: %.cu $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

Solver2D_CUDA_FIM.exe: $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS)

.PHONY : clean
clean :
	rm -f *.exe $(OBJ)

