CXX = nvcc -std=c++11
CXXFLAGS = -lnetcdf
OBJ  = PSPS_CUDA_ExtendedBF_3D_ComShader.o PSPS_CUDA_ExtendedBF_3D_Tools.o PSPS_CUDA_ExtendedBF_3D_Solver.o PSPS_CUDA_ExtendedBF_3D_Program.o PSPS_CUDA_ExtendedBF_3D_Main.o

%.o: %.cu $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

Solver3D_CUDA_ExtendedBF.exe: $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS)

.PHONY : clean
clean :
	rm -f *.exe $(OBJ)
