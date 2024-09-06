CXX = g++ -std=c++11
CXXFLAGS = -lnetcdf
OBJ  = PSPS_CPU_DK_2D_ComShader.o PSPS_CPU_DK_2D_Tools.o PSPS_CPU_DK_2D_Solver.o PSPS_CPU_DK_2D_Program.o PSPS_CPU_DK_2D_Main.o

%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

Solver2D_CPU_DK.exe: $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS)

.PHONY : clean
clean :
	rm -f *.exe $(OBJ)

