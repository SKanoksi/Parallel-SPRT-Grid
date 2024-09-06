CXX = g++
CXXFLAGS = -lnetcdf
OBJ  = Creator.o

%.o: %.c $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

Creator.exe: $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS)

.PHONY : clean
clean :
	rm -f *.exe $(OBJ)

