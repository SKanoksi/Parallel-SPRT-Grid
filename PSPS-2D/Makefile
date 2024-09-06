
all: cpu cuda opengl

cpu:
	cd ./Solver2D_CPU_DK && $(MAKE)
	cd ./Solver2D_CPU_FMM && $(MAKE)

cuda:
	cd ./Solver2D_CUDA_BF && $(MAKE)
	cd ./Solver2D_CUDA_ExtendedBF && $(MAKE)
	cd ./Solver2D_CUDA_FIM && $(MAKE)

opengl:
	cd ./Solver2D_OpenGL_BF && $(MAKE)
	cd ./Solver2D_OpenGL_ExtendedBF && $(MAKE)
	cd ./Solver2D_OpenGL_FIM && $(MAKE)

.PHONY : clean
clean :
	rm -f ./Solver2D*/*.exe ./Solver2D*/*.o

