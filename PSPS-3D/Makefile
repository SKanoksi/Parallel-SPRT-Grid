
all: cuda opengl

cuda:
	cd ./Solver3D_CUDA_BF && $(MAKE)
	cd ./Solver3D_CUDA_ExtendedBF && $(MAKE)
	cd ./Solver3D_CUDA_FIM && $(MAKE)

opengl:
	cd ./Solver3D_OpenGL_BF && $(MAKE)
	cd ./Solver3D_OpenGL_ExtendedBF && $(MAKE)
	cd ./Solver3D_OpenGL_FIM && $(MAKE)

.PHONY : clean
clean :
	rm -f ./Solver3D*/*.exe ./Solver3D*/*.o

