#!bin/bash
set -e # Exit if any error
# Need to copy all 3D programs (*.exe) and all required shaders (*.comSh), eg. PSPS_FIM_Solver3D.comsh,
# into this folder, before running this script !!!

Header="Z128"
TestModel="Zgrad3D_128.nc"
ResX="128"
ResY="128"
ResZ="128"
RadX="2"
RadY="2"
RadZ="2"
BlockX="4"
BlockY="4"
BlockZ="4"
SourceX="1"
SourceY="1"
SourceZ="1"
MaxLoop="29999"

# Select programs (Y = Yes, otherwise = No)
CUDA_BF="N"
CUDA_FIM="N"
CUDA_ExBF="N"
OpenGL_BF="Y"
OpenGL_FIM="Y"
OpenGL_ExBF="Y"
FIM_WonKi="N"

# ======================= Begin =======================

# ------------------- CUDA -------------------

# CUDA_BF
if [ "$CUDA_BF" == "Y" ]
then
	echo "Running : CUDA_BF"
		Model=$Header"_CUDA_BF.nc"
		Log=$Header"_CUDA_BF.txt"
		cp ./$TestModel $Model
		./Solver3D_CUDA_BF.exe $Model $RadX $RadY $RadZ $ResX $ResY $ResZ $SourceX $SourceY $SourceZ $MaxLoop > $Log
	echo "Finished: CUDA_BF"
fi

# CUDA_FIM
if [ "$CUDA_FIM" == "Y" ]
then
	echo "Running : CUDA_FIM"
		Model=$Header"_CUDA_FIM.nc"
		Log=$Header"_CUDA_FIM.txt"
		cp ./$TestModel $Model
		./Solver3D_CUDA_FIM.exe $Model $BlockX $BlockY $BlockZ $ResX $ResY $ResZ $SourceX $SourceY $SourceZ $MaxLoop > $Log
	echo "Finished: CUDA_FIM"
fi

# CUDA_ExBF
if [ "$CUDA_ExBF" == "Y" ]
then
	echo "Running : CUDA_ExBF"
		Model=$Header"_CUDA_ExBF.nc"
		Log=$Header"_CUDA_ExBF.txt"
		cp ./$TestModel $Model
		./Solver3D_CUDA_ExtendedBF.exe $Model $ResX $ResY $ResZ $SourceX $SourceY $SourceZ $MaxLoop > $Log
	echo "Finished: CUDA_ExBF"
fi

# ------------------- OpenGL -------------------

# Bellman-Ford versions

# OpenGL_BF
if [ "$OpenGL_BF" == "Y" ]
then
	echo "Running : OpenGL_BF"
		Model=$Header"_OpenGL_BF.nc"
		Log=$Header"_OpenGL_BF.txt"
		cp ./$TestModel $Model
		./Solver3D_OpenGL_BF.exe $Model $RadX $RadY $RadZ $ResX $ResY $ResZ $SourceX $SourceY $SourceZ $MaxLoop > $Log
	echo "Finished: OpenGL_BF"
fi

# Fast Iterative method

if [ "$OpenGL_FIM" == "Y" ]
then
	echo "Running : OpenGL_FIM"
		Model=$Header"_OpenGL_FIM.nc"
		Log=$Header"_OpenGL_FIM.txt"
		cp ./$TestModel $Model
		./Solver3D_OpenGL_FIM.exe $Model $BlockX $BlockY $BlockZ $ResX $ResY $ResZ $SourceX $SourceY $SourceZ $MaxLoop > $Log
	echo "Finished: OpenGL_FIM"
fi

# Extended Bellman-Ford versions

if [ "$OpenGL_ExBF" == "Y" ]
then
	echo "Running : OpenGL_ExBF"
		Model=$Header"_OpenGL_ExBF.nc"
		Log=$Header"_OpenGL_ExBF.txt"
		cp ./$TestModel $Model
		./Solver3D_OpenGL_ExtendedBF.exe $Model $ResX $ResY $ResZ $SourceX $SourceY $SourceZ $MaxLoop > $Log
	echo "Finished: OpenGL_ExBF"
fi

: <<'COMMENT'
# use dos2unix, if it is coded in window.
COMMENT



