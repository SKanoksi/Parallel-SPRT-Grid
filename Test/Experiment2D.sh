#!bin/bash
set -e # Exit if any error
# Need to copy all 2D programs (*.exe) and all required shaders (*.comSh), eg. PSPS_FIM_Solver2D.comsh,
# into this folder, before running this script !!!

Header="R2048"
TestModel="Rgrad2D_2048.nc"
ResX="2048"
ResY="2048"
RadX="4"
RadY="4"
BlockX="8"
BlockY="8"
SourceX="1024"
SourceY="1"
MaxLoop="499999"

# Select programs (Y = Yes, otherwise = No)
DK="Y"
FMM="Y"
CUDA_BF="N"
CUDA_FIM="N"
CUDA_ExBF="N"
OpenGL_BF="Y"
OpenGL_FIM="Y"
OpenGL_ExBF="Y"

# ======================= Begin =======================

# ------------------- CPU -------------------

# DK
if [ "$DK" == "Y" ]
then
	echo "Running : DK"
		Model=$Header"_CPU_DK.nc"
		Log=$Header"_CPU_DK.txt"
		cp ./$TestModel $Model
		./Solver2D_CPU_DK.exe $Model $RadX $RadY $SourceX $SourceY > $Log
	echo "Finished: DK"
fi

# FMM
if [ "$FMM" == "Y" ]
then
	echo "Running : FMM"
		Model=$Header"_CPU_FMM.nc"
		Log=$Header"_CPU_FMM.txt"
		cp ./$TestModel $Model
		./Solver2D_CPU_FMM.exe $Model $SourceX $SourceY > $Log
	echo "Finished: FMM"
fi


# ------------------- CUDA -------------------

# CUDA_BF
if [ "$CUDA_BF" == "Y" ]
then
	echo "Running : CUDA_BF"
		Model=$Header"_CUDA_BF.nc"
		Log=$Header"_CUDA_BF.txt"
		cp ./$TestModel $Model
		./Solver2D_CUDA_BF.exe $Model $RadX $RadY $ResX $ResY $SourceX $SourceY $MaxLoop > $Log
	echo "Finished: CUDA_BF"
fi

# CUDA_FIM
if [ "$CUDA_FIM" == "Y" ]
then
	echo "Running : CUDA_FIM"
		Model=$Header"_CUDA_FIM.nc"
		Log=$Header"_CUDA_FIM.txt"
		cp ./$TestModel $Model
		./Solver2D_CUDA_FIM.exe $Model $BlockX $BlockY $ResX $ResY $SourceX $SourceY $MaxLoop > $Log
	echo "Finished: CUDA_FIM"
fi

# CUDA_ExBF
if [ "$CUDA_ExBF" == "Y" ]
then
	echo "Running : CUDA_ExBF"
		Model=$Header"_CUDA_ExBF.nc"
		Log=$Header"_CUDA_ExBF.txt"
		cp ./$TestModel $Model
		./Solver2D_CUDA_ExtendedBF.exe $Model $ResX $ResY $SourceX $SourceY $MaxLoop > $Log
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
		./Solver2D_OpenGL_BF.exe $Model $RadX $RadY $ResX $ResY $SourceX $SourceY $MaxLoop > $Log
	echo "Finished: OpenGL_BF"
fi

# Fast Iterative method

if [ "$OpenGL_FIM" == "Y" ]
then
	echo "Running : OpenGL_FIM"
		Model=$Header"_OpenGL_FIM.nc"
		Log=$Header"_OpenGL_FIM.txt"
		cp ./$TestModel $Model
		./Solver2D_OpenGL_FIM.exe $Model $BlockX $BlockY $ResX $ResY $SourceX $SourceY $MaxLoop > $Log
	echo "Finished: OpenGL_FIM"
fi

# Extended Bellman-Ford versions

if [ "$OpenGL_ExBF" == "Y" ]
then
	echo "Running : OpenGL_ExBF"
		Model=$Header"_OpenGL_ExBF.nc"
		Log=$Header"_OpenGL_ExBF.txt"
		cp ./$TestModel $Model
		./Solver2D_OpenGL_ExtendedBF.exe $Model $ResX $ResY $SourceX $SourceY $MaxLoop > $Log
	echo "Finished: OpenGL_ExBF"
fi


: <<'COMMENT'
# use dos2unix, if it is coded in window.
COMMENT



