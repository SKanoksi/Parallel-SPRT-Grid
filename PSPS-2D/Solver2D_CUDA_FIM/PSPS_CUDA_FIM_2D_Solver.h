/*************************************
     Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#ifndef PSPS_CUDA_FIM_2D_SOLVER_H
#define PSPS_CUDA_FIM_2D_SOLVER_H

// Standard library
#include <limits>
#include <cmath>
#include <chrono>

#include "PSPS_CUDA_FIM_2D_ComShader.h"

#define DIM 2

class PSPS_CUDA_FIM_2D_Solver : public PSPS_CUDA_FIM_2D_ComShader
{
public:
    explicit PSPS_CUDA_FIM_2D_Solver();
    ~PSPS_CUDA_FIM_2D_Solver();

   bool Compute();
   
private:
   // Some useful numbers
   unsigned int length[DIM], front, side ;

	// Others
	void addDummyVertices();
	void removeDummyVertices();

};

#endif // PSPS_CUDA_FIM_2D_SOLVER_H
