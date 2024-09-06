/*************************************
     Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#ifndef PSPS_CPU_FMM_2D_SOLVER_H
#define PSPS_CPU_FMM_2D_SOLVER_H

// Standard library
#include <queue> // Prior queue
#include <limits>
#include <cmath>
#include <chrono>

#include "PSPS_CPU_FMM_2D_ComShader.h"

#define DIM 2

// Special containers
typedef std::pair<float,int> Vertex ;

struct LineSegment{
    std::vector<int> shift ;
    std::vector<float> length ;
};

class PSPS_CPU_FMM_2D_Solver : public PSPS_CPU_FMM_2D_ComShader
{
public:
    explicit PSPS_CPU_FMM_2D_Solver();
    ~PSPS_CPU_FMM_2D_Solver();

   bool Compute();

private:
   // Some useful numbers
   int length[DIM], radius[DIM], front, side ;

   // Some useful function
   inline bool isShortest(int x);
   inline float min(float a, float b);

   // All real main functions
   float invHxSquare, invHySquare, sumInv ;
   float newCandidateOfTraveltime(const int ptr);

	// Others
	void addDummyVertices();
	void removeDummyVertices();

};

#endif // PSPS_CPU_FMM_2D_SOLVER_H
