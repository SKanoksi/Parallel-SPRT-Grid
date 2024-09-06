/*************************************
    Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#ifndef PSPS_OPENGL_3D_SOLVER_H
#define PSPS_OPENGL_3D_SOLVER_H

// Standard library
#include <limits>
#include <cmath>
#include <chrono>

#include "PSPS_OpenGL_BF_3D_ComShader.h"

#define DIM 3

class PSPS_OpenGL_BF_3D_Solver : public PSPS_OpenGL_BF_3D_ComShader
{
public:
    explicit PSPS_OpenGL_BF_3D_Solver();
    ~PSPS_OpenGL_BF_3D_Solver();

   bool Init();
   bool Set();
   bool Compute();
   bool Retrieve();
   bool Clear();

private:
   // Compute shader
   GLuint Parallel_Solver ;
   // Maps
   GLuint Slowness_Map, Traveltime_Map[2], Raypath_Map, UpdateMap_Map[2] ;
   // Some useful numbers
   GLuint length[DIM], top, front, side ;

	// Others
	void addDummyVertices();
	void removeDummyVertices();

};

#endif // PSPS_OPENGL_3D_SOLVER_H
