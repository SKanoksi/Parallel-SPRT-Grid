/*************************************
    Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#ifndef PSPS_OPENGL_2D_PROGRAM_H
#define PSPS_OPENGL_2D_PROGRAM_H

// GLEW
#define GLEW_STATIC
#include <GL/glew.h>

// Standard library
#include <iostream>

#include "PSPS_OpenGL_ExtendedBF_2D_Solver.h"
#include "PSPS_OpenGL_ExtendedBF_2D_Tools.h"

class PSPS_OpenGL_ExtendedBF_2D_Program
{
public:
   explicit PSPS_OpenGL_ExtendedBF_2D_Program();
   ~PSPS_OpenGL_ExtendedBF_2D_Program();

   void runResumeProject(std::string &inFile, int maxloop);
   void runNewProject(std::string &inFile, std::vector<int> &param);

private:
   PSPS_OpenGL_ExtendedBF_2D_Solver* solver ;
   PSPS_OpenGL_ExtendedBF_2D_Tools* asis ;

};

#endif // PSPS_OPENGL_2D_PROGRAM_H

