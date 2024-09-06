/*************************************
    Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#ifndef PSPS_OPENGL_2D_TOOLS_H
#define PSPS_OPENGL_2D_TOOLS_H

// NetCDF
#include <netcdf.h>
// Standard library
#include <limits>
#include <cstring>
#include <cstdlib>
#include <cmath>

#include "PSPS_OpenGL_FIM_2D_ComShader.h"

#define DIM 2
#define TyPe GLfloat

class PSPS_OpenGL_FIM_2D_Tools : public PSPS_OpenGL_FIM_2D_ComShader
{
public:
   explicit PSPS_OpenGL_FIM_2D_Tools();
   ~PSPS_OpenGL_FIM_2D_Tools();

   /// INIT ///
   bool newProject(std::string &inFile, std::vector<int> &param);
   bool resumeProject(std::string &inFile, int maxloop);

   /// WRITE ///
   bool saveProject(std::string &inFile);

private:
   /// INIT ///
   bool readSlownessModel(const char* filePath);
   bool readResumeProject(const char* filePath);
   bool inputSpec(std::vector<int> &param);
   bool checkSpec();

   /// RESCALE ///
   bool rescaleSlowness(std::vector<GLuint> &res_old);

};

#endif // PSPS_OPENGL_2D_TOOLS_H

