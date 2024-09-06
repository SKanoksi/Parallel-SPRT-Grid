/*************************************
     Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#ifndef PSPS_CUDA_ExtendedBF_2D_TOOLS_H
#define PSPS_CUDA_ExtendedBF_2D_TOOLS_H

// NetCDF
#include <netcdf.h>
// Standard library
#include <limits>
#include <cstring>
#include <cstdlib>
#include <cmath>

#include "PSPS_CUDA_ExtendedBF_2D_ComShader.h"

#define DIM 2
#define TyPe float

class PSPS_CUDA_ExtendedBF_2D_Tools : public PSPS_CUDA_ExtendedBF_2D_ComShader
{
public:
   explicit PSPS_CUDA_ExtendedBF_2D_Tools();
   ~PSPS_CUDA_ExtendedBF_2D_Tools();

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
   bool rescaleSlowness(std::vector<unsigned int> &res_old);

};

#endif // PSPS_CUDA_ExtendedBF_2D_TOOLS_H

