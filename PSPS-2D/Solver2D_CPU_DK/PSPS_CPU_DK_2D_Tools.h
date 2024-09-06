/*************************************
     Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#ifndef PSPS_CPU_DK_2D_TOOLS_H
#define PSPS_CPU_DK_2D_TOOLS_H

// NetCDF
#include <netcdf.h>
// Standard library
#include <limits>
#include <cstring>
#include <cstdlib>

#include "PSPS_CPU_DK_2D_ComShader.h"

#define DIM 2
#define TyPe float

class PSPS_CPU_DK_2D_Tools : public PSPS_CPU_DK_2D_ComShader
{
public:
   explicit PSPS_CPU_DK_2D_Tools();
   ~PSPS_CPU_DK_2D_Tools();

   /// INIT ///
   bool newProject(std::string &inFile, std::vector<int> &param);

   /// WRITE ///
   bool saveProject(std::string &inFile);

private:
   /// INIT ///
   bool readSlownessModel(const char* filePath);
   bool inputSpec(std::vector<int> &param);
   bool checkSpec();

   /// RESCALE ///
   bool rescaleSlowness(std::vector<unsigned int> &res_old);

};

#endif // PSPS_CPU_DK_2D_TOOLS_H

