/*************************************
     Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#ifndef PSPS_CPU_FMM_2D_PROGRAM_H
#define PSPS_CPU_FMM_2D_PROGRAM_H

// NetCDF
#include <netcdf.h>
// Standard library
#include <iostream>

#include "PSPS_CPU_FMM_2D_Solver.h"
#include "PSPS_CPU_FMM_2D_Tools.h"

class PSPS_CPU_FMM_2D_Program
{
public:
   explicit PSPS_CPU_FMM_2D_Program();
   ~PSPS_CPU_FMM_2D_Program();

   void runNewProject(std::string &inFile, std::vector<int> &param);

private:
   PSPS_CPU_FMM_2D_Solver* solver ;
   PSPS_CPU_FMM_2D_Tools* asis ;

};

#endif // PSPS_CPU_FMM_2D_PROGRAM_H

