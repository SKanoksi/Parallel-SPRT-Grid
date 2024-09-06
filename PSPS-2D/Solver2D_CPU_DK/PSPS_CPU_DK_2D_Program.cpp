/*************************************
     Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#include "PSPS_CPU_DK_2D_Program.h"

#define DIM 2
#define TyPe float

PSPS_CPU_DK_2D_Program::PSPS_CPU_DK_2D_Program()
{
	solver = new PSPS_CPU_DK_2D_Solver ;
     asis = new PSPS_CPU_DK_2D_Tools ;
}

PSPS_CPU_DK_2D_Program::~PSPS_CPU_DK_2D_Program()
{
	delete solver ;
	delete  asis  ;
}

void PSPS_CPU_DK_2D_Program::runNewProject(std::string &inFile, std::vector<int> &param)
{
   /// Read ///
   if( !(asis->newProject(inFile,param)) ){ return; }

   /// COMPUTE ///
   if( !solver->Compute() ){
      std::cout << "Sorry, PSPS_CPU_DK = ERROR @.@.\n" ;
      return;
   }

   /// WRITE ///
   if( !asis->saveProject(inFile) ){
      std::cout << "Sorry, PSPS_CPU_DK cannot save the project.\n" ;
      return;
   }

return; }



