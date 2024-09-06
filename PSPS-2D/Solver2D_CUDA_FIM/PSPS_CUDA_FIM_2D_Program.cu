/*************************************
     Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#include "PSPS_CUDA_FIM_2D_Program.h"

#define DIM 2
#define TyPe float

PSPS_CUDA_FIM_2D_Program::PSPS_CUDA_FIM_2D_Program()
{
	solver = new PSPS_CUDA_FIM_2D_Solver ;
     asis = new PSPS_CUDA_FIM_2D_Tools ;
}

PSPS_CUDA_FIM_2D_Program::~PSPS_CUDA_FIM_2D_Program()
{
	delete solver ;
	delete  asis  ;
}

void PSPS_CUDA_FIM_2D_Program::runResumeProject(std::string &inFile, int maxloop)
{
   /// Read ///
   if( !(asis->resumeProject(inFile,maxloop)) ){ return; }

   /// COMPUTE ///
   if( !solver->Compute() ){
      std::cout << "Sorry, PSPS_CUDA_FIM = ERROR @.@.\n" ;
      return;
   }

   /// WRITE ///
   if( !asis->saveProject(inFile) ){
      std::cout << "Sorry, PSPS_CUDA_FIM cannot save the project.\n" ;
      return;
   }

return; }

void PSPS_CUDA_FIM_2D_Program::runNewProject(std::string &inFile, std::vector<int> &param)
{
   /// Read ///
   if( !(asis->newProject(inFile,param)) ){ return; }

   /// COMPUTE ///
   if( !solver->Compute() ){
      std::cout << "Sorry, PSPS_CUDA_FIM = ERROR @.@.\n" ;
      return;
   }

   /// WRITE ///
   // Untill we can really resize the slowness map.
   if( !asis->saveProject(inFile) ){
      std::cout << "Sorry, PSPS_CUDA_FIM cannot save the project.\n" ;
      return;
   }

return; }

