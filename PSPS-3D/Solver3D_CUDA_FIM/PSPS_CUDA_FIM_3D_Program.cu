/*************************************
    Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#include "PSPS_CUDA_FIM_3D_Program.h"

#define DIM 3
#define TyPe float

PSPS_CUDA_FIM_3D_Program::PSPS_CUDA_FIM_3D_Program()
{
	solver = new PSPS_CUDA_FIM_3D_Solver ;
     asis = new PSPS_CUDA_FIM_3D_Tools ;
}

PSPS_CUDA_FIM_3D_Program::~PSPS_CUDA_FIM_3D_Program()
{
	delete solver ;
	delete  asis  ;
}

void PSPS_CUDA_FIM_3D_Program::runResumeProject(std::string &inFile, int maxloop)
{
   /// Read ///
   if( !(asis->resumeProject(inFile,maxloop)) ){ return; }

   /// COMPUTE ///
   if( !solver->Compute() ){
      std::cout << "Sorry, PSPS_CUDA_FIM = ERROR @,@.\n" ;
      return;
   }

   /// WRITE ///
   if( !asis->saveProject(inFile) ){
      std::cout << "Sorry, PSPS_CUDA_FIM cannot save the project.\n" ;
      return;
   }

return; }

void PSPS_CUDA_FIM_3D_Program::runNewProject(std::string &inFile, std::vector<int> &param)
{
   /// Read ///
   if( !(asis->newProject(inFile,param)) ){ return; }

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



