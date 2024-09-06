/*************************************
    Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

// Real main program
#include "PSPS_CUDA_FIM_3D_Program.h"

int main(int argc, char* argv[])
{
   /*************************************/
   // Start the Program.
   /*************************************/

   // Read input
   std::vector<int> input ;
   std::string fileName ;
   // Rerun project
   if( argc==3 ) // program+filename+maxloop
   { 
       fileName = argv[1] ;
       std::cout << "************************************************\n" ;
       std::cout << "PSPS_3D_CUDA_FIM >> Select resume project.\n" ; 
       // Run
       PSPS_CUDA_FIM_3D_Program program ;
       program.runResumeProject(fileName, atoi(argv[2]) ); // Max loop
       std::cout << "PSPS_3D_CUDA_FIM >> End.\n" ;
       std::cout << "************************************************\n" ;
   
    // New project
   }else 
   if( argc==12 ) // program+filename+3*block+3*resolution+3*source+maxloop
   {  
       fileName = argv[1] ;
       input.push_back( atoi(argv[2]) ); // Block
       input.push_back( atoi(argv[3]) );
	   input.push_back( atoi(argv[4]) );
       input.push_back( atoi(argv[5]) ); // Resolution
       input.push_back( atoi(argv[6]) );
	   input.push_back( atoi(argv[7]) );
       input.push_back( atoi(argv[8])-1 ); // Source
       input.push_back( atoi(argv[9])-1 );
	   input.push_back( atoi(argv[10])-1 );
       input.push_back( atoi(argv[11]) ); // Max loop
       std::cout << "************************************************\n" ;
       std::cout << "PSPS_3D_CUDA_FIM >> Select new project.\n" ;   
       // Run
       PSPS_CUDA_FIM_3D_Program program ;
       program.runNewProject(fileName,input);
       std::cout << "PSPS_3D_CUDA_FIM >> End.\n" ;
       std::cout << "************************************************\n" ;
   
   }else{
       std::cout << "ERROR:: Too few or too many input arguments.\n" << std::endl; 
       std::cout << "Expect \"./Solver3D_CUDA_FIM.exe InputFile.nc BlockX BlockY BlockZ SizeX SizeY SizeZ SourceX SourceY SourceZ MaxIteration\".\n" << std::endl;
       return 2 ;
   }

   /*************************************/
   // End program
   /*************************************/

return 0; }


