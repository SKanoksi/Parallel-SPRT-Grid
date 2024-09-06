/*************************************
    Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

// Real main program
#include "PSPS_CUDA_ExtendedBF_3D_Program.h"

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
       std::cout << "PSPS_3D_CUDA_ExtendedBF >> Select resume project.\n" ; 
       // Run
       PSPS_CUDA_ExtendedBF_3D_Program program ;
       program.runResumeProject(fileName, atoi(argv[2]) ); // Max loop
       std::cout << "PSPS_3D_CUDA_ExtendedBF >> End.\n" ;
       std::cout << "************************************************\n" ;
   
    // New project
   }else 
   if( argc==9 ) // program+filename+3*resolution+3*source+maxloop
   {  
       fileName = argv[1] ;
       input.push_back( 2 ); // Radius
       input.push_back( 2 );
	   input.push_back( 2 );
       input.push_back( atoi(argv[2]) ); // Resolution
       input.push_back( atoi(argv[3]) );
	   input.push_back( atoi(argv[4]) );
       input.push_back( atoi(argv[5])-1 ); // Source
       input.push_back( atoi(argv[6])-1 );
	   input.push_back( atoi(argv[7])-1 );
       input.push_back( atoi(argv[8]) ); // Max loop
       std::cout << "************************************************\n" ;
       std::cout << "PSPS_3D_CUDA_ExtendedBF >> Select new project.\n" ;   
       // Run
       PSPS_CUDA_ExtendedBF_3D_Program program ;
       program.runNewProject(fileName,input);
       std::cout << "PSPS_3D_CUDA_ExtendedBF >> End.\n" ;
       std::cout << "************************************************\n" ;
   
   }else{
       std::cout << "ERROR:: Too few or too many input arguments.\n" << std::endl; 
       std::cout << "Expect \"./Solver3D_CUDA_ExtendedBF.exe InputFile.nc SizeX SizeY SizeZ SourceX SourceY SourceZ MaxIteration\".\n" << std::endl;
       return 2 ;
   }

   /*************************************/
   // End program
   /*************************************/

return 0; }


