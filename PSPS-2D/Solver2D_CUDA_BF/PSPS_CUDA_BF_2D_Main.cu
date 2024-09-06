/*************************************
     Parallel Shortest path solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

// Standard library
//#include <iostream>
//#include <cstdlib>

// Real main program
#include "PSPS_CUDA_BF_2D_Program.h"

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
       std::cout << "PSPS_2D_CUDA_BF >> Select resume project.\n" ; 
       // Run
       PSPS_CUDA_BF_2D_Program program ;
       program.runResumeProject(fileName, atoi(argv[2]) ); // Max loop
       std::cout << "PSPS_2D_CUDA_BF >> End.\n" ;
       std::cout << "************************************************\n" ;
   
    // New project
   }else 
   if( argc==9 ) // program+filename+2*radius+2*resolution+2*source+maxloop
   {  
       fileName = argv[1] ;
       input.push_back( atoi(argv[2]) ); // Radius
       input.push_back( atoi(argv[3]) );
       input.push_back( atoi(argv[4]) ); // Resolution
       input.push_back( atoi(argv[5]) );
       input.push_back( atoi(argv[6])-1 ); // Source
       input.push_back( atoi(argv[7])-1 );
       input.push_back( atoi(argv[8]) ); // Max loop
       std::cout << "************************************************\n" ;
       std::cout << "PSPS_2D_CUDA_BF >> Select new project.\n" ;   
       // Run
       PSPS_CUDA_BF_2D_Program program ;
       program.runNewProject(fileName,input);
       std::cout << "PSPS_2D_CUDA_BF >> End.\n" ;
       std::cout << "************************************************\n" ;
   
   }else{
       std::cout << "ERROR:: Too few or too many input arguments.\n" << std::endl; 
       std::cout << "Expect \"./Solver2D_CUDA_BF.exe InputFile.nc RadiusX RadiusY SizeX SizeY SourceX SourceY MaxIteration\".\n" << std::endl; 
       return 2 ;
   }

   /*************************************/
   // End program
   /*************************************/

return 0; }
