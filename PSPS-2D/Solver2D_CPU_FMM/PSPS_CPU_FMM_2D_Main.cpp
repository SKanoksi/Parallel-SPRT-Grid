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
#include "PSPS_CPU_FMM_2D_Program.h"

int main(int argc, char* argv[])
{
   /*************************************/
   // Start the Program.
   /*************************************/

   // Read input
   std::vector<int> input ;
   std::string fileName ;

   // New project
   if( argc==4 ) // program+filename+2*source
   {  
       fileName = argv[1] ;
       input.push_back( 1 ); // Radius
       input.push_back( 1 );
       input.push_back( atoi(argv[2])-1 ); // Source
       input.push_back( atoi(argv[3])-1 );
       std::cout << "************************************************\n" ;
       std::cout << "PSPS_2D_CPU_FMM >> Select new project.\n" ;   
       // Run
       PSPS_CPU_FMM_2D_Program program ;
       program.runNewProject(fileName,input);
       std::cout << "PSPS_2D_CPU_FMM >> End.\n" ;
       std::cout << "************************************************\n" ;
   
   }else{
       std::cout << "ERROR:: Too few or too many input arguments.\n" << std::endl; 
       std::cout << "Expect \"./Solver2D_CPU_FMM.exe InputFile.nc SourceX SourceY\".\n" << std::endl; 
       return 2 ;
   }

   /*************************************/
   // End program
   /*************************************/

return 0; }
