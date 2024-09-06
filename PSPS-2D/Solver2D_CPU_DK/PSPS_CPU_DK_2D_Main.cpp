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
#include "PSPS_CPU_DK_2D_Program.h"

int main(int argc, char* argv[])
{
   /*************************************/
   // Start the Program.
   /*************************************/

   // Read input
   std::vector<int> input ;
   std::string fileName ;

   // New project
   if( argc==6 ) // program+filename+2*radius+2*source
   {  
       fileName = argv[1] ;
       input.push_back( atoi(argv[2]) ); // Radius
       input.push_back( atoi(argv[3]) );
       input.push_back( atoi(argv[4])-1 ); // Source
       input.push_back( atoi(argv[5])-1 );
       std::cout << "************************************************\n" ;
       std::cout << "PSPS_2D_CPU_DK >> Select new project.\n" ;   
       // Run
       PSPS_CPU_DK_2D_Program program ;
       program.runNewProject(fileName,input);
       std::cout << "PSPS_2D_CPU_DK >> End.\n" ;
       std::cout << "************************************************\n" ;
   
   }else{
       std::cout << "ERROR:: Too few or too many input arguments.\n" << std::endl; 
       std::cout << "Expect \"./Solver2D_CPU_DK.exe InputFile.nc RadiusX RadiusY SourceX SourceY\".\n" << std::endl; 
       return 2 ;
   }

   /*************************************/
   // End program
   /*************************************/

return 0; }
