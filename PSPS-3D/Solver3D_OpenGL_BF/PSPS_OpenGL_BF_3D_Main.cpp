/*************************************
    Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

// GLEW
#define GLEW_STATIC
#include <GL/glew.h>
// GLFW
#include <GLFW/glfw3.h>
// Standard library
#include <iostream>
//#include <string>
//#include <vector>
#include <cstdlib>

// Real main program
#include "PSPS_OpenGL_BF_3D_Program.h"

static void error_callback(int error, const char *description);

int main(int argc, char* argv[])
{
   /*************************************/
   // Start the Program.
   /*************************************/

   // Set GLFW Error handler
	glfwSetErrorCallback(error_callback);

  	// Initialize GLFW --> to initialize Glew
  	if(!glfwInit()) {
    	std::cout << "ERROR: Failed to initialize GLFW3.\n" << std::endl;
      return -1;
  	}

	// Check request 4.3
   glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
   glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
   glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  	glfwWindowHint(GLFW_VISIBLE, GL_FALSE);
	GLFWwindow* window = glfwCreateWindow(1, 1, "-/-", NULL, NULL);
	if(!window){
      std::cout << "ERROR:: Could not open window with GLFW3.\n" << std::endl;
	   glfwTerminate();
	   return 1;
	}
	glfwMakeContextCurrent(window);

	// Initialize GLEW --> setup the OpenGL_BF Function pointers
   glewExperimental = GL_TRUE ;
	GLenum error = glewInit();
   if( error != GLEW_OK ){
      std::cout << "ERROR:: Failed to initialize GLEW.\n"
                << ":>> " << glewGetErrorString(error) << std::endl;
      return -1;
   }
   glGetError();

   /*************************************/
   // Main program
   /*************************************/

   // Read input
   std::vector<int> input ;
   std::string fileName ;
   // Rerun project
   if( argc==3 ) // program+filename+maxloop
   { 
       fileName = argv[1] ;
       std::cout << "************************************************\n" ;
       std::cout << "PSPS_3D_OpenGL_BF >> Select resume project.\n" ; 
       // Run
       PSPS_OpenGL_BF_3D_Program program ;
       program.runResumeProject(fileName, atoi(argv[2]) ); // Max loop
       std::cout << "PSPS_3D_OpenGL_BF >> End.\n" ;
       std::cout << "************************************************\n" ;
   
    // New project
   }else 
   if( argc==12 ) // program+filename+3*radius+3*resolution+3*source+maxloop
   {  
       fileName = argv[1] ;
       input.push_back( atoi(argv[2]) ); // Radius
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
       std::cout << "PSPS_3D_OpenGL_BF >> Select new project.\n" ;   
       // Run
       PSPS_OpenGL_BF_3D_Program program ;
       program.runNewProject(fileName,input);
       std::cout << "PSPS_3D_OpenGL_BF >> End.\n" ;
       std::cout << "************************************************\n" ;
   
   }else{
       std::cout << "ERROR:: Too few or too many input arguments.\n" << std::endl; 
       std::cout << "Expect \"./Solver3D_OpenGL_BF.exe InputFile.nc RadiusX RadiusY RadiusZ SizeX SizeY SizeZ SourceX SourceY SourceZ MaxIteration\".\n" << std::endl;
       return 2 ;
   }

   /*************************************/
   // End program
   /*************************************/
  	glfwTerminate();

return 0; }

static void error_callback(int error, const char *description)
{
    fprintf(stderr, "ERROR:: %s.\n", description);

return ;}
