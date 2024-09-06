/*************************************
    Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#include "PSPS_OpenGL_FIM_2D_Program.h"

#define DIM 2
#define TyPe GLfloat

PSPS_OpenGL_FIM_2D_Program::PSPS_OpenGL_FIM_2D_Program()
{
	solver = new PSPS_OpenGL_FIM_2D_Solver ;
     asis = new PSPS_OpenGL_FIM_2D_Tools ;
}

PSPS_OpenGL_FIM_2D_Program::~PSPS_OpenGL_FIM_2D_Program()
{
	delete solver ;
	delete  asis  ;
}

void PSPS_OpenGL_FIM_2D_Program::runResumeProject(std::string &inFile, int maxloop)
{
   /// Welcome ///
   const GLubyte* renderer = glGetString(GL_RENDERER);
   std::cout << "Renderer: " << renderer << " is applied.\n" ;

   /// Read ///
   if( !(asis->resumeProject(inFile,maxloop)) ){ return; }
   
   /// INIT ///
   std::cout << "\nPSPS_OpenGL_FIM >> Compiling compute shader.\n" ;
   if( !solver->Init() ){
      std::cout << "Sorry, PSPS_OpenGL_FIM cannot compile compute shaders.\n" ;
      return;
   }

   /// SET ///
   if( !solver->Set() ){
      std::cout << "Sorry, PSPS_OpenGL_FIM cannot set textures (Not enough Ram ?).\n" ;
      return;
   }

   /// COMPUTE ///
   if( !solver->Compute() ){
      std::cout << "Sorry, PSPS_OpenGL_FIM cannot apply compute shaders ??.\n" ;
      return;
   }

   /// RETRIEVE ///
   if( !solver->Retrieve() ){
      std::cout << "Sorry, PSPS_OpenGL_FIM cannot retrieve data from GPU (Not enough Ram ?).\n" ;
      return;
   }

   /// CLEAR ///
   if( !solver->Clear() ){
      std::cout << "Sorry, PSPS_OpenGL_FIM cannot clear buffers or compute shaders.\n" ;
      return;
   }

   /// WRITE ///
   if( !asis->saveProject(inFile) ){
      std::cout << "Sorry, PSPS_OpenGL_FIM cannot save the project.\n" ;
      return;
   }

return; }

void PSPS_OpenGL_FIM_2D_Program::runNewProject(std::string &inFile, std::vector<int> &param)
{
   /// Welcome ///
   const GLubyte* renderer = glGetString(GL_RENDERER);
   std::cout << "Renderer: " << renderer << " is applied.\n" ;

   /// Read ///
   if( !(asis->newProject(inFile,param)) ){ return; }
   
   /// INIT ///
   std::cout << "\nPSPS_OpenGL_FIM >> Compiling compute shader.\n" ;
   if( !solver->Init() ){
      std::cout << "Sorry, PSPS_OpenGL_FIM cannot compile compute shaders.\n" ;
      return;
   }

   /// SET ///
   if( !solver->Set() ){
      std::cout << "Sorry, PSPS_OpenGL_FIM cannot set textures (Not enough Ram ?).\n" ;
      return;
   }

   /// COMPUTE ///
   if( !solver->Compute() ){
      std::cout << "Sorry, PSPS_OpenGL_FIM cannot apply compute shaders ??.\n" ;
      return;
   }

   /// RETRIEVE ///
   if( !solver->Retrieve() ){
      std::cout << "Sorry, PSPS_OpenGL_FIM cannot retrieve data from GPU (Not enough Ram ?).\n" ;
      return;
   }

   /// CLEAR ///
   if( !solver->Clear() ){
      std::cout << "Sorry, PSPS_OpenGL_FIM cannot clear buffers or compute shaders.\n" ;
      return;
   }

   /// WRITE ///
   if( !asis->saveProject(inFile) ){
      std::cout << "Sorry, PSPS_OpenGL_FIM cannot save the project.\n" ;
      return;
   }

return; }



