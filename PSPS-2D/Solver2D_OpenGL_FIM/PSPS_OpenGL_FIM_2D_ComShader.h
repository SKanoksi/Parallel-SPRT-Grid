/*************************************
    Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#ifndef PSPS_OPENGL_2D_COMSHADER_H
#define PSPS_OPENGL_2D_COMSHADER_H

// GLEW
#define GLEW_STATIC
#include <GL/glew.h>

// Standard library
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#define TyPe GLfloat

class PSPS_OpenGL_FIM_2D_ComShader
{
public:
    PSPS_OpenGL_FIM_2D_ComShader();
    ~PSPS_OpenGL_FIM_2D_ComShader();

protected:
   /// SetUp ///
   static bool isFinish ;
   static GLuint nloop ;
   static int maxLoop ;
   static std::vector<GLuint> res, block, source, nWorkGroup ;
   static std::vector<TyPe> stride ;

   /// Data map ///
   static std::vector<TyPe> slowness, traveltime ;
   static std::vector<GLint> updateMap ;

   /// Compute shader ///
   bool readSource(const char *sourcePath, std::string &comCode);
   GLuint createComSh(const char *comPath, const GLint lengthCode);
   template <typename T>
   void ToString(T in, std::string &out)
   {
      std::stringstream stream ;
      stream << in ;
      out = stream.str() ;
   }

private:

};

#endif // PSPS_OPENGL_2D_COMSHADER_H
