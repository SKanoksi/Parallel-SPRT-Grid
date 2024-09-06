/*************************************
    Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#include "PSPS_OpenGL_BF_2D_ComShader.h"

#define DIM 2
#define TyPe GLfloat

PSPS_OpenGL_BF_2D_ComShader::PSPS_OpenGL_BF_2D_ComShader()
{

}

PSPS_OpenGL_BF_2D_ComShader::~PSPS_OpenGL_BF_2D_ComShader()
{

}

/// SetUp ///
bool PSPS_OpenGL_BF_2D_ComShader::isFinish = false;
GLuint PSPS_OpenGL_BF_2D_ComShader::nloop = 0;
int PSPS_OpenGL_BF_2D_ComShader::maxLoop = 0;

std::vector<GLuint> PSPS_OpenGL_BF_2D_ComShader::res(DIM,0) ;
std::vector<GLuint> PSPS_OpenGL_BF_2D_ComShader::block(DIM,0) ;
std::vector<GLuint> PSPS_OpenGL_BF_2D_ComShader::source(DIM,0) ;
std::vector<GLuint> PSPS_OpenGL_BF_2D_ComShader::nWorkGroup(DIM,0) ;
std::vector<TyPe> PSPS_OpenGL_BF_2D_ComShader::stride(DIM,0.0f) ;

/// Data Maps ///
std::vector<TyPe> PSPS_OpenGL_BF_2D_ComShader::slowness(DIM,0) ;
std::vector<TyPe> PSPS_OpenGL_BF_2D_ComShader::traveltime(DIM,1.0f/0.0f) ;
std::vector<GLint> PSPS_OpenGL_BF_2D_ComShader::raypath(DIM,0) ;
std::vector<GLint> PSPS_OpenGL_BF_2D_ComShader::updateMap(DIM,0) ;

bool PSPS_OpenGL_BF_2D_ComShader::readSource(const char *sourcePath, std::string &comCode)
{
    // 0. Read file
    std::ifstream file(sourcePath,std::ifstream::in);
    std::stringstream buffer;
    buffer << file.rdbuf();
    file.close();

    // 1. Retrieve the source code from filePath
    comCode = buffer.str();

return true; }

GLuint PSPS_OpenGL_BF_2D_ComShader::createComSh(const char *comShaderCode, const GLint lengthCode)
{
    // 2. Compile shaders
    GLint success; GLchar infoLog[512]; // Check
    GLuint comShader ;

    // Only Compute Shader (Optional)
    comShader = glCreateShader(GL_COMPUTE_SHADER);
    glShaderSource(comShader, 1, &comShaderCode, &lengthCode);
    glCompileShader(comShader);
    // Print compile errors if any
    glGetShaderiv(comShader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        glGetShaderInfoLog(comShader, 512, NULL, infoLog);
        std::cout << "ERROR::COMPUTE::SHADER::COMPILATION_FAILED\n" << infoLog << std::endl;
    }

    // 3.Create shader program
    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, comShader);

    glLinkProgram(shaderProgram);
    // Check for linking errors
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if(!success) {
        glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
    }
    // Should detach and delete shader after a successful link.
    glDetachShader(shaderProgram, comShader);
    glDeleteShader(comShader);

return shaderProgram ; }

