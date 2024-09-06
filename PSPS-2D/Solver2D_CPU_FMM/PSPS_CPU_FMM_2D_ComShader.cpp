/*************************************
     Parallel Shortest path solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#include "PSPS_CPU_FMM_2D_ComShader.h"

#define DIM 2
#define TyPe float

PSPS_CPU_FMM_2D_ComShader::PSPS_CPU_FMM_2D_ComShader()
{

}

PSPS_CPU_FMM_2D_ComShader::~PSPS_CPU_FMM_2D_ComShader()
{

}

/// SetUp ///S
std::vector<unsigned int> PSPS_CPU_FMM_2D_ComShader::res(DIM,0) ;
std::vector<unsigned int> PSPS_CPU_FMM_2D_ComShader::block(DIM,0) ;
std::vector<unsigned int> PSPS_CPU_FMM_2D_ComShader::source(DIM,0) ;
std::vector<TyPe> PSPS_CPU_FMM_2D_ComShader::stride(DIM,0.0f) ;

/// Data Maps ///
std::vector<TyPe> PSPS_CPU_FMM_2D_ComShader::slowness(DIM,0) ;
std::vector<TyPe> PSPS_CPU_FMM_2D_ComShader::traveltime(DIM,1.0f/0.0f) ;
std::vector<int> PSPS_CPU_FMM_2D_ComShader::raypath(DIM,0) ;


