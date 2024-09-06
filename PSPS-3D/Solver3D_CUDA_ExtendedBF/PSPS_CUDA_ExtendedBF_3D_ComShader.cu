/*************************************
    Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#include "PSPS_CUDA_ExtendedBF_3D_ComShader.h"

#define DIM 3
#define TyPe float

PSPS_CUDA_ExtendedBF_3D_ComShader::PSPS_CUDA_ExtendedBF_3D_ComShader()
{

}

PSPS_CUDA_ExtendedBF_3D_ComShader::~PSPS_CUDA_ExtendedBF_3D_ComShader()
{

}

/// SetUp ///
bool PSPS_CUDA_ExtendedBF_3D_ComShader::isFinish = false;
unsigned int PSPS_CUDA_ExtendedBF_3D_ComShader::nloop = 0;
int PSPS_CUDA_ExtendedBF_3D_ComShader::maxLoop = 0;

std::vector<unsigned int> PSPS_CUDA_ExtendedBF_3D_ComShader::res(DIM,0) ;
std::vector<unsigned int> PSPS_CUDA_ExtendedBF_3D_ComShader::block(DIM,0) ;
std::vector<unsigned int> PSPS_CUDA_ExtendedBF_3D_ComShader::source(DIM,0) ;
std::vector<unsigned int> PSPS_CUDA_ExtendedBF_3D_ComShader::nWorkGroup(DIM,0) ;
std::vector<TyPe> PSPS_CUDA_ExtendedBF_3D_ComShader::stride(DIM,0.0f) ;

/// Data Maps ///
std::vector<TyPe> PSPS_CUDA_ExtendedBF_3D_ComShader::slowness(DIM,0) ;
std::vector<TyPe> PSPS_CUDA_ExtendedBF_3D_ComShader::traveltime(DIM,1.0f/0.0f) ;
std::vector<int> PSPS_CUDA_ExtendedBF_3D_ComShader::raypath(DIM,0) ;
std::vector<int> PSPS_CUDA_ExtendedBF_3D_ComShader::updateMap(DIM,0) ;

