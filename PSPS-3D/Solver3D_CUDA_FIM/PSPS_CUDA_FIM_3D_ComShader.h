/*************************************
    Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#ifndef PSPS_CUDA_3D_COMSHADER_H
#define PSPS_CUDA_3D_COMSHADER_H

// Standard library
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#define TyPe float

class PSPS_CUDA_FIM_3D_ComShader
{
public:
    PSPS_CUDA_FIM_3D_ComShader();
    ~PSPS_CUDA_FIM_3D_ComShader();

protected:
   /// SetUp ///
   static bool isFinish ;
   static unsigned int nloop ;
   static int maxLoop ;
   static std::vector<unsigned int> res, block, source, nWorkGroup ;
   static std::vector<TyPe> stride ;

   /// Data map ///
   static std::vector<TyPe> slowness, traveltime ;
   static std::vector<int> updateMap ;

private:

};

#endif // PSPS_CUDA_3D_COMSHADER_H
