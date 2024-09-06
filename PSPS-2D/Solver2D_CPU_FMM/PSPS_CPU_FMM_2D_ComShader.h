/*************************************
     Parallel Shortest path solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#ifndef PSPS_CPU_FMM_2D_COMSHADER_H
#define PSPS_CPU_FMM_2D_COMSHADER_H

// Standard library
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#define TyPe float

class PSPS_CPU_FMM_2D_ComShader
{
public:
    PSPS_CPU_FMM_2D_ComShader();
    ~PSPS_CPU_FMM_2D_ComShader();

protected:
   /// SetUp ///
   static std::vector<unsigned int> res, block, source ;
   static std::vector<TyPe> stride ;

   /// Data map ///
   static std::vector<TyPe> slowness, traveltime ;
   static std::vector<int> raypath ;

};

#endif // PSPS_CPU_FMM_2D_COMSHADER_H
