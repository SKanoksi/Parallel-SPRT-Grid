/*************************************
    Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#define DIM 2
#define TyPe GLfloat

#include "PSPS_OpenGL_FIM_2D_Solver.h"

PSPS_OpenGL_FIM_2D_Solver::PSPS_OpenGL_FIM_2D_Solver()
{

}

PSPS_OpenGL_FIM_2D_Solver::~PSPS_OpenGL_FIM_2D_Solver()
{

}

bool PSPS_OpenGL_FIM_2D_Solver::Init()
{
   /// Create Compute Shader (Solver) ///
   // Compute some useful numbers.
   for(GLuint i=0 ; i<DIM ; ++i)
   {
      length[i] = res[i]+2 ;
   }
   side = 1 ;
   front = length[0] ;

   // Read code
   std::string comCode, temp ;
   if( !readSource("./PSPS_FIM_Solver2D.comSh",comCode) ){ // The code must be found ###
      std::cout << "ERROR:: PSPS_FIM_Solver2D.comSh is not found.\n" ;
      return false;
   }
   // Define some useful constant variables
     ToString( block[0], temp);
   comCode.replace( comCode.find("BLOCK_X"), 7,  temp.c_str() ); // The marks MUST be in the code ###
     ToString( block[1], temp);
   comCode.replace( comCode.find("BLOCK_Y"), 7, temp.c_str() );
     ToString( block[0]+2, temp);
   comCode.replace( comCode.find("SHARED_X"), 8,  temp.c_str() );
     ToString( block[1]+2, temp);
   comCode.replace( comCode.find("SHARED_Y"), 8, temp.c_str() );
     ToString( block[0], temp);
   comCode.replace( comCode.find("BLock_X"), 8,  temp.c_str() );
     ToString( block[1], temp);
   comCode.replace( comCode.find("BLock_Y"), 8, temp.c_str() );
     ToString( stride[0], temp);
   comCode.replace( comCode.find("STRIDE_X"), 8,  temp.c_str() );
     ToString( stride[1], temp);
   comCode.replace( comCode.find("STRIDE_Y"), 8, temp.c_str() );  
     ToString( block[0]+block[1], temp);
   comCode.replace( comCode.find("ITERPERBLOCK"), 12, temp.c_str() );
    
   // Compile code
   Parallel_Solver = createComSh(comCode.c_str(), comCode.length());

return true; }

bool PSPS_OpenGL_FIM_2D_Solver::Set()
{
   /// Add dummy vertices ONLY 2D ///
   addDummyVertices();

   /// Create all maps ///
   glGenTextures(1, &Slowness_Map);
   glGenTextures(2, &Traveltime_Map[0]);
   glGenTextures(2, &UpdateMap_Map[0]);

   /// Set all maps ///
   /// Set slowness
   glBindTexture(GL_TEXTURE_2D, Slowness_Map);
     glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, length[0], length[1] , 0, GL_RED, GL_FLOAT, (GLvoid*)&slowness[0]);
      // Needed by Nvidia, even if not used. !!!
     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      
      /// Set traveltime x2, updateMap x2
      // nloop%2==0 -> read
      glBindTexture(GL_TEXTURE_2D, Traveltime_Map[0]);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, length[0], length[1], 0, GL_RED, GL_FLOAT, (GLvoid*)&traveltime[0]);
        // Needed by Nvidia, even if not used. !!!
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
       glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glBindTexture(GL_TEXTURE_2D, UpdateMap_Map[0]);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R32I, nWorkGroup[0]+2, nWorkGroup[1]+2, 0, GL_RED_INTEGER, GL_INT, (GLvoid*)&updateMap[0]);
        // Needed by Nvidia, even if not used. !!!
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
       glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      // nloop%2==0 -> write
      glBindTexture(GL_TEXTURE_2D, Traveltime_Map[1]);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, length[0], length[1], 0, GL_RED, GL_FLOAT, (GLvoid*)&traveltime[0]);
        // Needed by Nvidia, even if not used. !!!
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glBindTexture(GL_TEXTURE_2D, UpdateMap_Map[1]);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R32I, nWorkGroup[0]+2, nWorkGroup[1]+2, 0, GL_RED_INTEGER, GL_INT, (GLvoid*)0);
        // Needed by Nvidia, even if not used. !!!
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
   glBindTexture(GL_TEXTURE_2D, 0);

return true; }

bool PSPS_OpenGL_FIM_2D_Solver::Compute()
{
   /// INITIAL
   // Open program
   glUseProgram(Parallel_Solver);
   // Set unswaped textures (Slowness)
   glBindImageTexture(0, Slowness_Map , 0, GL_FALSE, 0, GL_READ_ONLY , GL_R32F);

   // Create Finish_Buffer
   GLint running = 0 ;
   GLuint FinishSSBO ;
   glGenBuffers(1, &FinishSSBO);
   glBindBuffer(GL_SHADER_STORAGE_BUFFER, FinishSSBO);
      glBufferData(GL_SHADER_STORAGE_BUFFER, 1*sizeof(GLint), &running, GL_DYNAMIC_READ); // use other options ? ###
   // Set Finish_Buffer
   glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, FinishSSBO);

   // Create and Set Indirect_Buffer
   GLuint IndBO ; // dispatch buffer object
   glGenBuffers(1, &IndBO);
   glBindBuffer(GL_DISPATCH_INDIRECT_BUFFER, IndBO);
   static const struct {
      GLuint num_groups_x ;
      GLuint num_groups_y ;
      GLuint num_groups_z ;
   } indirect = { nWorkGroup[0], nWorkGroup[1], 1 };
   glBufferData(GL_DISPATCH_INDIRECT_BUFFER, sizeof(indirect), &indirect, GL_STATIC_DRAW);


   { ///***** =.\= Parallel Shortest Path Solver (OpenGL_FIM Compute Shader) =.\= *****///
    
      // Current Time
      time_t rawtime ;
      
      // Start
      std::cout << "PSPS_OpenGL_FIM >> Start Fast Iterative Method (Compute shader).\n" ;
      time(&rawtime);
      std::cout << "   START:: " << ctime(&rawtime) ;
      
      // Start timer
      auto start_time = std::chrono::high_resolution_clock::now();
      
      ///*** Main loop ***///
      GLint loop ;
      for(loop=0 ; loop<maxLoop ; ++loop)
      {
         // Set swaped textures (Traveltime & UpdateMap)
         GLuint swap[2] = {0,1} ;
         if( loop%2==1 ){ swap[0] = 1 ; swap[1] = 0 ;}
         glBindImageTexture( 1, Traveltime_Map[swap[0]], 0, GL_FALSE, 0, GL_READ_ONLY , GL_R32F);
         glBindImageTexture( 2,  UpdateMap_Map[swap[0]], 0, GL_FALSE, 0, GL_READ_ONLY , GL_R32I);
         glBindImageTexture( 3, Traveltime_Map[swap[1]], 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_R32F);
         glBindImageTexture( 4,  UpdateMap_Map[swap[1]], 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_R32I);

         /***( Solving =/.= )***/
         glDispatchComputeIndirect(0);

         // Barrier
         glMemoryBarrier(GL_ALL_BARRIER_BITS);
         // Check and Reset Finish_Buffer
         GLint *ptr = (GLint*) glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_WRITE);
         running = *ptr ;
         *ptr = 0 ;
         glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);
         //std::cout << "  Loop = " << loop << "  Running = " << running << " \n" ; // for checking
         if( running==0 ){ isFinish = true ; break; }

      }

      // Stop timer
      auto end_time = std::chrono::high_resolution_clock::now();
      
      // Finish
      time(&rawtime);
      std::cout << "  FINISH:: " << ctime(&rawtime) ;
      std::cout << "PSPS_OpenGL_FIM >> Finish Fast Iterative Method (Compute shader).\n\n" ;
        
      std::cout << " *** Runtime = " << std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count() << " seconds. ***\n";  
        
      // Save nloop
      nloop += loop ;
      // Display now state
      if( isFinish ){
         std::cout << " *** Shortest paths are all found, after total " << nloop << " iterations. ***\n\n" ;
      }else{
         std::cout << " *** Solver is not finish finding Shortest paths, after total " << nloop << " iterations. ***\n\n" ;
      }
   }


   /// CLEAR
   // Delete Indirect_Buffer
   glBindBuffer(GL_DISPATCH_INDIRECT_BUFFER, 0);
   glDeleteBuffers(1, &IndBO);
   // Delete Finish_Buffer
   glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
   glDeleteBuffers(1, &FinishSSBO);
   // Close program
   glUseProgram(0);


return true; }

bool PSPS_OpenGL_FIM_2D_Solver::Retrieve()
{
   /// Retrieve Maps ///
   // Traveltime x1
   glBindTexture(GL_TEXTURE_2D, Traveltime_Map[ (nloop%2==0) ? 1:0 ]);
     glGetTexImage(GL_TEXTURE_2D, 0, GL_RED, GL_FLOAT, &traveltime[0]);

   // if not finish -> UpdateMap x1
   if( !isFinish ){
      glBindTexture(GL_TEXTURE_2D, UpdateMap_Map[ (nloop%2==0) ? 1:0 ]);
        glGetTexImage(GL_TEXTURE_2D, 0, GL_RED_INTEGER, GL_INT, &updateMap[0]);
   }
   glBindTexture(GL_TEXTURE_2D, 0);

   /// Remove dummy vertices ONLY 2D ///
   removeDummyVertices();

return true; }

bool PSPS_OpenGL_FIM_2D_Solver::Clear()
{
   /// Clear Shader and Buffers ///
   glDeleteProgram(Parallel_Solver);
   glDeleteBuffers(1, &Slowness_Map);
   glDeleteBuffers(2, &Traveltime_Map[0]);
   glDeleteBuffers(2, &UpdateMap_Map[0]);

return true; }

/************************************* Private *************************************/

void PSPS_OpenGL_FIM_2D_Solver::addDummyVertices()
{
   GLuint at ;
   bool isNegative = false ;
   std::vector<TyPe> temp ;

   /// Slowness map.
   // Copy data
   temp.assign( slowness.begin(), slowness.end() ) ;
   // Reinput the data
   slowness.assign( length[0]*length[1], 1.0f/0.0f) ;
   at = front+side ;
   for(GLuint j=0 ; j<res[1] ; ++j)
   {
      for(GLuint i=0 ; i<res[0] ; ++i)
      {
         if( temp[res[0]*j+i] < 0 ){ isNegative = true ; } // Slowness must only be positive ? ###
         slowness[ at+i ] = std::abs( temp[ res[0]*j+i ]*stride[0] ) ;
      }
      at += res[0]+2*side ;
   }
   if( isNegative ){
      std::cout << "Warning !! Negative slowness is detected -> using absolute value.\n" ;
   }

   /// Traveltime map.
   // Copy data
   temp.assign( traveltime.begin(), traveltime.end() ) ;
   // Reinput the data
   traveltime.assign( length[0]*length[1], 1.0f/0.0f) ;
   at = front+side ;
   for(GLuint j=0 ; j<res[1] ; ++j)
   {
      for(GLuint i=0 ; i<res[0] ; ++i)
      {
         traveltime[ at+i ] = temp[ res[0]*j+i ] ;
      }
      at += res[0]+2*side ;
   }

   /// UpdateMap.
   // Copy data
   std::vector<GLint> tem ;
   tem.assign( updateMap.begin(), updateMap.end() ) ;
   // Reinput the data
   updateMap.assign( (nWorkGroup[0]+2)*(nWorkGroup[1]+2), 0) ;
   at = (nWorkGroup[0]+2)+1 ;
   for(GLuint j=0 ; j<nWorkGroup[1] ; ++j)
   {
      for(GLuint i=0 ; i<nWorkGroup[0] ; ++i)
      {
         updateMap[ at+i ] = tem[ nWorkGroup[0]*j+i ] ;
      }
      at += nWorkGroup[0]+2 ;
   }

return; }

void PSPS_OpenGL_FIM_2D_Solver::removeDummyVertices()
{
   GLuint at ;
   std::vector<TyPe> temp ;

   /// Slowness map
   // Copy data
   temp.assign( slowness.begin(), slowness.end() ) ;
   // Reinput the data
   slowness.resize( res[0]*res[1] ) ;
   at = front + side ;
   for(GLuint j=0 ; j<res[1] ; ++j)
   {
      for(GLuint i=0 ; i<res[0] ; ++i)
      {
         slowness[ res[0]*j+i ] = temp[ at+i ]/stride[0] ;
      }
      at += res[0]+2*side ;
   }

   /// Traveltime map
   // Copy data
   temp.assign( traveltime.begin(), traveltime.end() ) ;
   // Reinput the data
   traveltime.resize( res[0]*res[1] ) ;
   at = front + side ;
   for(GLuint j=0 ; j<res[1] ; ++j)
   {
      for(GLuint i=0 ; i<res[0] ; ++i)
      {
         traveltime[ res[0]*j+i ] = temp[ at+i ] ;
      }
      at += res[0]+2*side ;
   }

   /// UpdateMap.
   if( isFinish ){ return ; } //No more update map
   // Copy data
   std::vector<GLint> tem ;
   tem.assign( updateMap.begin(), updateMap.end() ) ;
   // Reinput the data
   updateMap.resize( nWorkGroup[0]*nWorkGroup[1] ) ;
   at = (nWorkGroup[0]+2)+1 ;
   for(GLuint j=0 ; j<nWorkGroup[1] ; ++j)
   {
      for(GLuint i=0 ; i<nWorkGroup[0] ; ++i)
      {
         updateMap[ nWorkGroup[0]*j+i ] = tem[ at+i ] ;
      }
      at += nWorkGroup[0]+2 ;
   }

return; }
