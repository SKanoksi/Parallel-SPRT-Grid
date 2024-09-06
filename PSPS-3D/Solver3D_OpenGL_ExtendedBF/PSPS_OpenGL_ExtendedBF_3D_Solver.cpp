/*************************************
    Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#define DIM 3
#define TyPe GLfloat

#include "PSPS_OpenGL_ExtendedBF_3D_Solver.h"

PSPS_OpenGL_ExtendedBF_3D_Solver::PSPS_OpenGL_ExtendedBF_3D_Solver()
{

}

PSPS_OpenGL_ExtendedBF_3D_Solver::~PSPS_OpenGL_ExtendedBF_3D_Solver()
{

}

bool PSPS_OpenGL_ExtendedBF_3D_Solver::Init()
{
   /// Create Compute Shader (Solver) ///
   // Compute some useful numbers.
   GLuint shared[DIM] ;
   for(GLuint i=0 ; i<DIM ; ++i)
   {
      length[i] = res[i]+block[i] ;
      shared[i] = 2*block[i] ;
   }
   side = block[0]/2 ;
   front = length[0]*block[1]/2 ;
   top = length[0]*length[1]*block[2]/2 ;

     // Read code
   std::string comCode, temp ;
   if( !readSource("./PSPS_ExtendedBF_Solver3D.comSh",comCode) ){ // The code must be found ###
      std::cout << "ERROR:: PSPS_ExtendedBF_Solver3D.comSh is not found (in the same folder as this program).\n" ;
      return false;
   }
   // Define some useful variables
     ToString(block[0],temp);
   comCode.replace( comCode.find("BLOCK_X"), 7,  temp.c_str() ); // The marks MUST be in the code ###
     ToString(block[1],temp);
   comCode.replace( comCode.find("BLOCK_Y"), 7, temp.c_str() );
     ToString(block[2],temp);
   comCode.replace( comCode.find("BLOCK_Z"), 7, temp.c_str() );
   ToString(shared[0],temp);
   comCode.replace( comCode.find("SHARED_X"), 8,  temp.c_str() );
     ToString(shared[1],temp);
   comCode.replace( comCode.find("SHARED_Y"), 8, temp.c_str() );
     ToString(shared[2],temp);
   comCode.replace( comCode.find("SHARED_Z"), 8, temp.c_str() );
     ToString(block[0]+1,temp);
   comCode.replace( comCode.find("GROUP_X"), 7,  temp.c_str() );
     ToString(block[1]+1,temp);
   comCode.replace( comCode.find("GROUP_Y"), 7, temp.c_str() );
     ToString(block[2]+1,temp);
   comCode.replace( comCode.find("GROUP_Z"), 7, temp.c_str() );
     ToString( GLuint(block[0]/2), temp);
   comCode.replace( comCode.find("RADIUS_X"), 8,  temp.c_str() );
     ToString( GLuint(block[1]/2), temp);
   comCode.replace( comCode.find("RADIUS_Y"), 8, temp.c_str() );
     ToString( GLuint(block[2]/2), temp);
   comCode.replace( comCode.find("RADIUS_Z"), 8, temp.c_str() );

   // Compile code
   Parallel_Solver = createComSh(comCode.c_str(), comCode.length());

return true; }

bool PSPS_OpenGL_ExtendedBF_3D_Solver::Set()
{
  /// Add dummy vertices ONLY 3D ///
   addDummyVertices();

   /// Create all maps ///
   glGenTextures(1, &Slowness_Map);
   glGenTextures(2, &Traveltime_Map[0]);
   glGenTextures(1, &Raypath_Map);
   glGenTextures(2, &UpdateMap_Map[0]);

   /// Set all maps ///
   /// Set slowness, raypath
   glBindTexture(GL_TEXTURE_3D, Slowness_Map);
     glTexImage3D(GL_TEXTURE_3D, 0, GL_R32F, length[0], length[1], length[2], 0, GL_RED, GL_FLOAT, (GLvoid*)&slowness[0]);
     // Needed by Nvidia, even if not used. !!!     
     glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
     glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
   glBindTexture(GL_TEXTURE_3D, Raypath_Map);
     glTexImage3D(GL_TEXTURE_3D, 0, GL_R32I, res[0], res[1], res[2], 0, GL_RED_INTEGER, GL_INT, (GLvoid*)&raypath[0]);
     // Needed by Nvidia, even if not used. !!!     
     glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
     glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

   /// Set traveltime x2, updateMap x2
      // nloop%2==0 -> read
      glBindTexture(GL_TEXTURE_3D, Traveltime_Map[0]);
        glTexImage3D(GL_TEXTURE_3D, 0, GL_R32F, length[0], length[1], length[2], 0, GL_RED, GL_FLOAT, (GLvoid*)&traveltime[0]);
        // Needed by Nvidia, even if not used. !!!     
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glBindTexture(GL_TEXTURE_3D, UpdateMap_Map[0]);
        glTexImage3D(GL_TEXTURE_3D, 0, GL_R32I, nWorkGroup[0]+2, nWorkGroup[1]+2, nWorkGroup[2]+2, 0, GL_RED_INTEGER, GL_INT, (GLvoid*)&updateMap[0]);
        // Needed by Nvidia, even if not used. !!!     
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      // nloop%2==0 -> write
      glBindTexture(GL_TEXTURE_3D, Traveltime_Map[1]);
        glTexImage3D(GL_TEXTURE_3D, 0, GL_R32F, length[0], length[1], length[2], 0, GL_RED, GL_FLOAT, (GLvoid*)&traveltime[0]);
        // Need by NVIDIA's texture !!!
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glBindTexture(GL_TEXTURE_3D, UpdateMap_Map[1]);
        glTexImage3D(GL_TEXTURE_3D, 0, GL_R32I, nWorkGroup[0]+2, nWorkGroup[1]+2, nWorkGroup[2]+2, 0, GL_RED_INTEGER, GL_INT, (GLvoid*)0);
        // Need for NVIDIA's texture !!!
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
   glBindTexture(GL_TEXTURE_3D, 0);

return true; }

bool PSPS_OpenGL_ExtendedBF_3D_Solver::Compute()
{
    /// INITIAL
   // Open program
   glUseProgram(Parallel_Solver);
   // Set unswaped textures (Slowness & Raypath)
   glBindImageTexture(0, Slowness_Map , 0, GL_TRUE, 0, GL_READ_ONLY , GL_R32F);
   glBindImageTexture(5, Raypath_Map  , 0, GL_TRUE, 0, GL_WRITE_ONLY, GL_R32I);

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
   } indirect = { nWorkGroup[0], nWorkGroup[1], nWorkGroup[2] };
   glBufferData(GL_DISPATCH_INDIRECT_BUFFER, sizeof(indirect), &indirect, GL_STATIC_DRAW);


   { ///***** =.\= Parallel Shortest Path Solver (OpenGL_ExtendedBF Compute Shader) =.\= *****///
    
      // Current Time
      time_t rawtime ;
      
      // Start
      std::cout << "PSPS_OpenGL_ExtendedBF >> Start \'Extended\' Bellman-Ford (Compute shader).\n" ;
      time(&rawtime);
      std::cout << "   START:: " << ctime(&rawtime) ;
      
      // Start timer
      auto start_time = std::chrono::high_resolution_clock::now();
      
      ///*** Main loop ***///
      GLint loop ;
      //GLuint totRunning = 0 ; 
      for(loop=0 ; loop<maxLoop ; ++loop)
      {
         // Set swaped textures (Traveltime & UpdateMap)
         GLuint swap[2] = {0,1} ;
         if( loop%2==1 ){ swap[0] = 1 ; swap[1] = 0 ;}
         glBindImageTexture( 1, Traveltime_Map[swap[0]], 0, GL_TRUE, 0, GL_READ_ONLY , GL_R32F);
         glBindImageTexture( 2,  UpdateMap_Map[swap[0]], 0, GL_TRUE, 0, GL_READ_ONLY , GL_R32I);
         glBindImageTexture( 3, Traveltime_Map[swap[1]], 0, GL_TRUE, 0, GL_WRITE_ONLY, GL_R32F);
         glBindImageTexture( 4,  UpdateMap_Map[swap[1]], 0, GL_TRUE, 0, GL_WRITE_ONLY, GL_R32I);

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
         //totRunning += running ;
         if( running==0 ){ isFinish = true ; break; }

      }
      //std::cout << "  Total Running = " << totRunning << "\n";
      
      // Stop timer
      auto end_time = std::chrono::high_resolution_clock::now();
      
      // Finish
      time(&rawtime);
      std::cout << "  FINISH:: " << ctime(&rawtime) ;
      std::cout << "PSPS_OpenGL_ExtendedBF >> Finish \'Extended\' Bellman-Ford (Compute shader).\n\n" ;
        
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

bool PSPS_OpenGL_ExtendedBF_3D_Solver::Retrieve()
{
   /// Retrieve Maps ///
   // Traveltime x1
   glBindTexture(GL_TEXTURE_3D, Traveltime_Map[ (nloop%2==0) ? 1:0 ]);
     glGetTexImage(GL_TEXTURE_3D, 0, GL_RED, GL_FLOAT, &traveltime[0]);
   // Raypath
   glBindTexture(GL_TEXTURE_3D, Raypath_Map);
     glGetTexImage(GL_TEXTURE_3D, 0, GL_RED_INTEGER, GL_INT, &raypath[0]);

   // if not finish -> UpdateMap x1
   if( !isFinish ){
      glBindTexture(GL_TEXTURE_3D, UpdateMap_Map[ (nloop%2==0) ? 1:0 ]);
        glGetTexImage(GL_TEXTURE_3D, 0, GL_RED_INTEGER, GL_INT, &updateMap[0]);
   }
   glBindTexture(GL_TEXTURE_3D, 0);

   /// Remove dummy vertices ONLY 3D ///
   removeDummyVertices();

return true; }

bool PSPS_OpenGL_ExtendedBF_3D_Solver::Clear()
{
   /// Clear Shader and Buffers ///
   glDeleteProgram(Parallel_Solver);
   glDeleteBuffers(1, &Slowness_Map);
   glDeleteBuffers(2, &Traveltime_Map[0]);
   glDeleteBuffers(1, &Raypath_Map);
   glDeleteBuffers(2, &UpdateMap_Map[0]);

return true; }

/************************************* Private *************************************/

void PSPS_OpenGL_ExtendedBF_3D_Solver::addDummyVertices()
{
   GLuint at ;
   bool isNegative = false ;
   std::vector<TyPe> temp ;

   /// Slowness map.
   // Copy data
   temp.assign( slowness.begin(), slowness.end() ) ;
   // Reinput the data
   slowness.assign( length[0]*length[1]*length[2], 1.0f/0.0f) ;
   at = top+front+side ;
   for(GLuint k=0 ; k<res[2] ; ++k)
   {
		for(GLuint j=0 ; j<res[1] ; ++j)
		{
			for(GLuint i=0 ; i<res[0] ; ++i)
			{
				if( temp[ (k*res[1]+j)*res[0]+i ] < 0 ){ isNegative = true ; } // Slowness must only be positive ? ###
				slowness[ at+i ] = std::abs( temp[ (k*res[1]+j)*res[0]+i ]*stride[0] ) ;
			}
			at += res[0]+2*side ;
		}
		at += 2*front ; 
   }
   if( isNegative ){
      std::cout << "Warning !! Negative slowness is detected -> using absolute value.\n" ;
   }

   /// Traveltime map.
   // Copy data
   temp.assign( traveltime.begin(), traveltime.end() ) ;
   // Reinput the data
   traveltime.assign( length[0]*length[1]*length[2], 1.0f/0.0f) ;
   at = top+front+side ;
   for(GLuint k=0 ; k<res[2] ; ++k)
   {
		for(GLuint j=0 ; j<res[1] ; ++j)
		{
			for(GLuint i=0 ; i<res[0] ; ++i)
			{
				traveltime[ at+i ] = temp[ (k*res[1]+j)*res[0]+i ] ;
			}
			at += res[0]+2*side ;
		}
		at += 2*front ; 
   }

   /// UpdateMap.
   // Copy data
   std::vector<GLint> tem ;
   tem.assign( updateMap.begin(), updateMap.end() ) ;
   // Reinput the data
   updateMap.assign( (nWorkGroup[0]+2)*(nWorkGroup[1]+2)*(nWorkGroup[2]+2), 0) ;
   at = (nWorkGroup[0]+2)*(nWorkGroup[1]+2)+(nWorkGroup[0]+2)+1 ;
   for(GLuint k=0 ; k<nWorkGroup[2] ; ++k)
   {
		for(GLuint j=0 ; j<nWorkGroup[1] ; ++j)
		{
			for(GLuint i=0 ; i<nWorkGroup[0] ; ++i)
			{
				updateMap[ at+i ] = tem[ (k*nWorkGroup[1]+j)*nWorkGroup[0]+i ] ;
			}
			at += nWorkGroup[0]+2 ;
		}
		at += 2*(nWorkGroup[0]+2) ;
   }

return; }

void PSPS_OpenGL_ExtendedBF_3D_Solver::removeDummyVertices()
{
   GLuint at ;
   std::vector<TyPe> temp ;

   /// Slowness map
   // Copy data
   temp.assign( slowness.begin(), slowness.end() ) ;
   // Reinput the data
   slowness.resize( res[0]*res[1]*res[2] ) ;
   at = top+front+side ;
   for(GLuint k=0 ; k<res[2] ; ++k)
   {
		for(GLuint j=0 ; j<res[1] ; ++j)
		{
			for(GLuint i=0 ; i<res[0] ; ++i)
			{
				slowness[ (k*res[1]+j)*res[0]+i ] = temp[ at+i ]/stride[0] ;
			}
			at += res[0]+2*side ;
		}
		at += 2*front ; 
   }
   
   /// Traveltime map
   // Copy data
   temp.assign( traveltime.begin(), traveltime.end() ) ;
   // Reinput the data
   traveltime.resize( res[0]*res[1]*res[2] ) ;
   at = top+front+side ;
   for(GLuint k=0 ; k<res[2] ; ++k)
   {
		for(GLuint j=0 ; j<res[1] ; ++j)
		{
			for(GLuint i=0 ; i<res[0] ; ++i)
			{
				traveltime[ (k*res[1]+j)*res[0]+i ] = temp[ at+i ] ;
			}
			at += res[0]+2*side ;
		}
		at += 2*front ; 
   }
   
   /// UpdateMap.
   if( isFinish ){ return ; } //No more update map
   // Copy data
   std::vector<GLint> tem ;
   tem.assign( updateMap.begin(), updateMap.end() ) ;
   // Reinput the data
   updateMap.resize( nWorkGroup[0]*nWorkGroup[1]*nWorkGroup[2] ) ;
   at = (nWorkGroup[0]+2)*(nWorkGroup[1]+2)+(nWorkGroup[0]+2)+1 ;
   for(GLuint k=0 ; k<nWorkGroup[2] ; ++k)
   {
		for(GLuint j=0 ; j<nWorkGroup[1] ; ++j)
		{
			for(GLuint i=0 ; i<nWorkGroup[0] ; ++i)
			{
				updateMap[ (k*nWorkGroup[1]+j)*nWorkGroup[0]+i ] = tem[ at+i ] ;
			}
			at += nWorkGroup[0]+2 ;
		}
		at += 2*(nWorkGroup[0]+2) ;
   }
   

return; }
