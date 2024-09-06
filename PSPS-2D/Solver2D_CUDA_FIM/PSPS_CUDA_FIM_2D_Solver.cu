/*************************************
     Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#define DIM 2
#define TyPe float

#include "PSPS_CUDA_FIM_2D_Solver.h"

PSPS_CUDA_FIM_2D_Solver::PSPS_CUDA_FIM_2D_Solver()
{

}

PSPS_CUDA_FIM_2D_Solver::~PSPS_CUDA_FIM_2D_Solver()
{

}

//**************** CUDA_FIM program ****************//

__constant__ int radius_X ;
__constant__ int radius_Y ;

surface<void, cudaSurfaceType2D>  SlownessMap, Traveltime_A, Traveltime_B ;
surface<void, cudaSurfaceType2D>  UpdateMap_A,  UpdateMap_B ;

//**************** AB functions ****************//

/// Calculate new traveltime 
__device__ float newCandidateOfTraveltime(const int2 ptrMap)
{
	float T_L, T_R, T_U, T_D, S ;
	surf2Dread( &T_L, Traveltime_A, ptrMap.x-sizeof(float), ptrMap.y);
	surf2Dread( &T_R, Traveltime_A, ptrMap.x+sizeof(float), ptrMap.y);
	surf2Dread( &T_D, Traveltime_A, ptrMap.x, ptrMap.y-1);
	surf2Dread( &T_U, Traveltime_A, ptrMap.x, ptrMap.y+1);

	surf2Dread( &S, SlownessMap, ptrMap.x, ptrMap.y);	
    float Tx = fminf( T_L, T_R), Ty = fminf( T_D, T_U);   
	
	// hx = hy Version, already included in S, (adding and removing dummy vertices)
	// Let Tx > Ty
	float temp ;
	if(Tx < Ty){ temp = Tx ; Tx = Ty ; Ty = temp ; }

	float TT = Ty + S ;
	if( 1.414213f*S > Tx-Ty ) 
	{	
		temp = 0.5f*( (Tx+Ty) + sqrtf( 2.0f*powf(S,2) - powf(Tx-Ty,2) ) ) ;
		if(temp > Tx){ TT = temp; }
	}
	

return TT ; }

/// Compare New and Old traveltime
__device__ void compareTraveltime(bool *isUpdated)
{
	int2 ptrMap = make_int2( sizeof(float)*(blockIdx.x*blockDim.x+threadIdx.x + radius_X), 
	                                        blockIdx.y*blockDim.y+threadIdx.y + radius_Y)  ; 
    float TT ;
	surf2Dread( &TT, Traveltime_A, ptrMap.x, ptrMap.y);

    	
	float newTT = newCandidateOfTraveltime(ptrMap) ;
	if( TT-newTT > 1e-10*newTT ) // Update, if less and relative change is more than 1e-10.
	{
	    TT = newTT ;
	    *isUpdated = true ;
	}
    
	// Write outTraveltime
	surf2Dwrite(TT, Traveltime_B, ptrMap.x, ptrMap.y);

}


//**************** Global CUDA ****************//

__global__ void cuda_PSPS_Solver(int *Running)
{
	
	// Define shared variables (1)
	__shared__ bool isAnyUpdated ;
	__shared__ bool needUpdate  ;


	/// Need update ?
	if( threadIdx.x==0 && threadIdx.y==0 ){
        needUpdate   = false ;
        isAnyUpdated = false ;
	}
	__syncthreads();
	if( threadIdx.x<3 && threadIdx.y<3 && (threadIdx.x==1 || threadIdx.y==1) ){
		int data ; 
		surf2Dread(&data, UpdateMap_A, sizeof(int)*(blockIdx.x+threadIdx.x), blockIdx.y+threadIdx.y);
	    if( data==1 ){
	        needUpdate = true ;
	    }
	}
	__syncthreads();
	if( !needUpdate ){ return; }
	
	
	/// Compare New and Old traveltime
	compareTraveltime(&isAnyUpdated);
	__syncthreads();
	
	
/*
    // For no update map scheme ### don't forget to set isAnyUpdated = false;
    if( isAnyUpdated ){    
        ++(*Running) ;
        atomicAdd(Running,1);
    }
*/


	/// IsUpdate ? -> update running and UpdateMap.
	if( threadIdx.x==0 && threadIdx.y==0 )
	{
        if( isAnyUpdated ){
            ++(*Running) ;
            //atomicAdd(Running,1);
        }
        surf2Dwrite( (isAnyUpdated) ? 1:0 , UpdateMap_B, sizeof(int)*(blockIdx.x+1), blockIdx.y+1); 
	}

}


//**************** End CUDA_FIM ****************//


bool PSPS_CUDA_FIM_2D_Solver::Compute()
{
   /// Create Compute Shader (Solver) ///
   // Compute some useful numbers.
   for(unsigned int i=0 ; i<DIM ; ++i)
   {
      length[i] = res[i]+block[i] ;
   }
   side = block[0]/2 ;
   front = length[0]*block[1]/2  ;
   
   //******************************************************************************************
   
   /// Add dummy vertices ONLY 2D ///
   addDummyVertices();

	// Set constant values	
	int tempRa[2] = {(int)block[0]/2,(int)block[1]/2} ;
    cudaMemcpyToSymbol(radius_X, &tempRa[0], sizeof(int));
	cudaMemcpyToSymbol(radius_Y, &tempRa[1], sizeof(int));

   // Set Maps
   cudaChannelFormatDesc channelDescFloat = cudaCreateChannelDesc(32,0,0,0, cudaChannelFormatKindFloat); 
   cudaChannelFormatDesc channelDescInt   = cudaCreateChannelDesc(32,0,0,0, cudaChannelFormatKindSigned); 
	   
   // Slowness
   cudaArray* cuArraySlowness ;
   cudaMallocArray(&cuArraySlowness, &channelDescFloat, length[0], length[1], cudaArraySurfaceLoadStore);
   cudaMemcpyToArray(cuArraySlowness, 0, 0, &slowness[0], length[0]*length[1]*sizeof(float), cudaMemcpyHostToDevice);
   cudaBindSurfaceToArray(SlownessMap, cuArraySlowness);
   // Traveltime x2
   cudaArray *cuArrayTraveltime_A, *cuArrayTraveltime_B ;
   cudaMallocArray(&cuArrayTraveltime_A, &channelDescFloat, length[0], length[1], cudaArraySurfaceLoadStore);
   cudaMallocArray(&cuArrayTraveltime_B, &channelDescFloat, length[0], length[1], cudaArraySurfaceLoadStore);
   cudaMemcpyToArray(cuArrayTraveltime_A, 0, 0, &traveltime[0], length[0]*length[1]*sizeof(float), cudaMemcpyHostToDevice);
   cudaMemcpyToArray(cuArrayTraveltime_B, 0, 0, &traveltime[0], length[0]*length[1]*sizeof(float), cudaMemcpyHostToDevice); 
   //cudaBindSurfaceToArray(Traveltime_A, cuArrayTraveltime_A);
   //cudaBindSurfaceToArray(Traveltime_B, cuArrayTraveltime_B);

   // UpdateMap x2
   cudaArray *cuArrayUpdateMap_A, *cuArrayUpdateMap_B ;
   cudaMallocArray(&cuArrayUpdateMap_A, &channelDescInt, nWorkGroup[0]+2, nWorkGroup[1]+2, cudaArraySurfaceLoadStore);
   cudaMallocArray(&cuArrayUpdateMap_B, &channelDescInt, nWorkGroup[0]+2, nWorkGroup[1]+2, cudaArraySurfaceLoadStore);
   cudaMemcpyToArray(cuArrayUpdateMap_A, 0, 0, &updateMap[0], (nWorkGroup[0]+2)*(nWorkGroup[1]+2)*sizeof(int), cudaMemcpyHostToDevice);
   //cudaBindSurfaceToArray(UpdateMap_A, cuArrayUpdateMap_A);
   //cudaBindSurfaceToArray(UpdateMap_B, cuArrayUpdateMap_B);

   //******************************************************************************************

   // Set status
   int running = 0 ; 
   int *cuRunning ;
   cudaMalloc(&cuRunning, sizeof(int));
   cudaMemcpy(cuRunning, &running, sizeof(int), cudaMemcpyHostToDevice);
   bool AB = true ;
   
   { ///***** =.\= Parallel Shortest Path Solver (CUDA_FIM) =.\= *****///

      // Current Time
      time_t rawtime ;
      
      // Start
      std::cout << "\nPSPS_CUDA_FIM >> Start Fast Iterative Method (CUDA_FIM).\n" ;
      time(&rawtime);
      std::cout << "   START:: " << ctime(&rawtime) ;

      // Start timer
      auto start_time = std::chrono::high_resolution_clock::now();
      
      ///*** Main loop ***///
      int loop ; 
	  dim3 dimBlock( block[0], block[1] ); 
      dim3 dimGroup( nWorkGroup[0], nWorkGroup[1] );
      
	  // swap==0  A-> B
      // swap==1  B-> A
	  for(loop=0 ; loop<maxLoop ; ++loop)
      {		   
         /// AB ///
		 cudaBindSurfaceToArray(Traveltime_A, cuArrayTraveltime_A);
         cudaBindSurfaceToArray(Traveltime_B, cuArrayTraveltime_B);
		 cudaBindSurfaceToArray(UpdateMap_A, cuArrayUpdateMap_A);
		 cudaBindSurfaceToArray(UpdateMap_B, cuArrayUpdateMap_B);
         /***( Solving =/.= )***/
		 cuda_PSPS_Solver<<<dimGroup,dimBlock>>>(cuRunning);

		 // Read status 
		 cudaMemcpy(&running, cuRunning, sizeof(int), cudaMemcpyDeviceToHost);
         //std::cout << "  Loop = " << loop << "  Running = " << running << " \n" ; 
         if( running==0 ){ isFinish = true ; AB = false ; break; }
         running = 0 ; 
		 cudaMemcpy(cuRunning, &running, sizeof(int), cudaMemcpyHostToDevice);
        
         ++loop;
         
         /// BA ///
		 cudaBindSurfaceToArray(Traveltime_A, cuArrayTraveltime_B);
         cudaBindSurfaceToArray(Traveltime_B, cuArrayTraveltime_A);
		 cudaBindSurfaceToArray(UpdateMap_A, cuArrayUpdateMap_B);
		 cudaBindSurfaceToArray(UpdateMap_B, cuArrayUpdateMap_A);
         /***( Solving =/.= )***/
		 cuda_PSPS_Solver<<<dimGroup,dimBlock>>>(cuRunning);

		 // Read status 
		 cudaMemcpy(&running, cuRunning, sizeof(int), cudaMemcpyDeviceToHost);
         //std::cout << "  Loop = " << loop << "  Running = " << running << " \n" ; 
         if( running==0 ){ isFinish = true ; AB = true ; break; }
         running = 0 ; 
		 cudaMemcpy(cuRunning, &running, sizeof(int), cudaMemcpyHostToDevice);
      }

      // Stop timer
      auto end_time = std::chrono::high_resolution_clock::now();
      
      // Finish
      time(&rawtime);
      std::cout << "  FINISH:: " << ctime(&rawtime) ;
      std::cout << "PSPS_CUDA_FIM >> Finish Fast Iterative Method (CUDA_FIM).\n\n" ;
              
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
   
   // Clear status
   cudaFree(cuRunning);   
   
   //******************************************************************************************

    /// Retrieve Maps ///
	if( !AB ){ 
		cudaMemcpyFromArray(&traveltime[0], cuArrayTraveltime_B, 0, 0, length[0]*length[1]*sizeof(float), cudaMemcpyDeviceToHost);
	}else{
		cudaMemcpyFromArray(&traveltime[0], cuArrayTraveltime_A, 0, 0, length[0]*length[1]*sizeof(float), cudaMemcpyDeviceToHost);
	}
	if( !isFinish )
	{
		if( !AB ){ 
			cudaMemcpyFromArray(&updateMap[0], cuArrayUpdateMap_B, 0, 0, (nWorkGroup[0]+2)*(nWorkGroup[1]+2)*sizeof(int), cudaMemcpyDeviceToHost);
		}else{
			cudaMemcpyFromArray(&updateMap[0], cuArrayUpdateMap_A, 0, 0, (nWorkGroup[0]+2)*(nWorkGroup[1]+2)*sizeof(int), cudaMemcpyDeviceToHost);
		}	
	}

   /// Remove dummy vertices ONLY 2D ///
   removeDummyVertices();

   /// Clear Shader and Buffers ///
   cudaFreeArray(cuArraySlowness); 
   cudaFreeArray(cuArrayTraveltime_A);
   cudaFreeArray(cuArrayTraveltime_B);
   cudaFreeArray(cuArrayUpdateMap_A);
   cudaFreeArray(cuArrayUpdateMap_B);
   
return true; }


/************************************* Private *************************************/

void PSPS_CUDA_FIM_2D_Solver::addDummyVertices()
{
   unsigned int at ;
   bool isNegative = false ;
   std::vector<TyPe> temp ;

   /// Slowness map.
   // Copy data
   temp.assign( slowness.begin(), slowness.end() ) ;
   // Reinput the data
   slowness.assign( length[0]*length[1], 1.0f/0.0f) ;
   at = front+side ;
   for(unsigned int j=0 ; j<res[1] ; ++j)
   {
      for(unsigned int i=0 ; i<res[0] ; ++i)
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
   for(unsigned int j=0 ; j<res[1] ; ++j)
   {
      for(unsigned int i=0 ; i<res[0] ; ++i)
      {
         traveltime[ at+i ] = temp[ res[0]*j+i ] ;
      }
      at += res[0]+2*side ;
   }

   /// UpdateMap.
   // Copy data
   std::vector<int> tem ;
   tem.assign( updateMap.begin(), updateMap.end() ) ;
   // Reinput the data
   updateMap.assign( (nWorkGroup[0]+2)*(nWorkGroup[1]+2), 0) ;
   at = (nWorkGroup[0]+2)+1 ;
   for(unsigned int j=0 ; j<nWorkGroup[1] ; ++j)
   {
      for(unsigned int i=0 ; i<nWorkGroup[0] ; ++i)
      {
         updateMap[ at+i ] = tem[ nWorkGroup[0]*j+i ] ;
      }
      at += nWorkGroup[0]+2 ;
   }

return; }

void PSPS_CUDA_FIM_2D_Solver::removeDummyVertices()
{
   unsigned int at ;
   std::vector<TyPe> temp ;

   /// Slowness map
   // Copy data
   temp.assign( slowness.begin(), slowness.end() ) ;
   // Reinput the data
   slowness.resize( res[0]*res[1] ) ;
   at = front + side ;
   for(unsigned int j=0 ; j<res[1] ; ++j)
   {
      for(unsigned int i=0 ; i<res[0] ; ++i)
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
   for(unsigned int j=0 ; j<res[1] ; ++j)
   {
      for(unsigned int i=0 ; i<res[0] ; ++i)
      {
         traveltime[ res[0]*j+i ] = temp[ at+i ] ;
      }
      at += res[0]+2*side ;
   }

   /// UpdateMap.
   if( isFinish ){ return ; } //No more update map
   // Copy data
   std::vector<int> tem ;
   tem.assign( updateMap.begin(), updateMap.end() ) ;
   // Reinput the data
   updateMap.resize( nWorkGroup[0]*nWorkGroup[1] ) ;
   at = (nWorkGroup[0]+2)+1 ;
   for(unsigned int j=0 ; j<nWorkGroup[1] ; ++j)
   {
      for(unsigned int i=0 ; i<nWorkGroup[0] ; ++i)
      {
         updateMap[ nWorkGroup[0]*j+i ] = tem[ at+i ] ;
      }
      at += nWorkGroup[0]+2 ;
   }

return; }
