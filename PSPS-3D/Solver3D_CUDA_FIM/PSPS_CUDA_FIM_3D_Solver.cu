/*************************************
    Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#define DIM 3
#define TyPe float

#include "PSPS_CUDA_FIM_3D_Solver.h"

PSPS_CUDA_FIM_3D_Solver::PSPS_CUDA_FIM_3D_Solver()
{

}

PSPS_CUDA_FIM_3D_Solver::~PSPS_CUDA_FIM_3D_Solver()
{

}

//**************** CUDA program ****************//

__constant__ int radius_X ;
__constant__ int radius_Y ;
__constant__ int radius_Z ;

surface<void, cudaSurfaceType3D>  SlownessMap, Traveltime_A, Traveltime_B ;
surface<void, cudaSurfaceType3D>  UpdateMap_A,  UpdateMap_B ;


/// Calculate new traveltime 
__device__ float newCandidateOfTraveltime(const int3 ptrMap)
{
	float T_L, T_R, T_U, T_D, T_F, T_B, S ;
	surf3Dread( &T_L, Traveltime_A, ptrMap.x-sizeof(float), ptrMap.y, ptrMap.z);
	surf3Dread( &T_R, Traveltime_A, ptrMap.x+sizeof(float), ptrMap.y, ptrMap.z);
	surf3Dread( &T_D, Traveltime_A, ptrMap.x, ptrMap.y-1, ptrMap.z);
	surf3Dread( &T_U, Traveltime_A, ptrMap.x, ptrMap.y+1, ptrMap.z);
	surf3Dread( &T_F, Traveltime_A, ptrMap.x, ptrMap.y, ptrMap.z-1);
	surf3Dread( &T_B, Traveltime_A, ptrMap.x, ptrMap.y, ptrMap.z+1);

	surf3Dread( &S, SlownessMap, ptrMap.x, ptrMap.y, ptrMap.z);	
    float Tx = fminf( T_L, T_R), Ty = fminf( T_D, T_U), Tz = fminf( T_F, T_B);   
	
	// hx = hy Version, already included in S (adding and removing dummy vertices)
	// Let Tx > Ty > Tz
	float temp ;
	if(Tx < Ty){ temp = Tx ; Tx = Ty ; Ty = temp ; }
	if(Ty < Tz){ temp = Ty ; Ty = Tz ; Tz = temp ; }
	if(Tx < Ty){ temp = Tx ; Tx = Ty ; Ty = temp ; }

	float TT = Tz + S ;
	if( 1.414213f*S > Ty-Tz) 
	{	
	    temp = 0.5f*( (Ty+Tz) + sqrt( 2.0f*pow(S,2) - pow(Ty-Tz,2) ) ) ;
	
		if(temp > Ty){ TT = temp; }else{ return TT ; }
		float tt ;
		if( 0 < ( tt=(3.0f*pow(S,2) -pow(Tx-Ty,2)-pow(Tx-Tz,2)-pow(Ty-Tz,2)) ) )
		{			
			temp = ( (Tx+Ty+Tz) + sqrt(tt) )/3.0f; 
			if(temp > Tx){ TT = temp ; }
		}
	}
	

return TT ; }

/// Compare New and Old traveltime
__device__ void compareTraveltime(bool *isUpdated)
{
	int3 ptrMap = make_int3( sizeof(float)*(blockIdx.x*blockDim.x+threadIdx.x + radius_X), 
	                                        blockIdx.y*blockDim.y+threadIdx.y + radius_Y, 
	                                        blockIdx.z*blockDim.z+threadIdx.z + radius_Z)  ; 
    float TT ;
    surf3Dread( &TT, Traveltime_A, ptrMap.x, ptrMap.y, ptrMap.z);
	
    	
	float newTT = newCandidateOfTraveltime(ptrMap) ;
	if( TT-newTT > 1e-10*newTT ) // Update, if less and relative change is more than 1e-10.
	{
	    TT = newTT ;
	    *isUpdated = true ;
	}
    
	// Write outTraveltime
	surf3Dwrite(TT, Traveltime_B, ptrMap.x, ptrMap.y, ptrMap.z);
	
	
}

//**************** Global CUDA ****************//
__global__ void cuda_PSPS_Solver(int *Running)
{
	
	// Define shared variables (1)
	__shared__ bool isAnyUpdated ;
	__shared__ bool needUpdate  ;


	/// Need update ?
	if( threadIdx.x==0 && threadIdx.y==0 && threadIdx.z==0 ){
        needUpdate   = false ;
        isAnyUpdated = false ;
	}
	__syncthreads();
	if( threadIdx.x<3 && threadIdx.y<3 && threadIdx.z<3 && 
	     ( ( threadIdx.x==1 && threadIdx.y==1 ) || 
	       ( threadIdx.x==1 && threadIdx.z==1 ) || 
	       ( threadIdx.y==1 && threadIdx.z==1 ) )
	){
		int data ; 
		surf3Dread(&data, UpdateMap_A, sizeof(int)*(blockIdx.x+threadIdx.x), blockIdx.y+threadIdx.y, blockIdx.z+threadIdx.z);
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
    // For no update map scheme  
    // ### don't forget to set isAnyUpdate = false and remove the __syncthreads() 
    if( isAnyUpdated ){    
//	    ++(*Running) ;
        atomicAdd(Running,1);
    }
*/


	/// IsUpdate ? -> update running and UpdateMap.
	if( threadIdx.x==0 && threadIdx.y==0 && threadIdx.z==0 )
	{
        if( isAnyUpdated ){
            ++(*Running) ;
            //atomicAdd(Running,1);
        }
	    surf3Dwrite( (isAnyUpdated) ? 1:0 , UpdateMap_B, sizeof(int)*(blockIdx.x+1), blockIdx.y+1, blockIdx.z+1);
	}


}


//**************** End CUDA ****************//

bool PSPS_CUDA_FIM_3D_Solver::Compute()
{
	/// Create Compute Shader (Solver) ///
   // Compute some useful numbers.
   for(unsigned int i=0 ; i<DIM ; ++i)
   {
      length[i] = res[i]+block[i] ;
   }
   side = block[0]/2 ;
   front = length[0]*block[1]/2 ;
   top = length[0]*length[1]*block[2]/2 ;
   
    //******************************************************************************************
    
   /// Add dummy vertices ONLY 3D ///
   addDummyVertices();

	// Set constant values	
	int tempRa[] = {(int)side, (int)block[1]/2, (int)block[2]/2} ;
	cudaMemcpyToSymbol(radius_X, &tempRa[0], sizeof(int));
	cudaMemcpyToSymbol(radius_Y, &tempRa[1], sizeof(int));
	cudaMemcpyToSymbol(radius_Z, &tempRa[2], sizeof(int));
	
	// Set Maps
   cudaChannelFormatDesc channelDescFloat = cudaCreateChannelDesc(32,0,0,0, cudaChannelFormatKindFloat); 
   cudaChannelFormatDesc channelDescInt   = cudaCreateChannelDesc(32,0,0,0, cudaChannelFormatKindSigned); 
	   
   // Slowness   
   cudaArray* cuArraySlowness ;
   cudaMalloc3DArray(&cuArraySlowness, &channelDescFloat, make_cudaExtent(length[0], length[1], length[2]), cudaArraySurfaceLoadStore);
   cudaMemcpy3DParms paramInS = {0} ;
	  paramInS.srcPtr = make_cudaPitchedPtr(&slowness[0], length[0]*sizeof(float), length[0], length[1]);
	  paramInS.dstArray = cuArraySlowness ;
      paramInS.kind =  cudaMemcpyHostToDevice ;
      paramInS.extent = make_cudaExtent( length[0], length[1], length[2]) ;
   cudaMemcpy3D(&paramInS);
   cudaBindSurfaceToArray(SlownessMap, cuArraySlowness);
   
   // Traveltime x2
   cudaArray *cuArrayTraveltime_A, *cuArrayTraveltime_B ;
   cudaMalloc3DArray(&cuArrayTraveltime_A, &channelDescFloat, make_cudaExtent( length[0], length[1], length[2]), cudaArraySurfaceLoadStore);
   cudaMalloc3DArray(&cuArrayTraveltime_B, &channelDescFloat, make_cudaExtent( length[0], length[1], length[2]), cudaArraySurfaceLoadStore);
   cudaMemcpy3DParms paramInT = {0} ;
	  paramInT.srcPtr = make_cudaPitchedPtr(&traveltime[0], length[0]*sizeof(float), length[0], length[1]) ;
	  paramInT.dstArray = cuArrayTraveltime_A ;
      paramInT.kind =  cudaMemcpyHostToDevice ;
      paramInT.extent = make_cudaExtent( length[0], length[1], length[2]) ;
   cudaMemcpy3D(&paramInT);
   	  paramInT.dstArray = cuArrayTraveltime_B ;
   cudaMemcpy3D(&paramInT);
   //cudaBindSurfaceToArray(Traveltime_A, cuArrayTraveltime_A);
   //cudaBindSurfaceToArray(Traveltime_B, cuArrayTraveltime_B);

   // UpdateMap x2
   cudaArray *cuArrayUpdateMap_A, *cuArrayUpdateMap_B ;
   cudaMalloc3DArray(&cuArrayUpdateMap_A, &channelDescInt, make_cudaExtent( nWorkGroup[0]+2, nWorkGroup[1]+2, nWorkGroup[2]+2), cudaArraySurfaceLoadStore);
   cudaMalloc3DArray(&cuArrayUpdateMap_B, &channelDescInt, make_cudaExtent( nWorkGroup[0]+2, nWorkGroup[1]+2, nWorkGroup[2]+2), cudaArraySurfaceLoadStore);
   cudaMemcpy3DParms paramInU = {0} ;
	  paramInU.srcPtr = make_cudaPitchedPtr(&updateMap[0], (nWorkGroup[0]+2)*sizeof(int), (nWorkGroup[0]+2), (nWorkGroup[1]+2)) ;
	  paramInU.dstArray = cuArrayUpdateMap_A ;
      paramInU.kind =  cudaMemcpyHostToDevice ;
      paramInU.extent = make_cudaExtent( nWorkGroup[0]+2, nWorkGroup[1]+2, nWorkGroup[2]+2) ;
   cudaMemcpy3D(&paramInU);
   //cudaBindSurfaceToArray(UpdateMap_A, cuArrayUpdateMap_A);
   //cudaBindSurfaceToArray(UpdateMap_B, cuArrayUpdateMap_B);

   //******************************************************************************************

    // Set status
   int running = 0 ; 
   int *cuRunning ;
   cudaMalloc(&cuRunning, sizeof(int));
   cudaMemcpy(cuRunning, &running, sizeof(int), cudaMemcpyHostToDevice);
   bool AB = true ;

   { ///***** =.\= Parallel Shortest Path Solver (CUDA) =.\= *****///

      ///*** Main loop ***///
      int loop ; 
	  dim3 dimBlock( block[0], block[1], block[2] ); 
      dim3 dimGroup( nWorkGroup[0], nWorkGroup[1], nWorkGroup[2] );

	  // Current Time
      time_t rawtime ;
      
      // Start
      std::cout << "\nPSPS_CUDA_FIM >> Start Fast iterative method (CUDA).\n" ;
      time(&rawtime);
      std::cout << "   START:: " << ctime(&rawtime) ;
      
      // Start timer
      auto start_time = std::chrono::high_resolution_clock::now();
	  
	  // swap==0  A-> B
      // swap==1  B-> A
	  for(loop=0 ; loop<maxLoop ; ++loop)
      {		   
         /// AB ///
		cudaBindSurfaceToArray(Traveltime_A, cuArrayTraveltime_A, channelDescFloat);
        cudaBindSurfaceToArray(Traveltime_B, cuArrayTraveltime_B, channelDescFloat);
        cudaBindSurfaceToArray(UpdateMap_A, cuArrayUpdateMap_A, channelDescInt);
        cudaBindSurfaceToArray(UpdateMap_B, cuArrayUpdateMap_B, channelDescInt);
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
        cudaBindSurfaceToArray(Traveltime_A, cuArrayTraveltime_B, channelDescFloat);
        cudaBindSurfaceToArray(Traveltime_B, cuArrayTraveltime_A, channelDescFloat);
        cudaBindSurfaceToArray(UpdateMap_A, cuArrayUpdateMap_B, channelDescInt);
        cudaBindSurfaceToArray(UpdateMap_B, cuArrayUpdateMap_A, channelDescInt);
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
      std::cout << "PSPS_CUDA_FIM >> Finish Fast iterative method (CUDA).\n\n" ;
        
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
	// Traveltime
    cudaMemcpy3DParms paramOutT = {0} ;
	  if( !AB ){
         paramOutT.srcArray = cuArrayTraveltime_B ; 
      }else{
		 paramOutT.srcArray = cuArrayTraveltime_A ;
	  }
	  paramOutT.dstPtr = make_cudaPitchedPtr(&traveltime[0], length[0]*sizeof(float), length[0], length[1]) ;
	  paramOutT.kind   = cudaMemcpyDeviceToHost ;
      paramOutT.extent = make_cudaExtent( length[0], length[1], length[2]) ;
    cudaMemcpy3D(&paramOutT);
	
	// UpdateMap
	if( !isFinish )
	{
		 cudaMemcpy3DParms paramOutU = {0} ;
			if( !AB ){
				paramOutU.srcArray = cuArrayUpdateMap_B ;
			}else{
				paramOutU.srcArray = cuArrayUpdateMap_A ;
			}
			paramOutU.dstPtr = make_cudaPitchedPtr(&updateMap[0], (nWorkGroup[0]+2)*sizeof(int), (nWorkGroup[0]+2), (nWorkGroup[1]+2)) ;
			paramOutU.kind =  cudaMemcpyDeviceToHost ;
			paramOutU.extent = make_cudaExtent( nWorkGroup[0]+2, nWorkGroup[1]+2, nWorkGroup[2]+2) ;
		cudaMemcpy3D(&paramOutU);
	}

   /// Remove dummy vertices ONLY 3D ///
   removeDummyVertices();

   /// Clear Shader and Buffers ///
   cudaFreeArray(cuArraySlowness); 
   cudaFreeArray(cuArrayTraveltime_A);
   cudaFreeArray(cuArrayTraveltime_B);
   cudaFreeArray(cuArrayUpdateMap_A);
   cudaFreeArray(cuArrayUpdateMap_B);

return true; }

/************************************* Private *************************************/

void PSPS_CUDA_FIM_3D_Solver::addDummyVertices()
{
   unsigned int at ;
   bool isNegative = false ;
   std::vector<TyPe> temp ;

   /// Slowness map.
   // Copy data
   temp.assign( slowness.begin(), slowness.end() ) ;
   // Reinput the data
   slowness.assign( length[0]*length[1]*length[2], 1.0f/0.0f) ;
   at = top+front+side ;
   for(unsigned int k=0 ; k<res[2] ; ++k)
   {
		for(unsigned int j=0 ; j<res[1] ; ++j)
		{
			for(unsigned int i=0 ; i<res[0] ; ++i)
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
   for(unsigned int k=0 ; k<res[2] ; ++k)
   {
		for(unsigned int j=0 ; j<res[1] ; ++j)
		{
			for(unsigned int i=0 ; i<res[0] ; ++i)
			{
				traveltime[ at+i ] = temp[ (k*res[1]+j)*res[0]+i ] ;
			}
			at += res[0]+2*side ;
		}
		at += 2*front ; 
   }

   /// UpdateMap.
   // Copy data
   std::vector<int> tem ;
   tem.assign( updateMap.begin(), updateMap.end() ) ;
   // Reinput the data
   updateMap.assign( (nWorkGroup[0]+2)*(nWorkGroup[1]+2)*(nWorkGroup[2]+2), 0) ;
   at = (nWorkGroup[0]+2)*(nWorkGroup[1]+2)+(nWorkGroup[0]+2)+1 ;
   for(unsigned int k=0 ; k<nWorkGroup[2] ; ++k)
   {
		for(unsigned int j=0 ; j<nWorkGroup[1] ; ++j)
		{
			for(unsigned int i=0 ; i<nWorkGroup[0] ; ++i)
			{
				updateMap[ at+i ] = tem[ (k*nWorkGroup[1]+j)*nWorkGroup[0]+i ] ;
			}
			at += nWorkGroup[0]+2 ;
		}
		at += 2*(nWorkGroup[0]+2) ;
   }

return; }

void PSPS_CUDA_FIM_3D_Solver::removeDummyVertices()
{
   unsigned int at ;
   std::vector<TyPe> temp ;

   /// Slowness map
   // Copy data
   temp.assign( slowness.begin(), slowness.end() ) ;
   // Reinput the data
   slowness.resize( res[0]*res[1]*res[2] ) ;
   at = top+front+side ;
   for(unsigned int k=0 ; k<res[2] ; ++k)
   {
		for(unsigned int j=0 ; j<res[1] ; ++j)
		{
			for(unsigned int i=0 ; i<res[0] ; ++i)
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
   for(unsigned int k=0 ; k<res[2] ; ++k)
   {
		for(unsigned int j=0 ; j<res[1] ; ++j)
		{
			for(unsigned int i=0 ; i<res[0] ; ++i)
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
   std::vector<int> tem ;
   tem.assign( updateMap.begin(), updateMap.end() ) ;
   // Reinput the data
   updateMap.resize( nWorkGroup[0]*nWorkGroup[1]*nWorkGroup[2] ) ;
   at = (nWorkGroup[0]+2)*(nWorkGroup[1]+2)+(nWorkGroup[0]+2)+1 ;
   for(unsigned int k=0 ; k<nWorkGroup[2] ; ++k)
   {
		for(unsigned int j=0 ; j<nWorkGroup[1] ; ++j)
		{
			for(unsigned int i=0 ; i<nWorkGroup[0] ; ++i)
			{
				updateMap[ (k*nWorkGroup[1]+j)*nWorkGroup[0]+i ] = tem[ at+i ] ;
			}
			at += nWorkGroup[0]+2 ;
		}
		at += 2*(nWorkGroup[0]+2) ;
   }
   

return; }
