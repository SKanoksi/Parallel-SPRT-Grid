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

#include "PSPS_CUDA_ExtendedBF_2D_Solver.h"

// CUDA runtime
#include <cuda_runtime.h>

PSPS_CUDA_ExtendedBF_2D_Solver::PSPS_CUDA_ExtendedBF_2D_Solver()
{

}

PSPS_CUDA_ExtendedBF_2D_Solver::~PSPS_CUDA_ExtendedBF_2D_Solver()
{

}

//**************** CUDA_ExtendedBF program ****************//

__constant__ int shared_X ;
__constant__ int shared_Y ;
__constant__ int group_X  ;
__constant__ int group_Y  ;
__constant__ int radius_X ;
__constant__ int radius_Y ;

surface<void, cudaSurfaceType2D>  SlownessMap, Traveltime_A, Traveltime_B ;
surface<void, cudaSurfaceType2D>  UpdateMap_A,  UpdateMap_B, RaypathMap ;

// Calculate Weight (Slowness*LineWeight)
__device__ float SPR(const int ix, const int iy, const int2 Source, const float *Slowness) // Vertex = relative directional vector wrt. source(posOnShmem)
{
	int nX = abs(ix), nY = abs(iy) ;
	const float2  fVec   = make_float2(ix,iy) ;
	const float2  stride = make_float2( 1.0f/nX, 1.0f/nY ) ;
	
	float   Weight = 0 ;
	float2  cross  = make_float2( 1.0f-0.5f*stride.x, 1.0f-0.5f*stride.y ) ;
	float crossOld = 1.0f ;

	// Main process
	while( 0 < nX || 0 < nY )
	{
		// Find the next largest crossing point wrt. crossOld
		if( cross.x<cross.y ){

			// Compute weight.
			int indexX = Source.x + __float2int_rn( 0.5f*(crossOld+cross.y)*fVec.x ) ;
			int indexY = Source.y + __float2int_rn( 0.5f*(crossOld+cross.y)*fVec.y ) ;
			Weight += Slowness[ indexY*shared_X + indexX ]*(crossOld-cross.y) ;

			// For next crossing point.
			crossOld = cross.y ;
			--nY ;
			cross.y -= stride.y ;
		
		}else{
			
			// Compute weight.
			int indexX = Source.x + __float2int_rn( 0.5f*(crossOld+cross.x)*fVec.x ) ;
			int indexY = Source.y + __float2int_rn( 0.5f*(crossOld+cross.x)*fVec.y ) ;
			Weight += Slowness[ indexY*shared_X + indexX ]*(crossOld-cross.x) ;

			// For next crossing point.
			crossOld = cross.x ;
			--nX ;
			cross.x -= stride.x ;
			
		}
	}
	// Weight at source
	Weight += Slowness[Source.y*shared_X + Source.x]*crossOld ;
    Weight *= sqrt( powf(__int2float_rn(ix),2) + powf(__int2float_rn(iy),2) )  ;

return Weight ;}


/// Calculate new traveltime 
__device__ float FIM(const float Tx, const float Ty, const float S)
{
return 0.5f*( Tx+Ty + sqrt(2.0f*powf(S,2) - powf(Tx-Ty,2)) ) ; } // If it is NaN, then newTt<Tt is false anyway.

// Upload Slowness x4 and Traveltime x4 to shared memory
__device__ void uploadShared(float *Slowness, float *Traveltime)
{
	float temp ;
	int2  Stride = make_int2( 2*radius_X, 2*radius_Y ) ;

	// (1,1)
	int2    ptrMap = make_int2( sizeof(float)*(blockIdx.x*blockDim.x + threadIdx.x), blockIdx.y*blockDim.y + threadIdx.y);
	int2 ptrShared = make_int2( threadIdx.x, threadIdx.y ) ;
	surf2Dread(&temp,SlownessMap, ptrMap.x, ptrMap.y);
	Slowness[ ptrShared.y*shared_X + ptrShared.x ] = temp ;
	surf2Dread(&temp,Traveltime_A, ptrMap.x, ptrMap.y);
	Traveltime[ ptrShared.y*shared_X + ptrShared.x ] = temp ;

	// (2,1)
       ptrMap.x += sizeof(float)*Stride.x ;
	ptrShared.x += Stride.x ;
	surf2Dread(&temp,SlownessMap, ptrMap.x, ptrMap.y);
	Slowness[ ptrShared.y*shared_X + ptrShared.x ] = temp ;
	surf2Dread(&temp,Traveltime_A, ptrMap.x, ptrMap.y);
	Traveltime[ ptrShared.y*shared_X + ptrShared.x ] = temp ;

	 // (2,2)
	   ptrMap.y += Stride.y ;
	ptrShared.y += Stride.y ;
	surf2Dread(&temp,SlownessMap, ptrMap.x, ptrMap.y);
	Slowness[ ptrShared.y*shared_X + ptrShared.x ] = temp ;
	surf2Dread(&temp, Traveltime_A, ptrMap.x, ptrMap.y);
	Traveltime[ ptrShared.y*shared_X + ptrShared.x ] = temp ;

	 // (1,2)
       ptrMap.x -= sizeof(float)*Stride.x ;
	ptrShared.x -= Stride.x ;
	surf2Dread(&temp,SlownessMap, ptrMap.x, ptrMap.y);
	Slowness[ ptrShared.y*shared_X + ptrShared.x ] = temp ;
	surf2Dread(&temp,Traveltime_A, ptrMap.x, ptrMap.y);
	Traveltime[ ptrShared.y*shared_X + ptrShared.x ] = temp ;

}

/// Compare New and Old traveltime
__device__ void compareTraveltime(const float *Traveltime, const float *Slowness, bool *isUpdated)
{
	const int2 posWorker =  make_int2( threadIdx.x+radius_X, threadIdx.y+radius_Y) ; // posWorker on shMem ###
	float Tt = Traveltime[posWorker.y*shared_X + posWorker.x] ;
	int rayPath = -100 ;
	const int stride = int(group_X) ;
	const int2 shift = make_int2(radius_X, radius_Y);


	// Start:: Comparing neighbor != (0,0) ***************************************
    
    float T1, T2, T3, T4, newTt ;
        
    // Node 1 (-2,-1)
    float Tt1 = SPR(-2,-1,posWorker,Slowness) ; // is also node 24 
    newTt = Tt1 + Traveltime[ posWorker.x-2 + (posWorker.y-1)*shared_X ] ; 
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y-1)+(shift.x-2) ;
	}
	
	// Node 2 (-1,-1)
	T1 = SPR(-1,-1,posWorker,Slowness) ;
	newTt = T1 + Traveltime[ posWorker.x-1 + (posWorker.y-1)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y-1)+(shift.x-1) ;
	}
	
	// Node 3 (-1,-2)
	T2 = SPR(-1,-2,posWorker,Slowness) ;
	newTt = T2 + Traveltime[ posWorker.x-1 + (posWorker.y-2)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y-2)+(shift.x-1) ;
	}
	
	// Node 4 (-2,-2)
	newTt = FIM(T2,Tt1, Slowness[ posWorker.x-2 + (posWorker.y-2)*shared_X ] )
	        + Traveltime[ posWorker.x-2 + (posWorker.y-2)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y-2)+(shift.x-2) ;
	}
	newTt = SPR(-2,-2,posWorker,Slowness)  + Traveltime[ posWorker.x-2 + (posWorker.y-2)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y-2)+(shift.x-2) ;
	}
	
	// Node 5 (0,-1)
	T3 = SPR(0,-1,posWorker,Slowness) ;
	newTt = T3 + Traveltime[ posWorker.x + (posWorker.y-1)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y-1)+(shift.x) ;
	}
	
	// Node 6 (1,-2)
	T1 = SPR(1,-2,posWorker,Slowness) ;
	newTt = T1 + Traveltime[ posWorker.x+1 + (posWorker.y-2)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y-2)+(shift.x+1) ;
	}
	
	// Node 7 (0,-2)
	newTt = FIM(min(T1,T2),T3, Slowness[ posWorker.x + (posWorker.y-2)*shared_X ] )
	        + Traveltime[ posWorker.x + (posWorker.y-2)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y-2)+(shift.x) ;
	}
	newTt = SPR(0,-2,posWorker,Slowness) + Traveltime[ posWorker.x + (posWorker.y-2)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y-2)+(shift.x) ;
	}
		
	// Node 8 (1,-1)
	T4 = SPR(1,-1,posWorker,Slowness) ;
	newTt = T4 + Traveltime[ posWorker.x+1 + (posWorker.y-1)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y-1)+(shift.x+1) ;
	}
	
	// Node 9 (2,-1)
	T3 = SPR(2,-1,posWorker,Slowness) ;
	newTt = T3 + Traveltime[ posWorker.x+2 + (posWorker.y-1)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y-1)+(shift.x+2) ;
	}
	
	// Node 10 (2,-2)
	newTt = FIM(T1,T3, Slowness[ posWorker.x+2 + (posWorker.y-2)*shared_X ] )
	        + Traveltime[ posWorker.x+2 + (posWorker.y-2)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y-2)+(shift.x+2) ;
	}
	newTt = SPR(2,-2,posWorker,Slowness)  + Traveltime[ posWorker.x+2 + (posWorker.y-2)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y-2)+(shift.x+2) ;
	}
	
	// Node 11 (1,0)
	T1 = SPR(1,0,posWorker,Slowness) ;
	newTt = T1 + Traveltime[ posWorker.x+1 + (posWorker.y)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y)+(shift.x+1) ;
	}
	
	// Node 12 (2,1)
	T4 = SPR(2,1,posWorker,Slowness) ;
	newTt = T4 + Traveltime[ posWorker.x+2 + (posWorker.y+1)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y+1)+(shift.x+2) ;
	}
	
	// Node 13 (2,0)
	newTt = FIM(T1,min(T3,T4), Slowness[ posWorker.x+2 + (posWorker.y)*shared_X ] )
	        + Traveltime[ posWorker.x+2 + (posWorker.y)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y)+(shift.x+2) ;
	}
	newTt = SPR(2,0,posWorker,Slowness)  + Traveltime[ posWorker.x+2 + (posWorker.y)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y)+(shift.x+2) ;
	}
	
	// Node 14 (1,1)
	T2 = SPR(1,1,posWorker,Slowness) ;
	newTt = T2 + Traveltime[ posWorker.x+1 + (posWorker.y+1)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y+1)+(shift.x+1) ;
	}
	
	// Node 15 (1,2)
	T2 = SPR(1,2,posWorker,Slowness) ;
	newTt = T2 + Traveltime[ posWorker.x+1 + (posWorker.y+2)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y+2)+(shift.x+1) ;
	}
	
	// Node 16 (2,2)
	newTt = FIM(T2,T4, Slowness[ posWorker.x+2 + (posWorker.y+2)*shared_X ] )
	        + Traveltime[ posWorker.x+2 + (posWorker.y+2)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y+2)+(shift.x+2) ;
	}
	newTt = SPR(2,2,posWorker,Slowness)  + Traveltime[ posWorker.x+2 + (posWorker.y+2)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y+2)+(shift.x+2) ;
	}
	
	// Node 17 (0,1)
	T3 = SPR(0,1,posWorker,Slowness) ;
	newTt = T3 + Traveltime[ posWorker.x + (posWorker.y+1)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y+1)+(shift.x) ;
	}
	
	// Node 18 (-1,2)
	T1 = SPR(-1,2,posWorker,Slowness) ;
	newTt = T1 + Traveltime[ posWorker.x-1 + (posWorker.y+2)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y+2)+(shift.x-1) ;
	}
	
	// Node 19 (0,2)
	newTt = FIM(min(T1,T2),T3, Slowness[ posWorker.x + (posWorker.y+2)*shared_X ] )
	        + Traveltime[ posWorker.x + (posWorker.y+2)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y+2)+(shift.x) ;
	}
	newTt = SPR(0,2,posWorker,Slowness)  + Traveltime[ posWorker.x + (posWorker.y+2)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y+2)+(shift.x) ;
	}
	
	// Node 20 (-1,1)
	T2 = SPR(-1,1,posWorker,Slowness) ;
	newTt = T2 + Traveltime[ posWorker.x-1 + (posWorker.y+1)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y+1)+(shift.x-1) ;
	}
	
	// Node 21 (-2,1)
	T3 = SPR(-2,1,posWorker,Slowness) ;
	newTt = T3 + Traveltime[ posWorker.x-2 + (posWorker.y+1)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y+1)+(shift.x-2) ;
	}
	
	// Node 22 (-2,2)
	newTt = FIM(T1,T3, Slowness[ posWorker.x-2 + (posWorker.y+2)*shared_X ] )
	        + Traveltime[ posWorker.x-2 + (posWorker.y+2)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y+2)+(shift.x-2) ;
	}
	newTt = SPR(-2,2,posWorker,Slowness)  + Traveltime[ posWorker.x-2 + (posWorker.y+2)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y+2)+(shift.x-2) ;
	}
	
	// Node 23 (-1,0)
	T1 = SPR(-1,0,posWorker,Slowness) ;
	newTt = T1 + Traveltime[ posWorker.x-1 + (posWorker.y)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y)+(shift.x-1) ;
	}
	
	// Node 24==1 (-2,-1)
	
	// Node 25 (-2,0)
	newTt = FIM(T1,min(T3,Tt1), Slowness[ posWorker.x-2 + (posWorker.y)*shared_X ] )
	        + Traveltime[ posWorker.x-2 + (posWorker.y)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y)+(shift.x-2) ;
	}
	newTt = SPR(-2,0,posWorker,Slowness) + Traveltime[ posWorker.x-2 + (posWorker.y)*shared_X ] ;
    if( newTt<Tt ){
		Tt = newTt ;
		rayPath = stride*(shift.y)+(shift.x-2) ;
	}

    // End:: Comparing neighbor != (0,0) ***************************************


	// Write outTraveltime
    int GlobalID_X = blockIdx.x*blockDim.x + threadIdx.x ;
    int GlobalID_Y = blockIdx.y*blockDim.y + threadIdx.y ;
    surf2Dwrite(Tt, Traveltime_B, sizeof(float)*(GlobalID_X+radius_X), GlobalID_Y+radius_Y);
		
	// if threadUpdate -> isUpdated[0] = true ;
	if( rayPath != -100 ){
      *isUpdated = true ;
	  surf2Dwrite(rayPath, RaypathMap, sizeof(int)*GlobalID_X, GlobalID_Y);
    }

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
	if( threadIdx.x<3 && threadIdx.y<3 ){
		int data ; 
		surf2Dread(&data, UpdateMap_A, sizeof(int)*(blockIdx.x+threadIdx.x), blockIdx.y+threadIdx.y);
	    if( data==1 ){
	        needUpdate = true ;
	    }
	}
	__syncthreads();
	if( !needUpdate ){ return; }



    // Define shared variables (2)
	extern __shared__ float sharedVar[] ;
	float *Slowness   = &sharedVar[0] ;  
	float *Traveltime = &sharedVar[shared_X*shared_Y] ; 
    //isAnyUpdated = false ; //%%%
	
    /// Upload Slowness x4 and Traveltime x4 to shared memory	
    uploadShared(Slowness,Traveltime); 
	__syncthreads();

	/// Compare New and Old traveltime
	compareTraveltime(Traveltime,Slowness,&isAnyUpdated);
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

//**************** End CUDA_ExtendedBF ****************//


bool PSPS_CUDA_ExtendedBF_2D_Solver::Compute()
{
   /// Create Compute Shader (Solver) ///
   // Compute some useful numbers.
   unsigned int shared[DIM] ;
   for(unsigned int i=0 ; i<DIM ; ++i)
   {
      length[i] = res[i]+block[i] ;
      shared[i] = 2*block[i] ;
   }
   side = block[0]/2 ;
   front = length[0]*block[1]/2  ;
   
   //******************************************************************************************
   
   /// Add dummy vertices ONLY 2D ///
   addDummyVertices();

	// Set constant values
	int tempSh[2] = {(int)shared[0],(int)shared[1]} ;
	cudaMemcpyToSymbol(shared_X, &tempSh[0], sizeof(int));
	cudaMemcpyToSymbol(shared_Y, &tempSh[1], sizeof(int));
	
	int tempGr[2] = {(int)block[0]+1,(int)block[1]+1} ;
	cudaMemcpyToSymbol(group_X, &tempGr[0], sizeof(int));
	cudaMemcpyToSymbol(group_Y, &tempGr[1], sizeof(int));
	
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
   // Raypath 
   cudaArray* cuArrayRaypath ;
   cudaMallocArray(&cuArrayRaypath, &channelDescInt, res[0], res[1], cudaArraySurfaceLoadStore);
   cudaMemcpyToArray(cuArrayRaypath, 0, 0, &raypath[0], res[0]*res[1]*sizeof(int), cudaMemcpyHostToDevice);
   cudaBindSurfaceToArray(RaypathMap, cuArrayRaypath);
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
   
   { ///***** =.\= Parallel Shortest Path Solver (CUDA_ExtendedBF) =.\= *****///

      // Current Time
      time_t rawtime ;
      
      // Start
      std::cout << "\nPSPS_CUDA_ExtendedBF >> Start \'Extended\' Bellman-Ford (CUDA_ExtendedBF).\n" ;
      time(&rawtime);
      std::cout << "   START:: " << ctime(&rawtime) ;

      // Start timer
      auto start_time = std::chrono::high_resolution_clock::now();
      
      ///*** Main loop ***///
      int loop ; 
	  dim3 dimBlock( block[0], block[1] ); 
      dim3 dimGroup( nWorkGroup[0], nWorkGroup[1] );
	  size_t nSharedVar = 2*shared[0]*shared[1]*sizeof(float) ;
      
	  // Running
	  for(loop=0 ; loop<maxLoop ; ++loop)
      {		   
         /// AB ///
	     cudaBindSurfaceToArray(Traveltime_A, cuArrayTraveltime_A);
         cudaBindSurfaceToArray(Traveltime_B, cuArrayTraveltime_B);
		 cudaBindSurfaceToArray(UpdateMap_A, cuArrayUpdateMap_A);
		 cudaBindSurfaceToArray(UpdateMap_B, cuArrayUpdateMap_B);
         /***( Solving =/.= )***/
		 cuda_PSPS_Solver<<<dimGroup,dimBlock,nSharedVar>>>(cuRunning);

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
		 cuda_PSPS_Solver<<<dimGroup,dimBlock,nSharedVar>>>(cuRunning);

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
      std::cout << "PSPS_CUDA_ExtendedBF >> Finish \'Extended\' Bellman-Ford (CUDA_ExtendedBF).\n\n" ;
      
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
    cudaMemcpyFromArray(&raypath[0], cuArrayRaypath, 0, 0, res[0]*res[1]*sizeof(float), cudaMemcpyDeviceToHost);
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
   cudaFreeArray(cuArrayRaypath);
   cudaFreeArray(cuArrayTraveltime_A);
   cudaFreeArray(cuArrayTraveltime_B);
   cudaFreeArray(cuArrayUpdateMap_A);
   cudaFreeArray(cuArrayUpdateMap_B);
   
return true; }


/************************************* Private *************************************/

void PSPS_CUDA_ExtendedBF_2D_Solver::addDummyVertices()
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

void PSPS_CUDA_ExtendedBF_2D_Solver::removeDummyVertices()
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
