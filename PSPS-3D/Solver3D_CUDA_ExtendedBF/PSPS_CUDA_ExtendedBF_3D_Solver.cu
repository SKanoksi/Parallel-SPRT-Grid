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

#include "PSPS_CUDA_ExtendedBF_3D_Solver.h"

PSPS_CUDA_ExtendedBF_3D_Solver::PSPS_CUDA_ExtendedBF_3D_Solver()
{

}

PSPS_CUDA_ExtendedBF_3D_Solver::~PSPS_CUDA_ExtendedBF_3D_Solver()
{

}

//**************** CUDA program ****************//

__constant__ int shared_X ;
__constant__ int shared_Y ;
__constant__ int shared_Z ;
__constant__ int group_X  ;
__constant__ int group_Y  ;
__constant__ int group_Z  ;
__constant__ int radius_X ;
__constant__ int radius_Y ;
__constant__ int radius_Z ;

surface<void, cudaSurfaceType3D>  SlownessMap, Traveltime_A, Traveltime_B ;
surface<void, cudaSurfaceType3D>  UpdateMap_A,  UpdateMap_B, RaypathMap ;

// Calculate Weight (Slowness*LineWeight): common function
__device__ float SPR(const float *Slowness, const int ix, const int iy, const int iz, const int3 Source) // Vertex = relative directional vector wrt. source(posOnShmem)
{
	      int3    nVec   = make_int3( abs(ix), abs(iy), abs(iz) ) ;
	const float3  fVec   = make_float3(ix, iy, iz) ;
	const float3  stride = make_float3( 1.0f/nVec.x, 1.0f/nVec.y, 1.0f/nVec.z) ;
	
    float   Weight = 0 ;
	float3  cross  = make_float3( 1.0f-0.5f*stride.x, 1.0f-0.5f*stride.y, 1.0f-0.5f*stride.z) ;
	float crossOld = 1.0f ;


	// Main process
	while( 0 < nVec.x || 0 < nVec.y || 0 < nVec.z )
	{
        int next = 2 ;
        if( cross.y<cross.x ){            
            if( cross.z<cross.x ){ next=0; }
        }else{
            if( cross.z<cross.y ){ next=1; }
        }

		// Find the next largest crossing point wrt. crossOld	 
		if( next==0 ){

			// Compute weight.
			int indexX = Source.x + __float2int_rn( 0.5f*(crossOld+cross.x)*fVec.x ) ;
			int indexY = Source.y + __float2int_rn( 0.5f*(crossOld+cross.x)*fVec.y ) ;
			int indexZ = Source.z + __float2int_rn( 0.5f*(crossOld+cross.x)*fVec.z ) ;
			Weight += Slowness[ (indexZ*shared_Y + indexY)*shared_X + indexX ]*(crossOld-cross.x) ;

			// For next crossing point.
			crossOld = cross.x ;
			--nVec.x ;
			cross.x -= stride.x ;
		
		}else{
		if( next==1 ){
				
			// Compute weight.
			int indexX = Source.x + __float2int_rn( 0.5f*(crossOld+cross.y)*fVec.x ) ;
			int indexY = Source.y + __float2int_rn( 0.5f*(crossOld+cross.y)*fVec.y ) ;
			int indexZ = Source.z + __float2int_rn( 0.5f*(crossOld+cross.y)*fVec.z ) ;
			Weight += Slowness[ (indexZ*shared_Y + indexY)*shared_X + indexX ]*(crossOld-cross.y) ;

			// For next crossing point.
			crossOld = cross.y ;
			--nVec.y ;
			cross.y -= stride.y ;
				
		}else{
        // if next == 2 
				
			// Compute weight.
			int indexX = Source.x + __float2int_rn( 0.5f*(crossOld+cross.z)*fVec.x ) ;
			int indexY = Source.y + __float2int_rn( 0.5f*(crossOld+cross.z)*fVec.y ) ;
			int indexZ = Source.z + __float2int_rn( 0.5f*(crossOld+cross.z)*fVec.z ) ;
			Weight += Slowness[ (indexZ*shared_Y + indexY)*shared_X + indexX ]*(crossOld-cross.z) ;

			// For next crossing point.
			crossOld = cross.z ;
			--nVec.z ;
			cross.z -= stride.z ;
				
		}
        }
			
	}
	
	// Weight at source
	Weight += Slowness[ (Source.z*shared_Y + Source.y)*shared_X + Source.x]*crossOld ;
    Weight *= sqrt( powf(__int2float_rn(ix), 2) + powf(__int2float_rn(iy), 2) + powf(__int2float_rn(iz), 2) ) ;

return Weight ;}

/// Calculate new traveltime 
__device__ float FIM(float Tx, float Ty, float Tz, const float S)
{     
    float TT = 3.0f*powf(S,2) -powf(Tx-Ty,2)-powf(Tx-Tz,2)-powf(Ty-Tz,2) ;
    if( TT < 0 )
    {
        // Let Tx = Largest (use TT as temp)
        if(Tx < Ty){ TT = Tx ; Tx = Ty ; Ty = TT ; } 
	    if(Tx < Tz){ TT = Tx ; Tx = Tz ; Tz = TT ; }	
	    TT = 0.5f*( Ty+Tz + sqrt(2.0f*powf(S,2) - powf(Ty-Tz,2)) ) ; 
	    
    }else{
 
        TT = ( Tx+Ty+Tz + sqrt(TT) )/3.0f ; 
 
    }
    
    	   
return TT ; } // If it is NaN, then newTt<Tt is false anyway.

// Greatest Common Divisor for edge reduction
__device__ int GCD(int a, int b)
{
    a = abs(a) ;
    b = abs(b) ;
    int c ;
    while( a!=0 )
    {
        c = a ;
        a = b%a ;
        b = c ;
    }

return b; }

__device__ void checkTraveltime(const float *Traveltime, float *Tt, int *rayPath, const int3 posWorker, int i, int j, int k, float weight)
{
    float newTt = weight + Traveltime[ posWorker.x+i +shared_X*(posWorker.y+j +shared_Y*(posWorker.z+k)) ] ; 
    if( newTt<*Tt ){
		*Tt = newTt ;
		*rayPath = ((radius_Z+k)*group_Y + radius_Y+j)*group_X + radius_X+i ;
	}

}

//**************** AB functions ****************//

// Upload Slowness x4 and Traveltime x4 to shared memory
__device__ void uploadShared(float *Slowness, float *Traveltime)
{
	float temp ;
	int3  Stride = make_int3( 2*radius_X, 2*radius_Y, 2*radius_Z ) ;
        
	// (1,1,1)
	int3    ptrMap = make_int3( sizeof(float)*(blockIdx.x*blockDim.x + threadIdx.x), blockIdx.y*blockDim.y + threadIdx.y,  blockIdx.z*blockDim.z + threadIdx.z);
	int3 ptrShared = make_int3( threadIdx.x, threadIdx.y, threadIdx.z ) ;
	surf3Dread(&temp,SlownessMap, ptrMap.x, ptrMap.y, ptrMap.z );
	Slowness[ (ptrShared.z*shared_Y + ptrShared.y)*shared_X + ptrShared.x ] = temp ;
	surf3Dread(&temp,Traveltime_A, ptrMap.x, ptrMap.y, ptrMap.z);
	Traveltime[ (ptrShared.z*shared_Y + ptrShared.y)*shared_X + ptrShared.x ] = temp ;

	// (2,1,1)
       ptrMap.x += sizeof(float)*Stride.x ;
	ptrShared.x += Stride.x ;
	surf3Dread(&temp,SlownessMap, ptrMap.x, ptrMap.y, ptrMap.z );
	Slowness[ (ptrShared.z*shared_Y + ptrShared.y)*shared_X + ptrShared.x ] = temp ;
	surf3Dread(&temp,Traveltime_A, ptrMap.x, ptrMap.y, ptrMap.z);
	Traveltime[ (ptrShared.z*shared_Y + ptrShared.y)*shared_X + ptrShared.x ] = temp ;

	 // (2,2,1)
	   ptrMap.y += Stride.y ;
	ptrShared.y += Stride.y ;
	surf3Dread(&temp,SlownessMap, ptrMap.x, ptrMap.y, ptrMap.z );
	Slowness[ (ptrShared.z*shared_Y + ptrShared.y)*shared_X + ptrShared.x ] = temp ;
	surf3Dread(&temp,Traveltime_A, ptrMap.x, ptrMap.y, ptrMap.z);
	Traveltime[ (ptrShared.z*shared_Y + ptrShared.y)*shared_X + ptrShared.x ] = temp ;

	 // (1,2,1)
       ptrMap.x -= sizeof(float)*Stride.x ;
	ptrShared.x -= Stride.x ;
	surf3Dread(&temp,SlownessMap, ptrMap.x, ptrMap.y, ptrMap.z );
	Slowness[ (ptrShared.z*shared_Y + ptrShared.y)*shared_X + ptrShared.x ] = temp ;
	surf3Dread(&temp,Traveltime_A, ptrMap.x, ptrMap.y, ptrMap.z);
	Traveltime[ (ptrShared.z*shared_Y + ptrShared.y)*shared_X + ptrShared.x ] = temp ;

	//********************************************************************
	
	// (1,2,2)
	   ptrMap.z += Stride.z ;
	ptrShared.z += Stride.z ;
	surf3Dread(&temp,SlownessMap, ptrMap.x, ptrMap.y, ptrMap.z );
	Slowness[ (ptrShared.z*shared_Y + ptrShared.y)*shared_X + ptrShared.x ] = temp ;
	surf3Dread(&temp,Traveltime_A, ptrMap.x, ptrMap.y, ptrMap.z);
	Traveltime[ (ptrShared.z*shared_Y + ptrShared.y)*shared_X + ptrShared.x ] = temp ;

	// (2,2,2)
       ptrMap.x += sizeof(float)*Stride.x ;
	ptrShared.x += Stride.x ;
	surf3Dread(&temp,SlownessMap, ptrMap.x, ptrMap.y, ptrMap.z );
	Slowness[ (ptrShared.z*shared_Y + ptrShared.y)*shared_X + ptrShared.x ] = temp ;
	surf3Dread(&temp,Traveltime_A, ptrMap.x, ptrMap.y, ptrMap.z);
	Traveltime[ (ptrShared.z*shared_Y + ptrShared.y)*shared_X + ptrShared.x ] = temp ;

	 // (2,1,2)
	   ptrMap.y -= Stride.y ;
	ptrShared.y -= Stride.y ;
	surf3Dread(&temp,SlownessMap, ptrMap.x, ptrMap.y, ptrMap.z );
	Slowness[ (ptrShared.z*shared_Y + ptrShared.y)*shared_X + ptrShared.x ] = temp ;
	surf3Dread(&temp,Traveltime_A, ptrMap.x, ptrMap.y, ptrMap.z);
	Traveltime[ (ptrShared.z*shared_Y + ptrShared.y)*shared_X + ptrShared.x ] = temp ;

	 // (1,1,2)
       ptrMap.x -= sizeof(float)*Stride.x ;
	ptrShared.x -= Stride.x ;
	surf3Dread(&temp,SlownessMap, ptrMap.x, ptrMap.y, ptrMap.z );
	Slowness[ (ptrShared.z*shared_Y + ptrShared.y)*shared_X + ptrShared.x ] = temp ;
	surf3Dread(&temp,Traveltime_A, ptrMap.x, ptrMap.y, ptrMap.z);
	Traveltime[ (ptrShared.z*shared_Y + ptrShared.y)*shared_X + ptrShared.x ] = temp ;


}


/// Compare New and Old traveltime
__device__ void compareTraveltime(const float *Traveltime, const float *Slowness, bool *isUpdated)
{
	int3 posWorker =  make_int3( threadIdx.x+radius_X, threadIdx.y+radius_Y,  threadIdx.z+radius_Z) ; // posWorker on shMem ###
	float Tt = Traveltime[ (posWorker.z*shared_Y + posWorker.y)*shared_X + posWorker.x ] ;
	int rayPath = -100 ;
	
	
	// Start:: Comparing neighbor != (0,0,0) ***************************************
    
	//-----------------------------------------------------------------------
	float  T1, T2,  T3,  T4, T5, Tt1 ;
	float Tz1, Tz2, Tz3, Tz4 ;
	float Tz5, Tz6, Tz7, Tz8 ;
	float Te1, Te2, Te3 ;
	
    // As 2D, START:: Level z = -2 ------------------------------------------

          Tt1 = SPR(Slowness,-2,-1,-2, posWorker); 
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,-1,-2, Tt1); 

	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-1,-1,-2, SPR(Slowness,-1,-1,-2, posWorker));	

	      T1  = SPR(Slowness,-1,-2,-2, posWorker); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-1,-2,-2,  T1);	

	      T5  = SPR(Slowness,-2,-2,-1, posWorker); 
	      Tz1 = T5 ;
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,-2,-1,  T5);

	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,-2,-2, FIM(T1,Tt1,T5, Slowness[posWorker.x-2 + shared_X*( posWorker.y-2 +shared_Y*( posWorker.z-2 )) ] ));
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,-2,-2, SPR(Slowness,-2,-2,-2, posWorker));
	
	// -----------------
	
	      T2  = SPR(Slowness,0,-1,-2, posWorker);
	      Te1 = T2 ;  
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,0,-1,-2,  T2);

	      T3  = SPR(Slowness,1,-2,-2, posWorker); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,1,-2,-2,  T3);

	      T5  = SPR(Slowness,0,-2,-1, posWorker);
	      Tz5 = T5 ; 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,0,-2,-1,  T5);

    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,0,-2,-2, FIM(min(T1,T3),T2,T5 , Slowness[posWorker.x + shared_X*( posWorker.y-2 +shared_Y*( posWorker.z-2 )) ] ));
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,0,-2,-2, SPR(Slowness,0,-2,-2, posWorker));

    // -----------------

	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,1,-1,-2, SPR(Slowness,1,-1,-2, posWorker));
	
	      T1  = SPR(Slowness,2,-1,-2, posWorker); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,-1,-2,  T1);
	
	      T5  = SPR(Slowness,2,-2,-1, posWorker); 
	      Tz2 = T5 ;
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,-2,-1,  T5);

    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,-2,-2, FIM(T3,T1,T5 , Slowness[posWorker.x+2 + shared_X*( posWorker.y-2 +shared_Y*( posWorker.z-2 )) ] ));   
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,-2,-2, SPR(Slowness,2,-2,-2, posWorker));
	
	// -----------------
	 
	      T2  = SPR(Slowness,1,0,-2, posWorker); 
	      Te2 = T2 ; 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,1,0,-2,  T2);
	
	      T3  = SPR(Slowness,2,1,-2, posWorker); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,1,-2,  T3);
	
	      T5  = SPR(Slowness,2,0,-1, posWorker); 
	      Tz6 = T5 ;
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,0,-1,  T5);

    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,0,-2, FIM(T2,min(T1,T3),T5 , Slowness[posWorker.x+2 + shared_X*( posWorker.y +shared_Y*( posWorker.z-2 )) ] ));   
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,0,-2, SPR(Slowness,2,0,-2, posWorker));
 
    // -----------------

	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,1,1,-2, SPR(Slowness,1,1,-2, posWorker));
	
	      T1  = SPR(Slowness,1,2,-2, posWorker); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,1,2,-2,  T1);
	
	      T5  = SPR(Slowness,2,2,-1, posWorker);
	      Tz3 = T5 ; 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,2,-1,  T5);

    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,2,-2, FIM(T1,T3,T5 , Slowness[posWorker.x+2 + shared_X*( posWorker.y+2 +shared_Y*( posWorker.z-2 )) ] ));   
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,2,-2, SPR(Slowness,2,2,-2, posWorker));
   
    // -----------------
   
	      T2  = SPR(Slowness, 0,1,-2, posWorker);
	      Te3 = T2 ; 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker, 0,1,-2, T2);
	
	      T3  = SPR(Slowness,-1,2,-2, posWorker); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-1,2,-2, T3);
	
          T5  = SPR(Slowness, 0,2,-1, posWorker); 
          Tz7 = T5 ;
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker, 0,2,-1, T5);

    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker, 0,2,-2, FIM(min(T1,T3),T2,T5 , Slowness[posWorker.x + shared_X*( posWorker.y+2 +shared_Y*( posWorker.z-2 )) ] ));   
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker, 0,2,-2, SPR(Slowness,0,2,-2, posWorker));
   
    // -----------------
	
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-1,1,-2, SPR(Slowness,-1,1,-2, posWorker));
	
	      T1  = SPR(Slowness,-2,1,-2, posWorker); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,1,-2, T1);
	
	      T5  = SPR(Slowness,-2,2,-1, posWorker); 
	      Tz4 = T5 ;
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,2,-1, T5);

    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,2,-2, FIM(T3,T1,T5 , Slowness[posWorker.x-2 + shared_X*( posWorker.y+2 +shared_Y*( posWorker.z-2 )) ] ));   
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,2,-2, SPR(Slowness,-2,2,-2, posWorker));
	
	// -----------------
	
	      T2  = SPR(Slowness,-1,0,-2, posWorker); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-1,0,-2, T2);
	
	      T5  = SPR(Slowness,-2,0,-1, posWorker); 
	      Tz8 = T5 ;
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,0,-1, T5);
	
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,0,-2, FIM(T2,min(T1,Tt1),T5 , Slowness[posWorker.x-2 + shared_X*( posWorker.y +shared_Y*( posWorker.z-2 )) ] ));   
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,0,-2, SPR(Slowness,-2,0,-2, posWorker));
    
    // -----------------
    // Special node (0,0)
    
    	  T4  = SPR(Slowness,0,0,-1, posWorker);
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,0,0,-1, T4);
    
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,0,0,-2, FIM(min(T2,Te2),min(Te1,Te3), T4, Slowness[posWorker.x + shared_X*( posWorker.y +shared_Y*( posWorker.z-2 )) ] )); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,0,0,-2, SPR(Slowness,0,0,-2, posWorker));

    // As 2D,   END:: Level z = -2 ------------------------------------------
   	// As 2D, START:: Level z = -1 ------------------------------------------

    // Comparing neighbor (With edge reduction)
	// For node (i,j,-1) except (0,0,-1) 
	for(int j = -radius_Y ; j<=radius_Y ; ++j )
	{
		for(int i = -radius_X ; i <= radius_X ; ++i )
		{
		    if( GCD(i,j)==1 )
		    {
                checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,i,j,-1, SPR(Slowness,i,j,-1, posWorker));
		    }
		}
	}
 
    // As 2D,   END:: Level z = -1 ------------------------------------------
    // As 2D, START:: Level z =  0 ------------------------------------------
 
 
          Tt1 = SPR(Slowness,-2,-1,0, posWorker); 
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,-1,0, Tt1); 

	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-1,-1,0, SPR(Slowness,-1,-1,0, posWorker));	

	      T1  = SPR(Slowness,-1,-2,0, posWorker); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-1,-2,0,  T1);	

          T2  = Tz1 ;
	      T5  = SPR(Slowness,-2,-2,1, posWorker); 
	      Tz1 = T5 ;
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,-2,1,  T5);

	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,-2,0, FIM(T1,Tt1,min(T2,T5), Slowness[posWorker.x-2 + shared_X*( posWorker.y-2 +shared_Y*( posWorker.z )) ] ));
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,-2,0, SPR(Slowness,-2,-2,0, posWorker));
	
	// -----------------
	
	      T3  = SPR(Slowness,0,-1,0, posWorker);
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,0,-1,0,  T3);

	      T4  = SPR(Slowness,1,-2,0, posWorker); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,1,-2,0,  T4 );

          T2  = Tz5 ; 
	      T5  = SPR(Slowness,0,-2,1, posWorker);
	      Tz5 = T5 ; 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,0,-2,1,  T5);
	
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,0,-2,0, FIM(min(T1,T4),T3,min(T2,T5), Slowness[posWorker.x + shared_X*( posWorker.y-2 +shared_Y*( posWorker.z )) ] )); 
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,0,-2,0, SPR(Slowness,0,-2,0, posWorker));

    // -----------------

	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,1,-1,0, SPR(Slowness,1,-1,0, posWorker));
	
	      T1  = SPR(Slowness,2,-1,0, posWorker); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,-1,0,  T1);
	
	      T2  = Tz2 ;
	      T5  = SPR(Slowness,2,-2,1, posWorker);
	      Tz2 = T5 ; 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,-2,1,  T5);

    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,-2,0, FIM(T4,T1,min(T2,T5) , Slowness[posWorker.x+2 + shared_X*( posWorker.y-2 +shared_Y*( posWorker.z )) ] ));   
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,-2,0, SPR(Slowness,2,-2,0, posWorker));
	
	// -----------------
	     
	      T3  = SPR(Slowness,1,0,0, posWorker);
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,1,0,0,  T3);
	
	      T4  = SPR(Slowness,2,1,0, posWorker); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,1,0,  T4);
 
          T2  = Tz6 ; 
	      T5  = SPR(Slowness,2,0,1, posWorker); 
	      Tz6 = T5 ;
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,0,1,  T5);
 
  	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,0,0, FIM(T3,min(T1,T4),min(T2,T5), Slowness[posWorker.x+2 + shared_X*( posWorker.y +shared_Y*( posWorker.z )) ] )); 
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,0,0, SPR(Slowness,2,0,0, posWorker));
 
    // -----------------

	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,1,1,0, SPR(Slowness,1,1,0, posWorker));
	
	      T1  = SPR(Slowness,1,2,0, posWorker); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,1,2,0,  T1);
	
	      T2  = Tz3 ;
	      T5  = SPR(Slowness,2,2,1, posWorker); 
	      Tz3 = T5 ;
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,2,1,  T5);

    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,2,0, FIM(T1,T4,min(T2,T5) , Slowness[posWorker.x+2 + shared_X*( posWorker.y+2 +shared_Y*( posWorker.z )) ] ));   
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,2,0, SPR(Slowness,2,2,0, posWorker));
   
    // -----------------
   
          T3  = SPR(Slowness, 0,1,0, posWorker);
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker, 0,1,0, T3);
	
	      T4  = SPR(Slowness,-1,2,0, posWorker); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-1,2,0, T4);
	
	      T2  = Tz7 ; 
	      T5  = SPR(Slowness, 0,2,1, posWorker); 
	      Tz7 = T5 ;
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker, 0,2,1,  T5);
 
 	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker, 0,2,0, FIM(min(T1,T4),T3,min(T2,T5), Slowness[posWorker.x + shared_X*( posWorker.y+2 +shared_Y*( posWorker.z )) ] )); 
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker, 0,2,0, SPR(Slowness,0,2,0, posWorker));
   
    // -----------------
	
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-1,1,0, SPR(Slowness,-1,1,0, posWorker));
	
	      T1  = SPR(Slowness,-2,1,0, posWorker); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,1,0, T1);
	
	      T2  = Tz4 ; 
	      T5  = SPR(Slowness,-2,2,1, posWorker); 
	      Tz4 = T5 ;
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,2,1, T5);

    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,2,0, FIM(T4,T1,min(T2,T5) , Slowness[posWorker.x-2 + shared_X*( posWorker.y+2 +shared_Y*( posWorker.z )) ] ));   
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,2,0, SPR(Slowness,-2,2,0, posWorker));
	
	// -----------------
	
	      T3  = SPR(Slowness,-1,0,0, posWorker);
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-1,0,0, T3);
	 
	      T2  = Tz8 ; 
	      T5  = SPR(Slowness,-2,0,1, posWorker); 
	      Tz8 = T5 ;
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,0,1,  T5);
	
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,0,0, FIM(T3,min(Tt1, T1),min(T2,T5), Slowness[posWorker.x-2 + shared_X*( posWorker.y +shared_Y*( posWorker.z )) ] )); 
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,0,0, SPR(Slowness,-2,0,0, posWorker));
    
 
    // As 2D,   END:: Level z =  0 ------------------------------------------
    // As 2D, START:: Level z = +1 ------------------------------------------
    
    // Comparing neighbor (With edge reduction)
	// For node (i,j,1) except (0,0,1) 
	for(int j = -radius_Y ; j<=radius_Y ; ++j )
	{
		for(int i = -radius_X ; i <= radius_X ; ++i )
		{
		    if( GCD(i,j)==1 )
		    {
                checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,i,j,1, SPR(Slowness,i,j,1, posWorker));
		    }
		}
	}
    
    // As 2D,   END:: Level z = +1 ------------------------------------------
    // As 2D, START:: Level z = +2 ------------------------------------------

          Tt1 = SPR(Slowness,-2,-1,2, posWorker); 
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,-1,2, Tt1); 

	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-1,-1,2, SPR(Slowness,-1,-1,2, posWorker));	

	      T1  = SPR(Slowness,-1,-2,2, posWorker); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-1,-2,2,  T1);	

	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,-2,2, FIM(T1,Tt1,Tz1, Slowness[posWorker.x-2 + shared_X*( posWorker.y-2 +shared_Y*( posWorker.z+2 )) ] ));
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,-2,2, SPR(Slowness,-2,-2,2, posWorker));
	
	// -----------------
	
	      T2  = SPR(Slowness,0,-1,2, posWorker); 
	      Te1 = T2 ; 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,0,-1,2,  T2);

	      T3  = SPR(Slowness,1,-2,2, posWorker); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,1,-2,2,  T3);

    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,0,-2,2, FIM(min(T1,T3),T2,Tz5 , Slowness[posWorker.x + shared_X*( posWorker.y-2 +shared_Y*( posWorker.z+2 )) ] ));
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,0,-2,2, SPR(Slowness,0,-2,2, posWorker));

    // -----------------

	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,1,-1,2, SPR(Slowness,1,-1,2, posWorker));
	
	      T1  = SPR(Slowness,2,-1,2, posWorker); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,-1,2,  T1);
	
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,-2,2, FIM(T3,T1,Tz2 , Slowness[posWorker.x+2 + shared_X*( posWorker.y-2 +shared_Y*( posWorker.z+2 )) ] ));   
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,-2,2, SPR(Slowness,2,-2,2, posWorker));
	
	// -----------------
	 
	      T2  = SPR(Slowness,1,0,2, posWorker); 
	      Te2 = T2 ; 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,1,0,2,  T2);
	
	      T3  = SPR(Slowness,2,1,2, posWorker); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,1,2,  T3);

    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,0,2, FIM(T2,min(T1,T3),Tz6 , Slowness[posWorker.x+2 + shared_X*( posWorker.y +shared_Y*( posWorker.z+2 )) ] ));   
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,0,2, SPR(Slowness,2,0,2, posWorker));
 
    // -----------------

	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,1,1,2, SPR(Slowness,1,1,2, posWorker));
	
	      T1  = SPR(Slowness,1,2,2, posWorker); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,1,2,2,  T1);
	

    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,2,2, FIM(T1,T3,Tz3 , Slowness[posWorker.x+2 + shared_X*( posWorker.y+2 +shared_Y*( posWorker.z+2 )) ] ));   
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,2,2,2, SPR(Slowness,2,2,2, posWorker));
   
    // -----------------
   
	      T2  = SPR(Slowness, 0,1,2, posWorker); 
	      Te3 = T2 ; 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker, 0,1,2, T2);
	
	      T3  = SPR(Slowness,-1,2,2, posWorker); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-1,2,2, T3);

    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker, 0,2,2, FIM(min(T1,T3),T2,Tz7 , Slowness[posWorker.x + shared_X*( posWorker.y+2 +shared_Y*( posWorker.z+2 )) ] ));   
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker, 0,2,2, SPR(Slowness,0,2,2, posWorker));
   
    // -----------------
	
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-1,1,2, SPR(Slowness,-1,1,2, posWorker));
	
	      T1  = SPR(Slowness,-2,1,2, posWorker); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,1,2, T1);

    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,2,2, FIM(T3,T1,Tz4 , Slowness[posWorker.x-2 + shared_X*( posWorker.y+2 +shared_Y*( posWorker.z+2 )) ] ));   
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,2,2, SPR(Slowness,-2,2,2, posWorker));
	
	// -----------------
	
	      T2  = SPR(Slowness,-1,0,2, posWorker); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-1,0,2, T2);
	
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,0,2, FIM(T2,min(T1,Tt1),Tz8 , Slowness[posWorker.x-2 + shared_X*( posWorker.y +shared_Y*( posWorker.z+2 )) ] ));   
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,-2,0,2, SPR(Slowness,-2,0,2, posWorker));
    
    // -----------------
    // Special node (0,0)	
   
          T4  = SPR(Slowness,0,0,1, posWorker);
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,0,0,1, T4);
    	
    checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,0,0,2, FIM(min(T2,Te2),min(Te1,Te3),T4 , Slowness[posWorker.x + shared_X*( posWorker.y +shared_Y*( posWorker.z+2 )) ] )); 
	checkTraveltime(Traveltime,&Tt,&rayPath,posWorker,0,0,2, SPR(Slowness,0,0,2, posWorker));


    // As 2D, END:: Level z = +2 --------------------------------------------
 
	//-----------------------------------------------------------------------
	
    // End:: Comparing neighbor != (0,0,0) ***************************************

	// Write outTraveltime
    int GlobalID_X = blockIdx.x*blockDim.x + threadIdx.x ;
    int GlobalID_Y = blockIdx.y*blockDim.y + threadIdx.y ;
	int GlobalID_Z = blockIdx.z*blockDim.z + threadIdx.z ;
    surf3Dwrite(Tt, Traveltime_B, sizeof(float)*(GlobalID_X+radius_X), GlobalID_Y+radius_Y,  GlobalID_Z+radius_Z);
		
	// if threadUpdate -> isUpdated[0] = true ;
	if( rayPath != -100 ){
      *isUpdated = true ;
	  surf3Dwrite(rayPath, RaypathMap, sizeof(int)*GlobalID_X, GlobalID_Y, GlobalID_Z);
    }

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
	if( threadIdx.x<3 && threadIdx.y<3 && threadIdx.z<3 ){
		int data ; 
		surf3Dread(&data, UpdateMap_A, sizeof(int)*(blockIdx.x+threadIdx.x), blockIdx.y+threadIdx.y, blockIdx.z+threadIdx.z);
	    if( data==1 ){
	        needUpdate = true ;
	    }
	}
	__syncthreads();
	if( !needUpdate ){ return; }



    // Define shared variables (2)
	extern __shared__ float sharedVar[] ;
	float *Slowness   = &sharedVar[0] ;  
	float *Traveltime = &sharedVar[ shared_X*shared_Y*shared_Z ] ; 
    //isAnyUpdated = false ; //%%%

	/// Upload Slowness x4 and Traveltime x4 to shared memory	
    uploadShared( Slowness, Traveltime); 
	__syncthreads();

	/// Compare New and Old traveltime
	compareTraveltime(Traveltime, Slowness, &isAnyUpdated);
	__syncthreads();

/*
    // For no update map scheme  
    // ### don't forget to set isUpdate[0] = false and remove the __syncthreads() 
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

bool PSPS_CUDA_ExtendedBF_3D_Solver::Compute()
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
   front = length[0]*block[1]/2 ;
   top = length[0]*length[1]*block[2]/2 ;
   
    //******************************************************************************************
    
   /// Add dummy vertices ONLY 3D ///
   addDummyVertices();

	// Set constant values
	int tempSh[] = { (int)shared[0], (int)shared[1], (int)shared[2]} ;
	cudaMemcpyToSymbol(shared_X, &tempSh[0], sizeof(int));
	cudaMemcpyToSymbol(shared_Y, &tempSh[1], sizeof(int));
	cudaMemcpyToSymbol(shared_Z, &tempSh[2], sizeof(int));
	
	int tempGr[] = { (int)block[0]+1, (int)block[1]+1, (int)block[2]+1} ;
	cudaMemcpyToSymbol(group_X, &tempGr[0], sizeof(int));
	cudaMemcpyToSymbol(group_Y, &tempGr[1], sizeof(int));
	cudaMemcpyToSymbol(group_Z, &tempGr[2], sizeof(int));
	
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
   // Raypath 
   cudaArray* cuArrayRaypath ;
   cudaMalloc3DArray(&cuArrayRaypath, &channelDescInt, make_cudaExtent(res[0], res[1], res[2]), cudaArraySurfaceLoadStore);
   cudaMemcpy3DParms paramInR = {0} ;
	  paramInR.srcPtr = make_cudaPitchedPtr(&raypath[0], res[0]*sizeof(int), res[0], res[1]) ;
	  paramInR.dstArray = cuArrayRaypath ;
      paramInR.kind =  cudaMemcpyHostToDevice ;
      paramInR.extent = make_cudaExtent( res[0], res[1], res[2]) ;
   cudaMemcpy3D(&paramInR);
   cudaBindSurfaceToArray(RaypathMap, cuArrayRaypath);
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
	  size_t nSharedVar = 2*shared[0]*shared[1]*shared[2]*sizeof(float) ;

	  // Current Time
      time_t rawtime ;
      
      // Start
      std::cout << "\nPSPS_CUDA_ExtendedBF >> Start \'Extended\' Bellman-Ford (CUDA).\n" ;
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
		 cuda_PSPS_Solver<<<dimGroup,dimBlock,nSharedVar>>>(cuRunning);

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
      std::cout << "PSPS_CUDA_ExtendedBF >> Finish \'Extended\' Bellman-Ford (CUDA).\n\n" ;
        
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

	// Raypath
	cudaMemcpy3DParms paramOutR = {0} ;
	  paramOutR.srcArray = cuArrayRaypath ;
	  paramOutR.dstPtr = make_cudaPitchedPtr(&raypath[0], res[0]*sizeof(int), res[0], res[1]) ;
      paramOutR.kind =  cudaMemcpyDeviceToHost ;
      paramOutR.extent = make_cudaExtent( res[0], res[1], res[2]) ;
    cudaMemcpy3D(&paramOutR);
	
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
   cudaFreeArray(cuArrayRaypath);
   cudaFreeArray(cuArrayTraveltime_A);
   cudaFreeArray(cuArrayTraveltime_B);
   cudaFreeArray(cuArrayUpdateMap_A);
   cudaFreeArray(cuArrayUpdateMap_B);

return true; }

/************************************* Private *************************************/

void PSPS_CUDA_ExtendedBF_3D_Solver::addDummyVertices()
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

void PSPS_CUDA_ExtendedBF_3D_Solver::removeDummyVertices()
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
