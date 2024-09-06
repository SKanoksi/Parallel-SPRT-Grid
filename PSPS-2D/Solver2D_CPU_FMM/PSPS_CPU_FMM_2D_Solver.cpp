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

#include "PSPS_CPU_FMM_2D_Solver.h"

PSPS_CPU_FMM_2D_Solver::PSPS_CPU_FMM_2D_Solver()
{

}

PSPS_CPU_FMM_2D_Solver::~PSPS_CPU_FMM_2D_Solver()
{

}

//**************** CPU_FMM program ****************//

// For simplicity when use them
inline bool PSPS_CPU_FMM_2D_Solver::isShortest(int x){ return (0<=x) ;}
inline float PSPS_CPU_FMM_2D_Solver::min(float a, float b){ if( a<b ){ return a; }else{ return b; } }

float PSPS_CPU_FMM_2D_Solver::newCandidateOfTraveltime(const int ptr)
{
	float Tx = min( traveltime[ ptr-1 ], traveltime[ ptr+1 ]);
	float Ty = min( traveltime[ ptr-length[0] ], traveltime[ ptr+length[0] ]);
	float S = slowness[ptr] ;
	
	// hx = hy Version, already included in S, (adding and removing dummy vertices)
	// Let Tx > Ty
	float temp ;
	if(Tx < Ty){ temp = Tx ; Tx = Ty ; Ty = temp ; }

	float TT = Ty + S ;
	if( 1.414213f*S > Tx-Ty ) 
	{	
		temp = 0.5f*( (Tx+Ty) + std::sqrt( 2.0f*std::pow(S,2) - std::pow(Tx-Ty,2) ) ) ;
		if(temp > Tx){ TT = temp; }
	}
	
	    
return TT ; }


bool PSPS_CPU_FMM_2D_Solver::Compute()
{
   /// Create Compute Shader (Solver) ///
   // Compute some useful numbers.
   for(unsigned int i=0 ; i<DIM ; ++i)
   {
      length[i] = int(res[i]+block[i]) ;
	  radius[i] = int(block[i]/2) ;
   }
   side  = radius[0] ;
   front = (length[0]+1)*radius[1] ; // = actual front + side
   int xGroup = block[0]+1 ;
   // Constant used for hx != hy version
   invHxSquare = 1.0f/stride[0]*stride[0] ;
   invHySquare = 1.0f/stride[1]*stride[1] ;
   sumInv = invHxSquare + invHySquare ;

   //******************************************************************************************

   /// Add dummy vertices ONLY 2D ///
   addDummyVertices();
   

   { ///***** =.\= Parallel Shortest Path Solver (CPU_FMM) =.\= *****///

    // Current Time
    time_t rawtime ;
      
    // Start
    std::cout << "\nPSPS_CPU_FMM >> Start Fast Marching Method (CPU_FMM).\n" ;
    time(&rawtime);
    std::cout << "   START:: " << ctime(&rawtime) ;
    
    // Start timer
    auto start_time = std::chrono::high_resolution_clock::now();

    // Declare Min queue
    std::priority_queue<Vertex, std::vector<Vertex>, std::greater<Vertex> > MyQ ;
    // PUSH source = Vertex(0.0f, sourceX, sourceY)
    MyQ.push( std::make_pair( 0.0f, length[0]*source[1] +source[0]+front) );


    ///*** Main loop ***///
    while( !MyQ.empty() )
    {
        // Extract the min TT of forward stars
		Vertex ptr = MyQ.top() ;
		MyQ.pop();

		// Check isShortest, if not check its neighbors
        if( !isShortest( raypath[ptr.second] ) )
		{
            // isShortest = false -> set isShortest = true 
            raypath[ptr.second] *= -1 ; 

			// Neighbor = Down
	      	int ptrMapNeighbor = ptr.second - length[0] ; 
		    if( !isShortest( raypath[ptrMapNeighbor] ) )
		    {
			    float newTT = newCandidateOfTraveltime(ptrMapNeighbor);
			    if( newTT<traveltime[ ptrMapNeighbor ] )
			    {
				    // Update traveltime
				    traveltime[ ptrMapNeighbor ] = newTT ;
				    raypath[ ptrMapNeighbor ] = -1 ; // isShortest(false), +1 is added.$$$
				    // Push it in Q
				    MyQ.push( std::make_pair(newTT, ptr.second-length[0] ) );
                }
		    }
		  
			// Neighbor = Left
	      	ptrMapNeighbor = ptr.second -1 ; 
		    if( !isShortest( raypath[ptrMapNeighbor] ) )
		    {
			    float newTT = newCandidateOfTraveltime(ptrMapNeighbor); 
			    if( newTT<traveltime[ ptrMapNeighbor ] )
			    {
				    // Update traveltime
				    traveltime[ ptrMapNeighbor ] = newTT ;
				    raypath[ ptrMapNeighbor ] = -1 ; // isShortest(false), +1 is added.$$$
				    // Push it in Q
				    MyQ.push( std::make_pair(newTT, ptr.second-1 ) );
                }
  	 	    }

			// Neighbor = Right
	      	ptrMapNeighbor = ptr.second+ 1 ; 
		    if( !isShortest( raypath[ptrMapNeighbor] ) )
		    {
			    float newTT = newCandidateOfTraveltime(ptrMapNeighbor);
			    if( newTT<traveltime[ ptrMapNeighbor ] )
			    {
				    // Update traveltime
				    traveltime[ ptrMapNeighbor ] = newTT ;
				    raypath[ ptrMapNeighbor ] = -1 ; // isShortest(false), +1 is added.$$$
				    // Push it in Q
				    MyQ.push( std::make_pair(newTT, ptr.second+1 ) );
                }
		    }

			// Neighbor = Up
	      	ptrMapNeighbor = ptr.second + length[0] ; 
		    if( !isShortest( raypath[ptrMapNeighbor] ) )
		    {
			    float newTT = newCandidateOfTraveltime(ptrMapNeighbor);
			    if( newTT<traveltime[ ptrMapNeighbor ] )
			    {
				    // Update traveltime
				    traveltime[ ptrMapNeighbor ] = newTT ;
				    raypath[ ptrMapNeighbor ] = -1 ; // isShortest(false), +1 is added.$$$
				    // Push it in Q
				    MyQ.push( std::make_pair(newTT, ptr.second+length[0] ) );
                }    
 		    }

    		// ### End checking neighbors ###
		}

    }

    // Stop timer
    auto end_time = std::chrono::high_resolution_clock::now();
      
    // Finish
    time(&rawtime);
    std::cout << "  FINISH:: " << ctime(&rawtime) ;
    std::cout << "PSPS_CPU_FMM >> Finish Fast Marching Method (CPU_FMM).\n\n" ;
    std::cout << " *** Runtime = " << std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count() << " seconds. ***\n\n";
          
    }


	//******************************************************************************************

    /// Remove dummy vertices ONLY 2D ///
    removeDummyVertices();
    
    
return true; }


/************************************* Private *************************************/

void PSPS_CPU_FMM_2D_Solver::addDummyVertices()
{
	// Not efficient but easy to code ^.^ (just copy the other version)
   unsigned int at ;
   bool isNegative = false ;
   std::vector<TyPe> temp ;

   /// Slowness map.
   // Copy data
   temp.assign( slowness.begin(), slowness.end() ) ;
   // Reinput the data
   slowness.assign( length[0]*length[1], 1.0f/0.0f) ;
   at = front ;
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
   at = front ;
   for(unsigned int j=0 ; j<res[1] ; ++j)
   {
      for(unsigned int i=0 ; i<res[0] ; ++i)
      {
         traveltime[ at+i ] = temp[ res[0]*j+i ] ;
      }
      at += res[0]+2*side ;
   }


   // Raypath = Checking shortest (Just for now)
   // Copy data
   std::vector<int> tem ;
   tem.assign( raypath.begin(), raypath.end() ) ;
   // Reinput the data
   raypath.assign( length[0]*length[1], 2) ; // Any positive number not only 2
   at = front ;
   for(unsigned int j=0 ; j<res[1] ; ++j)
   {
      for(unsigned int i=0 ; i<res[0] ; ++i)
      {
         raypath[ at+i ] = tem[ res[0]*j+i ] ;
      }
      at += res[0]+2*side ;
   }

return; }


void PSPS_CPU_FMM_2D_Solver::removeDummyVertices()
{
   // Not efficient but easy to code ^.^ (just copy the other version)
   unsigned int at ;
   std::vector<TyPe> temp ;
   
   /// Slowness map
   // Copy data
   temp.assign( slowness.begin(), slowness.end() ) ;
   // Reinput the data
   slowness.resize( res[0]*res[1] ) ;
   at = front ;
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
   at = front ;
   for(unsigned int j=0 ; j<res[1] ; ++j)
   {
      for(unsigned int i=0 ; i<res[0] ; ++i)
      {
         traveltime[ res[0]*j+i ] = temp[ at+i ] ;
      }
      at += res[0]+2*side ;
   }

   // Raypath = Checking shortest (Just for now)
   /// Copy data
   std::vector<int> tem ;
   tem.assign( raypath.begin(), raypath.end() ) ;
   // Reinput the data
   raypath.resize( res[0]*res[1] ) ;
   at = front ;
   for(unsigned int j=0 ; j<res[1] ; ++j)
   {
      for(unsigned int i=0 ; i<res[0] ; ++i)
      {
         raypath[ res[0]*j+i ] = tem[ at+i ]-1 ; // Remove +1 when write raypath. $$$
      }
      at += res[0]+2*side ;
   }

return; }
