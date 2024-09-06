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

#include "PSPS_CPU_DK_2D_Solver.h"

PSPS_CPU_DK_2D_Solver::PSPS_CPU_DK_2D_Solver()
{

}

PSPS_CPU_DK_2D_Solver::~PSPS_CPU_DK_2D_Solver()
{

}

//**************** CPU_DK program ****************//

// For simplicity when use them
inline bool PSPS_CPU_DK_2D_Solver::isShortest(int x){ return (0<=x) ;}

// Calculate "Line" Weight (Actual weight of each edge to neighbor = sum of Slowness*LineWeight)
void PSPS_CPU_DK_2D_Solver::setLineWeight()
{
	// Set LineWeight
	for(int j = -radius[1] ; j<=radius[1] ; ++j)
	for(int i = -radius[0] ; i<=radius[0] ; ++i)
	{
    	LineWeight.push_back( calLineSegment(i,j) ); // [0] = -radius[1]*radius[0]-radius[0] (the shift)
	}

return ; }

LineSegment  PSPS_CPU_DK_2D_Solver::calLineSegment(const int VertexX, const int VertexY)
{
	// Declare variable
	LineSegment save ;

	// Set variables
	int     nVec[2] = { std::abs(VertexX), std::abs(VertexY) };
	const float Stride[2] = { (float) 1.0f/nVec[0] , (float) 1.0f/nVec[1] } ;
	float  cross[2] = { (float) 1.0f-0.5f*Stride[0], (float) 1.0f-0.5f*Stride[1] } ;
	float  crossOld = 1.0f ;
	// Multiply by scale
	const float fVec[2] = { (float)VertexX, (float)VertexY } ;
	const float scale   = std::sqrt( fVec[0]*fVec[0]*stride[0]*stride[0] + fVec[1]*fVec[1]*stride[1]*stride[1] ) ;

	// Main process
	while( 0 < nVec[0] || 0 < nVec[1] )
	{
		// Find the next largest crossing point wrt. crossOld
		if( cross[0]<cross[1] ){

			// Compute weight.
			int indexX = int( roundf( 0.5f*(crossOld+cross[1])*fVec[0] ) ) ;
			int indexY = int( roundf( 0.5f*(crossOld+cross[1])*fVec[1] ) ) ;
			save.shift.push_back( length[0]*indexY+indexX );
			save.length.push_back( (crossOld-cross[1])*scale );

			// For next crossing point.
			crossOld = cross[1] ;
			--nVec[1] ;
			cross[1] -= Stride[1] ;

		}else{

			// Compute weight.
			int indexX = int( roundf( 0.5f*(crossOld+cross[0])*fVec[0] ) ) ;
			int indexY = int( roundf( 0.5f*(crossOld+cross[0])*fVec[1] ) ) ;
			save.shift.push_back( length[0]*indexY+indexX );
			save.length.push_back( (crossOld-cross[0])*scale );

			// For next crossing point.
			crossOld = cross[0] ;
			--nVec[0] ;
			cross[0] -= Stride[0] ;

		}
	}
	// Weight at source
	save.shift.push_back( 0 ) ; // = iMap(0,0) = length[0]*indexY+indexX
	save.length.push_back( crossOld*scale );

return save ;}

// Calculate Weight (Slowness*LineWeight)
float PSPS_CPU_DK_2D_Solver::findWeight(const int ptrMap, const int shiftMap)
{
	// find Weight
	float w = 0 ;
	for(unsigned int m = 0 ; m<LineWeight[shiftMap].shift.size() ; ++m )
	{
			w += slowness[ ptrMap+LineWeight[shiftMap].shift[m] ]*LineWeight[shiftMap].length[m] ;
	}

return w ; }

int PSPS_CPU_DK_2D_Solver::GCD(int a, int b)
{
    a = std::abs(a) ;
    b = std::abs(b) ;
    int c ;
    while( a!=0 )
    {
        c = a ;
        a = b%a ;
        b = c ;
    }

return b; }

bool PSPS_CPU_DK_2D_Solver::Compute()
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

   //******************************************************************************************

   /// Add dummy vertices ONLY 2D ///
   addDummyVertices();


   { ///***** =.\= Parallel Shortest Path Solver (CPU_DK) =.\= *****///

    // Current Time
    time_t rawtime ;
      
    // Start
    std::cout << "\nPSPS_CPU_DK >> Start Dijkstra's algorithm (CPU_DK).\n" ;
    time(&rawtime);
    std::cout << "   START:: " << ctime(&rawtime) ;

    // Start timer
    auto start_time = std::chrono::high_resolution_clock::now();
      
    // Initialization
    setLineWeight();
    // Declare Min queue
    std::priority_queue<Vertex, std::vector<Vertex>, std::greater<Vertex> > MyQ ;
    // PUSH source = Vertex(0.0f, indexOfSource)
    MyQ.push( std::make_pair( 0.0f, length[0]*int(source[1]) + int(source[0]) + front ) );


    ///*** Main loop ***///
    while( !MyQ.empty() )
    {
        // Extract the min TT of forward stars
		Vertex ptr = MyQ.top() ;
		MyQ.pop();

		// Check isShortest, if not check its neighbors
        int ptrMap = ptr.second ;
		if( !isShortest( raypath[ptr.second] ) )
		{
            // isShortest = false -> set isShortest = true 
            raypath[ptr.second] *= -1 ; 
 
            // ### Start checking neighbors WITH Edge reduction ### 
			for(int ry = -radius[1] ; ry <= radius[1] ; ++ry )
			for(int rx = -radius[0] ; rx <= radius[0] ; ++rx )
			{
			    if( GCD(rx,ry)!=1 ){ continue; } // Comment this line = No edge reduction *#*#*#*#*#*
                // Check neighbor = isShortest
				int ptrMapNeighbor = ptr.second + (length[0]*ry+rx) ; 
				if( !isShortest( raypath[ptrMapNeighbor] ) )
				{
					// Calculate new candidate of traveltime of the neighbor
                    int ray = (ry+radius[1])*xGroup + rx+radius[0] ;
					float newTT = ptr.first + findWeight(ptr.second,ray) ;// ptr.first = tt[ptrMap=ptr.second]
		
                    // If newTT smaller then save it and push in Q, else = it's already in Q.
					if( newTT<traveltime[ ptrMapNeighbor ] )
					{
						// Update traveltime
						traveltime[ ptrMapNeighbor ] = newTT ;
						raypath[ ptrMapNeighbor ] = -(ray+1) ; // isShortest(false), +1 is added.$$$
						// Push it in Q
						MyQ.push( std::make_pair(newTT, ptrMapNeighbor) );
                    }

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
    std::cout << "PSPS_CPU_DK >> Finish Dijkstra's algorithm (CPU_DK).\n\n" ;
    std::cout << " *** Runtime = " << std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count() << " seconds. ***\n\n";
    
    }


	//******************************************************************************************

    /// Remove dummy vertices ONLY 2D ///
    removeDummyVertices();

return true; }


/************************************* Private *************************************/

void PSPS_CPU_DK_2D_Solver::addDummyVertices()
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
         slowness[ at+i ] = std::abs( temp[ res[0]*j+i ] ) ;
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
   
   std::vector<int> tem ;
   // Copy data
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


void PSPS_CPU_DK_2D_Solver::removeDummyVertices()
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
         slowness[ res[0]*j+i ] = temp[ at+i ] ;
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
   
   std::vector<int> tem ;
   // Copy data
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
