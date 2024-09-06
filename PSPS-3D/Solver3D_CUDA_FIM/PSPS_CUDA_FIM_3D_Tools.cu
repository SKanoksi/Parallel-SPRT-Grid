/*************************************
     Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#include "PSPS_CUDA_FIM_3D_Tools.h"

#define TyPe float
#define NC_TyPe NC_FLOAT
#define DIM 3
#define NAME_DATA "Slowness"
#define GPU_MEMORY_LIMIT 0.90     // Will use at most 90% of total GPU memory. ###
#define SHARED_MEMORY_LIMIT 0.98  // Will use at most 98% of total GPU shared memory. ###

PSPS_CUDA_FIM_3D_Tools::PSPS_CUDA_FIM_3D_Tools()
{

}

PSPS_CUDA_FIM_3D_Tools::~PSPS_CUDA_FIM_3D_Tools()
{

}

bool PSPS_CUDA_FIM_3D_Tools::newProject(std::string &inFile, std::vector<int> &param)
{
   // Read slowness model.
   if( !readSlownessModel(inFile.c_str()) ){ return false; }

   // Input the spec.
   std::vector<unsigned int> res_old(res); // Save old resolution
   if( !inputSpec(param) ){ return false; }
   if( !checkSpec() ){ return false; }
   
   // Change the scale
   for(unsigned int i=0 ; i<DIM ; ++i)
   {
        stride[i] = (TyPe) (res_old[i]*stride[i])/res[i] ;
   }
   if( std::fabs(stride[0]-stride[1]) > 0.01*stride[1] || std::fabs(stride[2]-stride[1]) > 0.01*stride[1] )
   {
        std::cout << "\n!!! ERROR::FIM:: Difference of grid specing (hx,hy,hz) is more than 1 percent. !!!\n\n" ;
        return false;
   }

   // Rescale Slowness map
   for(unsigned int i=0 ; i<DIM ; ++i)
   {
      if( res[i]!=res_old[i] ){
         if( !rescaleSlowness(res_old) ){
            std::cout << "\nNewProject >> ERROR:: Rescaling Slowness model.\n" ;
            return false;
         }
         break;
      }
   }

   // Reallocate memory for other maps.
   unsigned int numData = 1 , numGroup = 1 ;
   for(unsigned int i=0 ; i<DIM ; ++i)
   {
       numData *= res[i] ;
      numGroup *= nWorkGroup[i] ;
   }
   traveltime.clear();
   traveltime.resize(  numData, 1.0f/0.0f);
    updateMap.clear();
    updateMap.resize( numGroup, 0);

   // Set source, Traveltime+updateMap Only 3D ###
   int atMap = 0, atGroup = 0 ;
   for(int i=DIM-1 ; 0<=i ; --i)
   {
        atMap *= res[i];
        atMap += source[i] ;
      atGroup *= nWorkGroup[i] ;
      atGroup += (int) source[i]/block[i] ;
   }
	//  atMap = res[0]*source[1] + source[0] ;
	//atGroup = nWorkGroup[0]*((int)source[1]/block[1]) + ((int)source[0]/block[0]) ;
    traveltime[atMap] = 0.0f ;
   updateMap[atGroup] = 1 ;

   // Set setup
      nloop = 0 ;
   isFinish = false ;

return true; }

bool PSPS_CUDA_FIM_3D_Tools::resumeProject(std::string &inFile, int maxloop)
{
   // Read the resume project.
   maxLoop = maxloop ;
   if( !readResumeProject(inFile.c_str()) ){ return false; }
   if( !checkSpec() ){ return false; }
   if( std::fabs(stride[0]-stride[1]) > 0.01*stride[1] || std::fabs(stride[2]-stride[1]) > 0.01*stride[1] )
   {
        std::cout << "\n!!! ERROR::FIM:: Difference of grid specing (hx,hy,hz) is more than 1 percent. !!!\n\n" ;
        return false;
   }
   
return true; }

bool PSPS_CUDA_FIM_3D_Tools::saveProject(std::string &inFile)
{
	// NetCDF error handling
	#ifndef nCDF
	  #define nCDF(Func) { int NC_Stat; if( NC_NOERR != (NC_Stat=Func) ){ fprintf(stderr, "ERROR::NETCDF:: %s\n", nc_strerror(NC_Stat)); std::exit(2);} }
	#else
	  #undef nCDF
	  #define nCDF(Func) { int NC_Stat; if( NC_NOERR != (NC_Stat=Func) ){ fprintf(stderr, "ERROR::NETCDF:: %s\n", nc_strerror(NC_Stat)); std::exit(2);} }
	#endif

   /// Read some data from input NetCDF file. ///
	// Open the input NetCDF file
	int ncid ;
	nCDF( nc_open(inFile.c_str(), NC_NOWRITE, &ncid) )

   // Get coordinates name
   std::string name[DIM] ;
   for(unsigned int i=0; i<DIM; ++i )
   {
      char coordName[NC_MAX_NAME+1];
      nCDF( nc_inq_dim(ncid,i+1,coordName,NULL) )
      name[i] = coordName ;
   }

   // Get max and min
   float min[DIM], max[DIM] ;
   int varid ;
   nCDF( nc_inq_varid(ncid, "minCoord", &varid) )
   nCDF( nc_get_var_float(ncid, varid, min) )
   nCDF( nc_inq_varid(ncid, "maxCoord", &varid) )
   nCDF( nc_get_var_float(ncid, varid, max) )
   // Get unit of Slowness
   char slownessUnit[NC_MAX_NAME+1] ;
   nCDF( nc_inq_varid(ncid, "Slowness", &varid) )
   nCDF( nc_get_att_text(ncid, varid, "Unit", slownessUnit) )

	// close the input NetCDF file
	nCDF( nc_close(ncid) )

   /**********************************/

   /// Write data to output NetCDF file.
   // Will use the old filename (Remove old file)
   if( remove(inFile.c_str())!=0 ){
      // Cannot remove the file -> Use new filename
      inFile = inFile + "_PSPS_CUDA_FIM.out" ;
   }
   // Create a NetCDF file
   nCDF( nc_create(inFile.c_str(), NC_CLOBBER, &ncid) ) // (may want to select cmode). ###

	// Define dimension and coordinates
	int dim_id, coord_id[DIM], group_id[DIM] ;
	nCDF( nc_def_dim(ncid, "Dimension", DIM, &dim_id) ) // Define dimension.
	for(unsigned int i=0 ; i<DIM ; ++i)
	{
		nCDF( nc_def_dim(ncid, name[i].c_str(), res[i], &coord_id[i]) ) // Define coordiantes
	}
	for(unsigned int i=0 ; i<DIM ; ++i)
	{
			if( !isFinish ){
         name[i] += "_nWorkGroup" ;
         nCDF( nc_def_dim(ncid, name[i].c_str(), nWorkGroup[i], &group_id[i]) ) // for UpdateMap
      }
	}

	// Define variables
    int min_id, max_id, stride_id ;
	int radius_id, source_id, state_id ;
	int slowness_id, traveltime_id, update_id ;
	nCDF( nc_def_var(ncid, "minCoord", NC_TyPe, 1, &dim_id, &min_id) ) // Define min[].
	nCDF( nc_def_var(ncid, "maxCoord", NC_TyPe, 1, &dim_id, &max_id) ) // Define max[].
	nCDF( nc_def_var(ncid, "strideCoord", NC_TyPe, 1, &dim_id, &stride_id) ) // Define stride[].

   nCDF( nc_def_var(ncid, "radius", NC_INT, 1, &dim_id, &radius_id) ) // Define radius[].
   nCDF( nc_def_var(ncid, "source", NC_INT, 1, &dim_id, &source_id) ) // Define radius[].
   nCDF( nc_def_var(ncid, "state", NC_INT, 1, &dim_id, &state_id) ) // Define state[] = isFinish+nloop.

	nCDF( nc_def_var(ncid, NAME_DATA, NC_TyPe, DIM, coord_id, &slowness_id) ) // Define Slowness map.
	nCDF( nc_def_var(ncid, "Traveltime", NC_TyPe, DIM, coord_id, &traveltime_id) ) // Define Traveltime map.
	if( !isFinish ){
      nCDF( nc_def_var(ncid, "UpdateMap", NC_INT, DIM, group_id, &update_id) ) // Define Update map.
	}
	nCDF( nc_put_att_text(ncid, slowness_id, "Unit", std::strlen(slownessUnit), slownessUnit) )
	nCDF( nc_enddef(ncid) ) // Finish all definitions.

	// Write scale and data
	nCDF( nc_put_var_float(ncid, min_id, &min[0]) )
	nCDF( nc_put_var_float(ncid, max_id, &max[0]) )
	nCDF( nc_put_var_float(ncid, stride_id, stride.data()) )

    // Uint is only avaliable after NetCDF 4, so it's a little bit clumsy. Only 3D ###
    int outputRadius[] = { int(block[0]/2), int(block[1]/2), int(block[2]/2) } ;
	nCDF( nc_put_var_int(ncid, radius_id, &outputRadius[0]) )
	int outputSource[] = { int(source[0]), int(source[1]), int(source[2]) } ;
	nCDF( nc_put_var_int(ncid, source_id, &outputSource[0]) )
	int outputState[] = { int(isFinish), int(nloop), 0 } ;
	nCDF( nc_put_var_int(ncid, state_id, &outputState[0]) )

	nCDF( nc_put_var_float(ncid, slowness_id, slowness.data() ) )
	nCDF( nc_put_var_float(ncid, traveltime_id, traveltime.data() ) )
	if( !isFinish ){
      nCDF( nc_put_var_int(ncid, update_id, updateMap.data() ) )
    }
	// close the NetCDF file
	nCDF( nc_close(ncid) )

	// Undefine NetCDF error handling
	#ifdef nCDF
	  #undef nCDF
	#endif

    // Clear all data
    slowness.clear();
    traveltime.clear();
    updateMap.clear();


return true ; }

/************************************* Private *************************************/

bool PSPS_CUDA_FIM_3D_Tools::readSlownessModel(const char* filePath) // Read res and slowness.
{
	// NetCDF error handling
	#ifndef nCDF
	  #define nCDF(Func) { int NC_Stat; if( NC_NOERR != (NC_Stat=Func) ){ fprintf(stderr, "ERROR::NETCDF:: %s\n", nc_strerror(NC_Stat)); std::exit(2);} }
	#else
	  #undef nCDF
	  #define nCDF(Func) { int NC_Stat; if( NC_NOERR != (NC_Stat=Func) ){ fprintf(stderr, "ERROR::NETCDF:: %s\n", nc_strerror(NC_Stat)); std::exit(2);} }
	#endif

	// Open NetCDF file
	int ncid ;
   nCDF( nc_open(filePath, NC_NOWRITE, &ncid) )

	// Check dimension
	int ndims ;
   nCDF( nc_inq(ncid,&ndims,NULL,NULL,NULL) )
   size_t tempDim ;
   nCDF( nc_inq_dim(ncid,0,NULL,&tempDim) )
   if( !(ndims-1==DIM || ndims-1==2*DIM) ){ // ndims = Dimension(1) + Coordinates(xDIM) (+ Group(xDIM) )
      std::cout << "NewProject >> Sorry, File format or Dimension are not valid in this program.\n" ;
      return false;
   }
   // Get coordinates
   for(unsigned int i=0; i<DIM; ++i )
   {
      char coordName[NC_MAX_NAME+1];
      size_t tempRes ;
      nCDF( nc_inq_dim(ncid,i+1,coordName,&tempRes) )
      res[i] = (unsigned int) tempRes ;
      std::cout << "NewProject >> Coordinate " << i+1
                << " (" << coordName << ") with "
                << res[i] << " data points.\n" ;
   }

   // Calculate number of vertices
    unsigned int numData = 1 ; 
    for(unsigned int i=0 ; i<DIM ; ++i)
	{
		numData *= res[i] ;
	}

   // Allocate memory
   slowness.resize( numData );

   // Get stride
   int varid ;
   nCDF( nc_inq_varid(ncid, "strideCoord", &varid) )
   nCDF( nc_get_var_float(ncid, varid, stride.data() ) )
   // Get Slowness map
   nCDF( nc_inq_varid(ncid, NAME_DATA, &varid) )
   nCDF( nc_get_var_float(ncid, varid, slowness.data() ) )

	// close the NetCDF file
	nCDF( nc_close(ncid) )

	// Undefine NetCDF error handling
	#ifdef nCDF
	  #undef nCDF
	#endif

return true; }

bool PSPS_CUDA_FIM_3D_Tools::inputSpec(std::vector<int> &param)
{
   block[0]  = std::abs( param[0] ) ;
   block[1]  = std::abs( param[1] ) ;
   block[2]  = std::abs( param[2] ) ;
   res[0]    = std::abs( param[3] ) ;
   res[1]    = std::abs( param[4] ) ;
   res[2]    = std::abs( param[5] ) ;
   source[0] = std::abs( param[6] ) ;
   source[1] = std::abs( param[7] ) ;
   source[2] = std::abs( param[8] ) ;
   maxLoop   = std::abs( param[9] ) ;
   
   for(unsigned int i=0 ; i<DIM ; ++i)
   {
      std::cout << "NewProject >> Coordinate " << i+1
                << " will have number of cells "
                << res[i] << ".\n" ;
      std::cout << "NewProject >> Coordinate " << i+1
                << " has \"Block\" = " << block[i]  << ".\n";
      // Check and Compute nWorkGroup
      if( block[i]%2 != 0 ){
          std::cout << "NewProject >> A \"Block\" is invalid, it MUST be divisible by 2 !!!\n" ;
          return false ;
      }
      if( res[i]%block[i] != 0 ){
          std::cout << "NewProject >> A number of cells is invalid, it MUST be divisible by Block !!!\n" ;
          return false ;
      }else{
          nWorkGroup[i] = (unsigned int) res[i]/block[i] ;
      }
   }
   std::cout << "NewProject >> The \"Source\" is at ("
             << source[0]+1 << "," << source[1]+1 << "," << source[2]+1 << ").\n" ; // In C [0,res-1] -> Display [1,res] 

return true ; }

bool PSPS_CUDA_FIM_3D_Tools::readResumeProject(const char* filePath) // Read all maps and Setup.
{
    bool flag = true ;

	// NetCDF error handling
	#ifndef nCDF
	  #define nCDF(Func) { int NC_Stat; if( NC_NOERR != (NC_Stat=Func) ){ fprintf(stderr, "ERROR::NETCDF:: %s\n", nc_strerror(NC_Stat)); std::exit(2);} }
	#else
	  #undef nCDF
	  #define nCDF(Func) { int NC_Stat; if( NC_NOERR != (NC_Stat=Func) ){ fprintf(stderr, "ERROR::NETCDF:: %s\n", nc_strerror(NC_Stat)); std::exit(2);} }
	#endif

	// Open NetCDF file
	int ncid ;
   nCDF( nc_open(filePath, NC_NOWRITE, &ncid) )

	// Check dimension
   int ndims ;
   nCDF( nc_inq(ncid,&ndims,NULL,NULL,NULL) )
   size_t tempDim ;
   nCDF( nc_inq_dim(ncid,0,NULL,&tempDim) )
   if( !(ndims-1==DIM || ndims-1==2*DIM) ){ // ndims = Dimension(1) + Coordinates(xDIM) (+ Group(xDIM) )
      std::cout << "ResumeProject >> Sorry, File format or Dimension are not valid in this program.\n" ;
      return false;
   }
   // Get coordinates
   for(unsigned int i=0; i<DIM; ++i )
   {
      char coordName[NC_MAX_NAME+1];
      size_t tempRes ;
      nCDF( nc_inq_dim(ncid,i+1,coordName,&tempRes) )
      res[i] = (unsigned int) tempRes ;
      std::cout << "ResumeProject >> Coordinate " << i+1
                << " (" << coordName << ") with number of cells "
                << res[i] << ".\n" ;
   }

   // Get radius, source and state of the project.
   int varid ;
   int inputRadius[DIM], inputSource[DIM], inputState[DIM] ;
   nCDF( nc_inq_varid(ncid, "radius", &varid) )
   nCDF( nc_get_var_int(ncid, varid, &inputRadius[0] ) )
   nCDF( nc_inq_varid(ncid, "source", &varid) )
   nCDF( nc_get_var_int(ncid, varid, &inputSource[0] ) )
   nCDF( nc_inq_varid(ncid, "state", &varid) )
   nCDF( nc_get_var_int(ncid, varid, &inputState[0] ) ) // [isFinish, nloop]

   // Uint is only avaliable after NetCDF 4, so it's a little bit clumsy here !!.
   if( inputState[0] != 0 ){
      std::cout << "ResumeProject >> This project is marked as \"Finish\".\n" ;
      return false;
   }else{
      isFinish = false ;
   }
   nloop = std::abs( inputState[1] ) ;
   for(unsigned int i=0 ; i<DIM ; ++i)
   {
       block[i] = 2*std::abs( inputRadius[i] ) ;
      source[i] = std::abs( inputSource[i] ) ;
      // Check and Compute nWorkGroup
      if( res[i]%block[i] != 0 ){
        std::cout << "ResumeProject >> A number of cells is invalid, it MUST be divided by 2*radius !!!\n" ;
        flag = false ;         
        break;
      }else{
         nWorkGroup[i] = (unsigned int) res[i]/block[i] ;
      }
      std::cout << "ResumeProject >> Coordinate " << i+1
                << " has \"Radius\" = " << block[i]/2  << ".\n";
   }
   std::cout << "ResumeProject >> And the \"Source\" is at ("
             << source[0]+1 << "," << source[1]+1 << "," << source[2]+1 << ").\n" ; // In C [0,res-1] -> Display [1,res]

   // Calculate number of vertices
   unsigned int numData = 1 , numGroup = 1 ;
   for(unsigned int i=0 ; i<DIM ; ++i)
   {
       numData *= res[i] ;
      numGroup *= nWorkGroup[i] ;
   }

   // Allocate memory
   slowness.resize( numData );
   traveltime.resize( numData );
   updateMap.resize( numGroup );

   // Get stride
   nCDF( nc_inq_varid(ncid, "strideCoord", &varid) )
   nCDF( nc_get_var_float(ncid, varid, stride.data() ) )
   // Get all maps.
   nCDF( nc_inq_varid(ncid, NAME_DATA, &varid) )
   nCDF( nc_get_var_float(ncid, varid, slowness.data() ) )
   nCDF( nc_inq_varid(ncid, "Traveltime", &varid) )
   nCDF( nc_get_var_float(ncid, varid, traveltime.data() ) )
   nCDF( nc_inq_varid(ncid, "UpdateMap", &varid) )
   nCDF( nc_get_var_int(ncid, varid, updateMap.data() ) )

	// close the NetCDF file
	nCDF( nc_close(ncid) )

	// Undefine NetCDF error handling
	#ifdef nCDF
	  #undef nCDF
	#endif

return flag; }

bool PSPS_CUDA_FIM_3D_Tools::checkSpec()
{
   bool flag = true ;

   /// Check limitation on this computer.
   // Check Radius in each direction.
   for(unsigned int i=0 ; i<DIM ; ++i)
   {
      if( source[i] >= res[i] ){
         std::cout << "CheckSpec >> The \"Source\" is outside the domain (resolution) !!!\n" ;
         flag = false ; 
         break; 
      }
      if( block[i] < 3 ){
        std::cout << "CheckSpec >> Update map will ERROR, if \"Block\" < 3  !!!\n" ;
        flag = false ;         
        break; 
      }
   }

return flag; }


texture<float, cudaTextureType3D, cudaReadModeElementType>  inMap ;
surface<void, cudaSurfaceType3D> outMap ;

__global__ void cudaRescaleSlowness()
{
	// Calculate normalized texture coordinates 
	unsigned int x = blockIdx.x * blockDim.x + threadIdx.x ;  
	unsigned int y = blockIdx.y * blockDim.y + threadIdx.y ; 
	unsigned int z = blockIdx.z * blockDim.z + threadIdx.z ; 
    int resX = gridDim.x*blockDim.x ;
    int resY = gridDim.y*blockDim.y ;
	int resZ = gridDim.z*blockDim.z ;
	float u = x/__int2float_rn(resX) ; 
	float v = y/__int2float_rn(resY) ;
	float r = z/__int2float_rn(resZ) ;
	 
	if( x<resX && y<resY && z<resZ  )
	{
		float data = tex3D(inMap,u,v,r);
		surf3Dwrite(data, outMap, sizeof(float)*x, y, z); // sizeof(float)* ??
	}
	
}


bool PSPS_CUDA_FIM_3D_Tools::rescaleSlowness(std::vector<unsigned int> &res_old)
{
   std::cout << "\nNewProject >> START rescaling the slowness model.\n" ;

   /// *** If we magnify -> the two outmost boundaries will be repeated (the same). !!! ***///

   /// Create and Set Textures ///
 
   // Input map
   cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32,0,0,0, cudaChannelFormatKindFloat); 
   cudaArray* cuArrayIn ;
   cudaMalloc3DArray(&cuArrayIn, &channelDesc, make_cudaExtent(res_old[0], res_old[1], res_old[2]) );
   cudaMemcpy3DParms paramIn = {0} ;
      paramIn.srcPtr = make_cudaPitchedPtr(slowness.data(), res_old[0]*sizeof(float), res_old[0], res_old[1]) ;
      paramIn.dstArray = cuArrayIn ;
      paramIn.kind =  cudaMemcpyHostToDevice ;
      paramIn.extent = make_cudaExtent(res_old[0], res_old[1], res_old[2]) ;
   cudaMemcpy3D(&paramIn);

   inMap.addressMode[0] = cudaAddressModeMirror ; 
   inMap.addressMode[1] = cudaAddressModeMirror ; 
   inMap.addressMode[2] = cudaAddressModeMirror ; 
   inMap.filterMode = cudaFilterModeLinear; 
   inMap.normalized = true ;
   cudaBindTextureToArray(inMap, cuArrayIn, channelDesc);

   // Output map
   cudaArray* cuArrayOut ;
   cudaMalloc3DArray(&cuArrayOut, &channelDesc, make_cudaExtent(res[0], res[1], res[2]), cudaArraySurfaceLoadStore);
   cudaBindSurfaceToArray(outMap, cuArrayOut);
     
   /// Run CUDA ///
   dim3 dimBlock( block[0], block[1], block[2] ); 
   dim3 dimGroup( nWorkGroup[0], nWorkGroup[1], nWorkGroup[2] );
   cudaRescaleSlowness<<<dimGroup,dimBlock>>>();
   
   /// Retrieve data ///
   slowness.clear(); // Needed ?
   slowness.resize( res[0]*res[1]*res[2] );
   cudaMemcpy3DParms paramOut = {0} ;
      paramOut.srcArray = cuArrayOut ;
      paramOut.dstPtr = make_cudaPitchedPtr(slowness.data(), res[0]*sizeof(float), res[0], res[1]) ;
      paramOut.kind =  cudaMemcpyDeviceToHost ;
      paramOut.extent = make_cudaExtent(res[0], res[1], res[2]) ;
   cudaMemcpy3D(&paramOut);

   /// Clear shader and buffers ///
   cudaUnbindTexture(inMap);
   cudaFreeArray(cuArrayIn); 
   cudaFreeArray(cuArrayOut);

   std::cout << "NewProject >> FINISH rescaling the slowness model.\n" ;

   
return true; }
