/*************************************
     Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#include "PSPS_OpenGL_BF_2D_Tools.h"

#define TyPe GLfloat
#define NC_TyPe NC_FLOAT
#define DIM 2
#define NAME_DATA "Slowness"
#define GPU_MEMORY_LIMIT 0.90     // Will use at most 90% of total GPU memory. ###
#define SHARED_MEMORY_LIMIT 0.98  // Will use at most 98% of total GPU shared memory. ###

PSPS_OpenGL_BF_2D_Tools::PSPS_OpenGL_BF_2D_Tools()
{

}

PSPS_OpenGL_BF_2D_Tools::~PSPS_OpenGL_BF_2D_Tools()
{

}

bool PSPS_OpenGL_BF_2D_Tools::newProject(std::string &inFile, std::vector<int> &param)
{
   // Read slowness model.
   if( !readSlownessModel(inFile.c_str()) ){ return false; }

   // Input the spec.
   std::vector<GLuint> res_old(res); // Save old resolution
   if( !inputSpec(param) ){ return false; }
   if( !checkSpec() ){ return false; }
   
   // Change the scale
   for(GLuint i=0 ; i<DIM ; ++i)
   {
        stride[i] = (TyPe) (res_old[i]*stride[i])/res[i] ;
   }

   // Rescale Slowness map
   for(GLuint i=0 ; i<DIM ; ++i)
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
   GLuint numData = 1 , numGroup = 1 ;
   for(GLuint i=0 ; i<DIM ; ++i)
   {
       numData *= res[i] ;
      numGroup *= nWorkGroup[i] ;
   }
   traveltime.clear();
   traveltime.resize(  numData, 1.0f/0.0f);
      raypath.clear();
      raypath.resize(  numData, 0);
    updateMap.clear();
    updateMap.resize( numGroup, 0);

   // Set source, Traveltime+updateMap Only 2D ###
   int atMap = 0, atGroup = 0 ;
   for(int i=DIM-1 ; 0<=i ; --i)
   {
        atMap *= res[i];
        atMap += source[i] ;
      atGroup *= nWorkGroup[i] ;
      atGroup += (int) source[i]/block[i] ;
   }
	//  atMap = res[0]*source[1] + source[0] ;
	//atGroup = nWorkGroup[0]*((GLint)source[1]/block[1]) + ((GLint)source[0]/block[0]) ;
    traveltime[atMap] = 0.0f ;
   updateMap[atGroup] = 1 ;

   // Set setup
      nloop = 0 ;
   isFinish = false ;

return true; }

bool PSPS_OpenGL_BF_2D_Tools::resumeProject(std::string &inFile, int maxloop)
{
   // Read the resume project.
   maxLoop = maxloop ;
   if( !readResumeProject(inFile.c_str()) ){ return false; }
   if( !checkSpec() ){ return false; }

return true; }

bool PSPS_OpenGL_BF_2D_Tools::saveProject(std::string &inFile)
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
   for(GLuint i=0; i<DIM; ++i )
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
      inFile = inFile + "_PSPS_OpenGL_BF.out" ;
   }
   // Create a NetCDF file
   nCDF( nc_create(inFile.c_str(), NC_CLOBBER, &ncid) ) // (may want to select cmode). ###

	// Define dimension and coordinates
	int dim_id, coord_id[DIM], group_id[DIM] ;
	nCDF( nc_def_dim(ncid, "Dimension", DIM, &dim_id) ) // Define dimension.
	for(GLuint i=0 ; i<DIM ; ++i)
	{
		nCDF( nc_def_dim(ncid, name[i].c_str(), res[i], &coord_id[i]) ) // Define coordiantes
	}
	for(GLuint i=0 ; i<DIM ; ++i)
	{
			if( !isFinish ){
         name[i] += "_nWorkGroup" ;
         nCDF( nc_def_dim(ncid, name[i].c_str(), nWorkGroup[i], &group_id[i]) ) // for UpdateMap
      }
	}

	// Define variables
   int min_id, max_id, stride_id ;
	int radius_id, source_id, state_id ;
	int slowness_id, traveltime_id, raypath_id, update_id ;
	nCDF( nc_def_var(ncid, "minCoord", NC_TyPe, 1, &dim_id, &min_id) ) // Define min[].
	nCDF( nc_def_var(ncid, "maxCoord", NC_TyPe, 1, &dim_id, &max_id) ) // Define max[].
	nCDF( nc_def_var(ncid, "strideCoord", NC_TyPe, 1, &dim_id, &stride_id) ) // Define stride[].

   nCDF( nc_def_var(ncid, "radius", NC_INT, 1, &dim_id, &radius_id) ) // Define radius[].
   nCDF( nc_def_var(ncid, "source", NC_INT, 1, &dim_id, &source_id) ) // Define radius[].
   nCDF( nc_def_var(ncid, "state", NC_INT, 1, &dim_id, &state_id) ) // Define state[] = isFinish+nloop.

	nCDF( nc_def_var(ncid, NAME_DATA, NC_TyPe, DIM, coord_id, &slowness_id) ) // Define Slowness map.
	nCDF( nc_def_var(ncid, "Traveltime", NC_TyPe, DIM, coord_id, &traveltime_id) ) // Define Traveltime map.
	nCDF( nc_def_var(ncid, "Raypath", NC_INT, DIM, coord_id, &raypath_id) ) // Define Raypath map.
	if( !isFinish ){
      nCDF( nc_def_var(ncid, "UpdateMap", NC_INT, DIM, group_id, &update_id) ) // Define Update map.
	}
	nCDF( nc_put_att_text(ncid, slowness_id, "Unit", std::strlen(slownessUnit), slownessUnit) )
	nCDF( nc_enddef(ncid) ) // Finish all definitions.

	// Write scale and data
	nCDF( nc_put_var_float(ncid, min_id, &min[0]) )
	nCDF( nc_put_var_float(ncid, max_id, &max[0]) )
	nCDF( nc_put_var_float(ncid, stride_id, stride.data()) )

    // Uint is only avaliable after NetCDF 4, so it's a little bit clumsy. Only 2D ###
    int outputRadius[] = { int(block[0]/2), int(block[1]/2) } ;
	nCDF( nc_put_var_int(ncid, radius_id, &outputRadius[0]) )
	int outputSource[] = { int(source[0]), int(source[1]) } ;
	nCDF( nc_put_var_int(ncid, source_id, &outputSource[0]) )
	int outputState[] = { int(isFinish), int(nloop) } ;
	nCDF( nc_put_var_int(ncid, state_id, &outputState[0]) )

	nCDF( nc_put_var_float(ncid, slowness_id, slowness.data() ) )
	nCDF( nc_put_var_float(ncid, traveltime_id, traveltime.data() ) )
	nCDF( nc_put_var_int(ncid, raypath_id, raypath.data() ) )
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
    raypath.clear();
    updateMap.clear();


return true ; }

/************************************* Private *************************************/

bool PSPS_OpenGL_BF_2D_Tools::readSlownessModel(const char* filePath) // Read res and slowness.
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
   for(GLuint i=0; i<DIM; ++i )
   {
      char coordName[NC_MAX_NAME+1];
      size_t tempRes ;
      nCDF( nc_inq_dim(ncid,i+1,coordName,&tempRes) )
      res[i] = (GLuint) tempRes ;
      std::cout << "NewProject >> Coordinate " << i+1
                << " (" << coordName << ") with "
                << res[i] << " data points.\n" ;
   }

   // Calculate number of vertices
    GLuint numData = 1 ; 
    for(GLuint i=0 ; i<DIM ; ++i)
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

bool PSPS_OpenGL_BF_2D_Tools::inputSpec(std::vector<int> &param)
{
   block[0]  = std::abs( 2*param[0] ) ;
   block[1]  = std::abs( 2*param[1] ) ;
   res[0]    = std::abs( param[2] ) ;
   res[1]    = std::abs( param[3] ) ;
   source[0] = std::abs( param[4] ) ;
   source[1] = std::abs( param[5] ) ;
   maxLoop   = std::abs( param[6] ) ;
   
   for(GLuint i=0 ; i<DIM ; ++i)
   {
      std::cout << "NewProject >> Coordinate " << i+1
                << " will have number of cells "
                << res[i] << ".\n" ;
      std::cout << "NewProject >> Coordinate " << i+1
                << " has \"Radius\" = " << block[i]/2  << ".\n";
      // Check and Compute nWorkGroup
      if( res[i]%block[i] != 0 ){
          std::cout << "NewProject >> A number of cells is invalid, it MUST be divisible by 2*radius !!!\n" ;
          return false ;
      }else{
          nWorkGroup[i] = (GLuint) res[i]/block[i] ;
      }
   }
   std::cout << "NewProject >> The \"Source\" is at ("
             << source[0]+1 << "," << source[1]+1 << ").\n" ; // In C [0,res-1] -> Display [1,res] 

return true ; }

bool PSPS_OpenGL_BF_2D_Tools::readResumeProject(const char* filePath) // Read all maps and Setup.
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
   for(GLuint i=0; i<DIM; ++i )
   {
      char coordName[NC_MAX_NAME+1];
      size_t tempRes ;
      nCDF( nc_inq_dim(ncid,i+1,coordName,&tempRes) )
      res[i] = (GLuint) tempRes ;
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
   for(GLuint i=0 ; i<DIM ; ++i)
   {
       block[i] = 2*std::abs( inputRadius[i] ) ;
      source[i] = std::abs( inputSource[i] ) ;
      // Check and Compute nWorkGroup
      if( res[i]%block[i] != 0 ){
        std::cout << "ResumeProject >> A number of cells is invalid, it MUST be divided by 2*radius !!!\n" ;
        flag = false ;         
        break;
      }else{
         nWorkGroup[i] = (GLuint) res[i]/block[i] ;
      }
      std::cout << "ResumeProject >> Coordinate " << i+1
                << " has \"Radius\" = " << block[i]/2  << ".\n";
   }
   std::cout << "ResumeProject >> And the \"Source\" is at ("
             << source[0]+1 << "," << source[1]+1 << ").\n" ; // In C [0,res-1] -> Display [1,res]

   // Calculate number of vertices
   GLuint numData = 1 , numGroup = 1 ;
   for(GLuint i=0 ; i<DIM ; ++i)
   {
       numData *= res[i] ;
      numGroup *= nWorkGroup[i] ;
   }

   // Allocate memory
   slowness.resize( numData );
   traveltime.resize( numData );
   raypath.resize( numData );
   updateMap.resize( numGroup );

   // Get stride
   nCDF( nc_inq_varid(ncid, "strideCoord", &varid) )
   nCDF( nc_get_var_float(ncid, varid, stride.data() ) )
   // Get all maps.
   nCDF( nc_inq_varid(ncid, NAME_DATA, &varid) )
   nCDF( nc_get_var_float(ncid, varid, slowness.data() ) )
   nCDF( nc_inq_varid(ncid, "Traveltime", &varid) )
   nCDF( nc_get_var_float(ncid, varid, traveltime.data() ) )
   nCDF( nc_inq_varid(ncid, "Raypath", &varid) )
   nCDF( nc_get_var_int(ncid, varid, raypath.data() ) )
   nCDF( nc_inq_varid(ncid, "UpdateMap", &varid) )
   nCDF( nc_get_var_int(ncid, varid, updateMap.data() ) )

	// close the NetCDF file
	nCDF( nc_close(ncid) )

	// Undefine NetCDF error handling
	#ifdef nCDF
	  #undef nCDF
	#endif

return flag; }

bool PSPS_OpenGL_BF_2D_Tools::checkSpec()
{
   bool flag = true ;

   /// Get Computer Spec (OpenGL_BF)
   GLint n[DIM], l[DIM], nInvo, shMem, gpuMem, sizeTex ;
   for(GLuint i=0 ; i<DIM ; ++i)
   {
       glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, i, &n[i]);
       glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_SIZE, i, &l[i]);
   }
   glGetIntegerv(GL_MAX_COMPUTE_WORK_GROUP_INVOCATIONS, &nInvo);
   glGetIntegerv(GL_MAX_COMPUTE_SHARED_MEMORY_SIZE, &shMem);
   glGetIntegerv(GL_MAX_SHADER_STORAGE_BLOCK_SIZE, &gpuMem);
   glGetIntegerv(GL_MAX_TEXTURE_SIZE, &sizeTex);
   shMem  = (GLint) shMem*SHARED_MEMORY_LIMIT ;
   gpuMem = (GLint) gpuMem*GPU_MEMORY_LIMIT   ;

   /// Check limitation on this computer.
   // Check Radius in each direction.
   for(GLuint i=0 ; i<DIM ; ++i)
   {
      if( block[i] > (GLuint)l[i] ){
         std::cout << "CheckSpec >> Cannot run on this Computer, a 2*\"Radius\" is larger than GL_MAX_COMPUTE_WORK_GROUP_SIZE !!!\n" ;
         flag = false ;         
         break; 
      }
      if( res[i] > (GLuint)sizeTex || res[i] > (GLuint)n[i] ){
         std::cout << "CheckSpec >> Cannot run on this Computer, a resolution value is larger than GL_MAX_TEXTURE_SIZE or GL_MAX_COMPUTE_WORK_GROUP_COUNT !!!\n" ;
         flag = false ;          
         break;
      }
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

   // Check limitation.
   GLuint WorkGroup_Size = 1, shGroup_Size = 1 ;
   for(GLuint i=0 ; i<DIM ; ++i)
   {
      WorkGroup_Size *= block[i] ;
      shGroup_Size *= 2*block[i] ;
   }
   if( WorkGroup_Size > (GLuint)nInvo ){
      std::cout << "CheckSpec >> Invalid \"Radius\" value, bigger than GL_MAX_COMPUTE_WORK_GROUP_INVOCATIONS !!!\n" ;
      return false ;
   }
   if( 2*sizeof(TyPe)*shGroup_Size >= (GLuint)shMem ){ // Slowness + Travel time
      std::cout << "CheckSpec >> Invalid \"Radius\" value, need shared memory larger than GL_MAX_COMPUTE_SHARED_MEMORY_SIZE !!!\n" ;
      return false ;
   }

   // Only for 2D ###
   GLuint needMem = 0 ;
   needMem += 3*sizeof(TyPe)*(res[0]+block[0])*(res[1]+block[1]) ; // ( Slowness+Traveltime x2 )*( res+dummy)
   needMem += sizeof(GLint)*res[0]*res[1] ; // Raypath
   needMem += 2*sizeof(GLint)*(nWorkGroup[0]+2)*(nWorkGroup[1]+2) ; // Update map x2 (but will use only 0,1)
   if( needMem >= (GLuint)gpuMem ){
      std::cout << "ResumeProject >> Your GPU memory limitation is met !!!\n" ;
      return false;
   }


return flag; }


bool PSPS_OpenGL_BF_2D_Tools::rescaleSlowness(std::vector<GLuint> &res_old)
{
   std::cout << "NewProject >> START rescaling the slowness model.\n" ;

   /// *** If we magnify -> the two outmost boundaries will be repeated (the same). !!! ***///

   /// Create Compute shader. ///
   // Set parameters in the ComSh code.
   std::string comCode, temp ;
   if( !readSource("./PSPS_Rescale2D.comSh",comCode) ){ // The code must be found ###
       std::cout << "ERROR:: PSPS_Rescale2D.comSh is not found in the same or at the parent folder.\n" ;
       return false;
   }
     ToString(block[0],temp);
   comCode.replace( comCode.find("BLOCK_X"), 7, temp.c_str() ); // The marks MUST be in the code ###
     ToString(block[1],temp);
   comCode.replace( comCode.find("BLOCK_Y"), 7, temp.c_str() );
     ToString(res[0],temp); temp += ".0f" ;
   comCode.replace( comCode.find("RES_X"), 5, temp.c_str() ); // The marks MUST be in the code ###
     ToString(res[1],temp); temp += ".0f" ;
   comCode.replace( comCode.find("RES_Y"), 5, temp.c_str() );
   // Compile code
   GLuint rescaleComSh = createComSh(comCode.c_str(), comCode.length());

   /// Create and Set Textures ///
   // Create input and output textures.
   GLuint inTexture, outTexture ;
   glGenTextures(1, &inTexture);
   glGenTextures(1, &outTexture);

   // Set input texture.
   glBindTexture(GL_TEXTURE_2D, inTexture);
      // Input image.
      glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, res_old[0], res_old[1], 0, GL_RED, GL_FLOAT, (GLvoid*)&slowness[0]);
      // Set texture.
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      //glGenerateMipmap(GL_TEXTURE_2D);
   //glBindTexture(GL_TEXTURE_2D, 0);

   // Set output texture.
   glBindTexture(GL_TEXTURE_2D, outTexture);
      glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, res[0], res[1], 0, GL_RED, GL_FLOAT, NULL);
      // Needed by Nvidia, even if not used. !!!
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
   glBindTexture(GL_TEXTURE_2D, 0);

   /// Run Compute shader ///
   // Active texture.
   glActiveTexture(GL_TEXTURE0);
   glBindTexture(GL_TEXTURE_2D, inTexture);
   // Run compute shader.
   glUseProgram(rescaleComSh);
      glUniform1i(glGetUniformLocation(rescaleComSh, "inMap"), 0);
      glBindImageTexture(0, outTexture, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_R32F);
         glDispatchCompute( nWorkGroup[0], nWorkGroup[1], 1);
      glMemoryBarrier(GL_ALL_BARRIER_BITS);
   glUseProgram(0);
   // Deactive texture.
   glBindTexture(GL_TEXTURE_2D, 0);

   /// Retrieve data ///
   glBindTexture(GL_TEXTURE_2D, outTexture);
     slowness.clear(); // Needed ?
     slowness.resize( res[0]*res[1] );
     glGetTexImage(GL_TEXTURE_2D, 0, GL_RED, GL_FLOAT, &slowness[0]);
   glBindTexture(GL_TEXTURE_2D, 0);

   /// Clear shader and buffers ///
   glDeleteProgram(rescaleComSh);
   glDeleteBuffers(1, &inTexture);
   glDeleteBuffers(1, &outTexture);

   std::cout << "NewProject >> FINISH rescaling the slowness model.\n" ;

   /// Check glGet "Error" -> return false; ###

return true; }
