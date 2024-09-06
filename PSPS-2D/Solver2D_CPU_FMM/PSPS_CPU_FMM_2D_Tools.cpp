/*************************************
     Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

#include "PSPS_CPU_FMM_2D_Tools.h"

#define TyPe float
#define NC_TyPe NC_FLOAT
#define DIM 2
#define NAME_DATA "Slowness"

PSPS_CPU_FMM_2D_Tools::PSPS_CPU_FMM_2D_Tools()
{

}

PSPS_CPU_FMM_2D_Tools::~PSPS_CPU_FMM_2D_Tools()
{

}

bool PSPS_CPU_FMM_2D_Tools::newProject(std::string &inFile, std::vector<int> &param)
{
   // Read slowness model.
   if( !readSlownessModel(inFile.c_str()) ){ return false; }

   // Input the spec.
   std::vector<unsigned int> res_old(res); // Save old resolution
   if( !inputSpec(param) ){ return false; }
   if( !checkSpec() ){ return false; }

   // Change the scale
   for(int i=0 ; i<DIM ; ++i)
   {
        stride[i] = (float) (res_old[i]*stride[i])/res[i] ;
   }
   if( std::fabs(stride[0]-stride[1])>0.01*stride[1] )
   {
        std::cout << "\n!!! Difference of grid specing (hx,hy) is more than 1 percent. !!!\n\n" ;
        return false;
   }
   
   // Reallocate memory for other maps.
   unsigned int numData = 1 ;
   for(unsigned int i=0 ; i<DIM ; ++i)
   {
       numData *= res[i] ;
   }
   traveltime.clear();
   traveltime.resize(  numData, 1.0f/0.0f);
      raypath.clear();
      raypath.resize(  numData, -1); // Any negative numbers.

   // Set source, Traveltime Only 2D ###
   int atMap = 0 ;
   for(int i=DIM-1 ; 0<=i ; --i)
   {
        atMap *= res[i];
        atMap += source[i] ;
   }
   // atMap = res[0]*source[1] + source[0] ;
    traveltime[atMap] = 0.0f ;


return true; }

bool PSPS_CPU_FMM_2D_Tools::saveProject(std::string &inFile)
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
      inFile = inFile + "_PSPS_CPU_FMM.out" ;
   }
   // Create a NetCDF file
   nCDF( nc_create(inFile.c_str(), NC_CLOBBER, &ncid) ) // (may want to select cmode). ###

	// Define dimension and coordinates
	int dim_id, coord_id[DIM] ;
	nCDF( nc_def_dim(ncid, "Dimension", DIM, &dim_id) ) // Define dimension.
	for(unsigned int i=0 ; i<DIM ; ++i)
	{
		nCDF( nc_def_dim(ncid, name[i].c_str(), res[i], &coord_id[i]) ) // Define coordiantes
	}

	// Define variables
   int min_id, max_id, stride_id ;
	int radius_id, source_id ;
	int slowness_id, traveltime_id, raypath_id ;
	nCDF( nc_def_var(ncid, "minCoord", NC_TyPe, 1, &dim_id, &min_id) ) // Define min[].
	nCDF( nc_def_var(ncid, "maxCoord", NC_TyPe, 1, &dim_id, &max_id) ) // Define max[].
	nCDF( nc_def_var(ncid, "strideCoord", NC_TyPe, 1, &dim_id, &stride_id) ) // Define stride[].

   nCDF( nc_def_var(ncid, "radius", NC_INT, 1, &dim_id, &radius_id) ) // Define radius[].
   nCDF( nc_def_var(ncid, "source", NC_INT, 1, &dim_id, &source_id) ) // Define radius[].

	nCDF( nc_def_var(ncid, NAME_DATA, NC_TyPe, DIM, coord_id, &slowness_id) ) // Define Slowness map.
	nCDF( nc_def_var(ncid, "Traveltime", NC_TyPe, DIM, coord_id, &traveltime_id) ) // Define Traveltime map.
	nCDF( nc_def_var(ncid, "Raypath", NC_INT, DIM, coord_id, &raypath_id) ) // Define Raypath map.
	nCDF( nc_put_att_text(ncid, slowness_id, "Unit", std::strlen(slownessUnit), slownessUnit) )
	nCDF( nc_enddef(ncid) ) // Finish all definitions.

	// Write scale and data
	nCDF( nc_put_var_float(ncid, min_id, &min[0]) )
	nCDF( nc_put_var_float(ncid, max_id, &max[0]) )
	nCDF( nc_put_var_float(ncid, stride_id, stride.data()) )

   // Uint is only avaliable after NetCDF 4, so it's a little bit clumsy.
   int outputRadius[] = { int(block[0]/2), int(block[1]/2) } ;
	nCDF( nc_put_var_int(ncid, radius_id, &outputRadius[0]) )
	int outputSource[] = { int(source[0]), int(source[1]) } ;
	nCDF( nc_put_var_int(ncid, source_id, &outputSource[0]) )

	nCDF( nc_put_var_float(ncid, slowness_id, slowness.data() ) )
	nCDF( nc_put_var_float(ncid, traveltime_id, traveltime.data() ) )
	nCDF( nc_put_var_int(ncid, raypath_id, raypath.data() ) )

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


return true ; }

/************************************* Private *************************************/

bool PSPS_CPU_FMM_2D_Tools::readSlownessModel(const char* filePath) // Read res and slowness.
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
      std::cout << "NewProject >> Coordinate: " << i+1
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


bool PSPS_CPU_FMM_2D_Tools::inputSpec(std::vector<int> &param)
{
   block[0]  = std::abs( 2*param[0] ) ;
   block[1]  = std::abs( 2*param[1] ) ;
   source[0] = std::abs( param[2] ) ;
   source[1] = std::abs( param[3] ) ;
   
   for(int i=0 ; i<DIM ; ++i)
   {
      std::cout << "NewProject >> Coordinate " << i+1
                << " will have number of cells "
                << res[i] << ".\n" ;
      std::cout << "NewProject >> Coordinate " << i+1
                << " has \"Radius\" = " << block[i]/2  << ".\n";      
   }
   std::cout << "NewProject >> The \"Source\" is at ("
             << source[0]+1 << "," << source[1]+1 << ").\n" ; // In C [0,res-1] -> Display [1,res] 

return true ; }

bool PSPS_CPU_FMM_2D_Tools::checkSpec()
{
   bool flag = true ;
   // Check Radius in each direction.
   for(unsigned int i=0 ; i<DIM ; ++i)
   {
      if( source[i] >= res[i] ){
         std::cout << "ResumeProject >> The \"Source\" is outside the domain (resolution) !!!\n" ;
        flag = false ;         
        break; 
      }
   }

return flag ; }



