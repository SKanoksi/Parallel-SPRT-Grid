/*************************************
    Parallel Shortest Path Solver
             developed by
 Somrath Kanoksirirath
 Fourth year physics student,
 Faculty of Science,
 Mahidol University, Thailand.
*************************************/

// NetCDF
#include <netcdf.h>
// Standard library
#include <vector>
#include <cstring>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>


/******************  Create an ideal case *********************/

#define FILE_NAME "Rgrad2D_2048.nc" // Define file name ###
#define NAME_DATA "Slowness"
#define UNIT_DATA "s/m"      // Set unit of the result.

/* Homogeneous */
//#define Function(x) 1.0f

/* 2D Z-Gradient */
//#define Function(x) 1.0f/x[1]

/* 2D R-Gradient */
#define Function(x) 1.0f/std::sqrt( x[0]*x[0]+x[1]*x[1] )

/* 3D Z-Gradient */
//#define Function(x) 1.0f/x[2]


// An alternative way to define the inputFunction.
/*
inline float Function(float *x)
{
	float output = 0 ;
   if( x[1]<300.0f ){
      output = 1.5f ;
   }else{
      output = 1.0f ;
   }

return output; }
*/

int main()
{
	// Define dimensions ###
	const unsigned int dim = 2 ; // Set the number of dimension. ###
	float min[dim] ;
	float max[dim] ;
	unsigned int res[dim] ;
	std::string name[dim] ;

	// Set the coordinates ### (In NetCDF: The last, do first.)
	min[0]  =     0.0f ;
	max[0]  =   500.0f ;
	res[0]  =   2048   ;
	name[0] =   "X-Meter" ;

	min[1]  =     0.0f ;
	max[1]  =   500.0f ;
	res[1]  =   2048   ;
	name[1] =   "Y-Meter" ;

	//min[2]  =  1000.0f ;
	//max[2]  =  4000.0f ;
	//res[2]  =   128    ;
	//name[2] =   "Z-Meter" ;


	/******************  Prepare data *********************/

	// Calculate stride and the number of data points
	float stride[dim] ;
	unsigned int numData = 1 ;
	for(unsigned int i=0 ; i<dim ; ++i)
	{
		stride[i] = (max[i]-min[i])/res[i] ;
		numData *= res[i] ;
	}

	// Allocate memory
	float *data ;
	try{
		data = new float[numData] ;
	}
	catch( std::bad_alloc& ba)
	{
		std::cerr << "Bad_alloc caught: " << ba.what() << std::endl;
	}

	// Calculate data using the input Function
	for(unsigned int num=0 ; num<numData ; ++num)
	{
		// Number to Index to Data
		int N = num ; // N should equal to 0, when finish the loop below.
		float x[dim] ;
		unsigned int index ;
		for(int i=0 ; i<dim ; ++i) // In NetCDF: the last, do first !!!
		{
			index = N % res[i] ;
			    N = (N-index)/res[i] ;
			 x[i] = min[i]+(index+0.5f)*stride[i] ;
		}
		data[num] = Function(x) ;

	}



	/******************  Write NetCDF *********************/

	// NetCDF error handling
	#ifndef nCDF
	  #define nCDF(Func) { int NC_Stat; if( NC_NOERR != (NC_Stat=Func) ){ fprintf(stderr, "ERROR::NETCDF:: %s\n", nc_strerror(NC_Stat)); std::exit(2);} }
	#else
	  #undef nCDF
	  #define nCDF(Func) { int NC_Stat; if( NC_NOERR != (NC_Stat=Func) ){ fprintf(stderr, "ERROR::NETCDF:: %s\n", nc_strerror(NC_Stat)); std::exit(2);} }
	#endif

	// Create NetCDF file
	int ncid ;
   nCDF( nc_create(FILE_NAME, NC_CLOBBER, &ncid) ) // (may want to select cmode). ###

	// Define dimension and coordinates
	int dim_id, coord_id[dim], min_id, max_id, stride_id, data_id ;
	nCDF( nc_def_dim(ncid, "Dimension", dim, &dim_id) ) // Define dimension.
	for(unsigned int i=0 ; i<dim ; ++i)
	{
		nCDF( nc_def_dim(ncid, name[i].c_str(), res[i], &coord_id[i]) ) // Define coordiantes.
	}
	// Define variables
	nCDF( nc_def_var(ncid, "minCoord", NC_FLOAT, 1, &dim_id, &min_id) ) // Define min[].
	nCDF( nc_def_var(ncid, "maxCoord", NC_FLOAT, 1, &dim_id, &max_id) ) // Define max[].
	nCDF( nc_def_var(ncid, "strideCoord", NC_FLOAT, 1, &dim_id, &stride_id) ) // Define stride[].
	nCDF( nc_def_var(ncid, NAME_DATA, NC_FLOAT, dim, coord_id, &data_id) ) // Define data.
	nCDF( nc_put_att_text(ncid, data_id, "Unit", std::strlen(UNIT_DATA), UNIT_DATA) )
	nCDF( nc_enddef(ncid) ) // Finish all definitions.

	// Write scale and data
	nCDF( nc_put_var_float(ncid, min_id, &min[0]) )
	nCDF( nc_put_var_float(ncid, max_id, &max[0]) )
	nCDF( nc_put_var_float(ncid, stride_id, &stride[0]) )
	nCDF( nc_put_var_float(ncid, data_id, &data[0]) )

	// close the NetCDF file
	nCDF( nc_close(ncid) )

	// Undefine NetCDF error handling
	#ifdef nCDF
	  #undef nCDF
	#endif

   std::cout << "SUCCESS : Writing an ideal case, " << FILE_NAME << ".\n" ;

return 0; }



