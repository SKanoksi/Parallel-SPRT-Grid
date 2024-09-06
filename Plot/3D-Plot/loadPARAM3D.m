function [res, min, stride, radius, source, name] = loadPARAM3D(FILE_PATH,FILE_NAME)
% /************************************
%      Parallel Shortest Path Solver
%           (loadPARAM3D.m)
% *************************************/

% Open NetCDF file
FILE_PATH = strcat(FILE_PATH,FILE_NAME);
ncid = netcdf.open(FILE_PATH,'NOWRITE'); % format of NetCDF ###

% Get dimension
[name0, dim] = netcdf.inqDim(ncid,0); % Dimension
if dim~=3
    fprintf(' ERROR : the file, %s, is not for Plot2D.\n',FILE_PATH);
    return;
end
[name1, res1] = netcdf.inqDim(ncid,1); % Coord 1
[name2, res2] = netcdf.inqDim(ncid,2); % Coord 2
[name3, res3] = netcdf.inqDim(ncid,3); % Coord 3
name = {name1, name2, name3};
res  = [res1, res2, res3];
clear dim res1 res2 res3 name0 name1 name2 name3

% Get max, min, stride and data
    min = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'minCoord') );
   %max = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'maxCoord') );
 stride = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'strideCoord') );
 radius = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'radius') );
 source = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'source') );
 source = source + 1 ; % C -> MatLab

% Close NetCDF file
netcdf.close(ncid);

end

