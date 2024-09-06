function [data,res] = loadPSPS3D(FILE_PATH,FILE_NAME,NAME_DATA)
% /************************************
%      Parallel Shortest Path Solver
%  loadPSPS3D.m = load data from PSPS.nc 
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
res  = [res1, res2, res3];
clear dim res1 res2 res3 name0 name1 name2 name3

% Get data
data_id = netcdf.inqVarID(ncid,NAME_DATA);
   data = netcdf.getVar(ncid, data_id);

% Close NetCDF file
netcdf.close(ncid);

% Rearrange data (Important)
data = reshape(data,res);

end

