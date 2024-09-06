function [rayX, rayY, rayZ] = raypath3D(data, radius, source, nx,ny,nz)
% /************************************
%      Parallel Shortest Path Solver
% (raypath3D.m used in Plot_3D.m)
% *************************************/

% Set initial and some numbers
res = size(data);
block = 2*radius+1 ;
rayX = nx ;
rayY = ny ;
rayZ = nz ;

% At most, passing all the vertices
for i = 1:res(1)*res(2)*res(3)
    
    % Recursive backtracking (tracing back)
    num = data(nx,ny,nz) ;         % = (nz*block(2)+ny)*block(1)+nx
    shiftX = mod(num,block(1)) ;
    num = (num-shiftX)/block(1) ;  % = (nz*block(2)+ny)
    shiftY = mod(num,block(2));
    shiftZ = (num-shiftY)/block(2);
    
    % Only for Parallel version ###
    nx = nx + shiftX - radius(1) ;
    ny = ny + shiftY - radius(2) ;
    nz = nz + shiftZ - radius(3) ;
    rayX = [rayX, nx] ;
    rayY = [rayY, ny] ;
    rayZ = [rayZ, nz] ;
    
    % Found the source
    if( nx==source(1) && ny==source(2) && nz==source(3) )
		break;
    end
    
end

end

