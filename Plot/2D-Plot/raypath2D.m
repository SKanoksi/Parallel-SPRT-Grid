function [rayX, rayY] = raypath2D(data,radius, source, nx,ny,isSerial,isFM)
% /************************************
%      Parallel Shortest Path Solver
%     (raypath2D.m used in Plot_2D.m)
% *************************************/

if( not(isFM) ) 

% Shortest path raytracing
% Set initial and some numbers
res = size(data) ;
block = 2.*radius+1 ;
rayX = nx ;
rayY = ny ;
% At most, passing all the vertices
for i = 1:res(1)*res(2)
    % Recursive backtracking
    num = data(nx,ny) ;
    shiftX = mod(num,block(1)) ;
    shiftY = (num-shiftX)/block(1) ;

    if(~isSerial)
      % For Parallel version ###
      nx = nx + shiftX - radius(1) ; 
      ny = ny + shiftY - radius(2) ;
    else
      % For Serial version ###
      nx = nx - shiftX + radius(1) ; 
      ny = ny - shiftY + radius(2) ; 
    end
    rayX = [rayX, nx] ;
    rayY = [rayY, ny] ;
    % Found the source
    if( nx==source(1) && ny==source(2) )
		break;
    end
    
end

else %
    
% Find raypath (FM)    
% data = -gradient(T)/abs(gradient(T)) ###
% Set initial and some numbers
dh = 1 ; % Set step size ###
res = size(data);
source = single(source);
rayX = nx ;
rayY = ny ;
% At most, passing all the vertices
for i = 1:res(1)*res(2)
    
    % Recursive backtracking (FM)
    at = [round(nx),round(ny)];
    nx = nx + dh*data(at(2),at(1),1) ;
    ny = ny + dh*data(at(2),at(1),2) ;
    
    % Check BC
    nx = max(nx,1.0);
    nx = min(nx,single(res(2))) ; 
    ny = max(ny,1.0);
    ny = min(ny,single(res(1))) ; 
    
    rayX = [rayX, nx] ;
    rayY = [rayY, ny] ;
    
    % Found the source
    dis = sqrt( (nx-source(1)).^2 + (ny-source(2)).^2 );
    if( dis < 2*dh )
		break;
    end
    
    % To check
    % plot(rayX,rayY,'-k','LineWidth',1.0);
    % drawnow ;
    
end

end % ###

end

