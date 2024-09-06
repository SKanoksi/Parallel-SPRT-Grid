function output = createAnalytic2D(TT,w,z0,x0,res,min,stride)
% /************************************
%      Parallel Shortest Path Solver
% createAnalytic = create analytical solution 2D
% - - - - - - - - - - - - - - - - - -
% TT = Analytical solution (Vector in X-Coord)
% (1)=Z-Coordinate , (2)=X-Coordinate
% *************************************/

% Allocate memory
output = zeros([res(2),res(1)]);
x = min(1) + stride(1)*((1:res(1))-0.5) ; % Vector ###

% Calculate
z = min(2) + 0.5*stride(2) ; % Scalar ###
for nz = 1:res(1)
    output(nz,:) = TT(w,z,z0,x,x0) ;
    z = z + stride(1) ;
end


end

