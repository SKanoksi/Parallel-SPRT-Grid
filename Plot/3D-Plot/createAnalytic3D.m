function output = createAnalytic3D(TT,w,x0,y0,z0,res,min,stride)
% /************************************
%      Parallel Shortest Path Solver
% createAnalytic = compute analytical solution 3D
% - - - - - - - - - - - - - - - - - -
% *************************************/

% Allocate memory
output = zeros(res);
x = min(1) + stride(1)*((1:res(1))-0.5) ; % Vector ###
y = min(2) + stride(2)*((1:res(2))-0.5) ; % Vector ###
z = min(3) + stride(3)*((1:res(3))-0.5) ; % Vector ###

% Calculate
for nz = 1:res(3)
    for ny = 1:res(2)
        output(:,ny,nz) = TT(w,x,x0,y(ny),y0,z(nz),z0) ;
    end
end

end

