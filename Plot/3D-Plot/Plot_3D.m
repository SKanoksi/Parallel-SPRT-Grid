% /************************************
%      Parallel Shortest Path Solver
%            (Plot_Map3D.m)
% *************************************/
clear;

%NAME_DATA = 'Slowness' ;
NAME_DATA = 'Traveltime' ;

% PSPS
PSPS_NAME = 'Z128_OpenGL_ExBF.nc' ;
PSPS_PATH = '/home/somrath/PSPSX/Code_Further/testrun_3d/' ;

% For raypath, see further below ###

%% /******************  Read NetCDF *********************/
% Numerical solution
[PSPS_data] = loadPSPS3D(PSPS_PATH,PSPS_NAME,NAME_DATA);
[PSPS_res, minCoord, stride, radius, source, name] = loadPARAM3D(PSPS_PATH,PSPS_NAME);
% % Reading from .nrrd (FIM-Original by Won-Ki Jeong)
% [PSPS_data, header] = nrrdread(strcat(PSPS_PATH,PSPS_NAME));
% PSPS_res = size(PSPS_data) ;
% % Need to swap x,y in the case of .nrrd ###
% minCoord = [0,0,1000] ;
% stride = 2*[1,1,1]*11.71875 ; % Manually set the stride ###
% %radius = [2,2,2] ; 
% source = [1,1,1] ; 


% Analytical solution

% Homogeneous ###
% Ana_func = @(w,x,x0,y,y0,z,z0) w*sqrt( (x-x0).^2+(y-y0).^2+(z-z0).^2 ); 

% Z-Gradient ### 
% Ana_func = @(w,x,x0,y,y0,z,z0) (1/w)*acosh( (z.^2+(x-x0).^2+(y-y0).^2+z0^2)./(2*z0*z) ); % Gradient #
%  
% w  = 1 ; % ###
% x0 = minCoord(1)+single(source(1)-0.5)*stride(1) ; % Source position
% y0 = minCoord(2)+single(source(2)-0.5)*stride(2) ;
% z0 = minCoord(3)+single(source(3)-0.5)*stride(3) ;
% Ana_data = createAnalytic3D(Ana_func,w,x0,y0,z0,PSPS_res,minCoord,stride);

% Reference solution --> Ana_data
% REF_NAME = 'SEG-Salt_Corner/OpenGL/SaltCorner_OpenGL_FIM.nc' ;
% [Ana_data] = loadPSPS3D(PSPS_PATH,REF_NAME,NAME_DATA);


%% /****************** Analyze data *********************/

% % Compute deviation
% Devia_data = abs( 100*(PSPS_data-Ana_data)./Ana_data) ;
% minError = min(min(min(Devia_data))) ;
% Devia_data(source(1),source(2),source(3)) = minError ; % Round-off error at source ###
% 
% % Calculate statistics
% maxError = max(max(max(Devia_data))) ;
% avgError = sum(sum(sum(Devia_data)))/PSPS_res(1)/PSPS_res(2)/PSPS_res(3) ;
% RMS = sqrt( sum(sum(sum( Devia_data.^2 )))/PSPS_res(1)/PSPS_res(2)/PSPS_res(3) )

% % Plot deviation from analytical solution (3D --> 2D)
 Plot_data = PSPS_data ;
% Level = 10 ; % Select 'Layer' ###
% Plot_data = squeeze(Plot_data(:,:,Level)); % Select 'Plane' ###

%% /****************** Display 'Plot_Data' *********************/
 Plot_res = size(Plot_data) ;
 x = minCoord(1)+(1:Plot_res(1))*stride(1);
 y = minCoord(2)+(1:Plot_res(2))*stride(2);
 z = minCoord(3)+(1:Plot_res(3))*stride(3);


 % "Color image"
% figure(1); clf;
% imagesc(Plot_data) % May use "mesh" and others ###
% colormap(jet)
% colorbar('eastoutside')
% xlabel('Number of cells'); ylabel('Number of cells');
% axis image


 % Plot isosurface
 % Prepare
 val = [0.7,1.2,1.7,2.2,2.7] ;
 [X,Y,Z] = meshgrid(x,y,z);
 Plot_data = flip(Plot_data,3); % Flip z
 % Start
 figure(1); clf;
 hold on
 for at = val 
     isosurface(X,Y,Z,Plot_data,at)
 end
 hold off 
 colorbar('eastoutside')
 % Reverse Z
 zt = get(gca, 'ZTick'); 
 set(gca, 'ZTick',zt, 'ZTickLabel',fliplr(zt))
 xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
 % Set lighting after shown.
 camlight ; 
% lighting phong;
 lightangle(-45,30);
 view(45,28)
 grid on
 axis image

 % Show Slice
% figure(1); clf;
% xslice = 1 ;
% yslice = 1 ;
% zslice = 1 ;
% slice(X,Y,Z,Plot_data,xslice,yslice,zslice);
% axis image

%% /****************** Display 'Raypath' *********************/
 % Read RayData
 [RayData] = loadPSPS3D(PSPS_PATH,PSPS_NAME,'Raypath');
 
 % Raypaths: set destinations
 % Only boundary
 layer = [1,128];
 nX = [     12*(1:10) , 12*(1:10)       , ...
         1*ones(1,10) , 128*ones(1,10)  ] ;
 nY = [  1*ones(1,10) , 128*ones(1,10)  , ...
            12*(1:10) , 12*(1:10)       ] ;
 % Only one depth layer surface
 nXX = nX ; nYY = nY ; nZ = layer(1)*ones(1,length(nX)) ; 
 for lz = layer(2:end) 
    nXX = [nXX, nX] ; 
    nYY = [nYY, nY] ;
    nZ = [nZ, lz*ones(1,length(nX))];
 end
 nX = nXX ; nY = nYY ;
 clear nXX nYY ;
 nX = round(nX); nY = round(nY) ; nZ = round(nZ) ;
 nRay = length(nZ) ;
 
 % Trace Back and Plot them.
 figure(2); clf;
 maxCoord = minCoord + stride.*PSPS_res ;
 hold on ;
 for i = 1:nRay
     [rayX,rayY,rayZ] = raypath3D(RayData,radius, source, nX(i),nY(i),nZ(i)); 
     rayX = minCoord(1) + stride(1)*( double(rayX) )-0.5*stride(1) ;
     rayY = minCoord(2) + stride(2)*( double(rayY) )-0.5*stride(2) ;
     rayZ = minCoord(3) + stride(3)*( double(rayZ) )-0.5*stride(3) ;
     plot3(maxCoord(1)-rayX,maxCoord(2)-rayY,maxCoord(3)-rayZ,'-k');
 end
 source =  maxCoord - (minCoord+stride.*double(source)-0.5*stride) ;
 plot3(source(1),source(2),source(3),'p',...
      'MarkerSize',16, ... 
      'MarkerEdgeColor',[1,1,1],...
      'MarkerFaceColor','r')
 hold off ;
 % % Reverse Z
 zt = get(gca, 'ZTick'); 
 set(gca, 'ZTick',zt, 'ZTickLabel',fliplr(zt))
 axis image ;
 xlabel('Distance(m)'); ylabel('Offset(m)'); zlabel('Depth(m)');
 view(51,19)
