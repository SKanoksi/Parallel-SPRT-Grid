% /************************************
%      Parallel Shortest Path Solver
%            (Plot_Map2D.m)
% *************************************/
clear;

%NAME_DATA = 'Slowness' ;
NAME_DATA = 'Traveltime' ;

% PSPS
isFM      = false ; % Compute rays from traveltime
isSerial  = false ;
PSPS_NAME = 'R2048_OpenGL_BF.nc' ;
PSPS_PATH = '/home/somrath/PSPSX/Code_Further/testrun_2d/' ;

% For raypath, see further below ###

%% /******************  Read NetCDF *********************/
% Numerical solution
[PSPS_data] = loadPSPS2D(PSPS_PATH,PSPS_NAME,NAME_DATA);
[PSPS_res, minCoord, stride, radius, source, name] = loadPARAM2D(PSPS_PATH,PSPS_NAME);
PSPS_data = PSPS_data' ; % Row (C++) <-> Column (MatLab)

% Analytical solution

% Homogeneous ###
% Ana_func = @(w,z,z0,x,x0) w*sqrt( (z-z0).^2+(x-x0).^2 ) ; 

% Z-Gradient ### 
% Ana_func = @(w,z,z0,x,x0) (1/w)*acosh( (z.^2 + (x-x0).^2 + z0^2)./(2*z0*z) ); % Gradient #

% R-Gradient ### 
%  Ana_func = @(w,z,z0,x,x0) (1/w)*sqrt( log( sqrt(x.^2+z.^2)/sqrt(x0.^2+z0.^2) ).^2 + ...
%             (atan2(-z,x)).^2 ); % R-Gradient #
% 
% w  = 1.0 ; % ###
% x0 = minCoord(1)+single(source(1)-0.5)*stride(1) ; % Source position
% z0 = minCoord(2)+single(source(2)-0.5)*stride(2) ;  
% Ana_data = createAnalytic2D(Ana_func,w,z0,x0,PSPS_res,minCoord,stride);


% Reference solution --> Ana_data
% REF_NAME = 'Marmousi_L/CPU/MarL_FMM.nc' ;
% [Ana_data] = loadPSPS2D(PSPS_PATH,REF_NAME,NAME_DATA);
% Ana_data = Ana_data' ; 


%% /****************** Analyze data *********************/
 
%  % Compute deviation
%  Devia_data = abs( 100*(PSPS_data-Ana_data)./Ana_data) ;
%  minError   = min(min(Devia_data)) ;
%  Devia_data(source(2),source(1)) = minError ; % Round-off error at source ###
%  Devia_data(isinf(Devia_data)) = minError ; % Round-off error around source ###
%  
%  % Calculate statistics
%  maxError = max(max(Devia_data)) ;
%  avgError = sum(sum(Devia_data))/PSPS_res(1)/PSPS_res(2) ;
%  RMS = sqrt( sum(sum( Devia_data.^2 ))/PSPS_res(1)/PSPS_res(2) )
 
 % Plot deviation from analytical solution
  Plot_data = PSPS_data ; % Row (C++) <-> Column (MatLab)
 
 %% /****************** Display 'Plot_Data' *********************/
  Plot_res = size(Plot_data) ;
  x = minCoord(1)+(1:Plot_res(1))*stride(1);
  y = minCoord(2)+(1:Plot_res(2))*stride(2);
 
  % "Color image"
  figure(1); clf;
  h = imagesc(Plot_data); % May use "mesh" and others ###
  colormap(jet) 
  hcb = colorbar('eastoutside');
  title(hcb,'second')
  %set(h, 'XData', y);
  %set(h, 'YData', x);
 
 
  % Set axis labels
  xticklabels = round( minCoord(1)+(0:576:Plot_res(1))*stride(1) ) ;
  xticks = linspace(1, Plot_res(1), numel(xticklabels));
  set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
  yticklabels = round( minCoord(2)+(0:188:Plot_res(2))*stride(2) ) ;
  yticks = linspace(1, Plot_res(2), numel(yticklabels));
  set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
 
  xlabel('Offset(m)'); ylabel('Depth(m)');
  xlabel('Number of cells'); ylabel('Number of cells');
 % title('Marmousi model: Traveltime');
  axis image
 
%% /****************** Display 'Raypath' *********************/
 % Raypaths: set destinations
 nX = [      80*(1:25) , 80*(1:25)       , ...
          1*ones(1,25) , 2048*ones(1,25)  ] ;
 nY = [   1*ones(1,25) , 2048*ones(1,25)  , ...
             80*(1:25) , 80*(1:25)       ] ;
 nX = round(nX); nY = round(nY) ;
 nRay = length(nX) ; 
 
 % Trace Back and Plot them. (SPR)
 if( not(isFM) )
     % Read RayData (SPR)
     [RayData] = loadPSPS2D(PSPS_PATH,PSPS_NAME,'Raypath');
 else
     % Find inverse ray direction from traveltime (FMM)
     [TT] = loadPSPS2D(PSPS_PATH,PSPS_NAME,'Traveltime');
     % Compute gradient (see 'edit gradient')
     [GX, GY] = gradient(TT'); % Row (C++) <-> Column (MatLab)
     % Normalized and negative gradient
     GL = sqrt(GX.^2 + GY.^2) ; 
     RayData(:,:,1) = -GX./GL ;
     RayData(:,:,2) = -GY./GL ;
     clear GX GY GL TT ;
 end
 figure(2); clf;
 hold on ;
 for i = 1:nRay
     [rayX,rayY] = raypath2D(RayData,radius, source, nX(i),nY(i),isSerial,isFM); 
     rayX = minCoord(1) + stride(1)*( double(rayX) )-0.5*stride(1) ;
     rayY = minCoord(2) + stride(2)*( double(rayY) )-0.5*stride(2) ;
     plot(rayX,rayY,'-k','LineWidth',1.0);
 end
 source = minCoord + stride.*double(source) - 0.5*stride ;
 plot(source(1),source(2),'p',...
     'MarkerSize',12, ... 
     'MarkerEdgeColor',[1,1,1],...
     'MarkerFaceColor','r')
 hold off ;
 axis image ;
 xlabel('Offset(m)'); ylabel('Depth(m)');
 title('Raypaths');
 view(0,-90) % Row (C++) <-> Column (MatLab)


