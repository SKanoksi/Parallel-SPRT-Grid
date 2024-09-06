% /************************************
%      Parallel Shortest path solver
%           (Plot_VelModel.m)
% *************************************/
clear;

FILE_PATH = '/home/somrath/Work/PSPSX/Programs/FIM-WongKi/StructuredEikonal-master/src/'  ;
FILE_NAME = 'test.nrrd' ;

 PNG_NAME = 'FIM_Test' ; 
NAME_DATA = 'VelModel' ;

%% /******************  Read Binary file *********************/
FILE_PATH = strcat(FILE_PATH,FILE_NAME);
[data, header] = nrrdread(FILE_PATH);
res = size(data) ;



%% Display the data

% Contour plot
% figure(1); clf;
% cLevels = 0:4:120 ;    % Contour levels
% cs = contour(1:res(2),1:res(1),data,cLevels); 
% %clabel(cs);
% view(90,90);

% Color image
figure(3); clf;
%clims = [0 10];
imagesc(data(:,:,15))%,clims)
colormap(jet)
colorbar
axis image
%view(90,90);

% Mesh plot
% figure(3); clf;
% mesh(1:res(2),1:res(1),data);
% view(43,38);

% Output PNGs
% print('-f1',[PNG_NAME,'_',NAME_DATA,'_Contour.png'],'-dpng')
% print('-f2',[PNG_NAME,'_',NAME_DATA,'_Map.png'],'-dpng')
% print('-f3',[PNG_NAME,'_',NAME_DATA,'_Mesh.png'],'-dpng')

