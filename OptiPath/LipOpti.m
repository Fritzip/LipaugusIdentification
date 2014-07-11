% This script demonstrates the use of the optimization toolbox to solve a
% simple pathfinding problem. We are given an airplane with a fixed
% airspeed and a vector field of wind in which to navigate. The objective
% is to find the path which minimizes the time taken to travel from the
% starting point to the finish point.
%

%% First make some random vector field of wind, and set parameters

clear; clf;

temp = load('mastertemplate.mat');
temp = temp.m;
m = load('Srecaf.mat');
m = m.Srecaf;
edge = load('Edge.mat');
edge = edge.m;

temp = (temp./max(temp(:)))+1;
m = (m./max(m(:)))+1;
edge = (edge./max(edge(:)))+1;

AirSpeed = 500;

[sizeY, sizeX] = size(m);

numWayPoints = sizeX;
rng(50);

W_x = log(m)-mean(log(m(:)));
W_y = W_x.*3;

[Xgrid,Ygrid] = meshgrid(1:sizeX,1:sizeY);
hq = quiver(Xgrid,Ygrid,W_x,W_y,'k');
hold on;
xlabel('Units = [px]');
axis equal tight
plot([1 sizeX],[sizeY sizeY]/2,'k.','markersize',16)

%% Add some color to make it more visible
L = (sqrt((Xgrid-sizeX).^2 + (Ygrid-sizeY/2).^2));
Favorability = ((sizeX-Xgrid).*W_x +  (sizeY/2-Ygrid).*W_y)./L; 
Favorability(~isfinite(Favorability)) = 0; 

hold on;
h_im = imagesc(Favorability); % This will be the background for the vector field

set(h_im,'Xdata',[1 sizeX],'Ydata',[1 sizeY]);
uistack(h_im,'bottom');

% Change the colormap...
colormap(interp1([0,1,2],[1 0 0; 1 1 1; 0 1 0],0:0.01:2));
caxis(max(abs(Favorability(:)))*[-1 1]);

h_colorbar = colorbar;
title(h_colorbar,'Tailwind (km/h)')
xlabel(h_colorbar,'Headwind')
%% Generate a default set of way points
xWayPoints = linspace(1,sizeX,numWayPoints+2)';
yWayPoints = sizeY/2 * ones(numWayPoints+2,1);

h_wp = plot(xWayPoints,yWayPoints,'color','k','linestyle','none','marker','.','markersize',16);

%% Generate a continuous path from the waypoints
PathPoints = WayPoints_To_Path([xWayPoints,yWayPoints],'linear',sizeX,sizeY,2*sizeX);
h_path = plot(PathPoints(:,1),PathPoints(:,2),'k','linewidth',2);

%% Calculate the time taken:
StraightLineTime = getScoreFromPath(PathPoints, m, edge, temp)

% Display in Hours/Minutes
%fprintf('Straight Line Travel Time: %d hours, %.1f minutes\n',floor(StraightLineTime),rem(StraightLineTime,1)*60);

%% For example, a different set:
% Randomly choose waypoints
xWayPoints = linspace(0,sizeX,numWayPoints+2)' .* (1+0.05*randn(size(xWayPoints)));
yWayPoints = sizeY/2 + (0.1*sizeY*randn(numWayPoints+2,1));
xWayPoints([1 end]) = [1 sizeX];
yWayPoints([1 end]) = sizeY/2;

%Plot them
delete([h_wp h_path]);
h_wp = plot(xWayPoints,yWayPoints,'color','k','linestyle','none','marker','.','markersize',16);

% Generate a continuous path from the waypoints
PathPoints = WayPoints_To_Path([xWayPoints,yWayPoints],'cubic',sizeX,sizeY,2*sizeX);
h_path = plot(PathPoints(:,1),PathPoints(:,2),'k','linewidth',2);

LineTime = getScoreFromPath(PathPoints, m, edge, temp)
%fprintf('Travel Time: %d hours, %.1f minutes\n',floor(LineTime),rem(LineTime,1)*60);

%% Find an optimal path using FMINCON
% Define Objective Function
objectiveFun = @(P) getScoreFromPath(P,m,edge,temp,sizeX,sizeY,'cubic');

% Set optimization options
opts = optimset('fmincon');
opts.Display = 'iter';
opts.Algorithm = 'active-set';
opts.MaxFunEvals = 2000;

% Initial Conditions
xWayPoints = linspace(0,sizeX,numWayPoints+2)';
yWayPoints = sizeY/2 * ones(numWayPoints+2,1);
ic = [xWayPoints(2:end-1)'; yWayPoints(2:end-1)'];
ic = ic(:);

% Bounds
lb = zeros(size(ic(:)));
ub = reshape([sizeX*ones(1,numWayPoints); sizeY*ones(1,numWayPoints)],[],1);

%Do the optimizaiton
optimalWayPoints = fmincon(objectiveFun, ic(:), [],[],[],[],lb,ub,[],opts)

%% Plot the optimal solution:
delete([h_wp h_path]);
optimalWayPoints = [1 sizeY/2; reshape(optimalWayPoints,2,[])'; sizeX sizeY/2];

xWayPoints = optimalWayPoints(:,1);
yWayPoints = optimalWayPoints(:,2);
h_wp = plot(xWayPoints,yWayPoints,'color','k','linestyle','none','marker','.','markersize',16);

PathPoints = WayPoints_To_Path([xWayPoints,yWayPoints],'cubic',sizeX,sizeY,101);
h_path = plot(PathPoints(:,1),PathPoints(:,2),'k','linewidth',2);
LineTime = getTimeFromPath(PathPoints,W_x,W_y,AirSpeed);
%fprintf('Optimal Travel Time: %d hours, %.1f minutes\n',floor(LineTime),rem(LineTime,1)*60);
