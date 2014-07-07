function sc = getScoreFromPath(PathPoints,m,edge,temp,sizeX,sizeY,METHOD)
% This is the main function that actually calculates the line integral
% along the path. This is the objective function for the optimizer.
%

% If we are called from the optimization routine (caller = 'optimizer')
% then we need to interpolate the fine path from the input control points.
if isvector(PathPoints)
    PathPoints = [1 sizeY/2; reshape(PathPoints,2,[])'; sizeX sizeY/2];
    PathPoints = WayPoints_To_Path(PathPoints,METHOD,sizeX,sizeY,101);
end

dP = diff(PathPoints);
slope = dP(:,2)./dP(:,1);
sc1 = sum(abs(slope))+1;

sc2 = [];
sc3 = [];
sc4 = [];
sc5 = [];
for i = 1:size(PathPoints,1)
    x = round(PathPoints(i,1));
    y = round(PathPoints(i,2));
    
    %sc2 = [sc2;  sum(m(y,x-5:x-1))/sum(m(y,x+1:x+5))];
    sc3 = [sc3 m(y,x)];
    sc4 = [sc4 edge(y,x)];
    sc5 = [sc5 temp(y,x)];
end
% disp(sc3')
% disp(sc4')
% disp(sc5')
sc3 = 1+(sc3/max(sc3));
sc4 = 1+(sc4/max(sc4));
sc5 = 1+(sc5/max(sc5));
sc = (sc1)/(mean(sc3)*mean(sc4)*mean(sc5));

% 
% % Interpolate the wind vector field at all the points in PathPoints.
% V_wind = interp2(W_x,PathPoints(1:end,1),PathPoints(1:end,2),'linear');
% 
% % Dot product the wind (Vwind_path) with the direction vector (dP) to get
% % the tailwind/headwind contribution
% V_add = (sum(V_wind.*dP,2))./sqrt(sum(dP.^2,2));
% dx = sqrt(sum(dP.^2,2))*100; %dx is the length of each subinterval in PathPoints
% dt = dx./(AirSpeed+V_add);  %dT = dP/dV
% TravelTime = sum(dt);