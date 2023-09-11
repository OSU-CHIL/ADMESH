function [h_lfs] = MedialAxisFunction(X,Y,D,gradD,R,hmin,hmax,Settings,sb)
% medial_axis - Compute the medial axis of the polygon in PTS
%
% Syntax:  [h_lfs,MAP,AOF] = medial_axis(X,Y,R,gradD,PTS,plots,D,method)
%
% Inputs:
%    PTS - Data structure with fields x & y
%    PTS(1).x = x-coordinates of first polygon
%    PTS(1).y = y-coordinates of first polygon
%    X - (nxm) x-coordinates to rectangular grid
%    Y - (n,m) y-coordinates to rectangular grid
%    gradD - gradient of the distance function
%    D - Distance Function
%
% Outputs:
%    h_lfs - mesh size based on local feature size
%    MAP
%    AOF
%
% Other m-files required: medial_distance_FMM, Runge_Kutta_MAD
% Subfunctions: none
% MAT-files required: none
%
% Author: Colton Conroy
% The Ohio State University
% email address: conroy.51@osu.edu
% August 2013; Last revision: 08-August-2013

%------------- BEGIN CODE --------------

%------------------------------------------------------------------------------
% LFS Status
%------------------------------------------------------------------------------
Status = Settings.R.Status;

%------------------------------------------------------------------------------
% Compute Local Feature Size, If On.
%------------------------------------------------------------------------------
if strcmp(Status,'On')
    
    sb.ProgressBar.setVisible(true)
    sb.ProgressBar.setIndeterminate(true);
    sb.setText('Locating medial axis points...')
    
    % Locate Medial Axis Points w/in the Boundary (negative distance values)
    % Medial Axis point = point w/large singularity
    threshold = 0.15; % threshold for medial axis
    
    [LY,LX] = size(D);
    
    MA = zeros(LY,LX); % Initialize Medial Axis
    
    j = 2:LY-1; i = 2:LX-1;
    
    %compute dot product
    inner_productx = ...
        (X(j,i-1)-X(j,i))./(sqrt((X(j,i-1)-X(j,i)).^2+(Y(j+1,i)-Y(j,i)).^2)).*gradD.x(j,i-1) + ...
        (X(j,i+1)-X(j,i))./(sqrt((X(j,i+1)-X(j,i)).^2+(Y(j+1,i)-Y(j,i)).^2)).*gradD.x(j,i+1) + ...
        (X(j,i+1)-X(j,i))./(sqrt((X(j,i+1)-X(j,i)).^2)).*gradD.x(j,i+1) + ...
        (X(j,i+1)-X(j,i))./(sqrt((X(j,i+1)-X(j,i)).^2+(Y(j-1,i)-Y(j,i)).^2)).*gradD.x(j,i+1) + ...
        (X(j,i-1)-X(j,i))./(sqrt((X(j,i-1)-X(j,i)).^2+(Y(j-1,i)-Y(j,i)).^2)).*gradD.x(j,i-1) + ...
        (X(j,i-1)-X(j,i))./(sqrt((X(j,i-1)-X(j,i)).^2)).*gradD.x(j,i-1);
    
    inner_producty = ...
        (Y(j+1,i)-Y(j,i))./(sqrt((X(j,i-1)-X(j,i)).^2+(Y(j+1,i)-Y(j,i)).^2)).*gradD.y(j+1,i) + ...
        (Y(j+1,i)-Y(j,i))./(sqrt((Y(j+1,i)-Y(j,i)).^2)).*gradD.y(j+1,i) + ...
        (Y(j+1,i)-Y(j,i))./(sqrt((X(j,i+1)-X(j,i)).^2+(Y(j+1,i)-Y(j,i)).^2)).*gradD.y(j+1,i) + ...
        (Y(j-1,i)-Y(j,i))./(sqrt((X(j,i+1)-X(j,i)).^2+(Y(j-1,i)-Y(j,i)).^2)).*gradD.y(j-1,i) + ...
        (Y(j-1,i)-Y(j,i))./(sqrt((Y(j-1,i)-Y(j,i)).^2)).*gradD.y(j-1,i) + ...
        (Y(j-1,i)-Y(j,i))./(sqrt((X(j,i-1)-X(j,i)).^2+(Y(j-1,i)-Y(j,i)).^2)).*gradD.y(j-1,i);
    
    inner_product = inner_productx + inner_producty;
    
    sb.setText('Computing mesh size from Medial Axis...')
    
    AOF = inner_product./8;  %compute average outward flux
    
    MA(j,i) = AOF > threshold;
    
    MA(D > 0) = 0; % Any points outside the domain make 0's
    
    MA = bwmorph(MA,'skel',inf'); % Thin out Medial Axis line
    
    MA = bwmorph(MA,'clean',inf'); % Remove isolated pixels

    MAD = double(bwdist(MA))*abs(Y(1)-Y(2));
    
    LFS = abs(D) + abs(MAD);
    h_lfs = LFS./R;
    
    % Enforce boundary conditions
    h_lfs(h_lfs < hmin) = hmin;
    h_lfs(h_lfs > hmax) = hmax;

    sb.ProgressBar.setIndeterminate(false);
    sb.ProgressBar.setVisible(false)
    
else % LFS Off
    
    h_lfs = ones(size(X)).*hmax;
    
end

end
