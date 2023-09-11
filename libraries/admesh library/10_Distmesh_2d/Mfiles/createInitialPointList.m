function p = createInitialPointList(PTS,hmin)
% createInitialPointList - Creates initial point set for distmesh2d
% function
%
% Syntax:  p = createInitialPointList(PTS,hmin)
%
% Inputs:
%    PTS    - ADMESH Edge Structure
%    hmin   - Minimum element size
%
% Outputs:
%    p - initial point list for distmesh2d
%
% Other m-files required: distmesh2d
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% Last revision: 25-April-2014

%---------------------------------------------------------------------
% Begin Code
%---------------------------------------------------------------------

% Find bounding box
xmin = min(vertcat(PTS.Poly(:).x)); xmax = max(vertcat(PTS.Poly(:).x));
ymin = min(vertcat(PTS.Poly(:).y)); ymax = max(vertcat(PTS.Poly(:).y));

% Create initial distribution in bounding box (equilateral triangles)
[x,y] = meshgrid(xmin:hmin:xmax,ymin:hmin*sqrt(3)/2:ymax);

% Shift even rows
x(2:2:end,:)=x(2:2:end,:)+hmin/2;

% Store coordinate list in p
p=[x(:),y(:)];

end