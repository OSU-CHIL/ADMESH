function p = ConstraintDensityControl(p,nfix,C,MeshFun)
% ConstraintDensityControl - Remove points too close to boundary
%
% Syntax:  p = ConstraintDensityControl(p,t,MeshFun)
%
% Inputs:
%    p          - Point List
%    t          - Connectivity
%    MeshFun    - Mesh Size Function
%
% Outputs:
%    p - Modified point list for distmesh2d
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

if ~isempty(C)

   % hmin = sqrt(3)/2*hmin;
    % Remove points in node list that will cause obscure triangulations
    % due to the constraints.
    n = (nfix+1:size(p,1))';
    I = cell(size(C,1),1);
    
    x1 = p(C(:,1),1); x2 = p(C(:,2),1);
    y1 = p(C(:,1),2); y2 = p(C(:,2),2);
    
    Poly.x = zeros(numel(x1),5);
    Poly.y = zeros(numel(x1),5);
    
    % Compute segment length
    dx = x1 - x2; dy = y1 - y2; eL = sqrt(dx.^2 + dy.^2);
    
    % Compute normal, nx, ny
    nx = dy./eL; ny = -dx./eL;
    
    % Compute midpoint of each segment
    xm = (x1 + x2)/2;
    ym = (y1 + y2)/2;      
  
    % Interpolate element size to mid points
    L = MeshFun(xm,ym);
    
    % Shorten L by sqrt(3)/4*L
    L = L.*sqrt(3)/8;

    Poly.x = [ x1 + nx.*L, x1 - nx.*L, x2 - nx.*L,  x2 + nx.*L, x1 + nx.*L ];
    Poly.y = [ y1 + ny.*L, y1 - ny.*L, y2 - ny.*L,  y2 + ny.*L, y1 + ny.*L ];
        
    % Remove nodes between constraints
    for i = 1:size(C,1)
        
        IN = PointInPolygon(p(n,1),p(n,2),Poly.x(i,:),Poly.y(i,:) ); % mex version
        
        I{i} = n(IN(:));
        
        drawnow; % for graphics
        
    end
    
    I = cell2mat(I);
    p(I,:) = [];

end

end