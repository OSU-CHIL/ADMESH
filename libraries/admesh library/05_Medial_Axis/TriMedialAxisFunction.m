function [h_lfs,MAP,MAD]   = TriMedialAxisFunction(PTS,X,Y,D,gradD,R,hmin,hmax,Settings,sb,delta)
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
    
    %----------------------------------------------------------------
    % Create a constraints delaunay triangulation
    %----------------------------------------------------------------
    
    % Total number of segments
    nsegments = numel(PTS.Poly);
    
    % Prepare inputs for constrained delaunay triangulation
    C  = cell(nsegments,1);     % Cell for constraints
    
    % Vertically concatenate all polygons
    p = [vertcat(PTS.Poly(:).x) vertcat(PTS.Poly(:).y)];
    
    % Create unique coordinate set with back pointers.
    [p,~,ic] = unique(p,'rows','stable');
    
    % Create constraints matrix
    l = 0; % Initialize cumulative sum variable
    e = 0; % Initialize end index in constraints
    for k = 1:nsegments
        l       = l + length(PTS.Poly(k).x); 
        s       = e + 1; e = l; % Increase index
        C{k}    = [ic(s:e-1) ic(s+1:e)]; % Store constraint matrix
    end
    
    % Convert cell to matrix.
    C = sortrows(sort(cell2mat(C),2));
    
    % Check for repeating edges. If we have repeating edges then that means one
    % polygon lies on another. Both repeating constraints MUST be removed from
    % the constraints matrix because the "isinterior" function will not work
    % properly.
    [~,idx1]    = unique(C,'rows','first'); % Return first occurence
    [~,idx2]    = unique(C,'rows','last');  % Return last occurence
    C           = C(intersect(idx1,idx2),:);% Keep intersections.
    
    % Perform delaunay triangulation
    dt      = delaunayTriangulation(p,C);
    in      = isInterior(dt);
    
    % Construct a triangulation to represent the domain triangles.
    tr = triangulation(dt(in, :), dt.Points);
    
    % Construct a set of edges that join the circumcenters of neighboring
    % triangles; the additional logic constructs a unique set of such edges.
    numt    = size(tr,1);
    T       = (1:numt)';
    neigh   = tr.neighbors();
    cc      = tr.circumcenter();
    xcc     = cc(:,1);
    ycc     = cc(:,2);
    idx1    = T < neigh(:,1);
    idx2    =  T < neigh(:,2);
    idx3    = T < neigh(:,3);
    neigh   = [T(idx1) neigh(idx1,1); T(idx2) neigh(idx2,2); T(idx3) neigh(idx3,3)]';
    
    % Find all points less than hmax
    ix = find(D <= hmax);
    
    % Space each polygon coordinate by delta and compile into data structure
    p = PTS2PointList(PTS,delta);
    
    % Find the nearest neighbor in p for each point in (X,Y)
    sb.setText('Computing Distance Function..')
    
    MAP.x = xcc(neigh);
    MAP.y = ycc(neigh);
    
    % Create nearest neighbor search object
    ns = KDTreeSearcher([MAP.x(:) MAP.y(:)],'distance','euclidean');
    
    MAD = hmax*ones(size(X));
    
    % Perform nearest neighbor search
    [~,MAD(ix)] = knnsearch(ns,[X(ix) Y(ix)]);
        
    h_lfs = ones(size(X)).*hmax;
    
    sb.ProgressBar.setIndeterminate(false);
    sb.ProgressBar.setVisible(false)
    
else % LFS Off
    
    h_lfs = ones(size(X)).*hmax;
    
end

end

