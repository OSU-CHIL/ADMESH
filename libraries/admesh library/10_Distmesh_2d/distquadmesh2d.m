function MESH = distquadmesh2d(PTS,DistFun,DistGxFun,DistGyFun,MeshFun,xyzFun,hmin,Settings)
% distmesh2d - Generates a mesh based on mesh size h
%
% Syntax:  [p,t] = distmesh2d(DistanceFun,MeshSizeFun,hmin,guiFig)
%
% Inputs:
%    DistanceFun - gridded interpolant of the Distance Function
%    MeshSizeFun - gridded interpolant of the Mesh Size Function
%    hmin        - minimum element size
%    guiFig - handle that identifies the figure
%
% Outputs:
%    p - points of delaunay triangulation
%    t - connectivity list
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

% DISTMESH2D 2-D Mesh Generator using Distance Functions.
% distmesh2d.m v1.1
% Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.
%
% Author of Adjustments: Dustin West, Colton Conroy, Ethan Kubatko
% The Ohio State University
% email address: dww.425@gmail.com
% August 2013; Last revision: 08-August-2013

%------------- BEGIN CODE --------------

%--------------------------------------------------------------------------
% Get plot window handle
%--------------------------------------------------------------------------
pH = findobj('Tag','ADMESH Plot Window');

%--------------------------------------------------------------------------
% GUI check - Check if user wants to view mesh generation
%--------------------------------------------------------------------------
viewMesh = strcmpi(Settings.View.Status,'On'); axes(pH);

%--------------------------------------------------------------------------
% Initialize variables
%--------------------------------------------------------------------------
% Specifies how far the points can move (relatively) before a retriangulation 
if viewMesh
    ttol=.5; % Better for visual
else
    ttol=4; % Speeds program up a lot
end
% "Internal pressure" of truss
Fscale=1.2;    
% Time step in Euler's method
deltat=.2;     
% Tolerance in the geometry evaluations.
geps = .001*hmin; 
% Frequency in which to check for nodes that are too close to one another
densityctrlfreq= 25;  
% Number of iterations
niter = 1000;  
% Initialize current positions
pold=inf; 
% Intitalize current element quality
qold = 0; 

uiStatusBar('Creating initial distribution of nodal points...')

%--------------------------------------------------------------------------
% Double the mesh size
%--------------------------------------------------------------------------
MeshFun.Values = MeshFun.Values.*100000;

%--------------------------------------------------------------------------
% Create initial distribution in bounding box (equilateral triangles)
%--------------------------------------------------------------------------
p = createInitialPointList(PTS,hmin);

%--------------------------------------------------------------------------
% Remove points outside the region, apply the rejection method
%--------------------------------------------------------------------------
p = rejectionMethod(p,DistFun,MeshFun,geps);

%--------------------------------------------------------------------------
% Apply mesh constraints and concatenate with p if constraints exist.
%--------------------------------------------------------------------------
[p,nC,C,MESH] = GetMeshConstraints(p,hmin,PTS);

N = size(p,1); % number of nodes

in = ((nC+1):N)'; % Vector of non-pfix indices

uiStatusBar('Generating mesh... ')

uiStatusBar(1/niter)

%--------------------------------------------------------------------------
% Remove previous plot before starting
%--------------------------------------------------------------------------
if viewMesh == 1; 
    
    h = findobj('tag', 'Edge Structure');
    h = [h;findobj('tag', 'Internal Constraint')];
    h = [h;findobj('tag', 'External Constraint')];
    h = [h;findobj('tag', 'Open Ocean')];
    h = [h;findobj('tag', 'Line Constraint')];
    
    if ~isempty(h); delete(h); end
    
    h = findobj('tag', 'Mesh');
    
    if ~isempty(h); delete(h); end 

end


%--------------------------------------------------------------------------
% Generate Triangular Mesh
%--------------------------------------------------------------------------
for k = 1:round(niter/2)
    
    % Retriangulation by the Delaunay algorithm
    if( max(sqrt(sum((p-pold).^2,2))/hmin) > ttol )  % Any large movement?
        
        % Perform Delaunay Triangulation
        if isempty(C)
            t = delaunayTriangulation(p);
        else
            t = delaunayTriangulation(p,C);
        end
        
        % Compute centroids & Interpolate distances & Keep interior triangles
        t=t(DistFun((p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3)<-geps,:);
                
        % Save current positions
        pold=p;
        
        % Describe each bar by a unique pair of nodes
        bars = unique(sort([t(:,[1,2]);t(:,[2,3]);t(:,[3,1])],2),'rows');
        
        % Graphical output of the current mesh
        if viewMesh
            
            % Delete current plot
            h = findobj('tag', 'Mesh'); if ~isempty(h); delete(h); end
            
            % Plot mesh
            patch('vertices',p,'faces',t,'edgecol','b','facecol','none','Tag','Mesh');

            % Plot Constraints
            if ~isempty(C)
     
                patch('vertices',p,'faces',C,'edgecol','r','facecol','none','Tag','Mesh');
                
            end
            
            drawnow
                        
        end
        
    end

    % Move mesh points based on bar lengths L and forces F
    barvec  = p(bars(:,1),:)-p(bars(:,2),:);              % List of bar vectors
    L       = sqrt(sum(barvec.^2,2));                     % L = Bar lengths
    hbars   = MeshFun((p(bars(:,1),:)+p(bars(:,2),:))/2); % hbar = mesh size
    L0      = hbars*Fscale*sqrt(sum(L.^2)/sum(hbars.^2)); % L0 = Desired lengths
    
    % Density control 
    if( (mod(k,densityctrlfreq) == 0) && (k < niter-5) )
        
        uiStatusBar(k/niter)
        
        if any(L0>2*L) % Remove points that are too close together
            
            p(setdiff(reshape(bars(L0>2*L,:),[],1),1:nC),:)=[];
            
            N=size(p,1); pold=inf;
            
            in = ((nC+1):N)';
            
            continue;
            
        elseif(k > niter/2)
            
            p = BoundaryDensityControl(p,t,C);
            
            p = ConstraintDensityControl(p,nC,C,MeshFun);

            N=size(p,1); pold=inf; in = ((nC+1):N)';
            
            continue;
            
        end
        
    end

    F            = max(L0-L,0);         % Bar forces (scalars)
    Fvec         = F./L*[1,1].*barvec;  % Bar forces (x,y components)
    Ftot         = full(sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],N,2));
    Ftot(1:nC,:) = 0;                   % Force = 0 at fixed points
    p            = p+deltat*Ftot;       % Update node positions
    
    % Bring outside points back to the boundary
    p(in,:) = projectBackToBoundary(DistFun,DistGxFun,DistGyFun,p(in,:));
        
    % Check element quality. Keep track of best triangulation
    if k > (round(niter/2)-50)
        [q, ~] = MeshQuality(p,t,0,'Triangle');
        if q > qold; P = p; T = t; qold = q; end
    end
        
end

uiStatusBar(niter/niter)

%--------------------------------------------------------------------------
% Transition to quad mesh
%--------------------------------------------------------------------------
uiStatusBar('Converting mesh to quadrilateral mesh... ')

T = BoundaryCleanUp(P,T,C); % Remove bad boundary elements
[p,t]=fixmesh(P,T);         % Fix mesh

[p,t] = tri2quad(p(:,[1,2]),t);

N = size(p,1); % number of nodes

in = ((nC+1):N)'; % Vector of non-pfix indices

% Describe each bar by a unique pair of nodes
bars = unique(sort([t(:,[1,2]);t(:,[2,3]);t(:,[3,4]);t(:,[4,1])],2),'rows');

%--------------------------------------------------------------------------
% Half the mesh size
%--------------------------------------------------------------------------
%MeshFun.Values = MeshFun.Values./100;

for j = k:niter
       
    if mod(j,densityctrlfreq) == 0
        uiStatusBar(j/niter)
    end
    
    % Graphical output of the current mesh
    if viewMesh
        
        % Delete current plot
        h = findobj('tag', 'Mesh'); if ~isempty(h); delete(h); end
        
        % Plot mesh
        patch('vertices',p,'faces',t,'edgecol','b','facecol','none','Tag','Mesh');
        
        % Plot Constraints
        if ~isempty(C)
            
            patch('vertices',p,'faces',C,'edgecol','r','facecol','none','Tag','Mesh');
            
        end
        
        drawnow
        
    end
        
    
    % Move mesh points based on bar lengths L and forces F
    barvec  = p(bars(:,1),:)-p(bars(:,2),:);              % List of bar vectors
    L       = sqrt(sum(barvec.^2,2));                     % L = Bar lengths
    hbars   = MeshFun((p(bars(:,1),:)+p(bars(:,2),:))/2); % hbar = mesh size
    L0      = hbars*Fscale*sqrt(sum(L.^2)/sum(hbars.^2)); % L0 = Desired lengths
        
    F            = max(L0-L,0);         % Bar forces (scalars)
    Fvec         = F./L*[1,1].*barvec;  % Bar forces (x,y components)
    Ftot         = full(sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],N,2));
    Ftot(1:nC,:) = 0;                   % Force = 0 at fixed points
    p            = p+deltat*Ftot;       % Update node positions
    
    % Bring outside points back to the boundary
    p(in,:) = projectBackToBoundary(DistFun,DistGxFun,DistGyFun,p(in,:));
    
end


uiStatusBar(niter/niter)


return

%--------------------------------------------------------------------------
% Clean up 
%--------------------------------------------------------------------------
uiStatusBar('Cleaning up final mesh... ')

T = BoundaryCleanUp(P,T,C); % Remove bad boundary elements

[p,t]=fixmesh(P,T);         % Fix mesh

%--------------------------------------------------------------------------
% store final mesh results
%--------------------------------------------------------------------------
MESH = createMeshStruct(t,p,MESH,PTS,xyzFun); 

end