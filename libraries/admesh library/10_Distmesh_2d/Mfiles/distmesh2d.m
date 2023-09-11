function MESH = distmesh2d(PTS,phi,MeshFun,xyzFun,hmin,Settings,sb,pH)
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
%--------------------------- BEGIN CODE -----------------------------------

%--------------------------------------------------------------------------
% GUI check - Check if user wants to view mesh generation
%--------------------------------------------------------------------------
viewMesh = strcmpi(Settings.View.Status,'On'); axes(pH);

%--------------------------------------------------------------------------
% Initialize variables
%--------------------------------------------------------------------------
% Specifies how far the points can move (relatively) before a retriangulation 
ttol = .1;
% "Internal pressure" of truss
Fscale  = 1.2;    
% Time step in Euler's method
deltat  = .2;     
% Tolerance in the geometry evaluations.
geps    = .001*hmin; 
% Frequency in which to check for nodes that are too close to one another
densityctrlfreq= 25;  
% Number of iterations
niter   = 500;  
% Initialize current positions
pold=inf; 
% Intitalize current element quality
qold    = 0; 

sb.setText('Creating initial distribution of nodal points...')

%--------------------------------------------------------------------------
% Create initial distribution in bounding box (equilateral triangles)
%--------------------------------------------------------------------------
p = createInitialPointList(PTS,hmin);

%--------------------------------------------------------------------------
% Remove points outside the region, apply the rejection method
%--------------------------------------------------------------------------
p = rejectionMethod(p,phi.f,MeshFun,geps);

%--------------------------------------------------------------------------
% Apply mesh constraints and concatenate with p if constraints exist.
%--------------------------------------------------------------------------
[p,nC,C,MESH] = GetMeshConstraints(p,hmin,PTS);

N = size(p,1); % number of nodes

in = ((nC+1):N)'; % Vector of non-pfix indices

sb.ProgressBar.setVisible(true)
set(sb.ProgressBar, 'Minimum',1, 'Maximum',niter, 'Value',1)
sb.setText(['Generating mesh... ' num2str(0,'%.2f') '%'])

%--------------------------------------------------------------------------
% Remove previous plot before starting
%--------------------------------------------------------------------------

if viewMesh == 1; 
    
    h = findobj(pH,'tag', 'Edge Structure');
    h = [h;findobj(pH,'tag', 'Internal Constraint')];
    h = [h;findobj(pH,'tag', 'External Constraint')];
    h = [h;findobj(pH,'tag', 'Open Ocean')];
    h = [h;findobj(pH,'tag', 'Line Constraint')];
    
    if ~isempty(h); delete(h); end
    
    h = findobj('tag', 'Mesh');
    
    if ~isempty(h); delete(h); end 

end

%--------------------------------------------------------------------------
% Generate Mesh
%--------------------------------------------------------------------------
for k = 1:niter

    % Retriangulation by the Delaunay algorithm    
    if( max(sqrt(sum((p-pold).^2,2))/hmin) > ttol )  % Any large movement?

        % Save current positions
        pold=p;
        
        % Perform Delaunay Triangulation
        if isempty(C)
            t = delaunayTriangulation(p);
        else
            t = delaunayTriangulation(p,C);
        end

        % Compute centroids & Interpolate distances & Keep interior triangles
        t=sort(t(phi.f((p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3)<-geps,:),2);

        % Describe each bar by a unique pair of nodes
        bars = unique([t(:,[1,2]);t(:,[1,3]);t(:,[2,3])],'rows');
                        
        % Graphical output of the current mesh
        if viewMesh
            
            % Delete current plot
            h = findobj('tag', 'Mesh'); if ~isempty(h); delete(h); end
            
            % Plot mesh
            patch(...
                'vertices',p,...
                'faces',t,...
                'edgecol',[0 .4 .8],...
                'facecol','none',...
                'Tag','Mesh');

            % Plot Constraints
            if ~isempty(C)
     
                patch('vertices',p,'faces',C,'edgecol','r','facecol','none','Tag','Mesh');
                
            end
            
            drawnow; pause(1e-5);

        end
        
    end
        
    % Move mesh points based on bar lengths L and forces F
    barvec  = p(bars(:,1),:)-p(bars(:,2),:);              % List of bar vectors
    L       = sqrt(sum(barvec.^2,2));                     % L = Bar lengths
    hbars   = MeshFun((p(bars(:,1),:)+p(bars(:,2),:))/2); % hbar = mesh size
    L0      = hbars*Fscale*sqrt(sum(L.^2)/sum(hbars.^2)); % L0 = Desired lengths
    
    % Density control
    if( (mod(k,densityctrlfreq) == 0) && (k < niter-5) )
        
        set(sb.ProgressBar,'Value',k)
        sb.setText(['Generating mesh... ' num2str((k/niter)*100,'%.0f') '%'])
        
        if any(L0>2*L) % Remove points that are too close together
            
            p(setdiff(reshape(bars(L0>2*L,:),[],1),1:nC),:)=[];
            
            N=size(p,1); 
            pold=inf;
            
            in = ((nC+1):N)';
            
            continue;
            
        elseif(k > niter/2)
            
            p = BoundaryDensityControl(p,t,C);
            
            p = ConstraintDensityControl(p,nC,C,MeshFun);
            
            N=size(p,1); 
            pold=inf; 
            in = ((nC+1):N)';
            
            continue;
            
        end
        
    end

    F            = max(L0-L,0);         % Bar forces (scalars)
    Fvec         = F./L*[1,1].*barvec;  % Bar forces (x,y components)
    Ftot         = full(sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],N,2));
    Ftot(1:nC,:) = 0;                   % Force = 0 at fixed points
    p            = p+deltat*Ftot;       % Update node positions
    
    % Bring outside points back to the boundary
    p(in,:) = projectBackToBoundary(phi,p(in,:));
    
    % Check element quality. Keep track of best triangulation
    if k > (niter-50)
        [q, ~] = MeshQuality(p,t,0,'Triangle');
        if q > qold; P = p; T = t; qold = q; end
    end
    
end

sb.ProgressBar.setVisible(false)

%--------------------------------------------------------------------------
% Clean up 
%--------------------------------------------------------------------------
sb.setText('Cleaning up final mesh...')

T = BoundaryCleanUp(P,T,C); % Remove bad boundary elements

[p,t]=fixmesh(P,T);         % Fix mesh

%p = P; t = T;

%--------------------------------------------------------------------------
% store final mesh results
%--------------------------------------------------------------------------
MESH = createMeshStruct(t,p,MESH,PTS,xyzFun); 

end
