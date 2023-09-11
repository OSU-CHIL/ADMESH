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
ttol = 1;
% "Internal pressure" of truss
Fscale  = 1.2;    
% Time step in Euler's method
deltat  = .2;     
% Tolerance in the geometry evaluations.
geps    = .001*hmin; 
% Frequency in which to check for nodes that are too close to one another
densityctrlfreq= 50;  
% Number of iterations
niter   = 5*1100;  
niter = 200;
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

%**************************************************************************
% Younghun temporarly added
%**************************************************************************

ConstraintsXY = [];
if ~isempty(PTS.Constraints)
for i = 1 : length(PTS.Constraints)
    ConstraintsXY = vertcat(ConstraintsXY, PTS.Constraints(i).xy, [nan nan]);
end
ConstraintsXY1 = vertcat(PTS.Constraints(:).xy);
end
% PTS.Constraints = []; 

%--------------------------------------------------------------------------
% Apply mesh constraints and concatenate with p if constraints exist.
%--------------------------------------------------------------------------
[p,nC,C,MESH] = GetMeshConstraints(p,hmin,PTS);
% nC = size(ConstraintsXY1,1); %....Younghun added
% p = [ConstraintsXY1; p]; %....Younghun added
N = size(p,1); % number of nodes

% %--------------------------------------------------------------------------
% % 2021-04-12 Younghun: Add channel branchs as fixed points (which is not
% % done by GetMeshConstraints function. 
% % Include this block only for using magnetic force approach
% %--------------------------------------------------------------------------
% temp = struct2table(PTS.Constraints);
% temp = temp.xy;
% pFix = [];
% for i = 1 : length(temp)
%     pFix = [pFix; temp{i}([1 end],:)];
% end
% nC = size(pFix,1);
% p = [pFix; p];
% N = size(p,1);

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
    
    drawnow; % for graphics

    % Retriangulation by the Delaunay algorithm 
    if( max(sqrt(sum((p-pold).^2,2))/hmin) > ttol )  % Any large movement?

        % Save current positions
        pold=p;
        
        %------------------------------------------------------------------
        % 2021-04-13 Younghun: reject nodes near constraint lines
        %------------------------------------------------------------------
        if k > (niter*0.9)
            id = find(-phi.f(p) < hmin*2);
            temp = struct2table(PTS.Constraints);
            temp = temp.xy;
            % temp = vertcat({[PTS.Poly.x,PTS.Poly.y]}, temp(:));
            
            [Vx,Vy,D] = Compute8SSED_v3(temp,p(id,1),p(id,2),hmin);
            id2 = find(abs(D) > 0 & abs(D) < hmin*0.4 & ~isnan(Vx) & ~isnan(Vy) & ~isinf(Vx) & ~isinf(Vy));
%             p(id(id2),:) = [];
            %--------------
            % 2021-04-25 Younghun: Apply projection instead of rejecting
            % the points
            %--------------
            id = id(id2);
            temp = knnsearch(p(1:nC,:),p(id,:),'k',2);
            for i = 1 : size(temp,1)
                j = find(ismember(C,temp(i,:),'rows'));
                if isempty(j)
                    j = find(ismember(C,temp(i,[2 1]),'rows'));
                end
                if isempty(j)
                    id(i) = 0; id2(i) = 0;
                    continue;
                end
%                 C(j,:) = [C(j,1) id(i)];
%                 C(end+1,:) = [id(i), C(j,2)];
            end
            id = nonzeros(id); id2 = nonzeros(id2);
            p(id,:) = p(id,:) + [Vx(id2), Vy(id2)];
            
            
            
            pold = p;
            N = size(p,1);
            in = ((nC+1):N)'; % Vector of non-pfix indices
        end
        
        % Perform Delaunay Triangulation
        if isempty(C)
            dt = delaunayTriangulation(p);
        else
            dt = delaunayTriangulation(p,C);
        end
        
        %------------------------------------------------------------------
        % 2021-03-18 Younghun: fix a issue (not sure when it happens) by a temporary solution
        %------------------------------------------------------------------
        if ~isequal(p,dt.Points)
            warning([...
                'The routine ''delaunayTriangulation'' results more/less points.',...
                'As a temporal solution, modify the variable dt so that match to original p.']);
            %                 id_p = (~ismember(dt.Points,p,'rows'));
            %                 dt.Points = p;
            
            %                 warning([...
            %                     'The routine ''delaunayTriangulation'' results more/less points.',...
            %                     'As a temporal solution, replace original p with Points from ''delaunayTriangulation''.']);
            p = dt.Points;
            pold = p;
            t = dt.ConnectivityList;
            N = size(p,1);
            in = ((nC+1):N)'; % Vector of non-pfix indices
            % id = find(ismember(dt.ConnectivityList(:,1),find(id_p)));
            % I = ind2sub(size(dt.ConnectivityList),id);
            % dt.ConnectivityList(I,:) = [];
            
        else

        % Compute centroids & Interpolate distances 
        ind = phi.f((p(dt(:,1),:)+p(dt(:,2),:)+p(dt(:,3),:))/3) < -geps;
        
        % Keep interior triangles
        t = sort(dt(ind,:),2);
        end
        

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
            if ~isempty(ConstraintsXY)
            plot(ConstraintsXY(:,1),ConstraintsXY(:,2),'r');
            drawnow; pause(1e-5);
            end
        end
        
    end
        
    % Move mesh points based on bar lengths L and forces F
    barvec  = p(bars(:,1),:)-p(bars(:,2),:);              % List of bar vectors
    L       = sqrt(sum(barvec.^2,2));                     % L = Bar lengths
    hbars   = MeshFun((p(bars(:,1),:)+p(bars(:,2),:))/2); % hbar = mesh size
    L0      = hbars*Fscale*sqrt(sum(L.^2)/sum(hbars.^2)); % L0 = Desired lengths
    
%     % Density control
%     if( (mod(k,densityctrlfreq) == 0) && (k < niter-5) )
%         
%         set(sb.ProgressBar,'Value',k)
%         sb.setText(['Generating mesh... ' num2str((k/niter)*100,'%.0f') '%'])
%         
%         if any(L0>2*L) % Remove points that are too close together
%             
%             p(setdiff(reshape(bars(L0>2*L,:),[],1),1:nC),:)=[];
%             
%             N=size(p,1); 
%             pold=inf;
%             
%             in = ((nC+1):N)';
%             
%             continue;
%             
%         elseif(k > niter/2)
%             
%             p = BoundaryDensityControl(p,t,C);
%             
%             p = ConstraintDensityControl(p,nC,C,MeshFun);
%             
%             N=size(p,1); 
%             pold=inf; 
%             in = ((nC+1):N)';
%             
%             continue;
%             
%         end
%         
%     end

    F            = max(L0-L,0);         % Bar forces (scalars)
    Fvec         = F./L*[1,1].*barvec;  % Bar forces (x,y components)
    Ftot         = full(sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],N,2));
    if k <= (niter*0.9)
    Ftot(1:nC,:) = 0;                   % Force = 0 at fixed points
    end
    p            = p+deltat*Ftot;       % Update node positions
    
    % Bring outside points back to the boundary
    p(in,:) = projectBackToBoundary(phi,p(in,:));
    
    %----------------------------------------------------------------------
    % 2021-04-08 Younghun: Apply magnetic-like force
    %----------------------------------------------------------------------
    if k > (niter*0.9)
    id = find(-phi.f(p) < hmin*2);
    temp = struct2table(PTS.Constraints);
    temp = temp.xy;
    temp = vertcat({[PTS.Poly.x,PTS.Poly.y]}, temp(:));
    
    [Vx,Vy,D] = Compute8SSED_v3(temp,p(id,1),p(id,2),hmin);
    id2 = find(abs(D) < hmin*0.4 & ~isnan(Vx) & ~isnan(Vy) & ~isinf(Vx) & ~isinf(Vy));
%     p(id(id2),:) = p(id(id2),:) + [Vx(id2),Vy(id2)];
    
    %--------------
    % 2021-04-25 Younghun: Apply projection instead of rejecting
    % the points
    %--------------
    id = id(id2);
    temp = knnsearch(p(1:nC,:),p(id,:),'k',2);
    for i = 1 : size(temp,1)
        j = find(ismember(C,temp(i,:),'rows'));
        if isempty(j)
            j = find(ismember(C,temp(i,[2 1]),'rows'));
        end
        if isempty(j)
            id(i) = 0; id2(i) = 0;
            continue;
        end
        %                 C(j,:) = [C(j,1) id(i)];
        %                 C(end+1,:) = [id(i), C(j,2)];
    end
    id = nonzeros(id); id2 = nonzeros(id2);
    p(id,:) = p(id,:) + [Vx(id2), Vy(id2)];
    
    end

   
    %------------------------------------------------------------------
    % 2021-04-15 Younghun: Perform triangulation once again because the updated
    % positions (p) possibly make elements crossing the constraints
    %------------------------------------------------------------------
    if k > (niter*0.9)
        % Perform Delaunay Triangulation
        if isempty(C)
            dt = delaunayTriangulation(p);
        else
            dt = delaunayTriangulation(p,C);
        end
        
        if ~isequal(p,dt.Points)
            warning([...
                'The routine ''delaunayTriangulation'' results more/less points.',...
                'As a temporal solution, modify the variable dt so that match to original p.']);
            %                 id_p = (~ismember(dt.Points,p,'rows'));
            %                 dt.Points = p;
            
            %                 warning([...
            %                     'The routine ''delaunayTriangulation'' results more/less points.',...
            %                     'As a temporal solution, replace original p with Points from ''delaunayTriangulation''.']);
            p = dt.Points;
            pold = p;
            t = dt.ConnectivityList;
            N = size(p,1);
            in = ((nC+1):N)'; % Vector of non-pfix indices
            
        end
            
        % Compute centroids & Interpolate distances
        ind = phi.f((p(dt(:,1),:)+p(dt(:,2),:)+p(dt(:,3),:))/3) < -geps;
        
        % Keep interior triangles
        t = sort(dt(ind,:),2);
    end
    
    % Check element quality. Keep track of best triangulation
    if k > (niter-50)
        [q, ~] = MeshQuality(p,t,0,'Triangle');
        if q > qold; P = p; T = t; qold = q; end
    end
    
    set(sb.ProgressBar,'Value',k)
    sb.setText(['Generating mesh... ' num2str((k/niter)*100,'%.0f') '%'])
    
end

sb.ProgressBar.setVisible(false)

% %--------------------------------------------------------------------------
% % 2021-04-13 Younghun: reject nodes near constraint lines
% %--------------------------------------------------------------------------
% id = find(-phi.f(P) < hmin*2);
% temp = struct2table(PTS.Constraints);
% temp = temp.xy;
% % temp = vertcat({[PTS.Poly.x,PTS.Poly.y]}, temp(:));
% 
% [Vx,Vy,D] = Compute8SSED_v3(temp,P(id,1),P(id,2),hmin);
% id2 = abs(D) > 0 & abs(D) < hmin*0.4 & ~isnan(Vx) & ~isnan(Vy) & ~isinf(Vx) & ~isinf(Vy);
% P(id(id2),:) = [];

% %--------------------------------------------------------------------------
% % 2021-03-18 Younghun: Perform triangulation once again because the updated
% % positions (p) possibly make elements crossing the constraints
% %--------------------------------------------------------------------------
% % [~,~,~,MESH] = GetMeshConstraints(p,hmin,PTS);
% if isempty(C)
%     dt = delaunayTriangulation(P);
% else
%     dt = delaunayTriangulation(P,C);
%     % 2021-03-18 Younghun added below to fix a issue (not sure when it happens) by a temporary solution
%     if ~isequal(P,dt.Points)
%         warning([...
%             'The routine ''delaunayTriangulation'' results more/less points.',...
%             'As a temporal solution, modify the variable dt so that match to original p.']);
% %         id_p = (~ismember(dt.Points,p,'rows'));
% % %         dt.Points = P;
%         P = dt.Points;
%         C = dt.Constraints;
%         % id = find(ismember(dt.ConnectivityList(:,1),find(id_p)));
%         % I = ind2sub(size(dt.ConnectivityList),id);
%         % dt.ConnectivityList(I,:) = [];
%         
%     end
% end
% % Compute centroids & Interpolate distances
% ind = phi.f((P(dt(:,1),:)+P(dt(:,2),:)+P(dt(:,3),:))/3) < -geps;
% 
% % Keep interior triangles
% T = sort(dt(ind,:),2);

% %--------------------------------------------------------------------------
% % 2021-04-12 Younghun: Apply magnetic force as a post-process
% %--------------------------------------------------------------------------
% id = find(-phi.f(p) < hmin);
% temp = struct2table(PTS.Constraints);
% temp = temp.xy;
% temp = vertcat({[PTS.Poly.x,PTS.Poly.y]}, temp(:));
% 
% [~,~,D] = Compute8SSED_v3(temp,p(id,1),p(id,2),hmin);
% id2 = D < 0 & abs(D) < hmin/2;
% 
% temp = cell2mat(temp);
% [id2B,dist2B] = knnsearch(temp,p(id,:));
% id2 = abs(D) < hmin/2 & dist2B < hmin;
% % id2 = dist2B < hmin*0.5;
% p1 = p;
% p1(id(id2),:) = temp(id2B(id2),:);


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
