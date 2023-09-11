function AssignPeriodicBC(varargin)
% Assign periodic boundary conditions

%--------------------------------------------------------------------------
% Get GUI data
%--------------------------------------------------------------------------
guiFig = findobj('Tag','ADmesh Figure'); guiH = % guidata(guiFig);

%--------------------------------------------------------------------------
% Check if there is a mesh currently in the plot window
%--------------------------------------------------------------------------
meshPatch = findobj('Tag','Mesh');
if isempty(meshPatch); return; end

%--------------------------------------------------------------------------
% Declare persisten variables
%--------------------------------------------------------------------------
persistent pH p t N trep fb h

%--------------------------------------------------------------------------
% Evaluate input
%--------------------------------------------------------------------------
if nargin == 3
    
    status = 'selecting edges';

else
    
    status = get(findobj('tag','Periodic BC'), 'Checked');
    
end


switch status
    
    case 'off' % Set up plot for assigning boundary conditions
        
        % Get plot window handle
        pH = findobj('Tag','ADMESH Plot Window');
        
        % Initialize plot handle
        h = [];
        
        % Get mesh info
        p = get(meshPatch(end),'vertices'); % Vertices
        t = get(meshPatch(end),'faces');    % ConnecivityList
        N = guiH.MESH.ElementConnectivity;  % ElementConnectivityList
        trep = triangulation(t,p(:,[1 2])); % Triangulation representation
        fb = sort(freeBoundary(trep),2);    % Free boundary
        
        set([pH ; meshPatch],'Buttondownfcn',{@AssignPeriodicBC,meshPatch(end)})
        
        set(findobj('tag','Periodic BC'), 'Checked', 'on')
        
        uiStatusBar('Select an element edge...')
        
    case 'on' % Save element connectivity
        
        set([pH ; meshPatch],'Buttondownfcn',{})
        
        if ~isempty(h)
            delete(h)
        end
        
        set(findobj('tag','Periodic BC'), 'Checked', 'off')
        
        if ~isempty(N)
            guiH.MESH.ElementConnectivity = N;
        end
        
        % guidata(guiFig,guiH)
        
    case 'selecting edges'
        
        % Get current point on the screen
        currPt=get(pH,'CurrentPoint');
        X = currPt(1,1); Y = currPt(1,2);    
        
        % Find nearest vertex to the point selected
        id = knnsearch(p(:,[1 2]),[X Y]);
        
        % Find all elements attached to the vertex
        ti = cell2mat(vertexAttachments(trep,id));
        
        % Find element with closest edge
        
        D = zeros(numel(ti),1); % Initialize cell
        
        for k = 1:numel(ti) % Loop over elements
  
            i = [t(ti(k),:), t(ti(k),1)]; % Create triangle
            
            % Compute distance
            D(k) = Point2EdgeDistance(X,Y,p(i,1),p(i,2));
                        
        end
        
        % Closest triangle
        [~,ix] = min(D);
        
        % Indices to triangle (4x2) matrix
        i = [t(ti(ix),:)', [t(ti(ix),2:end), t(ti(ix),1)]'];
        
        % Initialize 
        D = zeros(3,1);
        
        % Find edge selected
        for k = 1:3
            D(k) = Point2EdgeDistance(X,Y,p(i(k,:),1),p(i(k,:),2));
        end
        [~,id] = min(D);
        
        % Check selected edge with free boundary list
        lia = ismember(sort(i(id,:),2),fb,'rows');
        
        % Check to see if edge selected is a boundary edge
        if ~lia; return; end
        
        % Find column to assign neighboring element to
        ni = setdiff(t(ti(ix),:),i(id,:)) == t(ti(ix),:);
        
        % Plot edge selected
        hold on
        h = [h; plot(p(i(id,:),1),p(i(id,:),2),'r','linewidth',1.5)];
        
        uiStatusBar('Select neighboring element...')
        
        % Wait for user to select neighboring element
        IN = 0; 
        while IN == 0
            
            k = 1;
            while k
                k = waitforbuttonpress;
            end
            
            % Get current point on the screen
            currPt=get(pH,'CurrentPoint');
            X = currPt(1,1); Y = currPt(1,2);
            
            % Find nearest vertex
            id = knnsearch(p(:,[1 2]),[X Y]);
            
            % Find all elements attached to vertex
            tii = cell2mat(vertexAttachments(trep,id));
            
            % Find element
            for k = 1:numel(tii)
                
                i = [t(tii(k),:), t(tii(k),1)];
                
                IN = PointInPolygon(X,Y,p(i,1),p(i,2));
                
                if IN; id = tii(k); break; end
                
            end
            
        end
        
        % highlight triangle
        [x,y] = deal(p(:,1),p(:,2));
        
        patchinfo.xdata          = x(t(id,:))';
        patchinfo.ydata          = y(t(id,:))';
        patchinfo.FaceColor      = 'g';
        patchinfo.linewidth      = 1.0;
        pobj                     = patch(patchinfo);
        
        drawnow; pause(.005);
        
        uiStatusBar('Appending to element connectivity list...')
        
        pause(.5)
        
        delete(pobj)
        
        % Store neighbor id
        N(ti(ix),:)
        N(ti(ix),ni) = id;
        N(ti(ix),:)
        
end

end