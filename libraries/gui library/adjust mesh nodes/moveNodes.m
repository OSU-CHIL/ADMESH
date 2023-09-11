function moveNodes(varargin)
% moveNodes - Interactivly select and move nodal points
%
% Syntax:  moveNodes
%
% Inputs:
%
% Outputs:
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% August 2013; Last revision: 11-February-2013

%---------------------------- BEGIN CODE ----------------------------

%------------------------------------------------------------------------------
% Initialize persistent variables
%------------------------------------------------------------------------------
persistent meshPatch pH marker st cPatch idc

%-------------------------------------------------------------------------
% Check user input
%-------------------------------------------------------------------------
if nargin == 1
    action = varargin{1};
elseif nargin == 4
    [app,~,action,idx] = deal(varargin{:});
elseif nargin == 5
    [app,~,action,idx,~] = deal(varargin{:});
end

%-------------------------------------------------------------------------
% What action should we take?
%-------------------------------------------------------------------------
switch(action)
    
    %---------------------------------------------------------------------
    % User attempts to select node(s)
    %---------------------------------------------------------------------
    case 'select node'
                        
        % Get guidata
%         gui         = guidata(fig);
        
        % gui.sb.setText('Please wait...'); drawnow
        
%         set(gui.Window,'renderer','opengl'); drawnow
            
        % Get plot children handles
%         st          = gui.st;
        pH          = app.UIAxes;
        meshPatch   = findobj(pH,'Tag','Mesh');
        cPatch      = findobj(pH,'Tag','Mesh Constraint');
                             
        % Get current point on the screen
        currPt = get(pH,'CurrentPoint');
        x = currPt(1,1); y = currPt(1,2);

        % Compute distance & choose points within min(d)  
        [~,idx] = min(abs(x - meshPatch.Vertices(:,1)) + abs(y - meshPatch.Vertices(:,2)));
                
        % Check for user mess up
        if isempty(idx);
            
            moveNodes('quit')
            
        else
            
            % gui.sb.setText('Grabing node, don''t move the mouse yet...'); drawnow
            
            % Determine if we need to move constraints to
            if ~isempty(cPatch)
                                                
%                 gui.sb.ProgressBar.setVisible(true)
%                 set(gui.sb.ProgressBar, 'Minimum',1, 'Maximum',length(cPatch), 'Value',1)
                
                idc = zeros( length(cPatch),2);
                j   = 1;
                
                for k = 1: length(cPatch)
                    
%                     set(gui.sb.ProgressBar,'Value',k)
                    
                    ix = find(abs((meshPatch.Vertices(idx,1)-cPatch(k).XData)+abs(meshPatch.Vertices(idx,2)-cPatch(k).YData)) <= eps);
                    
                    if isempty(ix)
                        
                        continue;
                        
                    elseif length(ix) == 1
                        
                        idc(j,1) = ix;
                        idc(j,2) = k;
                        j = j + 1;
                        
                    elseif length(j) == 2
                        
                        idc(j,1) = ix(1);
                        idc(j,2) = k;
                        j = j + 1;
                        
                        idc(j,1) = ix(2);
                        idc(j,2) = k;
                        j = j + 1;
                        
                    end
                    
                end
                
%                 gui.sb.ProgressBar.setVisible(false); drawnow
                
                idc(idc(:,1) == 0,:) = [];
 
            end
            
            marker = plot(pH,...
                meshPatch.Vertices(idx,1),meshPatch.Vertices(idx,2),...
                'bo',...
                'MarkerFaceColor','b',...
                'MarkerSize',5);
            drawnow
            
            % If the user holds down the left mouse button, move some nodes
            % around
            idxSub = nan;
            
%             set(app.UIFigure,'WindowButtonMotionFcn',{@moveNodes, 'move nodes', idx,idxSub})
            moveNodes(app,[],'move nodes',idx,idxSub);
%             drawnow
            
            % When the user lets up on the mouse, quit.
%             set(app.UIFigure,'WindowButtonUpFcn',{@moveNodes, 'quit', idx})
%             moveNodes(app,[],'quit',idx);
            drawnow
            
            % Inform user
            app.ProgressBarButton.Text = 'Moving node...';
            drawnow;
            
        end
      
    %---------------------------------------------------------------------
    % User moves nodes around
    %---------------------------------------------------------------------
    case 'move nodes'
        
        % Get current point on the screen
%         currPt=get(pH,'CurrentPoint');
        currPt = app.UIAxes.CurrentPoint;
        X = currPt(1,1); Y = currPt(1,2);
        
        % Update new point on screen
        marker.XData = X; marker.YData = Y;
        
        % Update new point in node list
        meshPatch.Vertices(idx,[1 2]) = [X Y];
        
        % Update constraints
        for k = 1:size(idc,1)
            cPatch(idc(k,2)).XData(idc(k,1)) = X;
            cPatch(idc(k,2)).YData(idc(k,1)) = Y;
        end
        
        drawnow

        % Update plot and mesh info box
        %DisplayMeshInfo(fig,meshPatch.Vertices,meshPatch.Faces,'update')
        
        % Update coordinate display
%         st.setText(['( ' num2str(X,'%.3f') ' , ' num2str(Y,'%.3f') ' )']);

        drawnow
        
        %---------------------------------------------------------------------
        % User lets up on the mouse
        %---------------------------------------------------------------------
    case 'quit'
                
        % Get guidata
%         gui = guidata(app);
        
        % gui.sb.setText('Updating, please wait...');
                
        % Update elevation if it's available
%         currPt=get(pH,'CurrentPoint');
        currPt = app.UIAxes.CurrentPoint;
        X = currPt(1,1); Y = currPt(1,2);
        
        % Update new point on screen
        if isempty(app.xyzFun)
            meshPatch.Vertices(idx,[1 2]) = [X Y];
        else
            meshPatch.Vertices(idx,:) = [X Y app.xyzFun(X,Y)];
        end

        % Update constraints
        for k = 1:size(idc,1)
            cPatch(idc(k,2)).XData(idc(k,1)) = X;
            cPatch(idc(k,2)).YData(idc(k,1)) = Y;
        end

        % Remove temporary patch
        delete(marker); drawnow

        % Update plot and mesh info box
        DisplayMeshInfo(app,meshPatch.Vertices,meshPatch.Faces,'update')
        
%         gui = guidata(app);
        
        % Update mesh
        app.MESH.Points = meshPatch.Vertices;
        
        set(app,'WindowButtonMotionFcn',@CoordDisplay)
        
        set(app,'WindowButtonUpFcn','')
        
        set(app.UIFigure,'renderer','zbuffer')
        
        % guidata(fig, gui);
        
        % Clear persistend varables
        %Vars=whos;
        %PersistentVars=Vars([Vars.persistent]);
        %PersistentVarNames={PersistentVars.name};
        %clear(PersistentVarNames{:});
        
        % gui.sb.setText('Ready');
        
end

end
