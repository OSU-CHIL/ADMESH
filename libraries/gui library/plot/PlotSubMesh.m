function PlotSubMesh(gui,MESH)
% PlotMesh - Plots mesh
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% August 2013; Last revision: 12-August-2013

%--------------------------- BEGIN CODE ---------------------------------------


pH = gui.ViewAxes; % Plot Window handle

%------------------------------------------------------------------------------
% If there is something currently in the plot window, delete it
%------------------------------------------------------------------------------
% Delete current plot
h = findobj(pH,'tag', 'sub mesh'); if ~isempty(h); delete(h); end

colorbar('delete') % Delete current colorbar

%------------------------------------------------------------------------------
% Tell ADMESH where to plot
%------------------------------------------------------------------------------
axes(pH); hold on

%------------------------------------------------------------------------------
% Plot mesh
%------------------------------------------------------------------------------
patchinfo.Vertices          = [MESH.Points(:,1),MESH.Points(:,2)];
patchinfo.Faces             = MESH.ConnectivityList;
patchinfo.FaceColor         = 'none';
patchinfo.FaceVertexCData   = ones(size(MESH.Points,1),1)*[0 0 1];
patchinfo.EdgeColor         = 'b';
patchinfo.linewidth         = 1.0;
h                           = patch(patchinfo);

% % Set up right-click menu
% hcmenu = uicontextmenu;
% uimenu(hcmenu,'Label','Delete Element','Callback',@DeleteElement);
set(h,'tag','sub mesh')

%------------------------------------------------------------------------------
% Plot freeboundary
%------------------------------------------------------------------------------
[ff, fb] = freeBoundary(triangulation(MESH.ConnectivityList,MESH.Points(:,1),MESH.Points(:,2)));

x = reshape(fb(ff',1),2,numel(ff)/2);
y = reshape(fb(ff',2),2,numel(ff)/2);

plot(x,y,'k','linewidth',1.5)

%------------------------------------------------------------------------------
% Remove constraints
%------------------------------------------------------------------------------
delete(findobj('tag', 'Internal Constraint'))
delete(findobj('tag', 'External Constraint'))
delete(findobj('tag', 'Open Ocean'))
delete(findobj('tag', 'Line Constraint'))
delete(findobj('tag', 'Mesh','EdgeColor','g'))
delete(findobj('tag', 'Mesh','EdgeColor','k'))
h = findobj(pH,'tag', 'sub domain'); if ~isempty(h); delete(h); end

%------------------------------------------------------------------------------
% Plot boundaries
%------------------------------------------------------------------------------

for k = 1:numel(MESH.Constraints)
  
    % Determine line color
    if MESH.Constraints(k).num == -1 % Open Ocean
        
        nStr = MESH.Constraints(k).nodeStr;
        
        patchinfo.Faces             = [ nStr(1:end-1),nStr(2:end) ];
        patchinfo.EdgeColor         = 'b';
        patchinfo.linewidth         = 1.5;
        h                           = patch(patchinfo);
        set(h,'tag','sub mesh')
        
    elseif any(MESH.Constraints(k).num == [0 2 10 12 20 22 30]); % External Boundary
        
        nStr = MESH.Constraints(k).nodeStr;
        
        patchinfo.Faces             = [ nStr(1:end-1),nStr(2:end) ];
        patchinfo.EdgeColor         = 'k';
        patchinfo.linewidth         = 1.5;
        h                           = patch(patchinfo);
        set(h,'tag','sub mesh')
        
    elseif any(MESH.Constraints(k).num == [1 11 21]); % Internal Boundary
        
        nStr = MESH.Constraints(k).nodeStr;
        
        patchinfo.Faces             = [ nStr(1:end-1),nStr(2:end) ];
        patchinfo.EdgeColor         = 'g';
        patchinfo.linewidth         = 1.5;
        h                           = patch(patchinfo);
        set(h,'tag','sub mesh')
        
    elseif any(MESH.Constraints(k).num == [3 13 23]); % External Constraints
        
        nStr = MESH.Constraints(k).nodeStr;
        
        patchinfo.Faces             = [ nStr(1:end-1),nStr(2:end) ];
        patchinfo.EdgeColor         = [.7 .5 0];
        patchinfo.linewidth         = 1.5;
        h                           = patch(patchinfo);
        set(h,'tag','sub mesh')
        
    elseif any(MESH.Constraints(k).num == 18); % Channel Constraints 
        
        nStr = MESH.Constraints(k).nodeStr;
        
        patchinfo.Faces             = [ nStr(1:end-1),nStr(2:end) ];
        patchinfo.EdgeColor         = 'r';
        patchinfo.linewidth         = 1.5;
        h                           = patch(patchinfo);
        set(h,'tag','sub mesh')
        
    elseif any(MESH.Constraints(k).num == [4 5 24 25]); % Internal Constraints
        
        % Plot first side
        nStr = MESH.Constraints(k).nodeStr(:,1);
        
        patchinfo.Faces             = [ nStr(1:end-1),nStr(2:end) ];
        patchinfo.EdgeColor         = [.7 .5 0];
        patchinfo.linewidth         = 1.5;
        h                           = patch(patchinfo);
        set(h,'tag','sub mesh')
        
        % Plot second side
        nStr = MESH.Constraints(k).nodeStr(:,2);
        
        patchinfo.Faces             = [ nStr(1:end-1),nStr(2:end) ];
        patchinfo.EdgeColor         = [.7 .5 0];
        patchinfo.linewidth         = 1.5;
        h                           = patch(patchinfo);
        set(h,'tag','sub mesh')
        
    end
    
end

% %------------------------------------------------------------------------------
% % Show elevation if it is enabled
% %------------------------------------------------------------------------------
% PlotElevation; 

%------------------------------------------------------------------------------
% Display mesh info
%------------------------------------------------------------------------------
%DisplayMeshInfo(mesh.p,mesh.tri);

%------------------------------------------------------------------------------
% Set axis properties
%------------------------------------------------------------------------------
%[~,~] = AdjustAspectRatio(pH,per);

%guiH.per = per; % Save offset input in gui handle

%------------------------------------------------------------------------------
% Display Distance Scale
%------------------------------------------------------------------------------
% xmin = min(mesh.p(:,1));
% xmax = max(mesh.p(:,1));
% ymin = min(mesh.p(:,2));
% ymax = max(mesh.p(:,2));
% 
% DistanceScale(xmin,xmax,ymin,ymax)

drawnow

% Save gui data
% guidata(gui.Window,gui);



end