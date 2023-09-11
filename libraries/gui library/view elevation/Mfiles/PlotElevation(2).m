function PlotElevation
% PlotElevation - Plots elevation
%
% Syntax:
%
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

%------------------------------------------------------------------------------
% Get GUI data
%------------------------------------------------------------------------------
guiFig = findobj('Tag','ADmesh Figure'); guiH = % guidata(guiFig);

%------------------------------------------------------------------------------
% Plot handle
%-----------------------------------------------------------------------------
pH = guiH.PlotWindow;

%------------------------------------------------------------------------------
% If we're plotting, are we displaying elevation for edge structure or
% mesh?
%-----------------------------------------------------------------------------
obj1 = findobj(get(pH,'children'),'tag','Edge Structure');
obj2 = findobj(get(pH,'children'),'tag','Mesh');

if isempty(obj1) && isempty(obj2) || strcmp(guiH.cmap,'None') || isempty(guiH.xyzFun)
    
    return
    
end

% Are we plotting elevation for edge structure or mesh?
if ~isempty(obj1)
    
    dispType = 'Edge Structure Elevation';
    
elseif ~isempty(obj2)
    
    dispType = 'Mesh Elevation';
    
end

%------------------------------------------------------------------------------
% Tell ADMESH where to plot
%------------------------------------------------------------------------------
axes(pH); hold on

%------------------------------------------------------------------------------
% Plot elevation
%------------------------------------------------------------------------------
switch dispType
    
    case 'Edge Structure Elevation'
        
        % Get elevation colormapping
        uiStatusBar('Computing color mapping from elevation..')
        [xv,yv,FaceVertexCData,AlphaData,cmap] = GetElevationRGB(guiH.xyzFun,guiH.PTS,guiH.cmap);
        
        uiStatusBar('Displaying elevation...')
        h = imagesc(xv,yv,FaceVertexCData,'AlphaData',AlphaData);
        
        uiStatusBar('Ready')
        set(h,'tag','Elevation')
        uistack(h,'bottom')
        
        % Distribute xyz data
        z = guiH.xyzFun.Values;
        
        colorBar(xv,yv,z,cmap)
        
        obj = findobj(get(pH,'children'),'tag','DistanceScale');
        uistack(obj,'top')
        
    case 'Mesh Elevation'
        
        % Get elevation colormapping
        [FaceVertexCData,cmap] = GetElevationRGB(...
            guiH.MESH.Points(:,1),guiH.MESH.Points(:,2),guiH.xyzFun,guiH.cmap);
        
        set(obj2(end),'FaceVertexCData',FaceVertexCData,'FaceColor','interp','edgecol','k')
        
        colorBar(guiH.MESH.Points(:,1),guiH.MESH.Points(:,2), guiH.xyzFun.Values,cmap)
        
end


end