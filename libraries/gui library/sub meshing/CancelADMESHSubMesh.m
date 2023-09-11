function CancelADMESHSubMesh(varargin)

% Deal input
if nargin == 2
    [fig,eventData] = deal(varargin{:});
else
    [fig,~,eventData.Key] = deal(varargin{:});
end

% If c is pressed, cancel submeshing
if strcmpi(eventData.Key, 'c')
    
    % Get gui data
    gui = % guidata(fig);
    
    % Clear key press function
    set(gui.Window, 'windowkeypressfcn'    , '')
    
    % Reset admesh button callback
    set(gui.RunAdmeshButton, 'Callback'  , @ADmeshRoutine)
    
    % Clear any plots
    h = findobj(gui.ViewAxes,'Tag','sub domain');
    if ~isempty(h)
        delete(h);
    end
    h = findobj(gui.ViewAxes,'Tag','sub mesh');
    if ~isempty(h)
        delete(h);
    end
    
    delete(findobj('tag', 'Internal Constraint'))
    delete(findobj('tag', 'External Constraint'))
    delete(findobj('tag', 'Open Ocean'))
    delete(findobj('tag', 'Line Constraint'))
    delete(findobj('tag', 'Mesh','EdgeColor','g'))
    delete(findobj('tag', 'Mesh','EdgeColor','k'))
    
    % If there is a mesh in the window, restore it
    meshPatch       = findobj('Tag','Mesh');
    if ~isempty(meshPatch)
        t = gui.MESH.ConnectivityList;
        FaceVertexAlphaData = 100*ones(size(t,1),1);
        set(meshPatch(end) ,'faces',t,'FaceVertexAlphaData',FaceVertexAlphaData)
    end
    
    % guidata(fig,gui)
    
end

end