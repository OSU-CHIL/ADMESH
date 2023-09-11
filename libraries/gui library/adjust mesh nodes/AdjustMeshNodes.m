function AdjustMeshNodes(varargin)

%------------------------------------------------------------------------------
% Get GUI data
%------------------------------------------------------------------------------
[app,~,status] = deal(varargin{:});

%------------------------------------------------------------------------------
% Check if there is a mesh currently in the plot window
%------------------------------------------------------------------------------
meshPatch = findobj(app.UIAxes,'Tag','Mesh');

if isempty(meshPatch)
    
    editNodeH = findobj('Tag','Move Node Tag'); % Get tool handle
    set(editNodeH,'State','off')                % Turn off
    drawnow
    return
    
end

%------------------------------------------------------------------------------
% Or we moving nodes?
%------------------------------------------------------------------------------
switch lower(status)
    
    case 'on'
        
        % Turn off all other tools
        zoomH = findobj('Tag','Zoom Tag');
        panH  = findobj('Tag','Pan Tag');
        set([zoomH panH],'State','off')
        
        drawnow
        
        % Create button down function
%         set([app.UIAxes ; meshPatch],'Buttondownfcn',{@moveNodes,'select node',[]})
        set([app.UIFigure; meshPatch],'Buttondownfcn',{@moveNodes,'select node',[]})

        app.ProgressBarButton.Text = 'Select a node...'; drawnow
        
        drawnow
        
    case 'off'
        
        % Turn off
        editNodeH = findobj('Tag','Move Node Tag');
        set(editNodeH,'State','off')
        
        drawnow
        
        % Remove button down function
%         set([app.UIAxes ; meshPatch],'Buttondownfcn','')
        
        drawnow
        
end

end