function PanCallback(varargin)

% Load GUI data
gui = % guidata(varargin{1});

% Get axes handle
pH = gui.ViewAxes;

% Get zoom obj
panH = findobj('Tag', 'Pan Tag');

% Check if there is anything currently in the plot window
if isempty(get(pH, 'Children'))
    
    set(panH(1),'State','off')
    return
    
end

% Get status
status = get(panH(1),'state');

switch lower(status)
    
    case 'on'
        
        % Set all other toolbar tools off
        set(findobj('tag','Zoom Tag'),'State','off')
        
        % Turn zoom on
        h = pan(gui.Window);
        set(h,'Enable','on')
        
    case 'off'
        
        h = pan(gui.Window);
        set(h,'Enable','off')
        
end


end