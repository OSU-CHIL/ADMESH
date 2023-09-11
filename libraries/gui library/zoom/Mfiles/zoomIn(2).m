function zoomIn(varargin)


% Load GUI data
gui = % guidata(varargin{1});

% Get axes handle
pH = gui.ViewAxes;

% Get zoom obj
zoomH = findobj('Tag', 'Zoom Tag');

% Check if there is anything currently in the plot window
if isempty(get(pH, 'Children'))
    
    set(zoomH(1),'State','off')
    return
    
end

% Get status
status = get(zoomH(1),'state');

switch lower(status)
    
    case 'on'
        
        % Set all other toolbar tools off
        set(findobj('tag','Pan Tag'),'State','off')
        
        %interactivemouse(gui.Window,'on')
        
        % Turn zoom on
        h = zoom(gui.Window);
        set(h,'Direction','in');
        set(h,'RightClickAction','InverseZoom');
        set(h,'Enable','on')
        
    case 'off'
        
        %interactivemouse(gui.Window,'off')
        h = zoom(gui.Window);
        set(h,'Enable','off')
        
end


%     function zoompostcallback(fig,pH)
%         
%       gui = % guidata(fig);
%       
%       set(gui.ViewAxes,'DataAspectRatio',gui.daspect)
% 
%     end

end