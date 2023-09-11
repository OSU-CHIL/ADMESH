function CreateNodeString(varargin)
% CreateNodeString - Create Node String
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
% Determine if an edge structure is loaded
%-----------------------------------------------------------------------------
ES_Handle = findobj('Tag','Edge Structure');

if isempty(ES_Handle)
    
    toolH = findobj('Tag','create node string tool');
    
    set(toolH,'state','off')
    
    return

end 

%------------------------------------------------------------------------------
% Get status
%-----------------------------------------------------------------------------
[~,~,status] = deal(varargin{:});

if strcmp(status,'On') % If off
    
    zoomH = findobj('Tag','Zoom Tag');
    panH  = findobj('Tag','Pan Tag');
    
    % Set all other toolbar tools off
    set([zoomH panH],'State','off')
    
    % Change line style on the domain
    set(ES_Handle, 'marker','.','markersize',15)
    
    % Set the button down function
    set(ES_Handle, 'ButtonDownFcn',@SelectNodeString)
    
    % Perform set up
    SelectNodeString('Set Up');
    
elseif strcmp(status,'Off') % If off
    
    % Clear mouse click function
    set(ES_Handle, 'ButtonDownFcn','')
    
    % Remove markers
    set(ES_Handle, 'marker','none')
    
    % Turn off toggle tool
    toolH = findobj('Tag','create node string tool');
    
    set(toolH,'state','off')
    
    return
    
end

% %------------------------------------------------------------------------------
% % Get GUI data
% %------------------------------------------------------------------------------
% guiFig = findobj('Tag','ADmesh Figure'); guiH = % guidata(guiFig);
% 
% %------------------------------------------------------------------------------
% % Get Plot handle
% %-----------------------------------------------------------------------------
% pH = guiH.PlotWindow;


















