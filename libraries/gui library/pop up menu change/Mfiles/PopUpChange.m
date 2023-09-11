function PopUpChange(fig,~,Selection)
% RadioButtonCallback - GUI Callback that manipulates ui-controls on the
% GUI depending on parameter selection
%
% Syntax:  RadioButtonCallback(guiFig,~,Panel)
%
% Inputs:
%    guiFig - handle that identifies the figure
%    Panel - String specifying which panel to manipulate
%
% Outputs:
%    non
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% August 2013; Last revision: 08-August-2013

%------------- BEGIN CODE --------------

gui = % guidata(fig); % Get GUI data

switch Selection
    
    case 'Curvature'
        
        switch get(gui.CurvatureStatus,'value')
            case 1
                set(gui.CurvatureValue, 'Enable', 'On')
            case 2
                set(gui.CurvatureValue, 'Enable', 'Off')   
        end
                
    case 'LFS'
        
        switch get(gui.LFSStatus,'value')
            case 1
                set(gui.LFSValue, 'Enable', 'Off')
            case 2
                set(gui.LFSValue, 'Enable', 'On')   
        end
        
    case 'Elevation'
        
        switch get(gui.ElevStatus,'value')
            case 1
                set(gui.ElevValue, 'Enable', 'Off')
            case 2
                set(gui.ElevValue, 'Enable', 'On')   
        end
        
        
    case 'Tidal'
        
        switch get(gui.TidalStatus,'value')
            case 1
                set(gui.TidalValue, 'Enable', 'Off')
            otherwise
                set(gui.TidalValue, 'Enable', 'On')   
        end
        
    case 'Grading'
        
        switch get(gui.GradingStatus,'value')
            case 1
                set(gui.GradingValue, 'Enable', 'On')
            case 2
                set(gui.GradingValue, 'Enable', 'Off')   
        end
        
end
