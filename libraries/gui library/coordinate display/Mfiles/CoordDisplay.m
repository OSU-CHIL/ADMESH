function CoordDisplay(varargin)
% PlotHoverCallback - GUI Callback that displays coordinate
% file.
%
% Syntax:  PlotHoverCallback(varargin)
%
% Inputs:
%    guiFig - handle that identifies the figure
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

persistent pH st 

if isempty(pH) || ~ishandle(pH);
    
    gui = % guidata(varargin{1});
    
    if isfield(gui,'ViewAxes')
        
        pH = gui.ViewAxes;
        st = gui.st;
    end
    return
end

if ~ishandle(pH); return; end

%---------------------------------------------------------------------
% Check for data in plot window
%---------------------------------------------------------------------
if isempty(get(pH, 'Children'));
    st.setText('')
    return
end

%---------------------------------------------------------------------
% Is mouse in current window?
%---------------------------------------------------------------------
% Get current point
pt = get(pH, 'CurrentPoint');
xInd = pt(1, 1);
yInd = pt(1, 2);

% check if its within axes limits
xLim = get(pH, 'XLim');
yLim = get(pH, 'YLim');

if xInd < xLim(1) || xInd > xLim(2) || yInd < yLim(1) || yInd > yLim(2)
    st.setText('')
    drawnow
    return;
end

%---------------------------------------------------------------------
% Display current coordinates
%---------------------------------------------------------------------
st.setText(['( ' num2str(xInd,'%.3f') ' , ' num2str(yInd,'%.3f') ' )'])

drawnow

end