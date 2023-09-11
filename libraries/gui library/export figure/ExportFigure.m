function ExportFigure(varargin)
% ExportFigureCallback - GUI Callback that exports the current plot
% file.
%
% Syntax:  ExportFigureCallback(guiFig,~,~)
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

%----------------------- BEGIN CODE -----------------------------


%---------------------------------------------------------------------
% Get GUI Data
%---------------------------------------------------------------------
fig = varargin{1}; gui = % guidata(fig);

%------------------------------------------------------------------------------
% Use system background color for GUI components
%------------------------------------------------------------------------------
BGFigColor = get(0,'DefaultUicontrolBackgroundColor');

%------------------------------------------------------------------------------
% Create new figure window
%------------------------------------------------------------------------------
h = figure(...
    'Toolbar'       , 'figure'                                      ,...
    'NumberTitle'   , 'off'                                         ,...
    'Color'         , BGFigColor                                    ,...
    'Name'          , 'ADMESH',...
    'Renderer'      , 'zbuffer');

%------------------------------------------------------------------------------
% Copy plot
%------------------------------------------------------------------------------
hNew = copyobj(gui.ViewAxes,h);

%------------------------------------------------------------------------------
% Paste to new figure
%------------------------------------------------------------------------------
set(hNew, ...
    'units','normalized',...
    'Position', [.075,.125,.875,.8],...
    'box'    , 'on',...
    'XTickMode', 'auto',...
    'YTickMode', 'auto',...
    'ycolor' ,'k',...
    'xcolor' ,'k',...
    'Ticklength',[0 0],...
    'color','white')

% Re-name tags for axis and all of it's children
set(hNew,'tag','ADMESH Export')
set(get(hNew,'Children'),'tag','ADMESH Export')
%set(h,'DeleteFcn',[])

%------------------------------------------------------------------------------
% Set colormapping 
%------------------------------------------------------------------------------
cmapping = colormap(gui.ViewAxes); colormap(cmapping)

%------------------------------------------------------------------------------
% Create border
%------------------------------------------------------------------------------
xlim = get(hNew,'XLim');
ylim = get(hNew,'YLim');

% Compute border widths
bwy = diff(ylim)/100; 
bwx = bwy;

% transparent borders
patch([xlim(1)+bwx,xlim(2)-bwx,xlim(2)-bwx,xlim(1)+bwx],ylim(1) + bwy*[0,0,1,1],'k','FaceColor','none','clipping','off')
patch([xlim(1)+bwx,xlim(2)-bwx,xlim(2)-bwx,xlim(1)+bwx],ylim(2) - bwy*[0,0,1,1],'k','FaceColor','none','clipping','off')
patch(xlim(1) + bwx*[0,0,1,1],[ylim(1)+bwy,ylim(2)-bwy,ylim(2)-bwy,ylim(1)+bwy],'k','FaceColor','none','clipping','off')
patch(xlim(2) - bwx*[0,0,1,1],[ylim(1)+bwy,ylim(2)-bwy,ylim(2)-bwy,ylim(1)+bwy],'k','FaceColor','none','clipping','off')

% Number of bars
nB = 5;
xspace = linspace(xlim(1),xlim(2),nB*2);

for k = 2:2:nB*2
    
    xp = [xspace(k-1) xspace(k-1) xspace(k) xspace(k)];
    patch(xp,ylim(1) + bwy*[0,1,1,0],'k','FaceColor','k','clipping','off')
    patch(xp,ylim(2) - bwy*[0,1,1,0],'k','FaceColor','k','clipping','off')
    
end

% Number of bars
nB = 5;
yspace = linspace(ylim(1),ylim(2),nB*2);

for k = 2:2:nB*2
    
    yp = [yspace(k-1) yspace(k) yspace(k) yspace(k-1)];
    patch(xlim(1) + bwx*[0,0,1,1],yp,'k','FaceColor','k','clipping','off')
    patch(xlim(2) - bwx*[0,0,1,1],yp,'k','FaceColor','k','clipping','off')
    
end


xlabel('X (m)')
ylabel('Y (m)')

k = find(gui.FilePath == '\',1,'last');

title(gui.FilePath(k+1:end-4),'Interpreter','none')

daspect([1 1 1])

maxfig(h,1);


end