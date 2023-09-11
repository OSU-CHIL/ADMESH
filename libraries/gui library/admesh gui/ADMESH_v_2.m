function ADMESH_v_2(app)
%--------------------------------------------------------------------------
% ADmesh_v_1_1_12 - ADMESH Graphical User Interface Script
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% Last revision: 10-June-2014
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Check Matlab Version
%--------------------------------------------------------------------------
ver = version('-release'); % Check version

if str2double(ver(1:end-1)) < 2012
    errordlg('You need Matlab 2012 or newer in order to use ADMESH.','ADMESH')
end

%--------------------------------------------------------------------------
% Check for existing ADMESH program
%--------------------------------------------------------------------------
if ~isempty(findobj('Tag','ADMESH Figure'))
    return
end

%--------------------------------------------------------------------------
% Create the user interface for the application 
%--------------------------------------------------------------------------

gui = struct(); % Initialize gui struct

%   Set the figure window size values
scrsz = get(0, 'ScreenSize');
sx = scrsz(3); sy = scrsz(4); MaxSX = round(sx*.85); MaxSY = round(sy*.85);
XBorder = (sx-MaxSX)/2; YBorder = (sy-MaxSY)/2;

% Figure window------------------------------------------------------------
gui.Window = figure( ...
    'Name'                  , 'ADMESH - Advanced Mesh Generator | The C.H.I.L ~ The Ohio State University', ...
    'NumberTitle'           , 'off', ...
    'MenuBar'               , 'none', ...
    'Toolbar'               , 'none', ...
    'Visible'               , 'off' , ...
    'Renderer'              , 'opengl' ,...
    'Color'                 , get(0,'DefaultUicontrolBackgroundColor'),...
    'Tag'                   , 'ADMESH Figure'                       ,...
    'DockControls'          , 'off'                                 ,...
    'WindowButtonMotionFcn' , @CoordDisplay                    ,...
    'CloseRequestFcn'       , @CloseADMESH                          ,...
    'HandleVisibility'      , 'callback'                            ,...
    'Position'              ,[ XBorder, YBorder, MaxSX, MaxSY ]);
%--------------------------------------------------------------------------

% + File menu--------------------------------------------------------------
gui.FileMenu =uimenu( gui.Window, 'Label', 'File' );
%uimenu( gui.FileMenu, 'Label', 'Open File...', 'Callback', @OpenFile );     % Open File
%uimenu( gui.FileMenu, 'Label', 'Save As (.14)...', 'Callback', @SaveMesh14File );   % Save File
%uimenu( gui.FileMenu, 'Label', 'Save As (.mat)...', 'Callback', @SaveMatFile );   % Save File
%uimenu( gui.FileMenu, 'Label', 'Save As (.kml)...', 'Callback', @SaveKMLFilev3 );   % Save File
%--------------------------------------------------------------------------

% + Edit menu--------------------------------------------------------------
gui.EditMenu =uimenu( gui.Window, 'Label', 'Edit' );
%uimenu( gui.EditMenu, 'Label', 'Clear...', 'Callback', @ClearWindow ); % Clear
%--------------------------------------------------------------------------

% + Import menu------------------------------------------------------------
gui.ImportMenu =uimenu( gui.Window, 'Label', 'Import' );
%uimenu(gui.ImportMenu,'Label','Mesh','Callback', @ImportMesh);
%uimenu(gui.ImportMenu,'Label','NOAA Coastline File','Callback', @ImportNOAACoastline);
%uimenu(gui.ImportMenu,'Label','Elevation File','Callback', @ImportElevationData);
%uimenu(gui.ImportMenu,'Label','KML File','Callback', @ImportKMLFile);
%--------------------------------------------------------------------------

% + Tool menu--------------------------------------------------------------
gui.ToolMenu =uimenu( gui.Window, 'Label', 'Tools' );
%%uimenu(gui.ToolMenu,'Label','Export Figure','Callback', @ExportFigure);
%uimenu(gui.ToolMenu,'Label','Extract Sub Domain','Callback', @ExtractSubDomain);
%uimenu(gui.ToolMenu,'Label','Mesh Sub Domain','Callback', @MeshSubDomain);
%%uimenu(gui.ToolMenu,'Label','Convert Mesh','Callback', @convertTrimesh2Quadmesh);
%--------------------------------------------------------------------------

% + Display menu-----------------------------------------------------------
gui.DisplayMenu =uimenu( gui.Window, 'Label', 'View' );

subMenu1        =uimenu( gui.DisplayMenu, 'Label', 'Contours' );
gui.ContourStatus(1) =uimenu(subMenu1,'Label','On','Callback',{@SetContourStatus,'on'});
gui.ContourStatus(2) =uimenu(subMenu1,'Label','Off','Checked','On','Callback',{@SetContourStatus,'off'});

subMenu2        =uimenu( gui.DisplayMenu, 'Label', 'Color Map' );
gui.Cmap(1) =uimenu(subMenu2,'Label','Land & Sea','Checked','On','Callback',{@SetContourStatus,'Land & Sea'});
gui.Cmap(2) =uimenu(subMenu2,'Label','Jet'       ,'Checked','Off','Callback',{@SetContourStatus,'Jet'});
gui.Cmap(3) =uimenu(subMenu2,'Label','Parula'    ,'Checked','Off','Callback',{@SetContourStatus,'Parula'});

%%uimenu(gui.DisplayMenu,'Label','Sub Region','Callback', @ViewSubRegion);
%--------------------------------------------------------------------------

% Create Toolbar-----------------------------------------------------------
hToolbar = uitoolbar('Parent',gui.Window);


% + Zoom in Toggle Tool
uitoggletool(...
    'Parent'            , hToolbar,...
    'TooltipString'     , 'Zoom In',...
    'Separator'         , 'on',....
    'CData'             , GetToolBarIcon('zoom in'),...
    'HandleVisibility'  , 'Callback', ...
    'Tag'               , 'Zoom Tag',...
    'OnCallback'        , @zoomIn ,...
    'OffCallback'       , @zoomIn);

%  Reset to Original View Tool
uipushtool(...
    'Parent'            , hToolbar,...
    'TooltipString'     , 'Reset to Original View',...
    'Separator'         , 'on',....
    'CData'             , GetToolBarIcon('double tool arrow'),...
    'HandleVisibility'  , 'callback', ...
    'ClickedCallback'   , @resetView );

% Pan Toggle Tool
uitoggletool(...
    'Parent'            , hToolbar,...
    'TooltipString'     , 'Pan',...
    'Separator'         , 'on',....
    'CData'             , GetToolBarIcon('pan'),...
    'HandleVisibility'  , 'callback', ...
    'Tag'               , 'Pan Tag',...
    'OnCallback'        , @PanCallback ,...
    'OffCallback'       , @PanCallback);

% Select Node Toggle Tool
uitoggletool(...
    'Parent'            , hToolbar,...
    'Separator'         , 'on',....
    'TooltipString'     , 'Select/Adjust Nodes',...
    'CData'             , GetToolBarIcon('select node'),...
    'HandleVisibility'  , 'callback', ...
    'Tag'               , 'Move Node Tag',...
    'OnCallback'        , {@AdjustMeshNodes, 'On'} ,...
    'OffCallback'       , {@AdjustMeshNodes, 'Off'});

% Element Threshold Toggle Tool
uipushtool(...
    'Parent'            , hToolbar,...
    'Separator'         , 'on',....
    'TooltipString'     , 'Set Element Quality Threshold',...
    'CData'             , GetToolBarIcon('bad element threshold'),...
    'HandleVisibility'  , 'callback', ...
    'ClickedCallback'   , @DisplayBadElements );

% Zoom to Bad Element Toggle Tool
uipushtool(...
    'Parent'            , hToolbar,...
    'Separator'         , 'on',....
    'TooltipString'     , 'Zoom To Low Quality Element',...
    'CData'             , GetToolBarIcon('zoom to bad element'),...
    'HandleVisibility'  , 'callback', ...
    'ClickedCallback'   , @ZoomToBadElements );

uipushtool(...
    'Parent'            , hToolbar,...
    'Separator'         , 'on',....
    'TooltipString'     , 'Flag Location',...
    'CData'             , GetToolBarIcon('pin'),...
    'HandleVisibility'  , 'callback', ...
    'ClickedCallback'   , @PinLocation );
%--------------------------------------------------------------------------

% Arrange the main interface-----------------------------------------------
Layout = uix.VBox('Parent', gui.Window, 'Spacing', 5 );

% Top box
mainLayout = uix.HBoxFlex('Parent', Layout, 'Spacing', 5 );

% Bottom box
uix.Empty( 'Parent', Layout); % Fill space

% Set proportion ratio
set( Layout, 'Heights', [-40 -1] );

%--------------------------------------------------------------------------

% + Create the viewing panel-----------------------------------------------
% 2021-07-02 Younghun: Change from "uiextras" to "uix" to avoid error. Not
% sure how this change affects the codes
gui.ViewPanel = uix.Panel('Parent', mainLayout,'Padding',15);

% + Create the view for main plot
gui.ViewAxes = axes(...
    'Parent'    , gui.ViewPanel ,...
    'units'     , 'normalized',...
    'position'  , [0 0 1 1],...
    'color'     , 'none',...
    'Tag'       , 'ADMESH Plot Window',...
    'XTick'     , [],...
    'YTick'     , [],...
    'ZTick'     , [],...
    'box'       , 'off',...
    'PlotBoxAspectRatio', [1.5 1 1],...
    'DataAspectRatio', [1 1 1],...
    'ycolor'    , get(0,'DefaultUicontrolBackgroundColor'),...
    'xcolor'    , get(0,'DefaultUicontrolBackgroundColor'));

%--------------------------------------------------------------------------


% + Create the criteria panel
controlPanel = uix.Panel( 'Parent', mainLayout);

% + Adjust the main layout
set( mainLayout,'Widths', [-2,-.6] );

% + Create the criteria controls
controlLayout = uix.VBoxFlex( 'Parent', controlPanel,'Padding', 15, 'Spacing', 20 );

% + Create panel for settings
settingsPanel = uix.Panel( 'Parent', controlLayout,'Title', 'Mesh Criteria');

% + Create panel for display mesh info
resultsPanel = uix.Panel( 'Parent', controlLayout,'Title', 'Mesh Results');

% + Run Admesh Button
gui.RunAdmeshButton = uicontrol( 'Style', 'PushButton', ...
    'Parent', controlLayout, ...
    'String', 'Run ADMESH', ...
    'Callback', @ADmeshRoutine );

% Set proportion ratio
set( controlLayout, 'Heights', [-1 150 35] ); % Make the list fill the space

% Add uicontrol for mesh results
gui.resultsBox = uicontrol('style','text','string','','parent',resultsPanel,'FontSize',12,'HorizontalAlignment', 'left');

% Initialize grid for seperating criteria
g = uix.Grid( 'Parent', settingsPanel, 'Spacing', 10 ,'Padding',10);

% FIRST COLUMN
% uicontrol('style','text','string','Spatial Resolution','HorizontalAlignment', 'center','parent',g)
% uicontrol('style','text','string','Max Element Size (m):','HorizontalAlignment', 'center','parent',g)
% uicontrol('style','text','string','Min Element Size (m):','HorizontalAlignment', 'center','parent',g)
% uicontrol('style','text','string','Boundary Curvature:','HorizontalAlignment', 'center','parent',g)
% uicontrol('style','text','string','Local Feature Size:','HorizontalAlignment', 'center','parent',g)
% uicontrol('style','text','string','Elevation Gradients:','HorizontalAlignment', 'center','parent',g)
% uicontrol('style','text','string','Dominate Tide:','HorizontalAlignment', 'center','parent',g)
% uicontrol('style','text','string','Mesh Grading:','HorizontalAlignment', 'center','parent',g)
% uicontrol('style','text','string','View Mesh Generation:','HorizontalAlignment', 'center','parent',g)

% SECOND COLUMN
gui.ResolutionStatus    = uicontrol('style','popup','string',{'Low','Medium','High','High 2x','High 3x'},'parent',g,'backgroundcolor','w' );
gui.MaxElementSize      = uicontrol('style','edit','string','','HorizontalAlignment', 'center','parent',g,'backgroundcolor','w');
gui.MinElementSize      = uicontrol('style','edit','string','','HorizontalAlignment', 'center','parent',g,'backgroundcolor','w');
gui.CurvatureValue      = uicontrol('style','edit','string','20','HorizontalAlignment', 'center','parent',g,'backgroundcolor','w');
gui.LFSValue            = uicontrol('style','edit','string','1','HorizontalAlignment', 'center','parent',g,'backgroundcolor','w','enable','off');
gui.ElevValue           = uicontrol('style','edit','string','0.3','HorizontalAlignment', 'center','parent',g,'backgroundcolor','w','enable','off');
gui.TidalValue          = uicontrol('style','edit','string','25','HorizontalAlignment', 'center','parent',g,'backgroundcolor','w','enable','off');
gui.GradingValue        = uicontrol('style','edit','string','0.15','HorizontalAlignment', 'center','parent',g,'backgroundcolor','w');
gui.ViewStatus          = uicontrol('style','popup','string',{'Off','On'},'parent',g,'backgroundcolor','w' );

% THIRD COLUMN
uix.Empty( 'Parent', g); % Fill space
uix.Empty( 'Parent', g); % Fill space
uix.Empty( 'Parent', g); % Fill space
gui.CurvatureStatus = uicontrol('style','popup','string',{'On','Off'},'parent',g,'backgroundcolor','w','callback',{@PopUpChange, 'Curvature'} );
gui.LFSStatus       = uicontrol('style','popup','string',{'Off','On'},'parent',g,'backgroundcolor','w','callback',{@PopUpChange, 'LFS'} );
gui.ElevStatus      = uicontrol('style','popup','string',{'Off','On'},'parent',g,'backgroundcolor','w' ,'callback',{@PopUpChange, 'Elevation'});
gui.TidalStatus     = uicontrol('style','popup','string',{'Off','M2','S2','N2','K2'},'parent',g,'backgroundcolor','w','callback',{@PopUpChange, 'Tidal'} );
gui.GradingStatus   = uicontrol('style','popup','string',{'On','Off'},'parent',g,'backgroundcolor','w' ,'callback',{@PopUpChange, 'Grading'});
uix.Empty( 'Parent', g); % Fill space

% Set the sizes
set( g, 'Heights', [25 25 25 25 25 25 25 25 25]  );

% + Empty Admesh variables
gui.PTS                = [];       % Edge Structure
gui.xyzFun             = [];       % Elevation interpolant
gui.MESH               = [];       % Mesh data
gui.hmin               = [];       % Minimum element size
gui.MinEQ              = 0.3;      % Minimum element quality
gui.xLimits            = [0 1];    % x-axis limits
gui.yLimits            = [0 1];    % y-axis limits
gui.per                = 0;        % Percent offset
gui.DummyConstraint     = uicontrol('style','edit','string','0','HorizontalAlignment', 'center','parent',g,'backgroundcolor','w','enable','off');
% Make figure visible
set(gui.Window,'visible','On')
% drawnow; pause(.05)

% Maximize figure window
maxfig(gui.Window,1)

% Status bar
gui.sb = statusbar(gui.Window,'Ready.');

% Coordinate display (status text)
gui.st = javax.swing.JLabel('');
gui.sb.add(gui.st,'East');


% Save gui data
% guidata(gui.Window,gui)