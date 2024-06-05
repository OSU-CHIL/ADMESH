function [app,status] = ReadMat(file,app)
% OpenMeshCallback - GUI Callback that calls on Read_14 to read in a fort.14
% file.
%
% Syntax:  OpenMeshCallback(guiFig,~,~)
%
% Inputs:
%    guiFig - handle that identifies the figure
%
% Outputs:
%    non
%
% Other m-files required: READ_14
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
%---------------------------- BEGIN CODE --------------------------------------

%------------------------------------------------------------------------------
% Initialize output
%------------------------------------------------------------------------------
PTS = [];       % Edge Structure
xyzFun = [];    % Elevation
status = 1;     % Error flag
ElevationDataFilename = []; 
%------------------------------------------------------------------------------
% Check file
%------------------------------------------------------------------------------
progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',...
    'Checking file...','Indeterminate','on');

matData = struct2cell(whos('-file',file));

if ~any(strcmp(matData(1,:),'PTS'))
    
    msg = ['No edge structure exists in this file.' ...
        ' Make sure you are selecting the correct file for ADMESH.'];
    uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');
    
    app.ProgressBarButton.Text = 'Ready'; drawnow;

    status = 0;
    
    return
    
end

%------------------------------------------------------------------------------
% Erase any sub meshing stuff
%------------------------------------------------------------------------------
%CancelADMESHSubMesh('c')

%------------------------------------------------------------------------------
% Load edge structure 
%------------------------------------------------------------------------------
if any(strcmp(matData(1,:),'PTS'))
progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',...
    'Loading Edge Structure...','Indeterminate','on');
    
    % Initialize as 0
    PTS = 0;
    
    % Load PTS variable
    load(file, 'PTS')
    
    % Check fields
    if ~isfield(PTS,'Poly')
        
        msg = 'The edge structure has missing fields.';
        uiconfirm(app.UIFigure,msg,'ADMESH',...
            'Options',{'OK'},'DefaultOption',1,'Icon','Error');
        
        app.ProgressBarButton.Text = 'Ready'; drawnow;
        
        PTS = [];
        xyzFun = [];
        status = 0;
        return
        
    end

    if ~isfield(PTS,'Constraints')
        PTS.Constraints = [];
    end
    
    if ~isfield(PTS,'cpplon') || ~isfield(PTS,'cpplat')
        PTS.cpplon = []; PTS.cpplat = [];
    end
     
end

%------------------------------------------------------------------------------
% Load elevation data
%------------------------------------------------------------------------------
if any(strcmp(matData(1,:),'xyzFun'))
progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',...
    'Loading elevation data...','Indeterminate','on');
    
    % Initialize as 0
    xyzFun = [];
    
    % Load xyzFun variable
    load(file, 'xyzFun')
    
end

if any(strcmp(matData(1,:),'ElevationDataFilename'))
    
    % Load filename of elevation data
    load(file, 'ElevationDataFilename')
    
end


%------------------------------------------------------------------------------
% Load previous settings
%------------------------------------------------------------------------------
if any(strcmp(matData(1,:),'Settings'))
    
    progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',...
        'Loading settings...','Indeterminate','on');

    load(file, 'Settings')
        
    LoadSettings(Settings,app); 
    
end

app.PTS = PTS;
app.xyzFun = xyzFun;
app.ElevationDataFilename = ElevationDataFilename;

close(progdlg);

end