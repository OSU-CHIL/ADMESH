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
app.ProgressBarButton.Text = 'Checking file...'; drawnow;

matData = struct2cell(whos('-file',file));

if ~any(strcmp(matData(1,:),'PTS'))
    
    warndlg(['No edge structure exists in this file.' ...
        ' Make sure you are selecting the correct file for ADMESH.'],'Error');
    
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
app.ProgressBarButton.Text = 'Loading Edge Structure...'; drawnow;

if any(strcmp(matData(1,:),'PTS'))
    
    % Initialize as 0
    PTS = 0;
    
    % Load PTS variable
    load(file, 'PTS')
    
    % Check fields
    if ~isfield(PTS,'Poly')
        
        warndlg('The edge structure has missing fields.','Error');
        
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
app.ProgressBarButton.Text = 'Loading elevation data...'; drawnow;

if any(strcmp(matData(1,:),'xyzFun'))
    
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
    
    app.ProgressBarButton.Text = 'Loading previous settings...'; drawnow;
    
    load(file, 'Settings')
        
    LoadSettings(Settings,app); 
    
end

app.PTS = PTS;
app.xyzFun = xyzFun;
app.ElevationDataFilename = ElevationDataFilename;

app.ProgressBarButton.Text = 'Ready'; drawnow;

end