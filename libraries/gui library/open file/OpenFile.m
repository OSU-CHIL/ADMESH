function OpenFile(varargin)
% OpenFile - GUI Callback that loads a mesh or edge structure file
% file.
%
% Syntax:  OpenMatCallback(varargin)
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

%---------------------------- BEGIN CODE ---------------------------------

%------------------------------------------------------------------------------
% Get GUI data
%------------------------------------------------------------------------------
% gui = % guidata(varargin{1});
app = varargin{1};

%------------------------------------------------------------------------------
% Get filename & location from user
%------------------------------------------------------------------------------
app.ProgressBarButton.Text = 'Select a file...';
app.ProgressBarButton.Icon = '';

% Ask the user to select a file
[filename, pathname] = uigetfile(...
    {'*.mat;*.14;*.grd;*.2dm;*.shp','Files (*.mat,*.14,*.grd,*.2dm,*.shp)'},'Select a file');

% If user cancels
if filename == 0
    app.ProgressBarButton.Text = 'Ready.';
    return
end

%------------------------------------------------------------------------------
% Determine the type of file we're reading in
%------------------------------------------------------------------------------
[~,~,ext] = fileparts(filename);

% Turn off colormap
SetContourStatus(app,'off');

switch ext
    
    case '.mat'
        
        % Reset domain variables
        app.FilePath   = [pathname filename];
        app.PTS        = [];
        app.xyzFun     = [];
        app.MESH       = [];
        
        % Read in mat file
        [app,status] = ReadMat([pathname filename],app);
        
        if status == 0
            app.ProgressBarButton.Text = 'Ready.'; drawnow;
            errordlg('There was an error reading in the file.','ADMESH')
            return
        end
        
        % Plot Boundary
        PlotEdgeStructure(app,.1);
        
        
    case '.14'
        
        % Reset domain variables
        app.FilePath   = [pathname filename];
        app.PTS        = [];
        app.xyzFun     = [];
        app.MESH       = [];
        
        [app.MESH,app.xyzFun,status] = Read14File([pathname filename],app);
        
        if status == 0
            app.ProgressBarButton.Text = 'Ready.'; drawnow;
            errordlg('There was an error reading in the file.','ADMESH')
            return
        end
                
        % Plot mesh
        PlotMesh(app,.1);
        
    case '.grd'
        
        % Reset domain variables
        app.FilePath   = [pathname filename];
        app.PTS        = [];
        app.xyzFun     = [];
        app.MESH       = [];
        
        [app.MESH,app.xyzFun,status] = Read14File([pathname filename],app);
        
        if status == 0
            app.ProgressBarButton.Text = 'Ready.'; drawnow;
            errordlg('There was an error reading in the file.','ADMESH')
            return
        end
        
        % Plot mesh
        PlotMesh(app,.1);

    case '.2dm'
        
        % Reset domain variables
        app.FilePath   = [pathname filename];
        app.PTS        = [];
        app.xyzFun     = [];
        app.MESH       = [];
        
        [app.MESH,app.xyzFun,status] = Read2DM([pathname filename]);
        
        if status == 0
            app.ProgressBarButton.Text = 'Ready.'; drawnow;
            errordlg('There was an error reading in the file.','ADMESH')
            return
        end

        app.ProgressBarButton.Text = 'Ready.'; drawnow;
        
    case '.shp'
        
        choice = questdlg(...
            'Do you want to save settings?', ...
            'ADmesh', ...
            'Yes','No','Yes');
    
        if strcmpi(choice,'yes')
            [file, path] = uiputfile(...
                {'*.mat','Files (*.mat)'},'Select a file to save settings');
            if ~file % if abort selecting a file
                path = [];
                file = [];
            end
        else
            path = [];
            file = [];
        end
        
        % Reset domain variables
        app.FilePath   = [path file];
        app.PTS        = [];
        app.xyzFun     = [];
        app.MESH       = [];
        
        [app,status] = ReadShapefile([pathname filename],app);
        
        if ~isempty(app.FilePath)
            PTS = app.PTS;
            Settings = SaveSettings(app);
            save(app.FilePath,'PTS','Settings');
        end
        
        if status == 0
            app.ProgressBarButton.Text = 'Ready.'; drawnow;
            errordlg('There was an error reading in the file.','ADMESH')
            return
        end
        
        % Plot Boundary
        PlotEdgeStructure(app,.1);
        
end

end