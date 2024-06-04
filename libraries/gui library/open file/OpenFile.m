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
SetContourStatus(app,'Off');

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
        
        % Convert to cartesian coordinates if needed.
        [app.PTS,status] = CoordinateConversion(app,app.PTS,'auto');

        if status
            app.xyzFun = CoordinateConversion(app,app.xyzFun,'forward',app.PTS.cpplon,app.PTS.cpplat);
            PlotEdgeStructure(app,.1);
        end

    case {'.14','.grd'}
        
        % Reset domain variables
        app.FilePath   = [pathname filename];
        app.PTS        = [];
        app.xyzFun     = [];
        app.MESH       = [];
        
        [app.MESH,~,status] = Read14File([pathname filename],app);
        
        if status == 0
            app.ProgressBarButton.Text = 'Ready.'; drawnow;
            errordlg('There was an error reading in the file.','ADMESH')
            return
        end

        % Plot mesh
        PlotMesh(app,.1);
        
        % Convert to cartesian coordinates if needed.
        [app.MESH,status] = CoordinateConversion(app,app.MESH,'auto');

        if status
            PlotMesh(app,.1);
        end

    case '.2dm'
        
        % Reset domain variables
        app.FilePath   = [pathname filename];
        app.PTS        = [];
        app.xyzFun     = [];
        app.MESH       = [];
        
        [app.MESH,~,status] = Read_2DM([pathname filename]);
        
        if status == 0
            app.ProgressBarButton.Text = 'Ready.'; drawnow;
            errordlg('There was an error reading in the file.','ADMESH')
            return
        end

        % Convert to cartesian coordinates if needed.
        app.MESH   = CoordinateConversion(app,app.MESH,'auto');

        app.ProgressBarButton.Text = 'Ready.'; drawnow;
        
    case '.shp'
        
        msg = 'Do you want to create a file to save settings?';
        choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
            'Options',{'Yes','No'},'DefaultOption',1,'Icon','Warning');
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
        
        % Convert to cartesian coordinates if needed.
        [app.PTS,status] = CoordinateConversion(app,app.PTS,'auto');
        
        if status
            app.xyzFun = CoordinateConversion(app,app.xyzFun,'forward',app.PTS.cpplon,app.PTS.cpplat);
            PlotEdgeStructure(app,.1);
        end

end

end