function ImportElevationData(varargin)
% LoadBathyFile - GUI Callback that loads a bathymetry scatter set
%
% Syntax:  LoadBathyFile(guiFig,~,~)
%
% Inputs:
%    guiFig - handle that identifies the figure
%
% Outputs:
%    non
%
% Other m-files required: non
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% August 2013; Last revision: 08-August-2013

%------------- BEGIN CODE --------------

%--------------------------------------------------------------------------
% Get GUI data
%--------------------------------------------------------------------------
app = varargin{1};

% Turn off colormap
SetContourStatus(app,'Off');

%--------------------------------------------------------------------------
% Check to see if edge structure is loaded
%--------------------------------------------------------------------------
PTS = app.PTS; xyzFun = app.xyzFun;

%--------------------------------------------------------------------------
% Check for edge structure
%--------------------------------------------------------------------------
if isempty(PTS)
    
%     warndlg('You must first load an edge structure file (.mat or .14) file.'...
%         ,'Error');
%     return
    
end

%--------------------------------------------------------------------------
% Check for existing bathymetry data
%--------------------------------------------------------------------------
if ~isempty(xyzFun) % Bathymetry data is loaded
    
    msg = ['There is already existing bathymetry. Would you like '...
        'to replace the existing bathymetry or use this ' ...
        'bathymetry temporarily?.'];
    choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'Replace','Use Temporarily','Cancel'},'DefaultOption',1,'Icon','Warning');
    drawnow; pause(0.05);  % this innocent line prevents the Matlab hang
    
    if strcmp(choice,'Cancel')
        return
    end
    
else
    
    choice = 'Replace';
    
end

%--------------------------------------------------------------------------
% Get the file name from the user
%--------------------------------------------------------------------------
app.ProgressBarButton.Text = 'Select a (.xyz) or (.asc) file....'; drawnow;

[file, path] = uigetfile({'*.xyz;*.asc;*.tiff;*.tif','Elevation Files (*.xyz,*.asc,*.tiff,*.tif)'},...
    'Select a file',cd);

%--------------------------------------------------------------------------
% Did the user make a selection?
%--------------------------------------------------------------------------
if ~file
    app.ProgressBarButton.Text = 'Ready'; drawnow;
    return
end

%--------------------------------------------------------------------------
% Determine the file type
%--------------------------------------------------------------------------
app.ProgressBarButton.Text = 'Reading in data....'; drawnow;

% Get file extension
[~,~,ext] = fileparts([path file]);

%--------------------------------------------------------------------------
% Read in XYZ file
%--------------------------------------------------------------------------
if strcmp(ext,'.xyz')
    
    fid=fopen([path,file], 'r');    % Open file
    
    try
        
        app.ProgressBarButton.Text = 'Reading in elevation data...'; drawnow;
        xyz = textscan(fid, '%f %f %f'); % Read in data
        
        fclose(fid); % Close file
    catch
        
        fclose(fid);
        
        % Assign string to variable
        ERRORSTRING = ['There was a problem reading in this file. '...
            'Check to make sure the file contains only three columns '...
            'of numerical data'];
        
        % Assign window name
        DLGNAME = 'ESRI ASCII Grid Read-in Error';
        
        % Display error message
        errordlg(ERRORSTRING,DLGNAME)
        
        app.ProgressBarButton.Text = 'Ready'; drawnow;
        
        return
        
    end
    
    % Convert cell array to matrix
    xyz = cell2mat(xyz);
    
    % Convert coordinates
    if isfield(PTS,'cpplon') && ~isempty(PTS.cpplon) % Younghun added
        [xyz(:,1),xyz(:,2)] = Geo2Meters(xyz(:,1),xyz(:,2),PTS.cpplon,PTS.cpplat);
    end
    
    app.ProgressBarButton.Text = 'Creating elevation interpolant function...'; drawnow;
    xyzFun = scatteredInterpolant(xyz(:,1),xyz(:,2),xyz(:,3),'linear','nearest');
end

%--------------------------------------------------------------------------
% Read in ASC file
%--------------------------------------------------------------------------
if strcmp(ext,'.asc')
    
    try        

        % Open file and read in header information
        fid=fopen([path,file],'r');
        
        % Read in header
        app.ProgressBarButton.Text = 'Reading in elevation data...'; drawnow;
        textHeader = textscan(fid, '%s %f', 6);
        
        % Check for optional no data value
        match = strcmpi('nodata_value', textHeader{1}(end));
        
        if match
            nheader = 6;
        else
            nheader = 5;
        end
        
        % Set the file position indicator to the beginning of the file
        frewind(fid); 
        
        % Read in header
        textHeader = textscan(fid, '%s %f', nheader);
        
        % Get the number of rows
        nrows = textHeader{2}(strcmpi('nrows', textHeader{1}),1);
        
        % Get the number of columns
        ncols = textHeader{2}(strcmpi('ncols', textHeader{1}),1);
        
        % Read in bathymetry grid
        z = fscanf(fid,'%lg',[ncols nrows]);
        
        % Close the file
        fclose(fid);
        
        % If z is empty then there was an error reading the header
        if isempty(z)
            
            % Assign string to variable
            ERRORSTRING = ['There was a problem reading in this file. '...
                'Check to make sure the file contains only three columns '...
                'of numerical data'];
            
            % Assign window name
            DLGNAME = 'ESRI ASCII Grid Read-in Error';
            
            % Display error message
            errordlg(ERRORSTRING,DLGNAME)
            
            % gui.sb.setText('Ready')
            return
        end
        
        % If no data value exist, substitute NaN's
        if nheader == 6
            
            nodata = textHeader{2}(strcmpi('nodata_value', textHeader{1}),1);
            
            z(z==nodata) = NaN;
            
        end
                
        % Fill in NaN's and transpose/flip
        if any(isnan(z(:)))
            msg = 'NaNs are found in the elevation data. How do you proceed?';
            choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
                'Options',{'In-paint NaNs','Ignore NaNs'},'DefaultOption',2,'Icon','Warning');
            switch choice
                case 'In-paint NaNs'
                    z = flipud(inpaint_nans(z',4));
                case 'Ignore NaNs'
                    z = flipud(z');
            end
        end       
                
        % Create X & Y with meshgrid
        
        yllc = textHeader{2}(strcmpi('yllcorner', textHeader{1}),1);
        
        xllc = textHeader{2}(strcmpi('xllcorner', textHeader{1}),1); 
        
        cs = textHeader{2}(strcmpi('cellsize', textHeader{1}),1); 
        
        y = (yllc+(0.5*cs))+((nrows-1)*cs):-cs:yllc+(0.5*cs);
        
        x = xllc+(0.5*cs):cs:(xllc+(0.5*cs))+((ncols-1)*cs);
        
        [x,y] = meshgrid(x,fliplr(y));
               
    catch
        app.ProgressBarButton.Text = 'Ready'; drawnow;
        return
    end
    
    % Convert coordinates
    if isfield(PTS,'cpplon') && ~isempty(PTS.cpplon)
        app.ProgressBarButton.Text = 'Converting to cartesian coordinate system...'; drawnow;
        [x,y] = Geo2Meters(x,y,PTS.cpplon,PTS.cpplat);
    end
    
    %-------------------------------------------------------------------------
    % Generate scattered interpolant function for bathymetry
    %--------------------------------------------------------------------------
    app.ProgressBarButton.Text = 'Creating an interpolant function for the elevation data set...'; drawnow;

    switch choice
        case 'In-paint NaNs'
            xyzFun = griddedInterpolant(x',y',z','linear','nearest');
        case 'Ignore NaNs'
            I = ~isnan(z);
            xyzFun = scatteredInterpolant(x(I),y(I),z(I),'nearest','nearest');
        otherwise
            xyzFun = griddedInterpolant(x',y',z','linear','nearest');
    end
        
end

if any(strcmpi(ext,{'.tiff','.tif'}))
    
    app.ProgressBarButton.Text = 'Reading in elevation data...'; drawnow;
    
    Z = imread([path,file]);
    Z = double(Z);
    Z = flipud(Z);
    Z(Z < -1e20) = nan;
    Z(Z==0) = nan;
    gtinfo = geotiffinfo([path,file]);
    
    ax = gtinfo.BoundingBox;
    x = linspace(ax(1,1),ax(2,1),size(Z,2));
    y = linspace(ax(1,2),ax(2,2),size(Z,1));
    [X,Y] = meshgrid(x,y);
    
    app.ProgressBarButton.Text = 'Creating an interpolant function for the elevation data set...'; drawnow;
    I = ~isnan(Z);
    xyzFun = scatteredInterpolant(X(I),Y(I),Z(I),'nearest','nearest');
    
    ElevationDataFilename = [path,file];
    app.ElevationDataFilename = ElevationDataFilename;

end

xyzFun.Values = -xyzFun.Values;
xyzFun = CoordinateConversion(app,xyzFun,'auto',PTS.cpplon,PTS.cpplat);
app.xyzFun = xyzFun;

% Update GUI data
% guidata(fig,gui);
pause(.01)

% Replace data in file by user request
if strcmp(choice,'Replace') && ~isempty(app.FilePath)
    if any(strcmpi(ext,{'.tiff','.tif'}))
        try
            save(app.FilePath,'ElevationDataFilename','xyzFun','-append');
        end
    else
        try
            save(app.FilePath,'xyzFun','-append')
        end
    end
end

app.ProgressBarButton.Text = 'Ready'; drawnow;


end