function mesh = READ_2DM(file)
% READ_14 - Reads in domain and it's attributes
%
% Syntax:  [PTS, TRI, bathy, IBtype, CPPLAT, CPPLON] = READ_2DM(varargin)
%
% Inputs:
%    varargin - cell containing the GUI handles
%
% Outputs:
%    PTS - Data structure with fields x & y
%    PTS(1).x = x-coordinates of first polygon
%    PTS(1).y = y-coordinates of first polygon
%    xyz - coordinates and elevation at each node in tri - used for
%    elevation in ADmesh
%    IBtype - Cell Structure defining boundary conditions in PTS
%    CPPLAT, CPPLON - lat lon references for conversion from geographic
%    coordinates to cartesian.
%
% See also: Open14Callback
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
%------------- BEGIN CODE --------------

%------------------------------------------------------------------------------
% Ask the user if there is a bathymetry file to be included
%------------------------------------------------------------------------------

% Construct a questdlg with two options
msg = 'Is there a bathymetry file to go along with this?';
choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'Yes','No'},'DefaultOption',2,'Icon','Warning');

% Handle response
switch choice
    case 'Yes'
        % Ask the user to select a file
        [bathyfilename, bathypathname] = uigetfile({'*.dat'},'Select a file');
        bathymetry_included = 'Yes';
    case 'Nope'
        bathymetry_included = 'No';
end

pause(.01)

%------------------------------------------------------------------------------
% Read the mesh from the 2dm file
%------------------------------------------------------------------------------
fid=fopen(file, 'r');  % Open file

uiStatusBar('Scanning mesh file...')

% Read in file into cell g
g = textscan(fid,'%s','delimiter','\n'); g = g{1}; % Unpack cell

fclose(fid); % Close the file

% Find the total number of header lines
k = strfind(g, 'MESH'); nHeader = sum(~cellfun('isempty',k));

% Find the total number of all element connectivities
k = strfind(g, 'E3T '); nElems = sum(~cellfun('isempty',k));

% Find the total number of all element vertices
k = strfind(g, 'ND '); nVerts = sum(~cellfun('isempty',k));

% Find the total number of all element strings
k = strfind(g, 'NS '); nStrings = sum(~cellfun('isempty',k));

% Re-open file and read in data
fid=fopen(file, 'r');

%------------------------------------------------------------------------------
% Read in element connectivity
%------------------------------------------------------------------------------
uiStatusBar('Reading in element connectivity...')
C = textscan(fid, '%*s %*f %f %f %f %*d %*[^\n]', nElems,'Headerlines',nHeader);
mesh.tri = cell2mat(C); clear C

%------------------------------------------------------------------------------
% Read in vertices
%------------------------------------------------------------------------------
uiStatusBar('Reading in vertices...')
C = textscan(fid, '%*s %*d %f %f %f %*[^\n]', nVerts);
[x,y,z] = deal(C{:}); clear C

if size(x,1) == 1; x = x'; end
if size(y,1) == 1; y = y'; end
if size(z,1) == 1; z = z'; end

% Read in Open Ocean Boundary Node Strings
uiStatusBar('Reading in nodal strings...')
C = textscan(fid, '%*s %d %d %d %d %d %d %d %d %d %d %*[^\n]', nStrings);
C = cell2mat(C);

fclose(fid); % Close file

% Convert to cartesian coordinates if needed.
if ~sum(~(x < 180)) && ~sum(~(x > -180)) && ...
        ~sum(~(y <  90)) && ~sum(~(y) >  -90)
    
    uiStatusBar('Converting to cartesian coordinate system...')
    
    % Convert to XY (meters)
    [x,y,mesh.cpplon,mesh.cpplat] = Geo2Meters(x,y);

end

% Assign nodal coordinates
mesh.p = [x y z];

C = reshape(C',1,numel(C));

% Removes zeros in node string matrix
C(C == 0) = []; C = C';

nNodeStrings = sum(C<0);

% Intitialize PTS Constraint structure
mesh.Constraints(nNodeStrings,1) = struct('num',[],'type',[],'nodeStr',[],'data',[]);
starting = 1;
ending = find(C<0,1,'first')-1;

for i = 1:nNodeStrings
    
    mesh.Constraints(i).type = 'Open Ocean';
    
    mesh.Constraints(i).num  = -1;
    
    mesh.Constraints(i).nodeStr = C(starting:ending);
    
    %mesh.Constraints(i).xy = [x(C(starting:ending)), y(C(starting:ending))];
    
    C(ending + 1) = abs(C(ending + 1));
    
    starting = ending+2;
    ending = find(C<0,1,'first')-1;
    
end

%------------------------------------------------------------------------------
% Determine node strings
%------------------------------------------------------------------------------
uiStatusBar('Detecting boundary segments...')

% Get Mainland and island boundaries
trep = triangulation(mesh.tri, mesh.p(:,1), mesh.p(:,2)); tf = freeBoundary(trep)';

% Assign (x,y) edge segments values
xBound = x(tf); yBound = y(tf);

% Initialize counter, k
k = 1;

% Initialize edge 1
Poly(k).x(1,1) = xBound(1,1); Poly(k).y(1,1) = yBound(1,1);
Poly(k).x(2,1) = xBound(2,1); Poly(k).y(2,1) = yBound(2,1);

% Initialize node string 1
node(k).Str(1,1) = tf(1,1);
node(k).Str(2,1) = tf(2,1);

ptr = 2; % Pointer to current end point

% Mark recorded points as nan's
xBound(1,1) = nan; yBound(1,1) = nan;
xBound(2,1) = nan; yBound(2,1) = nan;

% Initialize searching logical value
searching = 1;

while searching
    
    % Find the next connecting edge with same end point
    ind = find(Poly(k).x(ptr,1) == xBound & Poly(k).y(ptr,1) == yBound);
    
    % If we're empty, are we done or moving on to the next boundary?
    if isempty(ind) 
        
        % Close current boundary
        Poly(k).x(end+1,1) = Poly(k).x(1,1);        % Polygon x value
        Poly(k).y(end+1,1) = Poly(k).y(1,1);        % Polygon y value
        node(k).Str(end+1,1) = node(k).Str(1,1);    % Node 
        PArea(k) = polyarea(Poly(k).x, Poly(k).y);  %#ok<AGROW> % Area 
        
        % Check to see if we're done
        if all(isnan(xBound)); searching = 0; continue; end
        
        k = k+1; % Start new structure for next boundary
        
        ind = find(~isnan(xBound(1,:)),1,'first'); % Find a point to start
        
        % Initialize edge 1 of polygon k
        Poly(k).x(1,1) = xBound(1,ind); Poly(k).y(1,1) = yBound(1,ind);
        Poly(k).x(2,1) = xBound(2,ind); Poly(k).y(2,1) = yBound(2,ind);
        
        % Initialize node string 1 of polygon k
        node(k).Str(1,1) = tf(1,ind);
        node(k).Str(2,1) = tf(2,ind);
        
        % Mark recorded points as nan's
        xBound(1,ind) = nan; yBound(1,ind) = nan;
        xBound(2,ind) = nan; yBound(2,ind) = nan;
        
        ptr = 2; % Pointer to current end point
        
        continue; % Continue with loop
    end
    
    % convert linear index to subscript
    [r,c] = ind2sub(size(xBound),ind);
    
    if r == 1 % if r == 1, we want to grab the opposite end point
        
        % increment pointer
        ptr = ptr+1;
        
        % Store polygon values
        Poly(k).x(ptr,1) = xBound(2,c); Poly(k).y(ptr,1) = yBound(2,c);
        
        % Store node value
        node(k).Str(ptr,1) = tf(2,c);
        
    else
        
        % increment pointer
        ptr = ptr+1;
        
        % Store polygon values
        Poly(k).x(ptr,1) = xBound(1,c); Poly(k).y(ptr,1) = yBound(1,c);
        
        % Store node value
        node(k).Str(ptr,1) = tf(1,c);
        
    end
    
    
    % Mark recorded points as nan's
    xBound(1,c) = nan; yBound(1,c) = nan;
    xBound(2,c) = nan; yBound(2,c) = nan;
    
end

clear Poly

% Order boundaries based on area (largest area considered mainland)
[~,I] = sort(PArea,'descend');

ns = nNodeStrings;

p = 1; % Counter

for k = 1:numel(I)
    
    % Remove open ocean boundaries from mainland node string
    if k == 1 && ns ~= 0
        
        % Assign mainland boundary node string
        string = node(k).Str(1:end-1); % Not supposes to be closed
                
        % Remove open ocean boundary nodes 
        for i = 1:ns
            
            % Keep end nodes for connectivity
            ix = ismember(string,mesh.Constraints(i).nodeStr(2:end-1)); 

            % Seperate by a single nan
            string(ix) = nan;
            
        end
        
        % Extract segments from NaN seperated vector
        pos = [true, isnan(string(:, 1)).', true];
        ini = strfind(pos, [true, false]);
        fin = strfind(pos, [false, true]) - 1;
        C   = cell(1, length(ini));
        for iC = 1:length(ini)
            C{iC} = string(ini(iC):fin(iC), :);
        end
        
        for i = 1:numel(C)
            
            % Add mesh constraint attributes for external boundary
            mesh.Constraints(p+ns).type = 'External Boundary';
            mesh.Constraints(p+ns).num  = 0;
            mesh.Constraints(p+ns).nodeStr = C{i};
            
            %is = ix(i)+1;
            
            p = p + 1;
            
        end
        
    elseif k == 1
        
        mesh.Constraints(p+ns).type = 'External Boundary';
        mesh.Constraints(p+ns).num  = 0;
        mesh.Constraints(p+ns).nodeStr = node(k).Str;
        
        p = p + 1;
        
    else
        
        mesh.Constraints(p+ns).type = 'Internal Boundary';
        mesh.Constraints(p+ns).num  = 1;
        mesh.Constraints(p+ns).nodeStr = node(k).Str;
        
        p = p + 1;
        
        
    end

end

if all(mesh.p(:,3) == 1) || all(mesh.p(:,3) == 0)
    mesh.p = mesh.p(:,[1 2]);
end

%------------------------------------------------------------------------------
% Read bathymetry data, if available
%------------------------------------------------------------------------------

if strcmp(bathymetry_included, 'Yes')
    
    fid=fopen([bathypathname,bathyfilename], 'r');  % Open file
    
    % Check to make sure file is still there
    if(fid < 0)
        
        warndlg('Did you move this file? It''s no longer there','Error');
        
        uiStatusBar('Ready')

        mesh = [];
        return
    end
    
    %statusbar(guiFig, 'Check consistency...'); pause(.0001)
    
    nVerts_bathy = textscan(fid, '%*s %f %*[^\n]', 1, 'Headerlines', 3);
    nElems_bathy = textscan(fid, '%*s %f %*[^\n]', 1);
    
    if( nVerts_bathy{1}-nVerts ~= 0 || nElems_bathy{1}-nElems ~= 0 )
        
        warndlg('The dimensions of bathymetry file do not match the 2dm file.','Error');
        
        uiStatusBar('Ready')

        mesh = [];
        return
    end
    
    uiStatusBar('Reading in bathymetry data...')
    
    % Read in bathymetry
    try
        C = textscan(fid, '%f', nVerts_bathy{1},'Headerlines', 4);
        
        fclose(fid); % Close file
        
        %------------------------------------------------------------------------------
        % Generate scattered interpolant function for bathymetry
        %------------------------------------------------------------------------------
        uiStatusBar('Creating an interpolant function for the elevation data set...')
        
        %z = zeros(length(x),1);
        z(:,1) = C{1};

        mesh.p = [mesh.p z];
        
    catch
        
        warndlg('There is a problem reading this file.','Error');
        uiStatusBar('Ready')
        mesh = [];
        return

    end
       
end

end