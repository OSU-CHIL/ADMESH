function PTS = Constraints2PTS(MESH)
% Constraints2PTS - Extracts mainland and islands from
% boundary conditions and stores in a data structure
%
% Syntax:  PTS = Constraints2PTS(PTS)
%
% Inputs:
%   PTS - structure array containing the follwing fields 
%         containing boundary conditions
%   	Constraints(n)  - field containing the following information
%                       num  - boundary type number (ADCIRC numbering system)
%                       type - string description of boundary type
%                       xy   - x & y coordinates of constraint
%                       data - additional data for constraint
%                   
%
% Outputs:
%   PTS - structure array containing the follwing fields 
%         containing boundary conditions
%       Poly(n)         - field containing the following information
%                       x   - x-coordinates of polygon
%                       y   - y-coordinates of first polygon
%   	Constraints(n)  - field containing the following information
%                       num  - boundary type number (ADCIRC numbering system)
%                       type - string description of boundary type
%                       xy   - x & y coordinates of constraint
%                       data - additional data for constraint
%
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% Last revision: 01-March-2013

%---------------------------------BEGIN CODE-----------------------------------

%------------------------------------------------------------------------------
% Assign mesh constraints to PTS data structure
%------------------------------------------------------------------------------
% Number of constraints
nc = numel(MESH.Constraints);

% Copy mesh constraints
PTS.Constraints = MESH.Constraints;

% Remove node string field
PTS.Constraints = rmfield(PTS.Constraints,'nodeStr');

% Make nodeStr xy coordinates
for k = 1:nc
    
    % Handle internal barrier constraints differently
    if any(PTS.Constraints(k).num == [4 5 24 25])
        
        ns = MESH.Constraints(k).nodeStr;
        
        PTS.Constraints(k).xy = [MESH.Points(ns(:,1),[1 2]);
                         flipud( MESH.Points(ns(:,2),[1 2]));
                                 MESH.Points(ns(1,1),[1 2])];
        
    else
        
        ns = MESH.Constraints(k).nodeStr;
        
        PTS.Constraints(k).xy = MESH.Points(ns(:,1),[1 2]);
        
    end
    
end

% Assign cpplon & cpplat
PTS.cpplon = MESH.cpplon;
PTS.cpplat = MESH.cpplat;

clear MESH

%------------------------------------------------------------------------------
% Get External boundary. Find which indices contain external boundaries
%------------------------------------------------------------------------------

% Find boundary indices in BC cell structure (external, open ocean,
% barrier)
I = find(ismember([PTS.Constraints.num],[0 2 10 12 20 22 30 -1 3 13 23]));

% Input all points in a 2xN matrix in xBound and yBound
Cx = cell(1,length(I)); Cy = cell(1,length(I));

for i = 1:length(I)
    Cx{i} = [ PTS.Constraints(I(i)).xy(1:end-1,1)' ;...
        PTS.Constraints(I(i)).xy(2:end,1)' ];
    Cy{i} = [ PTS.Constraints(I(i)).xy(1:end-1,2)' ;...
        PTS.Constraints(I(i)).xy(2:end,2)' ];
end

xBound = cell2mat(Cx); yBound = cell2mat(Cy); clear Cx Cy I
%------------------------------------------------------------------------------
% Order external boundary points
%------------------------------------------------------------------------------

% Initialize edge 1
PTS.Poly(1).x(1,1) = xBound(1,1); PTS.Poly(1).y(1,1) = yBound(1,1);
PTS.Poly(1).x(2,1) = xBound(2,1); PTS.Poly(1).y(2,1) = yBound(2,1);

ptr = 2; % Pointer to current end point

% Mark recorded points as nan's
xBound(1,1) = nan; yBound(1,1) = nan;
xBound(2,1) = nan; yBound(2,1) = nan;

Looping = 1;

while Looping
    
    % Find the next connecting edge with same end point
    ind = find(PTS.Poly(1).x(ptr) == xBound & PTS.Poly(1).y(ptr) == yBound);
    
    if isempty(ind) % If we're empty, we are done!?
        
        % Remove double points
        xy = unique([PTS.Poly(1).x PTS.Poly(1).y],'rows','stable');
        
        PTS.Poly(1).x = xy(:,1);
        PTS.Poly(1).y = xy(:,2);
        
        % Close current boundary
        PTS.Poly(1).x(end+1) = PTS.Poly(1).x(1);
        PTS.Poly(1).y(end+1) = PTS.Poly(1).y(1);
        
        Looping = 0; continue;
        
    end
    
    % convert linear index to subscript
    [r,c] = ind2sub(size(xBound),ind);
    
    if r == 1 % if r == 1, we want to grab the opposite end point
        ptr = ptr+1;
        PTS.Poly(1,1).x(ptr) = xBound(2,c); PTS.Poly(1,1).y(ptr) = yBound(2,c);
    else
        ptr = ptr+1;
        PTS.Poly(1,1).x(ptr) = xBound(1,c); PTS.Poly(1,1).y(ptr) = yBound(1,c);
    end
    
    % Mark recorded points as nan's
    xBound(1,c) = nan; yBound(1,c) = nan;
    xBound(2,c) = nan; yBound(2,c) = nan;
    
end

%------------------------------------------------------------------------------
% Remove External Boundary Data from Constraints
%------------------------------------------------------------------------------
I = ismember([PTS.Constraints.num],[0 2 10 12 20 22 30]);
PTS.Constraints(I) = [];

%------------------------------------------------------------------------------
% Get internal boundary data. Find which cells contain internal boundaries
%------------------------------------------------------------------------------

% Find boundary indices in BC cell structure
%I = find(ismember([PTS.Constraints.num],[1 11 21 ]));
I = find(ismember([PTS.Constraints.num],[1 11 21 4 5 24 25]));

j = 2; k = 1;

% Input all points in a 2xN matrix in xBound and yBound
Cx = cell(1,length(I)); Cy = cell(1,length(I));

for i = 1:length(I)
    
    % Check if boundary closes in on itself
    if isequal(PTS.Constraints(I(i)).xy(1,1:2),PTS.Constraints(I(i)).xy(end,1:2))
        
        xy = unique(PTS.Constraints(I(i)).xy,'rows','stable');
        
        PTS.Poly(j).x = [xy(:,1); xy(1,1)];
        PTS.Poly(j).y = [xy(:,2); xy(1,2)];
        
        j = j+1;
    else
        
        Cx{i} = [ PTS.Constraints(I(i)).xy(1:end-1,1)' ;...
            PTS.Constraints(I(i)).xy(2:end,1)' ];
        Cy{i} = [ PTS.Constraints(I(i)).xy(1:end-1,2)' ;...
            PTS.Constraints(I(i)).xy(2:end,2)' ];
        k = k+1;
    end
    
end

% Check to see if there are any internal boundaries or barriers that did
% not close
if all(cellfun(@isempty,Cx))
    return
else
    %Cx(~cellfun(@isempty,Cx)) = [];
end

xBound = cell2mat(Cx);
yBound = cell2mat(Cy); clear Cx Cy 

%------------------------------------------------------------------------------
% Order internal boundary points
%------------------------------------------------------------------------------

% Initialize edge 1
PTS.Poly(j).x(1,1) = xBound(1,1); PTS.Poly(j).y(1,1) = yBound(1,1);
PTS.Poly(j).x(2,1) = xBound(2,1); PTS.Poly(j).y(2,1) = yBound(2,1);

ptr = 2; % Pointer to current end point

% Mark recorded points as nan's
xBound(1,1) = nan; yBound(1,1) = nan;
xBound(2,1) = nan; yBound(2,1) = nan;

Looping = 1;

while Looping
    
    % Find the next connecting edge with same end point
    ind = find(PTS.Poly(j).x(ptr,1) == xBound & PTS.Poly(j).y(ptr,1) == yBound);
    
    if isempty(ind) % If we're empty, are we done or moving on to the next boundary?
        
        % Remove double points
        xy = unique([PTS.Poly(j).x PTS.Poly(j).y],'rows','stable');
        
        PTS.Poly(j).x = xy(:,1);
        PTS.Poly(j).y = xy(:,2);
        
        % Close current boundary
        PTS.Poly(j).x(end+1,1) = PTS.Poly(j).x(1,1);
        PTS.Poly(j).y(end+1,1) = PTS.Poly(j).y(1,1);
        
        % Check to see if we're done
        ind = ~isnan(xBound);
        
        if sum(ind(:)) == 0; Looping = 0; continue; end
        
        j = j+1; % Start new structure for next boundary
        
        ind = find(~isnan(xBound(1,:)),1,'first'); % Find a point to start
        
        PTS.Poly(j).x(1,1) = xBound(1,ind); PTS.Poly(j).y(1,1) = yBound(1,ind);
        PTS.Poly(j).x(2,1) = xBound(2,ind); PTS.Poly(j).y(2,1) = yBound(2,ind);
        
        % Mark recorded points as nan's
        xBound(1,ind) = nan; yBound(1,ind) = nan;
        xBound(2,ind) = nan; yBound(2,ind) = nan;
        
        ptr = 2; % Pointer to current end point
        
        continue; % Continue with loop
    end
    
    % convert linear index to subscript
    [r,c] = ind2sub(size(xBound),ind);
    
    if r == 1 % if r == 1, we want to grab the opposite end point
        
        ptr = ptr+1;
        PTS.Poly(j).x(ptr,1) = xBound(2,c); PTS.Poly(j).y(ptr,1) = yBound(2,c);
        
    else
        
        ptr = ptr+1;
        PTS.Poly(j).x(ptr,1) = xBound(1,c); PTS.Poly(j).y(ptr,1) = yBound(1,c);
        
    end
    
    % Mark recorded points as nan's
    xBound(1,c) = nan; yBound(1,c) = nan;
    xBound(2,c) = nan; yBound(2,c) = nan;
    
end

%------------------------------------------------------------------------------
% Remove Internal Boundary Data from Constraints
%------------------------------------------------------------------------------
PTS.Constraints(I) = [];

end
