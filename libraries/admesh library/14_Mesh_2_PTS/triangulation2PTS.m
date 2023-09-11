function PTS = triangulation2PTS(MESH)

% triangulation2PTS - Extracts mainland and islands from
% triangulation and stores in a data structure
%
% Syntax:  PTS = triangulation2PTS(tri,xyz)
%
% Inputs:
%   tri - connectivity
%   xyz - x,y and z coordinates
%
% Outputs:
%    PTS - Data structure with fields x & y
%    PTS(1).x = x-coordinates of first polygon
%    PTS(1).y = y-coordinates of first polygon
%
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% August 2013; Last revision: 11-October-2013

%------------- BEGIN CODE --------------

% Assign cpplon & cpplat
PTS.cpplon = MESH.cpplon;
PTS.cpplat = MESH.cpplat;

%tri = MESH.ConnectivityList;
x   = MESH.Points(:,1);
y   = MESH.Points(:,2);

%------------------------------------------------------------------------------
% Get Mainland and island boundaries
%------------------------------------------------------------------------------
trep = triangulation(MESH.ConnectivityList, x, y);

fe = freeBoundary(trep)';

% [fe,~,c] = unique(...
%     sort([MESH.ConnectivityList(:,[1 2]); ...
%     MESH.ConnectivityList(:,[2 3]); ...
%     MESH.ConnectivityList(:,[3 1])],2),'rows');
% 
% fe = sort(fe(histc(c,unique(c)) == 3,:),2)'

xBound = x(fe); yBound = y(fe);

%------------------------------------------------------------------------------
% Order boundary points
%------------------------------------------------------------------------------
j = 1;
% initialize edge 1
% Edge 1
Temp_PTS(j).x(1) = xBound(1,1); Temp_PTS(j).y(1) = yBound(1,1);
Temp_PTS(j).x(2) = xBound(2,1); Temp_PTS(j).y(2) = yBound(2,1);

ptr = 2; % Pointer to current end point

% Mark recorded points as nan's
xBound(1,1) = nan; yBound(1,1) = nan;
xBound(2,1) = nan; yBound(2,1) = nan;

Looping = 1;

while Looping
    
    % Find the next connecting edge with same end point
    ind = find(Temp_PTS(j).x(ptr) == xBound & Temp_PTS(j).y(ptr) == yBound,1,'first');
    
    if isempty(ind) % If we're empty, are we done or moving on to the next boundary?
        
        % Close current boundary
        Temp_PTS(j).x(end+1) = Temp_PTS(j).x(1);
        Temp_PTS(j).y(end+1) = Temp_PTS(j).y(1);
        
        % Compute boundary area
        PArea(j) = polyarea(Temp_PTS(j).x, Temp_PTS(j).y);   %#ok<*AGROW>
                
        % Check to see if we're done
        ind = ~isnan(xBound);
        
        if sum(ind(:)) == 0; Looping = 0; continue; end
        
        j = j+1; % Start new structure for next boundary
        
        ind = find(~isnan(xBound(1,:)),1,'first'); % Find a point to start
        
        Temp_PTS(j).x(1) = xBound(1,ind); Temp_PTS(j).y(1) = yBound(1,ind);
        Temp_PTS(j).x(2) = xBound(2,ind); Temp_PTS(j).y(2) = yBound(2,ind);
        
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
        Temp_PTS(j).x(ptr) = xBound(2,c); Temp_PTS(j).y(ptr) = yBound(2,c);
        
    else
        
        ptr = ptr+1;
        Temp_PTS(j).x(ptr) = xBound(1,c); Temp_PTS(j).y(ptr) = yBound(1,c);
        
    end
    
    % Mark recorded points as nan's
    xBound(1,c) = nan; yBound(1,c) = nan;
    xBound(2,c) = nan; yBound(2,c) = nan;
    
end

%------------------------------------------------------------------------------
% Order boundaries based on area
%------------------------------------------------------------------------------
[~,I] = sort(PArea,'descend');

PTS.Poly(numel(I),1) = struct('x',[],'y',[]);

for i = 1:numel(I)
    PTS.Poly(i).x = Temp_PTS(I(i)).x'; PTS.Poly(i).y = Temp_PTS(I(i)).y';
end

PTS.Constraints = [];

clear Temp_PTS

end