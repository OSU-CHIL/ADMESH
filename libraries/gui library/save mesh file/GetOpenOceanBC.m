function [NOPE,NETA,NVDLL,NBDV] = GetOpenOceanBC(t,p,PTS)

%figure
%triplot(t, p(:,1), p(:,2))

%------------------------------------------------------------------------------
% If there are no open ocean boundaries, return
%------------------------------------------------------------------------------
if ~any(ismember([PTS.Constraints.num],-1))
    NOPE = nan; NETA = nan; NVDLL = nan; NBDV = nan;
    return
end

%------------------------------------------------------------------------------
% Find which cells contain open ocean boundaries
%------------------------------------------------------------------------------
ix = find(ismember([PTS.Constraints.num],-1));

% Assign NOPE
NOPE = length(ix);

% Initialize NETA
NETA = 0;

% Initialize NVDLL
NVDLL = zeros(size(ix));

NBDV = cell(size(ix));

%------------------------------------------------------------------------------
% For each node string in IBtype, define the new node string in the new mesh
%------------------------------------------------------------------------------
for i = 1:length(ix)
    
    drawnow; % for graphics
    
    % Find starting node & ending node in the exisiting IBtype
    Snode = PTS.Constraints(ix(i)).xy(1,:); 
    Enode = PTS.Constraints(ix(i)).xy(end,:);
    
    % Use Snode & Enode to determine the new starting and ending node
    % for the open ocean boundary in p.
    [~,loca] = min(sum([abs(Snode(1,1) - p(:,1)).^2  abs(Snode(1,2) - p(:,2)).^2],2));
    [~,locb] = min(sum([abs(Enode(1,1) - p(:,1)).^2  abs(Enode(1,2) - p(:,2)).^2],2));
    
    %------------------------------------------------------------------------------
    % Determine if original open ocean node list is going cw or ccw
    %------------------------------------------------------------------------------
    
    % Get first edge in IBtype
    edgex = [PTS.Constraints(ix(i)).xy(1,1) PTS.Constraints(ix(i)).xy(2,1)];
    edgey = [PTS.Constraints(ix(i)).xy(1,2) PTS.Constraints(ix(i)).xy(2,2)];
        
    % Compute Mid-point on edge
    Midp = [mean(edgex) , mean(edgey)]; % MID point
    
    % Compute normals
    L = sqrt(diff(edgex).^2 + diff(edgey).^2);
    nx =  diff(edgey)./L;
    ny = -diff(edgex)./L;
    
    % Determine which way to project the mid-point by guessing an
    % initial projection and testing if that projection is inside (1)
    % or outside (0) the domain PTS.
    in = PointInPolygon(Midp(1) - nx*(L/2), Midp(2) - ny*(L/2), PTS.Poly(1).x, PTS.Poly(1).y);
    
    if in; sgn = -1; elseif in == 0; sgn = +1; end
    
    %------------------------------------------------------------------------------
    % Now we know which direction to project the midpoint, use
    % this test on the starting point of the node string in p
    % to see which end point should be the the next node and so on.
    %------------------------------------------------------------------------------
    
    %------------------------------------------------------------------------------
    % Get boundary edge segments of triangulation in the form xb (2xm) and
    % yb (2xm). Each row corresponds the end points of each edge segment.
    %------------------------------------------------------------------------------
    trep = triangulation(t, p(:,1), p(:,2)); % Re-triangulate
    fb = freeBoundary(trep); % Get boundary edge node pairs list
    xb = p(:,1); yb = p(:,2); xb = xb(fb)'; yb = yb(fb)'; % Extract coordinates
    
    % Determine which edge is going the same direction as the original
    
    % Find the starting location in [xb,yb], test the first edge associated
    % with the starting location
    k = find( ((xb(1,:) == p(loca,1)) & (yb(1,:) == p(loca,2))) );
    
    edgex = [xb(1,k) xb(2,k)];
    edgey = [yb(1,k) yb(2,k)];
    
    % Compute Mid-point on edge
    Midp = [mean(edgex) , mean(edgey)]; % MID point
    
    % Compute normals
    L = sqrt(diff(edgex).^2 + diff(edgey).^2);
    nx =  diff(edgey)./L;
    ny = -diff(edgex)./L;
    
    % Test projection obtained from above.
    in = PointInPolygon(Midp(1) + sgn*nx*(L/2), Midp(2) + sgn*ny*(L/2), PTS.Poly(1).x, PTS.Poly(1).y);
    
    % If test fails (in == 0), then choose the other edge in the xb, yb
    % edge list.
    if in
        
        % Passed test
        % Remove points from edge list, xb & yb
        xb(:,k) = nan;
        yb(:,k) = nan;
        
    else
        
        k = find( ((xb(2,:) == p(loca,1)) & (yb(2,:) == p(loca,2))) );
        edgex = [xb(2,k) xb(1,k)];
        edgey = [yb(2,k) yb(1,k)];
        % Remove points from edge list, xb & yb
        xb(:,k) = nan;
        yb(:,k) = nan;
        
    end
    
    % Enter loop, connecting end points until end point is reached
    Looping = 1;
    ptr = 2; % Pointer to current end point
    
    while Looping
        
        % Find the next connecting edge with same end point
        ind = find(edgex(ptr) == xb & edgey(ptr) == yb,1,'first');
        
        % convert linear index to subscript
        [r,c] = ind2sub(size(xb),ind);
        
        if r == 1 % if r == 1, we want to grab the opposite end point
            
            ptr = ptr+1;
            edgex(ptr) = xb(2,c); edgey(ptr) = yb(2,c);
            
            % Mark recorded points as nan's
            xb(:,c) = nan; yb(:,c) = nan;
            
        else
            
            ptr = ptr+1;
            edgex(ptr) = xb(1,c); edgey(ptr) = yb(1,c);
            
            % Mark recorded points as nan's
            xb(:,c) = nan; yb(:,c) = nan;
            
        end
        
        % Mark recorded points as nan's
        %xb(:,c) = nan; yb(:,c) = nan;
        
        % Check if we've reached end point
        if ( edgex(ptr) == p(locb,1) && edgey(ptr) == p(locb,2) )
            Looping = 0;
            continue;
        end
        
        
    end
    
    % Transpose vectors
    edgex = edgex';
    edgey = edgey';
    
    % Get running total for NETA
    NETA = length(edgex) + NETA;
    
    % Store NVDLL
    NVDLL(i) = length(edgex);
    
    % Define nodal points
    for n = 1:length(edgex)
        
        loc = find( edgex(n) == p(:,1) & edgey(n) == p(:,2) );
        
        NBDV{i}(n) = loc;
        
    end

end

end