function [sFB,nodestr] = compileNodeStrings(sb,fb)
% compileNodeStrings - Takes edges of sub domain (sb) and full domain
% (fb) and determines (1) the interior node strings for constraining in
% admesh (nodestr) and (2) the sub domain boundary node string.
%
% Syntax:  [sFB,nodestr] = compileNodeStrings(sb,fb
%
% Inputs:
%
% Outputs:
%
% See also: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% October 2013; Last revision: 21-October-2013

%------------------------------- BEGIN CODE -----------------------------------

% Transpose inputs
sb      = sb';    % Sub domain free boundary edges
tedges  = sb ;    % Copy of sub domain
fb      = fb';    % Full domain free boundary edges

% Order sub domain polygon
j = 1; sFB(j).str = sb(:,1); % Initialize starting point

% Remove edge from search
tedges(:,1) = nan;

% Initalize pointer to node string
ptr = 2;

while 1

    % Find the next connecting edge with same end point
    ind = find(sFB(j).str(ptr) == tedges,1,'first');
    
    if isempty(ind) % If we're empty, are we done or moving on to the next boundary?
        
        % Check to see if we're done
        if all(isnan(tedges)); break; end
        
        j = j+1; % Start new structure for next boundary
        
        ind = find(~isnan(tedges(1,:)),1,'first'); % Find a point to start
        
        % Assign first edge
        sFB(j).str = sb(:,ind);
        
        % Mark recorded points as nan's
        tedges(:,1) = nan;
        
        ptr = 2; % Pointer to current end point
        
        continue; % Continue with loop
    end
    
    % convert linear index to subscript
    [r,c] = ind2sub(size(sb),ind);
    
    if r == 1 % if r == 1, we want to grab the opposite end point
        
        ptr = ptr+1;
        sFB(j).str(ptr) = sb(2,c);
        
    else
        
        ptr = ptr+1;
        sFB(j).str(ptr) = sb(1,c);
        
    end
    
    % Mark recorded points as nan's
    tedges(:,c) = nan;
    
end

% Find edges of interior node string by keeping edges in sb that are not in
% fb. 
Ec = setdiff(sb',fb','rows')';% Interior element edges to constrain

nodestr(length(sFB)).str = []; % Initialize

for k = 1:length(sFB)
    
    % find all nodes in sFB(k).str that have Ec
    ni = ismember(sFB(k).str, Ec(:));
    
    % Create a column vector representing indices in x & y
    ix = (1:length(sFB(k).str))';
    
    % Set indices that do not below to zero
    ix(ni == 0) = 0;
    
    % Shift indices so we can extract one column vector of points
    shiftN = abs(find(ix == 0,1,'last') - length(ix));
    
    if isempty(shiftN); shiftN = 0; end
    
    % Shift
    ix = circshift(ix,shiftN);
    
    % Remove zeros
    ix(ix == 0) = [];

    % Assign node string
    nodestr(k).str = sFB(k).str(ix);

end



end