function t = BoundaryCleanUp(p,t,C)
% BoundaryCleanUp - Remove elements with elements quality below .1 around
% the boundary
%
% Syntax:  t = BoundaryCleanUp(p,t,C)
%
% Inputs:
%    p          - Point List
%    t          - Connectivity
%    C          - Constraints matrix
%
% Outputs:
%    t          - Modified connectivity list 
%
% Other m-files required: distmesh2d
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% Last revision: 25-April-2014

%---------------------------------------------------------------------
% Begin Code
%---------------------------------------------------------------------

% Turn off warning ID for triangulation. We get the warning because 
% there are points in p that are not referenced by t.
warning('off','MATLAB:triangulation:PtsNotInTriWarnId')

% Get triangulation representation
trep = triangulation(t, p);

% Get free boundary facets of the triangulation
ff = freeBoundary(trep);

% Get triangles attached to specified edge
ti = cell2mat(edgeAttachments(trep,ff)')';

% Calculate element edge lengths
a = sqrt((p(t(ti,2),1)-p(t(ti,1),1)).^2+(p(t(ti,2),2)-p(t(ti,1),2)).^2);
b = sqrt((p(t(ti,3),1)-p(t(ti,2),1)).^2+(p(t(ti,3),2)-p(t(ti,2),2)).^2);
c = sqrt((p(t(ti,1),1)-p(t(ti,3),1)).^2+(p(t(ti,1),2)-p(t(ti,3),2)).^2);

% Compute mesh quality q for each boundary element
q = ((b+c-a).*(c+a-b).*(a+b-c))./(a.*b.*c);

% Find elements with quality less than tolerance (.1). 
badT = ti(q < .15);

% We want to make sure we don't remove elements with constrained edges
if ~isempty(C)
    
    % Get triangles attached to constrainted edges
    tc = cell2mat(edgeAttachments(trep,C)')';
    
    % Find common members and remove
    badT(ismember(badT,tc)) = [];
    
end

% Remove bad boundary elements
t(badT,:) = [];

end