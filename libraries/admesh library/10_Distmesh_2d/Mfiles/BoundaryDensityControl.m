function p = BoundaryDensityControl(p,t,C)
% BoundaryDensityContorl - Remove points too close to boundary
%
% Syntax:  p = BoundaryDensityContorl(p,t,MeshFun)
%
% Inputs:
%    p          - Point List
%    t          - Connectivity
%    MeshFun    - Mesh Size Function
%
% Outputs:
%    p - Modified point list for distmesh2d
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

%warning('on','MATLAB:triangulation:PtsNotInTriWarnId')
% Get free boundary facets of the triangulation
ff = freeBoundary(triangulation(t, p));

% Get trangles attached to specified edge
ti = cell2mat(edgeAttachments(trep,ff)')';

% Seperate interior node 
iNodes = sum([t(ti,:) -[ff zeros(size(ff,1),1)]],2);

% Calculate element edge lengths
a = sqrt((p(t(ti,2),1)-p(t(ti,1),1)).^2+(p(t(ti,2),2)-p(t(ti,1),2)).^2);
b = sqrt((p(t(ti,3),1)-p(t(ti,2),1)).^2+(p(t(ti,3),2)-p(t(ti,2),2)).^2);
c = sqrt((p(t(ti,1),1)-p(t(ti,3),1)).^2+(p(t(ti,1),2)-p(t(ti,3),2)).^2);

% Compute mesh quality q for each boundary element
q = ((b+c-a).*(c+a-b).*(a+b-c))./(a.*b.*c);

% Find elements with quality less than tolerance (.3). Project interior
% node to boundary. 
badQ = iNodes(q < .2);

if isempty(badQ); return; end

% Make sure we don't remove points we're constraining
if ~isempty(C)
    badQ = setdiff(badQ,C(:));
end

% Remove node
p(badQ,:) = [];

end