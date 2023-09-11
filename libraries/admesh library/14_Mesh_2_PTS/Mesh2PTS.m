function PTS = Mesh2PTS(MESH)
% Mesh2PTS - Reads in domain and it's attributes
%
% Syntax:  [PTS,xyzFun,xyzConn] = Mesh2PTS(mesh)
%
% Inputs:
%    mesh
%
% Outputs:
%
% See also: Open14Callback
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% October 2013; Last revision: 21-October-2013

%------------------------------- BEGIN CODE -----------------------------------

%------------------------------------------------------------------------------
% Create PTS polygon structure
%------------------------------------------------------------------------------
if ~isempty(MESH.Constraints)
    PTS = Constraints2PTS(MESH);
else
    PTS = triangulation2PTS(MESH);
end

%------------------------------------------------------------------------------
% Create scattered interpolant data set
%------------------------------------------------------------------------------
% uiStatusBar('Creating elevation interpolant function...')
% 
% if nargout == 2
%     
%     if all(MESH.Points(:,3) == 1) || all(MESH.Points(:,3) == 0)
%         
%         xyzFun = [];
%     
%     else
%         
%         uiStatusBar(1/5)
%         
%         xmin = min(vertcat(PTS.Poly(:).x));
%         xmax = max(vertcat(PTS.Poly(:).x));
%         ymin = min(vertcat(PTS.Poly(:).y));
%         ymax = max(vertcat(PTS.Poly(:).y));
%         
%         Offsetx = (xmax-xmin)*.05;
%         Offsety = (ymax-ymin)*.05;
%         
%         % Generate unique x and y grid vectors
%         xRange = unique(MESH.Points(:,1)');
%         yRange = unique(MESH.Points(:,2)');
%         
%         % Produce a grid that won't run out of memory
%         maxVecLength = 1000;
%         
%         % Space vectors accordingly
%         if numel(xRange) > maxVecLength
%             spacing = floor(numel(xRange)/maxVecLength);
%             xRange = xRange(1:spacing:end);
%         end
%         
%         if numel(yRange) > maxVecLength
%             spacing = floor(numel(yRange)/maxVecLength);
%             yRange = yRange(1:spacing:end);
%         end
%         
%         xRange = [xmin-Offsetx, xRange, xmax+Offsetx];
%         yRange = [ymin-Offsety, yRange, ymax+Offsety];
%         
%         % Create grid
%         [xq,yq] = meshgrid(xRange,yRange);
%         
%         uiStatusBar(2/5)
%         % Interpolate scatter set
%         zq = griddata(MESH.Points(:,1),MESH.Points(:,2),MESH.Points(:,3),xq,yq);
%         
%         uiStatusBar(3/5)
%         % All elevation points outside domain are NaN's
%         IN = PointsInDomain(xq,yq,PTS);
%         zq(~IN) = nan;
%         
%         uiStatusBar(4/5)
%         % Interpolate NaN values via inpaint method
%         zq = inpaint_nans(zq,4);
%         
%         % Create gridded interpolant
%         xyzFun = griddedInterpolant(xq',yq',zq','linear','nearest');
%         
%         uiStatusBar(5/5)
% 
%     end
%     
% end
% 
% if nargout == 1
%     
%     varargout{1} = PTS;
%     
% else
%     
%     varargout{1} = PTS;
%     varargout{2} = xyzFun;
%     
% end

end