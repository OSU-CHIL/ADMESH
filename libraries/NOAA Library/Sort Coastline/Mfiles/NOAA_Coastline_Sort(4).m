function [PTS] = NOAA_Coastline_Sort(new_coast)
% NOAA_Coastline_Sort - Sorts out coastline data and puts into PTS data
% structure
%
% Syntax:  [PTS] = NOAA_Coastline_Sort(new_coast)
%
% Inputs:  
%    new_coast - matrix of coastline data
%
% Outputs:
%    PTS - Data structure with fields x & y
%           PTS(1).x = x-coordinates of first polygon
%           PTS(1).y = y-coordinates of first polygon
%
% Other m-files required: none
% Subfunctions: 
% MAT-files required: none

% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% August 2013; Last revision: 08-August-2013

%------------- BEGIN CODE --------------

%Check for double points in new_coast before sorting into edge structure
new_coast = unique(new_coast,'rows','stable');

iseg=find(isnan(new_coast(:,1)));

% initialize PTS
PTS.Poly = repmat(struct('x',1,'y',1), length(iseg)-1, 1 );

for i = 1:length(iseg)-1
    
    if i == 1
        PTS.Poly(i).x = new_coast(iseg(i)+1:iseg(i+1)-1,1);
        PTS.Poly(i).y = new_coast(iseg(i)+1:iseg(i+1)-1,2);
    else
        PTS.Poly(i).x = [new_coast(iseg(i)+1:iseg(i+1)-1,1);new_coast(iseg(i)+1,1)];
        PTS.Poly(i).y = [new_coast(iseg(i)+1:iseg(i+1)-1,2);new_coast(iseg(i)+1,2)];
    end
end

end