function [Trial] = min_sort(Trial)
% min_sort- Medial Axis Subroutine
%
% Syntax:  [Trial] = min_sort(Trial)
%
% Inputs:
%
% Outputs:
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Colton Conroy
% The Ohio State University
% email address: conroy.51@osu.edu
% August 2013; Last revision: 08-August-2013

%------------- BEGIN CODE --------------
u_min = min(Trial(:,2,2));
k = find(Trial(:,2,2) == u_min,1);
uu = Trial(1,:,:);
Trial(1,:,:) = Trial(k,:,:);
Trial(k,:,:) = uu;
end