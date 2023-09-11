function [h_lfs] = medial_distance_FMM(MAP,X,Y,D,R)
% medial_distance_FMM - Medial Distance Fast Marching Method
%
% Syntax:  [h_lfs] = medial_distance_FMM(MAP,X,Y,D,R)
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

% Medial Axis Distance Subprogram
debug = 0;
deltax  = X(1,2)-X(1,1);
deltay  = Y(2,1)-Y(1,1);
LX = length(X(1,:));
LY = length(Y(:,1));
MAD = zeros(size(X));
% Set up heap structures
[Accept,Trial] = heap(MAP,X,Y);
% calculate u values for nodes in Trial
LT = length(Trial(:,1,1));
for k = 1:LT % Initial Trial values will only be delta(x/y) away from MA
    if Trial(k,2,2) == 1
        Trial(k,2,2) = deltax; % replace backpointer w/value for u
    elseif Trial(k,2,2) == -1
        Trial(k,2,2) = deltax;
    elseif Trial(k,2,2) == 2
        Trial(k,2,2) = deltay;
    elseif Trial(k,2,2) == -2
        Trial(k,2,2) = deltay;
    end
end
% Main Loop
if debug == 1
    keyboard
end
i = 0;
while min(Trial(:,2,2)) ~= inf
    i = i+1;
    Trial = min_sort(Trial);
    jj = Trial(1,1,1); ii = Trial(1,2,1);
    MAD(jj,ii) = Trial(1,2,2);
    Accept(end+1,:,:) = Trial(1,:,:);
    Trial(1,2,2) = inf;
    % Look at neighbors
    kk = Trial(1,1,2);
    % right neighbor
    duplicate1 = find(Accept(:,1,2) == kk+LY);
    duplicate2 = find(Trial(:,1,2) == kk+LY);
    if numel(duplicate1) == 0
        if numel(duplicate2) == 0
            last_node = LX*LY;
            if kk+LY > last_node
                at_boundary = 1; %#ok<NASGU>
            else
                Trial(end+1,1,1) = jj;
                Trial(end,2,1) = ii+1;
                Trial(end,1,2) = kk+LY;
                Trial(end,2,2) = deltax+MAD(jj,ii);
            end
        elseif numel(duplicate2) ~= 0
            u_new = deltax+MAD(jj,ii);
            if u_new < Trial(duplicate2,2,2)
                Trial(duplicate2,2,2) = u_new;
            end
        end
    end
    % left neighbor
    duplicate1 = find(Accept(:,1,2) == kk-LY);
    duplicate2 = find(Trial(:,1,2) == kk-LY);
    if numel(duplicate1) == 0
        if numel(duplicate2) == 0
            if ii-1 ~= 0
                Trial(end+1,1,1) = jj;
                Trial(end,2,1) = ii-1;
                Trial(end,1,2) = kk-LY;
                Trial(end,2,2) = deltax+MAD(jj,ii);
            else
                at_boundary = 2; %#ok<NASGU>
            end
        elseif numel(duplicate2) ~= 0
            u_new = deltax+MAD(jj,ii);
            if u_new < Trial(duplicate2,2,2)
                Trial(duplicate2,2,2) = u_new;
            end
        end
    end
    % Top neighbor
    duplicate1 = find(Accept(:,1,2) == kk+1);
    duplicate2 = find(Trial(:,1,2) == kk+1);
    if numel(duplicate1) == 0
        if numel(duplicate2) == 0
            if jj+1 > LY
                at_boundary = 3; %#ok<NASGU>
            else
                Trial(end+1,1,1) = jj+1;
                Trial(end,2,1) = ii;
                Trial(end,1,2) = kk+1;
                Trial(end,2,2) = deltax+MAD(jj,ii);
            end
        elseif numel(duplicate2) ~= 0
            u_new = deltax+MAD(jj,ii);
            if u_new < Trial(duplicate2,2,2)
                Trial(duplicate2,2,2) = u_new;
            end
        end
    end
    % Bottom neighbor
    duplicate1 = find(Accept(:,1,2) == kk-1);
    duplicate2 = find(Trial(:,1,2) == kk-1);
    if numel(duplicate1) == 0
        if numel(duplicate2) == 0
            if jj-1 ~= 0
                Trial(end+1,1,1) = jj-1;
                Trial(end,2,1) = ii;
                Trial(end,1,2) = kk-1;
                Trial(end,2,2) = deltax+MAD(jj,ii);
            else
                at_boundary = 4; %#ok<NASGU>
            end
        elseif numel(duplicate2) ~= 0
            u_new = deltax+MAD(jj,ii);
            if u_new < Trial(duplicate2,2,2)
                Trial(duplicate2,2,2) = u_new;
            end
        end
    end
    
end
% Local feature size
LFS = abs(D) + abs(MAD);
h_lfs = LFS./R;
end