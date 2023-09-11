function [Accept,Trial] = heap(MAP,X,Y)
% heap - Medial Axis Subroutine
%
% Syntax:  [Accept,Trial] = heap(MAP,X,Y)
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

debug = '0';
% Initialize variables
Accept = zeros(1,1,1);
k = 0;
kk = 0;
LX = length(X(1,:));
LY = length(Y(:,1));
Trial(1,1,2) = 0;
if debug == '1'
    keyboard
end
% Accept heap
for i = 1:LX
    for j = 1:LY
        k = k+1;
        if MAP(j,i) == 1
            kk = kk + 1;
            Accept(kk,1,1) = j;
            Accept(kk,2,1) = i;
            Accept(kk,1,2) = k; % node #
            Accept(kk,2,2) = 0; % b/c accepted
        end
    end
end

% trial heap
j = 0;
for i = 1:length(Accept(:,1,1))
    % x-dir of back. mesh
    jj = Accept(i,1,1);
    ii = Accept(i,2,1)+1;
    kk = Accept(i,1,2);
    if MAP(jj,ii) ~= 1
        duplicate = find(Trial(:,1,2) == kk+LY);
        if numel(duplicate) == 0
            j = j+1;
            Trial(j,1,1) = jj;
            Trial(j,2,1) = ii;
            Trial(j,1,2) = kk+LY; % node #
            Trial(j,2,2) = 1; % 1 = neighbor to right is not accepted
        end
    end
    clear ii;
    ii = Accept(i,2,1)-1;
    if MAP(jj,ii) ~= 1
        duplicate = find(Trial(:,1,2) == kk-LY);
        if numel(duplicate) == 0
            j = j+1;
            Trial(j,1,1) = jj;
            Trial(j,2,1) = ii;
            Trial(j,1,2) = kk-LY; % node #
            Trial(j,2,2) = -1; % -1 = neighbor to left is not accepted
        end
    end
    % y-dir of back. mesh
    clear jj; clear ii;
    jj = Accept(i,1,1)+1;
    ii = Accept(i,2,1);
    if MAP(jj,ii) ~= 1
        duplicate = find(Trial(:,1,2) == kk+1);
        if numel(duplicate) == 0
            j = j+1;
            Trial(j,1,1) = jj;
            Trial(j,2,1) = ii;
            Trial(j,1,2) = kk+1; % node %
            Trial(j,2,2) = 2; % 2 = neighbor above is not accepted
        end
    end
    clear jj;
    jj = Accept(i,1,1)-1;
    if MAP(jj,ii) ~= 1
        duplicate = find(Trial(:,1,2) == kk-1);
        if numel(duplicate) == 0
            j = j+1;
            Trial(j,1,1) = jj;
            Trial(j,2,1) = ii;
            Trial(j,1,2) = kk-1; % node #
            Trial(j,2,2) = -2; % neighbor below is not accepted
        end
    end
end
end