function fb = dtFreeBoundary(varargin)

% Two methods.

% (1) return free boundary for entire domain
if nargin == 1
    
    MESH = varargin{1};
    
    % Get triangulation representation
    trep = triangulation(MESH.ConnectivityList,MESH.Points(:,[1 2]));
    
    fb = sort(freeBoundary(trep),2);
    
    %[fb,~,c] = unique(...
    %    sort([MESH.ConnectivityList(:,[1 2]); ...
    %    MESH.ConnectivityList(:,[2 3]); ...
    %    MESH.ConnectivityList(:,[3 1])],2),'rows');
    
    %fb = sort(fb(histc(c,unique(c)) == 3,:),2);
    
else % (2) return free boundary from sub domain
    
    % For some reason freeBoundary will not work with sub domains
    
    [MESH,ti] = deal(varargin{:});
    
    [fb,~,c] = unique(...
        sort([MESH.ConnectivityList(ti,[1 2]); ...
        MESH.ConnectivityList(ti,[2 3]); ...
        MESH.ConnectivityList(ti,[3 1])],2),'rows');
    
    ix = histc(c,unique(c)) == 3;
    
    %if all(~ix); ix = histc(c,unique(c)) == 9; end
    
    fb = sort(fb(ix,:),2);
    
end




end