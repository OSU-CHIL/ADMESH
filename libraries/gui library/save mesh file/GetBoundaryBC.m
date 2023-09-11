function varargout = GetBoundaryBC(plist,p,PTS,C,Ctype,Cdata)

%------------------------------------------------------------------------------
% Check for contraints
%------------------------------------------------------------------------------
if isempty(Ctype)
    
    NBOU = 0;
    NVEL    = 0;
    NVELL   = [];
    IBTYPE  = [];
    DATA    = [];
    c = 1;
    
else
    
    %------------------------------------------------------------------------------
    % Constraints exist, let's store these guys for writing out to our file
    %------------------------------------------------------------------------------
    nsegments = 0;              % number of actual contraint segments
    Cc = cell(size(Ctype));     % Cell of constraints
    ic = find(isnan(C(:,1)));   % Locations of NaN seperations
    s = 1;                      % Starting index
    ir = [];                    % Index in Ctype to remove
    for i = 1:length(Ctype)
        
        if any([Ctype{i}] == [3 13 23 18]) % Interior Barrier & Channel line
            
            nsegments = nsegments + 1;
            
        elseif any([Ctype{i}] == [4  24 5 25]) % Internal Barrier
            
            % There are actual 2 segments back to back that belong to one
            % They are combined for purpose of ADMESH
            nsegments = nsegments + .5;
            
            ir = [ir i]; 
            
        end
        
        Cc{i} = unique(C(s:ic(i)-1,:),'stable'); s = ic(i) +1; % Store constraint in cell
        
    end
    
    % Remove extra IBtype values for internal barriers
    ir(2:2:end) = []; Ctype(ir) = [];
    
    % Initialize variables
    NBOU    = nsegments;
    NVEL    = 0;
    NVELL   = cell(nsegments,1);
    IBTYPE  = cell(nsegments,1);
    DATA    = cell(nsegments,1);
    c = 1; % Index for variables above
    n = 1; % Index for Cc
    
    for i = 1:nsegments
        
        %------------------------------------------------------------------------------
        % External Barrier
        %------------------------------------------------------------------------------
        if any([Ctype{i}] == [3 13 23])
            
            % Update total number of nodes in all boundary segments
            NVEL = length(Cc{n}) + NVEL;
            
            % Update # of nodes in each segment
            NVELL{c} = length(Cc{n});
            
            % Update type of boundary
            IBTYPE{c} = Ctype{i};
            
            % Save data for boundary segment
            DATA{c} = [Cc{n} Cdata{n}];
            
            % Remove open boundary segment from plist
            plist = RemoveFromPlist(p,plist,Cc{n});
            
            %------------------------------------------------------------------------------
            % Channel Line
            %------------------------------------------------------------------------------
        elseif any([Ctype{i}] == 18)
            
            % Update total number of nodes in all boundary segments
            NVEL = length(Cc{n}) + NVEL;
            
            % Update # of nodes in each segment
            NVELL{c} = length(Cc{n});
            
            % Update type of boundary
            IBTYPE{c} = Ctype{i};

            % Save data for boundary segment
            DATA{c} = [Cc{n} Cdata{n}];
            
            % Remove open boundary segment from plist
            plist = RemoveFromPlist(p,plist,Cc{n});
            
            %--------------------------------------------------------------------------
            % % Internal barrier
            %--------------------------------------------------------------------------
        elseif any([Ctype{i}] == [4 24 5 25])
            
            
            % Update total number of nodes in all boundary segments
            NVEL = length(Cc{n}(:,1))*2 + NVEL;
            
            % Update # of nodes in each segment
            NVELL{c} = length(Cc{n}(:,1));
            
            % Update type of boundary
            IBTYPE{c} = Ctype{i};
            
            DATA{c} = [Cc{n} Cc{n+1} Cdata{n}];
            
            % Remove open boundary segment from plist
            plist = RemoveFromPlist(p,plist,Cc{n});
            
            plist = RemoveFromPlist(p,plist,Cc{n+1});
            
            % If internal barrier end points (end edges) touch the external or
            % internal boundary then include these edge points as such.
            
            % Edge 1
            xi = [p(Cc{n}(1),1); p(Cc{n+1}(1),1)];
            yi = [p(Cc{n}(1),2); p(Cc{n+1}(1),2)];
            
            d = zeros(2,length(PTS.Poly));
            
            % Compute distance to main boundary
            for k = 1:length(PTS.Poly)
                d(:,k) = Point2EdgeDistance(xi,yi,PTS.Poly(k).x,PTS.Poly(k).y);
            end
            
            % Store the minimum distance value and the index in PTS
            val = min(min(d));
            
            if val <= eps % Edge is on a boundary
                
                % Append connecting edges of internal barriers to plist.
                plist = [plist ;...
                    [nan nan]; ...
                    [xi yi] ;...
                    [nan nan]]; %#ok<*AGROW> % 1st edge
                
            else
                
                % Remove edge from plist
                for q = 1:2
                    
                    id2 = (xi(q) == plist(:,1) & yi(q) == plist(:,2));
                    
                    plist(id2,:) = nan;
                    
                end
                
            end
            
            % Edge 2
            xi = [p(Cc{n}(end),1); p(Cc{n+1}(end),1)];
            yi = [p(Cc{n}(end),2); p(Cc{n+1}(end),2)];
            
            d = zeros(2,length(PTS.Poly));
            
            % Compute distance to main boundary
            for k = 1:length(PTS)
                d(:,k) = Point2EdgeDistance(xi,yi,PTS.Poly(k).x,PTS.Poly(k).y);
            end
            
            % Store the minimum distance value and the index in PTS
            val = min(min(d));
            
            if val <= eps % Edge is on a boundary
                
                % Append connecting edges of internal barriers to plist.
                plist = [plist ;...
                    [nan nan]; ...
                    [xi yi] ;...
                    [nan nan]]; %#ok<*AGROW> % 1st edge
                
            else
                
                % Remove edge from plist
                for q = 1:2
                    
                    id2 = (xi(q) == plist(:,1) & yi(q) == plist(:,2));
                    
                    plist(id2,:) = nan;
                    
                end
                
            end
            
            
            
            n = n+1;
            
        end
        
        % Update counter
        c = c+1;
        
        n = n+1;
        
    end
    
end

%------------------------------------------------------------------------------
% Get internal and external boundary info
%------------------------------------------------------------------------------

% Join all edge segments possible
[plist]=join_cst(plist,eps);

% Compile segements into cell and designate a external or internal boundary
% condition
i1 = all(~isnan(plist),2); i2 = i1(:)';
idx = [strfind([~i2(1),i2],[0 1]); strfind([i2, ~i2(end)],[1 0])];
seg = mat2cell(plist(i1,:),diff(idx)+1,size(plist,2));

% Update output variables
NBOU = length(seg) + NBOU;
DATA    = [DATA   ; cell(length(seg),1)];
IBTYPE  = [IBTYPE ; cell(length(seg),1)];
NVELL   = [NVELL  ; cell(length(seg),1)];

% Go through each segment. Determine if segment is internal or external
for i = 1:length(seg)
    
    % Compute distance to main boundary
    d = Point2EdgeDistance(seg{i}(:,1),seg{i}(:,2),PTS.Poly(1).x,PTS.Poly(1).y);
    
    val = min(d);
    
    if val <= eps % External boundary
        
        IBTYPE{c} = 0;
        
        % Remove double points
        if ( seg{i}(1,1) ==  seg{i}(end,1) && seg{i}(1,2) == seg{i}(end,2) )
            
            seg{i} = [unique(seg{i},'rows','stable') ; seg{i}(1,:)];
        else
            seg{i} = unique(seg{i},'rows','stable');
            
        end
        
        
    else % Maybe an Internal boundary
        
        if ( seg{i}(1,1) ==  seg{i}(end,1) && seg{i}(1,2) == seg{i}(end,2) )
            
            IBTYPE{c} = 1;
            
            % Remove double points
            seg{i} = [unique(seg{i},'rows','stable') ; seg{i}(1,:)];
            
        else
            
            IBTYPE{c} = 0;
            
            % Remove double points
            seg{i} = unique(seg{i},'rows','stable');
            
        end
        
    end
    
    NVEL = length(seg{i}(:,1)) + NVEL;
    
    NVELL{c} = length(seg{i}(:,1));
    
    % Locate nodal values
    for k = 1:length(seg{i}(:,1))
        
        DATA{c}(k,1) = find((seg{i}(k,1) == p(:,1) & seg{i}(k,2) == p(:,2)));
        
    end
    
    
    c = c + 1;
    
end

if NBOU ~= 0 && NVEL ~= 0
    
    % Update outputs
    varargout{1} = NBOU;
    varargout{2} = NVEL;
    varargout{3} = NVELL;
    varargout{4} = IBTYPE;
    varargout{5} = DATA;
    
else
    
    % Update outputs
    varargout{1} = nan;
    varargout{2} = nan;
    varargout{3} = nan;
    varargout{4} = nan;
    varargout{5} = nan;
    
end
