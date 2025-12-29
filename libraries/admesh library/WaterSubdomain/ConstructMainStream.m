function MA_new = ConstructMainStream(MA,UIFigure)

msg = 'Constructing main stream...';
progdlg = uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg);

%% ========================================================================
% Construct main stream
%==========================================================================
Size = MA.Size;
BranchNodes = MA.BranchNodes;
BranchNodes1 = MA.BranchNodes;
[temp1,temp2] = cellfun(@(x) ind2sub(MA.Size,x),BranchNodes,'UniformOutput',0);
BranchXY = cellfun(@(x,y) [x(:),y(:)],temp1,temp2,'UniformOutput',0);

BranchNodeEnds = cellfun(@(x) [x(1) x(end)],MA.BranchNodes,'UniformOutput',0);
BranchNodeEnds = cell2mat(BranchNodeEnds(:));

NodeEndList = unique(BranchNodeEnds(:));
MergeID = zeros(length(NodeEndList),2);
for k = 1 : length(NodeEndList)
    I1 = find(BranchNodeEnds(:,1) == NodeEndList(k));
    I2 = find(BranchNodeEnds(:,2) == NodeEndList(k));
    I = [I1; I2];
    if length(I) == 1
        continue;
    end
    if ~isempty(I1) || ~isempty(I2)
        x = cell(length(I1)+length(I2),1);
        y = cell(length(I1)+length(I2),1);
        for i = 1 : length(I1)
            x{i} = BranchXY{I1(i)}(:,1);
            y{i} = BranchXY{I1(i)}(:,2);
        end
        for j = 1 : length(I2)
            x{length(I1) + j} = flipud(BranchXY{I2(j)}(:,1));
            y{length(I1) + j} = flipud(BranchXY{I2(j)}(:,2));
        end
    end

    x1 = x;
    y1 = y;
    
    SmoothingRMSE = km2deg(.5e-3);
    K = inf(length(x1));
    for i = 1 : length(x1)
        for j = i+1 : length(x1)
            x2 = [flip(x1{j}); x1{i}(2:end)];
            y2 = [flip(y1{j}); y1{i}(2:end)];
           
            FixedPoints = [x2([1,end]),y2([1,end]); x1{i}(1), y1{i}(1)];
            PI = ComputePathCurvature([x2,y2],FixedPoints,SmoothingRMSE);
            id = length(x1{j});
            K(i,j) = abs(PI.k(PI.p(id)));
        end
    end
    [i,j] = find(K == min(K(:)),1);
    MergeID(k,:) = I([i,j]);
    
    progdlg.Value = k/length(NodeEndList);
end
close(progdlg);

MergeID1 = MergeID;
BranchNodes1 = BranchNodes;
for k = 1 : size(MergeID1,1)
    I1 = MergeID1(k,1);
    I2 = MergeID1(k,2);
    if I1 == I2
        continue;
    end
    
    nodes1 = BranchNodes1{I1};
    nodes2 = BranchNodes1{I2};
    nodes1 = nodes1(:);
    nodes2 = nodes2(:);
    
    if nodes1(end) == nodes2(1)
        nodes = [nodes1; nodes2(2:end)];
    elseif nodes1(end) == nodes2(end)
        nodes = [nodes1; flipud(nodes2(1:end-1))];
    elseif nodes1(1) == nodes2(1)
        nodes = [flipud(nodes1); nodes2(2:end)];
    elseif nodes1(1) == nodes2(end)
        nodes = [flipud(nodes1); flipud(nodes2(1:end-1))];
    else
        % This case happens when a loop
        continue;
    end
    
    BranchNodes1{I1} = nodes;
    BranchNodes1{I2} = [];
    
    J1 = MergeID1 == I2;
    MergeID1(J1) = I1;
end

I = cellfun(@(x) isempty(x),BranchNodes1);
BranchNodes1(I) = [];
BranchNodes1 = BranchNodes1(:);


MA_new.Size = MA.Size;
MA_new.BranchNodes = BranchNodes1;
% MA_new.BranchLength = MainStreamLength;

return;

figure; hold on; axis equal;
for i = 1 : length(BranchXY1)
%     plot(BranchXY{i}(:,1),BranchXY{i}(:,2));
    plot(BranchXY1{i}(:,1),BranchXY1{i}(:,2));
end


MainStreamBranch = {};
MainStreamLength = [];
wbar = waitbar(0);
for i = 1 : length(BranchConnectivity)
    
    % Case 1: Isolated branch
    if isempty(BranchConnectivity{i})
        MainStreamBranch{end+1} = i;
        MainStreamLength(end+1) = BranchLength(i);
        continue;
    end
    
    x1 = BranchXY{i}(:,1);
    y1 = BranchXY{i}(:,2);
    t1 = 1 : length(x1);
    
    x1 = csaps(t1,x1,1e-2,t1);
    y1 = csaps(t1,y1,1e-2,t1);
    
    for j = 1 : length(BranchConnectivity{i})
        k = BranchConnectivity{i}(j);
        x2 = BranchXY{k}(:,1);
        y2 = BranchXY{k}(:,2);
        t2 = 1 : length(x2);
    
        x2 = csaps(t2,x2,1e-2,t2);
        y2 = csaps(t2,y2,1e-2,t2);
        
    end
    
    % Case 2: Not a level-1 branch
    if ~all(BranchConnectivity{i} > 0) && ~all(BranchConnectivity{i} < 0)
        continue;
    end
    
    % Case 3: Level 1 branch
    BranchChunk = i;
    while 1
        ChunkAdd = [];
        for j = 1 : length(BranchChunk)
            k = BranchChunk(j);
            ChunkAdd = [ChunkAdd; setdiff(abs(BranchConnectivity{k}),BranchChunk)];
        end
        if isempty(ChunkAdd)
            break;
        end
        BranchChunk = [BranchChunk; ChunkAdd];
    end
    
    for j = 1 : length(BranchChunk)
        k = BranchChunk(j);
        if all(BranchConnectivity{k} > 0) || all(BranchConnectivity{k} < 0)
            BranchChunk(j) = -BranchChunk(j);
        end
    end
    BranchChunk = sort(BranchChunk);
    
    
    while 1
        Path2 = [];
        for j = 1 : nnz(BranchChunk < 0)
            k = -BranchChunk(j);
            Path = {k};
            while 1
                PathNew = [];
                for kk = 1 : length(Path)
                    PathAdd = setdiff(abs(BranchConnectivity{Path{kk}(end)}),unique(vertcat(Path{:})));
                    PathAdd = setdiff(PathAdd,vertcat(MainStreamBranch{:}));
                    %             PathAdd = setdiff(PathAdd,unique([Path{:}]));
                    if isempty(PathAdd)
                        PathNew{end+1} = Path{kk}(:);
                    else
                        for kkk = 1 : length(PathAdd)
                            PathNew{end+1} = [Path{kk}(:); PathAdd(kkk)];
                        end
                    end
                end
                if isequal(size(Path),size(PathNew)) && all(cellfun(@(x,y) all(ismember(x,y)),Path,PathNew))
                    break;
                end
                Path = PathNew;
            end
            Path2 = [Path2, Path];
        end
        
        for j = 1 : length(Path2)
            Path2{j} = setdiff(Path2{j},vertcat(MainStreamBranch{:}),'stable');
        end
        
        PathLength = zeros(length(Path2),1);
        for j = 1 : length(Path2)
            PathLength(j) = sum(BranchLength(Path2{j}));
        end
        [maxlength,imax] = max(PathLength);
        
        if isempty(Path2{imax})
            break;
        end
        MainStreamBranch{end+1} = Path2{imax};
        MainStreamLength(end+1) = maxlength;
    end
    waitbar(i/length(BranchConnectivity),wbar);
end
delete(wbar);

%--------------------------------------------------------------------------
% Handle not included branches... especially short bridges
%--------------------------------------------------------------------------
id = find(~ismember(1:length(BranchConnectivity),vertcat(MainStreamBranch{:})));
temp_Length = BranchLength(id);
id = mat2cell(id(:),ones(length(id),1),1);

MainStreamBranch = [MainStreamBranch(:); id(:)];
MainStreamLength = [MainStreamLength(:); temp_Length(:)];


BranchNodes_new = {};
MainStreamNodes = cell(length(MainStreamBranch),1);

wbar = waitbar(0);
for i = 1 : length(MainStreamBranch)
    if length(MainStreamBranch{i}) == 1
        k = MainStreamBranch{i};
        MainStreamNodes{i} = BranchNodes{k};
        JointConnectivity_new(i,:) = [BranchNodes{k}(1),BranchNodes{k}(end)];
        continue;
    end
    
    k = MainStreamBranch{i}(1);
    k1 = MainStreamBranch{i}(2);
    if any(BranchConnectivity{k} == k1)
        iMainStreamNodes = (BranchNodes{k});
    elseif any(BranchConnectivity{k} == -k1)
        iMainStreamNodes = flip(BranchNodes{k});
    end
        
    for j = 2 : length(MainStreamBranch{i})
        k1 = MainStreamBranch{i}(j-1);
        k = MainStreamBranch{i}(j);
        if any(BranchConnectivity{k} == -k1)
            iMainStreamNodes = [iMainStreamNodes; (BranchNodes{k})];
        elseif any(BranchConnectivity{k} == k1)
            iMainStreamNodes = [iMainStreamNodes; flip(BranchNodes{k})];
        end
    end
    
    iMainStreamNodes = unique(iMainStreamNodes,'rows','stable');
    MainStreamNodes{i} = iMainStreamNodes;
    JointConnectivity_new(i,:) = [iMainStreamNodes(1),iMainStreamNodes(end)];
    waitbar(i/length(MainStreamBranch),wbar);
end
delete(wbar);

MA_new.Size = MA.Size;
MA_new.nBranch = length(MainStreamNodes);
MA_new.BranchNodes = MainStreamNodes;
MA_new.JointConnectivity = JointConnectivity_new;
MA_new.BranchLength = MainStreamLength;


