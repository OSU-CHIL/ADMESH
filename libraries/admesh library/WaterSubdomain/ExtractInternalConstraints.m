function [Constraints,pgon_old,MA_1D_mainstream,MA_1D_mainstreamLand,W2D_boundaryID] = ExtractInternalConstraints(ShorelineFile,A_min,dx,delta_cw,delta_lfs,delta_A,pruning_rho,pruning_theta,length_min,Nx,Ny)

deg2m = deg2km(1)*1e3;
Scale = deg2m;
%--------------------------------------------------------------------------
% Load example file (variable "Points")
%--------------------------------------------------------------------------
Shoreline = shaperead(ShorelineFile);

for i = 1 : length(Shoreline)
    Shoreline(i).np = length(Shoreline(i).X);
end

%% ========================================================================
% Parameter setting
%==========================================================================
%%========================================================================
% Set boundary points from shorelines and set up parameters
%==========================================================================
%--------------------------------------------------------------------------
% Construct polyshape with filtering small islands
%--------------------------------------------------------------------------
Shoreline1 = struct2table(Shoreline);
x = Shoreline1.X(:);
y = Shoreline1.Y(:);
x = cellfun(@(x) x(~isnan(x)),x,'UniformOutput',0);
y = cellfun(@(x) x(~isnan(x)),y,'UniformOutput',0);
A = cellfun(@(x,y) polyarea(x,y),x,y,'UniformOutput',0);
A = cell2mat(A);
I = A > A_min;
pgon = polyshape(x(I),y(I));
pgon_removed = polyshape(x(~I),y(~I));

%--------------------------------------------------------------------------
% Load polyshape data (if desired)
%--------------------------------------------------------------------------
pgon_noholes = rmholes(pgon);

I = find(isnan(pgon.Vertices(:,1)));
I = [0; I(:); size(pgon.Vertices,1)+1];
XY_b = {pgon_noholes.Vertices}; % boundary including the dd-box
k1 = 1;
    
for i = 1 : length(I)-1
    i1 = I(i) + 1;
    i2 = I(i+1) - 1;
    J = [i1:i2, i1];
    
    if nnz(~ismember(XY_b{1},pgon.Vertices(J,:))) > 0
        k1 = k1 + 1;
        XY_b{k1} = pgon.Vertices(J,:);
    end
end

%% ========================================================================
% Decompose domain
%==========================================================================
[dd_pgons,dd_Boxes,dd_ID,xg,yg] = DecomposePolyshape(pgon,Nx,Ny,dx,delta_cw*4);

%% ========================================================================
% Compute VDT in DDM
%==========================================================================
[Vxg,Vyg,MaskW,MaskL] = ComputeVDT_DDM(pgon,dd_ID,dd_Boxes,xg,yg,dx);


%% ========================================================================
% Compute Medial Axis
%==========================================================================
MaskMA_raw = ComputeMA_DDM(Vxg,Vyg,dd_ID);

%% ========================================================================
% Thinning MA
%%=========================================================================
MaskMA_thinned = ThinningMA(MaskMA_raw);
fprintf('Thinning is done\n');

%% ========================================================================
% Construct branches
%==========================================================================
MA_branch = ConstructMedialAxisBranch(MaskMA_thinned);
MA_branch = RemoveJointDuplicates(MA_branch);
MA = MA_branch;
fprintf('Branch construction is done\n');

%% ========================================================================
% Pruning corners
%==========================================================================
MA_pruned = SkeletonPruning_Choi_level1branch_v2(Vxg,Vyg,MA,pruning_rho,pruning_theta,dx);
fprintf('Pruning corner branches is done\n');

%% ========================================================================
% Compute distance map to medial axis
%==========================================================================
MA = MA_pruned;
D2MA = DistanceToMA(MA,{MaskL,MaskW},dd_ID,dx);

%% ========================================================================
% Compute lfs
%==========================================================================
lfs = sparse(size(MaskW,1),size(MaskW,2));
for iDD = 1 : size(dd_ID,1)
    I = dd_ID{iDD,1};
    J = dd_ID{iDD,2};
  
    lfs(I,J) = sqrt(Vxg(I,J).^2 + Vyg(I,J).^2) + abs(D2MA(I,J));
    
    fprintf('Compute width function (%d/%d)\n',iDD,size(dd_ID,1));
end
fprintf('Computing local feature size is done\n');
clear Dg Vxg Vyg;
% clear D2MA;
%% ========================================================================
% Split 1D and 2D masks based on lfs
%==========================================================================
[Water2D,Water1D] = SplitMask(MaskW,lfs,delta_lfs,delta_A,MA_pruned);
[Land2D,Land1D] = SplitMask(MaskL,lfs,delta_lfs,delta_A,MA_pruned);

%% ========================================================================
% Fill 2D mask regions
%==========================================================================
delta_filling = delta_lfs*sqrt(2);
[Mask2DCell,Mask1DCell] = FillMask2D({Water2D,Land2D},{Water1D,Land1D},MA_pruned,lfs,delta_filling,delta_A,xg,yg,XY_b,dx);

Water2D = Mask2DCell{1};
Water1D = Mask1DCell{1};
Land2D = Mask2DCell{2};
Land1D = Mask1DCell{2};

clear Mask2DCell Mask1DCell;


%% ========================================================================
% Transfer 1D masks to 2D masks based on connectivity.
%==========================================================================
Ratio = .7;

BranchNodesID = vertcat(MA_pruned.BranchNodes{:});
BranchNodesID = unique(BranchNodesID);

CC = bwconncomp(Land1D);
if CC.NumObjects > 0
%--------------------------------------------------------------------------
% Compute PA ratio
%--------------------------------------------------------------------------
CCRP = regionprops(CC,'MaxFeretProperties','MinFeretProperties',...
    'Perimeter','Area','ConvexArea','BoundingBox','SubarrayIdx',...
    'MajorAxisLength','MinorAxisLength','Centroid');
CCRP = struct2table(CCRP,'AsArray',1);
CC.PAR = CCRP.Perimeter./CCRP.Area; % Isoperimetric ratio
CC.IPR = (CCRP.Perimeter).^2./CCRP.Area; % Isoperimetric ratio
CC.FDR = CCRP.MaxFeretDiameter./CCRP.MinFeretDiameter; % Feret diameter ratio
CC.CAR = CCRP.ConvexArea./CCRP.Area;
CC.AR = CCRP.MajorAxisLength./CCRP.MinorAxisLength;

%--------------------------------------------------------------------------
% Use Euclidean buffer with size MinorAxisLength and check ratio of water2d and land2d
%--------------------------------------------------------------------------
Buffer = cell(CC.NumObjects,1);
for i = 1 : CC.NumObjects
    id = CC.PixelIdxList{i};
    [temp_I,temp_J] = ind2sub(MA_pruned.Size,id);
    
    Box = CCRP.BoundingBox(i,:);
%     BufferSize = round(CCRP.MaxFeretDiameter(i)/2);
%     BufferSize = round(max(Box(3:4))/2);
    BufferSize = round(CCRP.MinorAxisLength(i)/2);
    
    minI = max(1,floor(Box(2))-BufferSize);
    maxI = min(MA_pruned.Size(1),floor(Box(2)+Box(4))+BufferSize);
    Buffer_I = minI : maxI;
    
    minJ = max(1,floor(Box(1))-BufferSize);
    maxJ = min(MA_pruned.Size(2),floor(Box(1)+Box(3))+BufferSize);
    Buffer_J = minJ : maxJ;
    
    [Buffer_I,Buffer_J] = meshgrid(Buffer_I,Buffer_J);
    Buffer_I = Buffer_I(:);
    Buffer_J = Buffer_J(:);
    
    [~,dist] = knnsearch([temp_I,temp_J],[Buffer_I(:),Buffer_J(:)]);
    
    I = dist < BufferSize & dist > 0;
    if nnz(I) == 0
        0;
    end
    Buffer{i} = sub2ind(MA_pruned.Size,Buffer_I(I),Buffer_J(I));
end
temp_WLratio = cellfun(@(x) nnz(Water2D(x))/nnz(Land2D(x)),Buffer);
I = cellfun(@(x) isempty(x),Buffer);
temp_WLratio(I) = 10;
temp_WLratio(isnan(temp_WLratio)) = 0;

temp_K = temp_WLratio > 2 & CC.IPR > 30;
temp_K = vertcat(CC.PixelIdxList{temp_K});
L1DtoW2D = false(CC.ImageSize(1),CC.ImageSize(2));
L1DtoW2D(temp_K) = 1;

else
    L1DtoW2D = false(CC.ImageSize(1),CC.ImageSize(2));
end
    
%% ========================================================================
% Construct new water mask
%%=========================================================================
%--------------------------------------------------------------------------
% Filter out isolated regions
%--------------------------------------------------------------------------
NewW2D = Water2D | L1DtoW2D;
CC = bwconncomp(NewW2D);
for i = 1 : CC.NumObjects
    id = CC.PixelIdxList{i};
    if ~any(Water2D(id))
        L1DtoW2D(id) = 0;
    end
end
clear NewW2D;

%--------------------------------------------------------------------------
% Identify branch nodes indexes in L1DtoW2D mask
%--------------------------------------------------------------------------
BranchNodesID = vertcat(MA_pruned.BranchNodes{:});
BranchNodesID = unique(BranchNodesID);
K = BranchNodesID(L1DtoW2D(BranchNodesID));

L1DtoW2Dskel = false(size(L1DtoW2D));
L1DtoW2Dskel(K) = 1;

MA_LandBarrier = ConstructMedialAxisBranch(L1DtoW2Dskel);
MA_LandBarrier = RemoveJointDuplicates(MA_LandBarrier);

%% ========================================================================
% Patch Water 2D mask with L1DtoW2D mask
%==========================================================================
L1DtoW2D_XY = Mask2XY(L1DtoW2D,xg,yg);
if ~isempty(L1DtoW2D_XY)
pgon_x = cellfun(@(x) x(:,1),L1DtoW2D_XY,'UniformOutput',0);
pgon_y = cellfun(@(x) x(:,2),L1DtoW2D_XY,'UniformOutput',0);
pgon_x = cellfun(@(x) [x; nan],pgon_x,'UniformOutput',0);
pgon_y = cellfun(@(x) [x; nan],pgon_y,'UniformOutput',0);

MiterLimit = 2; % This value is tentative. Experiments would be need
pgon_add = polyshape(pgon_x,pgon_y);
pgon_add = polybuffer(pgon_add,dx*sqrt(2),'JointType','miter','MiterLimit',MiterLimit);

pgon_new = union(pgon,pgon_add);
else
    pgon_new = pgon;
end
% Fill small holes
pgon_new_holes = holes(pgon_new);
pgon_new_holes_area = area(pgon_new_holes);
I = pgon_new_holes_area/(dx^2) < 10;
if nnz(I) > 0
pgon_new_holes = union(pgon_new_holes(I));
pgon_new = union(pgon_new,pgon_new_holes);
end
% Update polygon (backup old one)
pgon_old = pgon;
pgon = pgon_new;

%% ========================================================================
% Begin again with updated polygon
%==========================================================================
pgon_noholes = rmholes(pgon);

I = find(isnan(pgon.Vertices(:,1)));
I = [0; I(:); size(pgon.Vertices,1)+1];
XY_b = {pgon_noholes.Vertices}; % boundary including the dd-box
k1 = 1;
    
for i = 1 : length(I)-1
    i1 = I(i) + 1;
    i2 = I(i+1) - 1;
    J = [i1:i2, i1];
    
    if nnz(~ismember(XY_b{1},pgon.Vertices(J,:))) > 0
        k1 = k1 + 1;
        XY_b{k1} = pgon.Vertices(J,:);
    end
end


%% Decompose domain
[dd_pgons,dd_Boxes,dd_ID,xg,yg] = DecomposePolyshape(pgon,Nx,Ny,dx,delta_cw*4);

%% ========================================================================
% Compute VDT in DDM
%==========================================================================
[Vxg,Vyg,MaskW,MaskL] = ComputeVDT_DDM(pgon,dd_ID,dd_Boxes,xg,yg,dx);
% Land mask is not required at this step
Vxg(MaskL) = 0;
Vyg(MaskL) = 0;

%% Compute Medial Axis
MaskMA_raw = ComputeMA_DDM(Vxg,Vyg,dd_ID);

%% ========================================================================
% Thinning MA
%%=========================================================================
MaskMA_thinned = ThinningMA(MaskMA_raw);
fprintf('Thinning is done\n');

%% ========================================================================
% Construct branches (with redundant joints)
%==========================================================================
MA_branch = ConstructMedialAxisBranch(MaskMA_thinned);
MA_branch = RemoveJointDuplicates(MA_branch);
MA = MA_branch;
fprintf('Branch construction is done\n');

%% ========================================================================
% Pruning corners
%==========================================================================
MA_pruned = SkeletonPruning_Choi_level1branch_v2(Vxg,Vyg,MA,pruning_rho,pruning_theta,dx);
fprintf('Pruning corner branches is done\n');

%% ========================================================================
% Compute distance map to medial axis
%==========================================================================
MA = MA_pruned;
D2MA = DistanceToMA(MA,{MaskW},dd_ID,dx);

%% ========================================================================
% Compute lfs
%==========================================================================
lfs = sparse(size(MaskW,1),size(MaskW,2));
for iDD = 1 : size(dd_ID,1)
    I = dd_ID{iDD,1};
    J = dd_ID{iDD,2};
  
    lfs(I,J) = sqrt(Vxg(I,J).^2 + Vyg(I,J).^2) + abs(D2MA(I,J));

    fprintf('Compute lfs (%d/%d)',iDD,size(dd_ID,1));
end
fprintf('Computing local feature size is done\n');
clear Dg D2MA Vxg Vyg;

%% ========================================================================
% Find water masks based on lfs (filter with connectivity to MA)
% (This is a part of script below, but written in function)
%==========================================================================
[Water2D,Water1D] = SplitMask(MaskW,lfs,delta_lfs,delta_A,MA_pruned);

%% ========================================================================
% Fill 2D mask regions
%==========================================================================
tempW2D = Water2D;
tempW2D(L1DtoW2D) = 1;
tempW1D = Water1D;
tempW1D(L1DtoW2D) = 1;

temp_MA = MA_pruned;
temp_MA.BranchNodes = vertcat(temp_MA.BranchNodes(:),MA_LandBarrier.BranchNodes(:));

delta_filling = delta_lfs*sqrt(2);
[~,~,M1DtoM2DCell] = FillMask2D({tempW2D},{tempW1D},temp_MA,lfs,delta_filling,delta_A,xg,yg,XY_b,dx);

M1DtoM2D = M1DtoM2DCell{1};

tempW2D = Water2D;
tempW2D(M1DtoM2D) = 1;
tempW1D = Water1D;
tempW1D(M1DtoM2D) = 0;

BranchNodesID = vertcat(temp_MA.BranchNodes{:});
BranchNodesID = unique(BranchNodesID);

[I,J] = ind2sub(size(tempW2D),BranchNodesID);
K = tempW1D(BranchNodesID);
MA1D = sparse(I,J,K,size(tempW2D,1),size(tempW2D,2));

[I,J] = ind2sub(size(tempW2D),BranchNodesID);
K = tempW2D(BranchNodesID);
MA2D = sparse(I,J,K,size(tempW2D,1),size(tempW2D,2));

CC = bwconncomp(full(tempW1D));
for i = 1 : CC.NumObjects
    id = CC.PixelIdxList{i};
    if ~any(MA1D(id))
        tempW2D(id) = 1;
        tempW1D(id) = 0;
    end
    fprintf('Fill 2D mask regions (%d/%d)',i,CC.NumObjects);
end

%----------------------------------------------------------------------
% Remove 2D mask not containing any MA points (it may need only when 
% the filling results isolated mask, which should not happen)
%----------------------------------------------------------------------
CC = bwconncomp(full(tempW2D));
for i = 1 : CC.NumObjects
    id = CC.PixelIdxList{i};
    if ~any(MA2D(id)) || numel(id) < delta_A
        tempW2D(id) = 0;
        tempW1D(id) = 1;
    end
end

Water2D = tempW2D;
Water1D = tempW1D;

%% ========================================================================
% Construct MA1D
%==========================================================================
MA_1D = ConstructMA1D(MA_pruned,Water1D,Water2D);
MA_1D = RemoveJointDuplicates(MA_1D);

%% ========================================================================
% Extract Water2D mask boundary
%==========================================================================
[~,W2D_boundaryID] = Mask2XY(Water2D,xg,yg);

temp_x = cellfun(@(x) xg(x(:,2)),W2D_boundaryID,'UniformOutput',0);
temp_y = cellfun(@(x) yg(x(:,1)),W2D_boundaryID,'UniformOutput',0);
pgon_W2D = polyshape(temp_x,temp_y);
%% ========================================================================
% Connect branches to Water2D boundary
%==========================================================================
MA = MA_1D;
MA_connected = ConnectMA1Dto2DArea_v2(MA,W2D_boundaryID,2);
MA_connected = PruneLevel1BranchByLength(dx,MA_connected,length_min,MA_connected.FlagConnected2D);
MA_connected = RemoveJointDuplicates(MA_connected);

MA = MA_LandBarrier;
if ~isempty(MA_LandBarrier.BranchNodes)
MA_connectedLand = ConnectMA1Dto2DArea_v2(MA,W2D_boundaryID,2);
MA_connectedLand = PruneLevel1BranchByLength(dx,MA_connectedLand,length_min,MA_connectedLand.FlagConnected2D);
MA_connectedLand = RemoveJointDuplicates(MA_connectedLand);
else
    MA_connectedLand = MA_LandBarrier;
end

%% ========================================================================
% Construct mainstreams
%==========================================================================
MA = MA_connected;
MA_1D_mainstream = ConstructMainStream_v3(MA);
MA_1D_mainstream = RemoveJointDuplicates(MA_1D_mainstream);

MA = MA_connectedLand;
MA_1D_mainstreamLand = ConstructMainStream_v3(MA);
MA_1D_mainstreamLand = RemoveJointDuplicates(MA_1D_mainstreamLand);

%% ========================================================================
% Construct output constraints
%==========================================================================
f = @(x) ind2xy(x.BranchNodes,x.Size,xg,yg);
ConstraintsCell{1} = f(MA_1D_mainstream);
ConstraintsCell{2} = f(MA_1D_mainstreamLand);
ConstraintsCell{3} = sub2xy(W2D_boundaryID,xg,yg);

k = 0;
ConstraintNum = -[18 17 19];
Constraints = [];
for i = 1 : length(ConstraintsCell)
    for j = 1 : length(ConstraintsCell{i})
        k = k + 1;
        Constraints(k).num = ConstraintNum(i);
        Constraints(k).xy = ConstraintsCell{i}{j};
        Constraints(k).type = 'line';
    end
end









