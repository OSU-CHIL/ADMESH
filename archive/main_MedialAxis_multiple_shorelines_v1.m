%==========================================================================
% main_MedialAxis_multiple_shorelines_v1
% - This function does this, this, and this...
% 
% Update history
% 2021-03-16 (v1): - written by Younghun Kang
% 
%==========================================================================

%% ========================================================================
% Load datasets
%==========================================================================
% clear all; clc; close all;
addpath(genpath('admesh-lib-part'));
%--------------------------------------------------------------------------
% Load example file (variable "Points")
%--------------------------------------------------------------------------
ShorelineFile = ['D:\Research\PR-PREEVENTS\Data\lower_Neches\ArcGIS\',...
                'ClintExtent/ShorelineCopy.shp'];
Shoreline     = shaperead(ShorelineFile);

for i = 1 : length(Shoreline)
    Shoreline(i).np = length(Shoreline(i).X);
end
% [~,index] = sortrows([Shoreline.np].'); Shoreline = Shoreline(index(end:-1:1)); clear index;
%% ========================================================================
% Parameter setting
%==========================================================================
close all;
DEG2KM = 6378*pi/180;
%%========================================================================
% Set boundary points from shorelines and set up parameters
%==========================================================================
fid = [1571 702 703 718 727 796 797 798 819 825 826 828 829 832 836];
% fid = 1571;
id = fid+1;
XY_b = cell(length(id),1);
x_b = []; y_b = [];
for i = 1 : length(id)
    j = id(i);
    XY_b{i} = [Shoreline(j).X; Shoreline(j).Y]';
    x_b = [x_b(:); XY_b{i}(:,1); nan]; 
    y_b = [y_b(:); XY_b{i}(:,2); nan]; 
end
% x_b = XY_b(:,1); y_b = XY_b(:,2);
% delta = 1e-5;
delta = 1/(DEG2KM*1e3);
hmax = 2*delta;
dx = 1*delta;

delta_lfs = 10/(DEG2KM*1e3); % meter to degree
delta_A = 500/(dx*DEG2KM*1e3)^2; % square meter to number of pixel
rho_pruning = 4*delta_lfs;
rho = 4*dx;
dtheta = pi*.9;
SmoothingParam = 0.01;

%% ========================================================================
% Compute medial axis
%==========================================================================
%--------------------------------------------------------------------------
% Construct background grid
%--------------------------------------------------------------------------
x = min(x_b) : dx : max(x_b);
y = min(y_b) : dx : max(y_b);
[X,Y] = meshgrid(x,y);

%--------------------------------------------------------------------------
% Compute distance from shoreline and gradient of it
%--------------------------------------------------------------------------
PTS = [];
for i = 1 : length(XY_b)
    PTS.Poly(i).x = XY_b{i}(:,1);
    PTS.Poly(i).y = XY_b{i}(:,2);
end
[D,gradD] = SignedDistanceFunction_v2(PTS,X,Y,delta,hmax);
% D = round(D,16);
id_out = D > 0;
% D(D > 0) = nan;
fprintf('Computing distance map is done\n');

% GradD = sqrt(gradD.x.^2 + gradD.y.^2);
[GradDx,GradDy] = gradient(D,delta);
GradD = sqrt(GradDx.^2 + GradDy.^2);

%--------------------------------------------------------------------------
% Compute 8SSED map
%--------------------------------------------------------------------------
[Vx,Vy] = Compute8SSED(X,Y,x_b,y_b);
Vx(id_out) = 0; Vy(id_out) = 0;
Vmag = sqrt(Vx.^2 + Vy.^2);
fprintf('Computing 8SSED map is done\n');

%--------------------------------------------------------------------------
% Find medial axis
%--------------------------------------------------------------------------
% id_MA = GradD < 0.9 & D <= 0;
% id_MA = Skeletonization_Choi_v2(BX,BY,Vx,Vy,2*dx,rho_pruning,dtheta);
id_MA_original = Skeletonization_DivSSED_v2(X,Y,Vx,Vy);
id_MA_original(id_out) = 0;
fprintf('Skeletonization is done\n');

%% ========================================================================
% Thinning MA
%%=========================================================================
bw = id_MA_original;

%--------------------------------------------------------------------------
% Fill one pixel holes (this part is temporarily added to improve the 
% quality of MA branches, which is constructed later)
%--------------------------------------------------------------------------
bw = bwmorph(bw,'fill');

%--------------------------------------------------------------------------
% (???) no idea why.. but bwskel looks better than bwmorph...
%--------------------------------------------------------------------------
% id_MA_thinned = bwmorph(bw,'thin');
id_MA_thinned = bwskel(bw);
% id_MA_thinned = bwmorph(bw,'skel',inf);

%--------------------------------------------------------------------------
% Remove isolated pixels
%--------------------------------------------------------------------------
id_MA_thinned = bwmorph(id_MA_thinned,'clean',inf);

% %--------------------------------------------------------------------------
% % I didn't mean it but this part fills holes
% % (it's better to merge joints latter...)
% % 2021-03-15 This part is commented out as it behaves weird
% %--------------------------------------------------------------------------
% [B_MA_thinned,L_MA_thinned] = bwboundaries(id_MA_thinned);
% nL = zeros(max(L_MA_thinned(:)),1);
% for i = 1 : max(L_MA_thinned(:))
%     nL(i) = nnz(L_MA_thinned == i);
% end
% A_MA_thinned = nL;
% id_MA_thinned = logical(L_MA_thinned > 0);
fprintf('Thinning is done\n');

%% ========================================================================
% Construct branches (with redundant joints)
%==========================================================================
MA_branch = ConstructMedialAxisBranch(id_MA_thinned);
MA = MA_branch;
fprintf('Branch construction is done\n'); 



%% ========================================================================
% Merge redundant joints
%==========================================================================
% MA_merged = MergeJoints2ClosestGridPoint(MA);
% MA = MA_merged;
% id_MA_joint = MA_merged.ID_Joint;
% fprintf('Merging joints is done\n');


%% ========================================================================
% Pruning corners
%==========================================================================
MA_pruned = SkeletonPruning_Choi_level1branch(X,Y,Vx,Vy,D,MA,rho_pruning,dtheta);
MA = MA_pruned;
fprintf('Pruning corner branches is done\n');

% PLOT = 'pruned'; debug_plots;
%% ========================================================================
% Separate 1D & 2D area based on local feature size
%==========================================================================
id_MA_pruned = vertcat(MA.ID{:});

bw_MA_pruned = zeros(size(X));
bw_MA_pruned(id_MA_pruned) = 1;
dMA_pruned = bwdist(bw_MA_pruned)*dx;
dMA_pruned(id_out) = nan;

lfs = (abs(dMA_pruned) + abs(D));
%--------------------------------------------------------------------------
% Truncate based on local feature size
%--------------------------------------------------------------------------
lfs_2D = lfs;
lfs_2D(id_out) = nan;
id_L1D = lfs_2D < delta_lfs;
lfs_2D(id_L1D) = nan;

%--------------------------------------------------------------------------
% Get boundary
%--------------------------------------------------------------------------
bw_2D = lfs_2D;
bw_2D(~isnan(bw_2D)) = 1;
bw_2D(isnan(bw_2D)) = 0;
bw_2D = bw_2D';
[B,L2D] = bwboundaries(bw_2D);
L2D = L2D';

%--------------------------------------------------------------------------
% Compute area of each labeled region
%--------------------------------------------------------------------------
nL = zeros(max(L2D(:)),1);
for i = 1 : max(L2D(:))
    nL(i) = nnz(L2D == i);
end
A = nL;

%--------------------------------------------------------------------------
% Remove small regions
%--------------------------------------------------------------------------
id = find(A < delta_A);
B(id) = [];
for i = 1 : length(id)
    L2D(L2D == id(i)) = 0;
end
%--------------------------------------------------------------------------
% Re-label 2D layers
%--------------------------------------------------------------------------
% [L2D_boundary,L2D] = bwboundaries(L2D);
idL2D = unique(L2D);
k = 0;
for i = 1 : length(idL2D)
    id = find(L2D == idL2D(i));
    if ~isempty(id)
        L2D(id) = k;
        k = k + 1;
    end
end

%--------------------------------------------------------------------------
% Make 1D layer
%--------------------------------------------------------------------------
L0 = double(~id_out);
L1D = L0;
L1D(L2D > 0) = 0;
[L1D_boundary,L1D] = bwboundaries(L1D);

%--------------------------------------------------------------------------
% Get indices for MA1D and MA2D
%--------------------------------------------------------------------------
id_MA1D = intersect(id_MA_pruned,find(L1D > 0));
id_MA2D = intersect(id_MA_pruned,find(L2D > 0));
XY_MA1D = [X(id_MA1D),Y(id_MA1D)];
XY_MA2D = [X(id_MA2D),Y(id_MA2D)];

fprintf('1D/2D area seperation is done\n');

%% ========================================================================
% Fill 2D area
%==========================================================================
L2D_Filled = L2D;
L1D_Filled = L1D;
for i = 1 : length(id_MA2D)
    k = id_MA2D(i);
    iL2D = L2D(k);
    
    x = X(k); y = Y(k);
    iVmag = Vmag(k);
    
    num_nghb_buffer = ceil(iVmag/dx)+1;
    
    [I,J] = ind2sub(size(X),k);
    
    I1 = max(I - num_nghb_buffer,1);
    I2 = min(I + num_nghb_buffer,size(X,1));
    
    J1 = max(J - num_nghb_buffer,1);
    J2 = min(J + num_nghb_buffer,size(X,2));
    
    I3 = I1:I2;
    J3 = J1:J2;
    [I3,J3] = meshgrid(I3,J3);
    I3 = I3(:);
    J3 = J3(:);
    
    k_nghb = sub2ind(size(X),I3,J3);
    
%     k_nghb = k + [-N-1, -N, -N+1, -1, 1, N-1, N, N+1];
    
    tempx = X(k_nghb);
    tempy = Y(k_nghb);
    
    dist = sqrt( (tempx - x).^2 + (tempy - y).^2);
    
    id = (dist < iVmag);
    id = k_nghb(id);
    id(id_out(id)) = [];
%     I4 = I3(id);
%     J4 = J3(id);
    
    L2D_Filled(id) = iL2D;
    L1D_Filled(id) = 0;
end
%--------------------------------------------------------------------------
% Re-label layers as some layers can be merged
%--------------------------------------------------------------------------
L2D_Filled(L2D_Filled > 0) = 1;
L1D_Filled(L1D_Filled > 0) = 1;
[~,L2D_Filled] = bwboundaries(L2D_Filled);
[~,L1D_Filled] = bwboundaries(L1D_Filled);

L2D_Filled(id_out) = 0;
L1D_Filled(id_out) = 0;

id_MA1D_Filled = intersect(id_MA_pruned, find(L1D_Filled > 0));
id_MA2D_Filled = intersect(id_MA_pruned, find(L2D_Filled > 0));

fprintf('Filling 2D area is done\n');

%% ========================================================================
% Remove MA points in 2D area
%==========================================================================
nBranch = 0;
ID = [];
Conn = [];
XY = [];
for i = 1 : length(MA.ID)
    id = MA.ID{i};
    is_MA1D = L1D_Filled(id) > 0;
    if any(is_MA1D)
        nBranch = nBranch + 1;
        id = id(is_MA1D);
        ID{nBranch} = id;
        Conn(nBranch,:) = [id(1), id(end)];
        XY{nBranch} = [X(id),Y(id)];
    end
end

MA_1D.nBranch = nBranch;
MA_1D.Size = MA.Size;
MA_1D.ID = ID;
MA_1D.Conn = Conn;
MA_1D.XY = XY;

MA = MA_1D;
%%
MA_1D_mainstream = ConstructMainStream(MA);
MA = MA_1D_mainstream;

%% ========================================================================
% Get boundary of L2D
%==========================================================================
L2D_new_boundary = GetLayerBoundary(L2D_Filled,X,Y,x_b,y_b,'Smoothing');

fprintf('Get boundary of 2D layer is done\n');
%% ========================================================================
% Make MA1D is connected to 2D area
%==========================================================================
XY_L2D_boundary = vertcat(L2D_new_boundary{:});

temp_MA = MA;
nBranch = temp_MA.nBranch;
ID = temp_MA.ID;
Conn = temp_MA.Conn;
XY = cell(length(ID),1);

for i = 1 : length(ID)
    id = ID{i};
    k = L1D_Filled(id) > 0;
    id = id(k);
    if isempty(id)
        continue;
    end
    x = X(id);
    y = Y(id);
    
    dist1 = sqrt((XY_L2D_boundary(:,1) - x(1)).^2 + (XY_L2D_boundary(:,2) - y(1)).^2);
    dist2 = sqrt((XY_L2D_boundary(:,1) - x(end)).^2 + (XY_L2D_boundary(:,2) - y(end)).^2);
    
    if any(dist1 < (2)*dx)
        [~,iend] = min(dist1);
%         MA_connected.ID{i} = vertcat(iend,MA_connected.ID{i});
        x = vertcat(XY_L2D_boundary(iend,1),x);
        y = vertcat(XY_L2D_boundary(iend,2),y);
        Conn(i,1) = -Conn(i,1);
    end
    
    if any(dist2 < (2)*dx)
        [~,iend] = min(dist2);
%         MA_connected.ID{i} = vertcat(MA_connected.ID{i},iend);
         x = vertcat(x,XY_L2D_boundary(iend,1));
         y = vertcat(y,XY_L2D_boundary(iend,2));
         Conn(i,2) = -Conn(i,2);
    end
    XY{i} = [x,y];    
end
[~,id] = unique(Conn,'rows');

MA_connected.nBranch = length(id);
MA_connected.ID = ID(id);
MA_connected.Conn = Conn(id,:);
MA_connected.XY = XY(id);

MA = MA_connected;
% MA_connected.nBranch = nBranch;
% MA_connected.ID = ID;
% MA_connected.Conn = Conn;
% MA_connected.XY = XY;

fprintf('Connecting MA1D to 2D area is done\n');

% %--------------------------------------------------------------------------
% % Plot result..
% %--------------------------------------------------------------------------
% figure; hold on;
% plot(x_b,y_b,'k');
% for i = 1 : length(L2D_new_boundary)
%     myPlotOnTop(L2D_new_boundary{i}(:,1),L2D_new_boundary{i}(:,2),'b','linewidth',1.5); 
% end
% for i = 1 : length(MA_connected.ID)
%     x = MA_connected.XY{i}(:,1);
%     y = MA_connected.XY{i}(:,2);
%     myPlotOnTop(x,y,'r','linewidth',1.5);
% end
% axis equal; box on;

%% ========================================================================
% Smoothing MA branches
%==========================================================================
MA = MA_connected;
MA_smoothed = SmoothingMA(MA,0.1);
MA = MA_smoothed;
fprintf('Smoothing is done\n');

figure; hold on;
plot(x_b,y_b,'k');
axis equal; box on;
for i = 1 : length(L2D_new_boundary)
    myPlotOnTop(L2D_new_boundary{i}(:,1),L2D_new_boundary{i}(:,2),'b','linewidth',1.5); 
end
for i = 1 : length(MA_connected.ID)
    x = MA_smoothed.XY{i}(:,1);
    y = MA_smoothed.XY{i}(:,2);
    myPlotOnTop(x,y,'r','linewidth',1.5);
end

return;

%% ========================================================================
% Construct main stream
%==========================================================================
ID = MA_branch.ID;
ID1 = MA_branch.ID;
Conn = MA_branch.Conn;
Conn1 = MA_branch.Conn;
IC = [];
ICONN = [];

ID_new = ID;
Conn_new = [];
i = 0;
while 1
    if i > size(Conn,1)
        break;
    end
    i = i + 1;
% for i = 1 : length(MA_connected.ID)
%     Conn = Conn1;
    nc = 0;
    iconn = Conn(i,:);
    id = find(Conn == iconn(1));
    id(id == i) = [];
    IC = [];
    if iconn(1) == 0 || nnz(Conn == iconn(1)) == 1
        j = find(Conn == iconn(2));
        j(j == iconn(2)) = [];
        if iconn(2) == 0 || isempty(j)
            continue;
        end
        nc = nc + 1;
        IC{nc} = [iconn(1) iconn(2)];
%         ICONN{nc} = i;
        Conn(i,:) = [];
        ID(i) = [];
        k = 0;
        while 1
            k = k + 1;
            if k > length(IC)
                break;
            end
            iid = IC{k};
%             if Conn(iid(end)) == 0
%                 continue;
%             end
%             iiconn = ICONN{k};
            j = find(Conn == iid(end));
            
            for jj = 1 : length(j)
                jjj = mod(j(jj) + size(Conn,1),2*size(Conn,1) -1)+1;
                if Conn(jjj) > 0
                IC{end+1} = [iid,Conn(jjj)];
%                 ICONN{end+1} = [iiconn,mod(j(jj)-1,size(Conn,1))+1];
                else
                    IC{end+1} = [iid,Conn(jjj),-1];
                end
            end
            j = mod(j-1,size(Conn,1))+1;
            Conn(j,:) = [];
            ID(j) = [];
        end        
    elseif iconn(2) == 0 || nnz(Conn == iconn(2)) == 1
        iconn = flip(iconn);
        j = find(Conn == iconn(1));
        j(j == iconn(2)) = [];
        if iconn(2) == 0 || isempty(j)
            continue;
        end
        nc = nc + 1;
        IC{nc} = [iconn(1) iconn(2)];
%         ICONN{nc} = i;
        Conn(i,:) = [];
        ID(i) = [];
        k = 0;
        while 1
            k = k + 1;
            if k > length(IC)
                break;
            end
            iid = IC{k};
%             if Conn(iid(end)) == 0
%                 continue;
%             end
%             iiconn = ICONN{k};
            j = find(Conn == iid(end));
            
            for jj = 1 : length(j)
                jjj = mod(j(jj) + size(Conn,1),2*size(Conn,1));
                if Conn(jjj) > 0
                IC{end+1} = [iid,Conn(jjj)];
%                 ICONN{end+1} = [iiconn,mod(j(jj)-1,size(Conn,1))+1];
                else
                    IC{end+1} = [iid,Conn(jjj),-1];
                end
            end
            j = mod(j-1,size(Conn,1))+1;
            Conn(j,:) = [];
            ID(j) = [];
        end          
    end
    
    if isempty(IC)
        continue;
    end
    IC = flip(IC);
%     ICONN = flip(ICONN);
    iremove = [];
    for j = 1 : length(IC)
        for k = j+1 : length(IC)
            if all(ismember(IC{k},IC{j}))
                IC{k} = [];
                iremove = [iremove,k];
            end
        end
    end
    IC(iremove) = [];

%     Conn = MA_connected.Conn;
    LENGTH = zeros(length(IC),1);
    ICONN = cell(length(IC),1);
    for j = 1 : length(IC)
        id = IC{j};
        id(id < 0) = [];
        ICONN{j} = [];
        for k = 1 : length(id)-1
            kk = find(ismember(Conn1,id([k k+1]),'rows'));
            if isempty(kk)
                kk = -find(ismember(Conn1,id([k+1 k]),'rows'));
            end
            ICONN{j} = [ICONN{j}, kk];
            LENGTH(j) = LENGTH(j) + length(ID1{abs(kk)});
        end
    end
    
    [~,id_main] = max(LENGTH);
    
    id_new = [];
    id_main = ICONN{id_main};
    for j = 1 : length(id_main)
        if id_main(j) > 0
            id_add = ID1{id_main(j)};
        else
            id_main(j) = -id_main(j);
            id_add = flip(ID1{id_main(j)});
        end
        id_new = vertcat(id_new,id_add);
        ID_new{id_main(j)} = [];
    end
    if ~isempty(id_new)
        ID_new{end+1} = id_new;
    end
    
%     ID_new = vertcat(ID{ICONN{id_main}});
    
%     ID{ICONN{id_main}(1)} = ID_new;
%     Conn1(ICONN{id_main}(1),:) = [ID_new(1), ID_new(end)];
    
%     Conn1(ICONN{id_main}(2:end),:) = 0;
    
    0;
    
    
end

id_remove = [];
for i = 1 : length(ID_new)
    if isempty(ID_new{i})
        id_remove = [id_remove,i];
    end
end
ID_new(id_remove) = [];

for i = 1 : length(ID_new)
    Conn_new(i,:) = [ID_new{i}(1),ID_new{i}(end)];
end


figure; hold on;
plot(x_b,y_b,'k');
axis equal; box on;
for ii = 1 : length(L2D_new_boundary)
    plot(L2D_new_boundary{ii}(:,1),L2D_new_boundary{ii}(:,2),'b'); 
end
for ii = 1 : length(ID_new)
    plot(X(ID_new{ii}),Y(ID_new{ii}),'linewidth',1.5);
end
axis(ax);

% figure; hold on;
% plot(SPx,SPy,'k');
% temp = IC{end};
% plot(BX(temp),BY(temp),'linewidth',1.5);
% axis equal; axis(ax);
%     
% for ii = 1 : length(temp)
%     j = ID1{temp(ii)};
%     plot(BX(j),BY(j),'linewidth',1.5);
% end
% axis equal;


%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Stops here...
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%% ========================================================================
% Fill 2D area
%==========================================================================
% L2D_Fill = zeros(size(BX));
L2D_Filled = L2D;
L1D_Filled = L1D;
dx = mean(mean(diff(X,1,2)));
dy = mean(mean(diff(Y,1,1)));
dx = max(dx,dy);
for i = 1 : max(L2D(:))
    id = (L2D(id_MA2D) == i);
    id_MA2D_i = id_MA2D(id);
    for j = 1 : length(id_MA2D_i)
        k = id_MA2D_i(j);
        x = X(k); y = Y(k);
        iVmag = Vmag(k);
        
        num_nghb_buffer = ceil(iVmag/dx)+1;
        
        [I,J] = ind2sub(size(X),k);
        
        I1 = max(I - num_nghb_buffer,1);
        I2 = min(I + num_nghb_buffer,size(X,1));
        
        J1 = max(J - num_nghb_buffer,1);
        J2 = min(J + num_nghb_buffer,size(X,2));
        
        I3 = I1:I2;
        J3 = J1:J2;
        [I3,J3] = meshgrid(I3,J3);
        I3 = I3(:);
        J3 = J3(:);
        
        k_nghb = sub2ind(size(X),I3,J3);
        
        %     k_nghb = k + [-N-1, -N, -N+1, -1, 1, N-1, N, N+1];
        
        tempx = X(k_nghb);
        tempy = Y(k_nghb);
        
        dist = sqrt( (tempx - x).^2 + (tempy - y).^2);
        
        id = (dist < iVmag);
        id = k_nghb(id);
        id(id_out(id)) = [];
        %     I4 = I3(id);
        %     J4 = J3(id);
        
        L2D_Filled(id) = i;
        L1D_Filled(id) = 0;
    end
end


L2D_Filled = L2D;
L1D_Filled = L1D;
for i = 1 : length(id_MA2D)
    k = id_MA2D(i);
    iL2D = L2D(k);
    
    x = X(k); y = Y(k);
    iVmag = Vmag(k);
    
    num_nghb_buffer = ceil(iVmag/dx)+1;
    
    [I,J] = ind2sub(size(X),k);
    
    I1 = max(I - num_nghb_buffer,1);
    I2 = min(I + num_nghb_buffer,size(X,1));
    
    J1 = max(J - num_nghb_buffer,1);
    J2 = min(J + num_nghb_buffer,size(X,2));
    
    I3 = I1:I2;
    J3 = J1:J2;
    [I3,J3] = meshgrid(I3,J3);
    I3 = I3(:);
    J3 = J3(:);
    
    k_nghb = sub2ind(size(X),I3,J3);
    
%     k_nghb = k + [-N-1, -N, -N+1, -1, 1, N-1, N, N+1];
    
    tempx = X(k_nghb);
    tempy = Y(k_nghb);
    
    dist = sqrt( (tempx - x).^2 + (tempy - y).^2);
    
    id = (dist < iVmag);
    id = k_nghb(id);
    id(id_out(id)) = [];
%     I4 = I3(id);
%     J4 = J3(id);
    
    L2D_Filled(id) = iL2D;
    L1D_Filled(id) = 0;
end

temp = nan(size(X));
temp(L1D > 0) = 1; temp(L2D > 0) = 2;
figure; hold on;
surf(X,Y,temp,'EdgeColor','none','FaceColor','interp');
myPlotOnTop(x_b,y_b,'k');
axis equal; box on; 
mycmap('gray',[0 3]);



temp = nan(size(X));
temp(L1D_Filled > 0) = 1; temp(L2D_Filled > 0) = 2;
figure; hold on;
surf(X,Y,L2D_Filled,'EdgeColor','none','FaceColor','interp');
myPlotOnTop(x_b,y_b,'k');
axis equal; box on;  axis(ax);
mycmap('gray',[0 3]);

%% Find boundary
[I,J] = ind2sub(size(X),find(id_2D_fill));
temp1 = sub2ind(size(X),min(I+1,size(X,1)),J);
temp2 = sub2ind(size(X),max(I-1,1),J);
temp3 = sub2ind(size(X),I,min(J+1,size(X,2)));
temp4 = sub2ind(size(X),I,max(J-1,1));

temp = vertcat(temp1,temp2,temp3,temp4);
id_2D_expand = false(size(X));
id_2D_expand(temp) = true;
figure; plot(X(id_2D_expand),Y(id_2D_expand),'.');

id = inpolygon(x_b,y_b,X(temp),Y(temp));
figure; plot(x_b(id),y_b(id),'.r');

% figure; hold on;
% plot(SPx,SPy,'k');
% plot(BX(L2D_new>0),BY(L2D_new>0),'sb','MarkerFaceColor','b');
% plot(BX(L2D>0),BY(L2D>0),'sk','MarkerFaceColor','k');
% view(2); axis equal; mycmap('binary',[0 1]);

% figure; hold on;
% plot(SPx,SPy,'k');
% surf(BX,BY,L2D,'edgecolor','none','facecolor','interp');
% view(2); axis equal; mycmap('binary',[0 1]);

% figure; hold on;
% plot(SPx,SPy,'k');
% surf(BX,BY,L2D_new,'edgecolor','none','facecolor','interp');
% view(2); axis equal; mycmap('binary',[0 1]);








%% ========================================================================
% Remove short branch !!!experimental!!!
%==========================================================================
delta_dist = 30;
iBranchLong = true(nBranch,1);
for i = 1 : nBranch
    x = MA.XY{i}(:,1);
    y = MA.XY{i}(:,2);
    
    dist = ComputeFlowlineDistance(x,y);
    dist = dist(end)*deg2km(1)*1e3;
    
    if dist < delta_dist && any(num_nghb_MA(MA.Conn(i,:)) == 1)
        iBranchLong(i) = 0;
    end
end
% iBranchLong = logical(iBranchLong);
MA.XY_new = MA.XY(iBranchLong);
% MA1D_branch.XY_new(iBranchLong) = [];

figure; hold on;
plot(x_b,y_b,'k');
for i = 1:nnz(iBranchLong)
    plot(MA.XY_new{i}(:,1),MA.XY_new{i}(:,2),'linewidth',1.5);
end
myPlotOnTop(X(id_MA_joint),Y(id_MA_joint),'ob','linewidth',2);
axis equal tight;



%% ========================================================================
% Generate 1D mesh
%==========================================================================
addpath('master-MeshGeneration1D');
Setting.h_min          = 10*km2deg(1e-3);
Setting.h_max          = 100*km2deg(1e-3);
Setting.h0             = 1*Setting.h_min;
Setting.K              = 20;
Setting.g              = 0.5;
Setting.Smoothing      = true;
Setting.SmoothingParam = 1e-2;

Points = MA.XY;

Mesh1D = cell(length(Points),1);
for iFL = 1 : length(Points)
    x = Points{iFL}(:,1);
    y = Points{iFL}(:,2);
    
    %----------------------------------------------------------------------
    % Add fixed points for junctions
    %----------------------------------------------------------------------
    fixedPoints = [];
    for i = 1 : length(Points)
        fixedPoints(i,:) = [Points{i}(end,1),Points{i}(end,2)];
    end
    id = ismember(fixedPoints,[x,y],'rows');
    fixedPoints = fixedPoints(id,:);
    fixedPoints = unique(fixedPoints,'rows');
    
    %----------------------------------------------------------------------
    % Generate 1D mesh
    %----------------------------------------------------------------------
    Mesh1D{iFL} = MeshGeneration1D([x,y],fixedPoints,Setting);
end


figure; hold on;
plot(x_b,y_b,'k');
% plot(BP(:,1),BP(:,2),'k','linewidth',2);
% plot(FL,'b','linewidth',2);
for i = 1 : length(Points)
%     plot(Points{i}(:,1),Points{i}(:,2),'k');
%     plot(sx{i},sy{i},'k');
    if ~isempty(Mesh1D{i})
        plot(Mesh1D{i}.X,Mesh1D{i}.Y,'-*b');
    end
end
daspect([1 1 1]);


% figure; plot(x,y);
% hold on; plot(fixedPoints(:,1),fixedPoints(:,2),'o');
% axis equal;

