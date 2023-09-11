%% ========================================================================
% Load datasets
%==========================================================================
% clear all; clc; close all;
addpath(genpath('admesh-lib-part'));
%--------------------------------------------------------------------------
% Load example file (variable "Points")
%--------------------------------------------------------------------------
ShorelineFile = ['D:\Research\PR-NSF-PREEVENTS\Data\Lower Neches\ArcGIS\',...
                'ClintExtent/ShorelineCopy.shp'];
Shoreline         = shaperead(ShorelineFile);

for i = 1 : length(Shoreline)
    Shoreline(i).np = length(Shoreline(i).X);
end
[~,index] = sortrows([Shoreline.np].'); Shoreline = Shoreline(index(end:-1:1)); clear index;
%% ========================================================================
% Parameter setting
%==========================================================================
close all;
DEG2KM = 6378*pi/180;
%%========================================================================
% Pick a shoreline and set up parameters
%==========================================================================
i = 8;
SP = [Shoreline(i).X; Shoreline(i).Y]';

% delta = 1e-5;
delta = 3/(DEG2KM*1e3);
hmax = 2*delta;
dx = 1*delta;

delta_lfs = 100/(DEG2KM*1e3); % meter to degree
delta_A = 500/(dx*DEG2KM*1e3)^2; % square meter to number of pixel
rho_pruning = 4*delta_lfs;
rho = 4*dx;
dtheta = pi*.9;

figure; plot(SP(:,1),SP(:,2),'k'); axis equal; 
%%========================================================================
% Simple test case
%==========================================================================
% delta = 0.5e-2;
% hmax = 2*delta;
% dx = 1*delta;
% 
% delta_lfs = 0.2;
% deltaA = 10;
% rho_pruning = 2*delta_lfs;
% rho = .5*delta_lfs;
% dtheta = pi*.6;


%%%%%% Rectangular channel
% SP = [0 1; 0 -1; 3 -1; 3 -delta_lfs/2; 4 -delta_lfs/2; 4 delta_lfs/2; 3 delta_lfs/2; 3 1; 0 1];

%%%%%% Curved channel
% t1 = linspace(0,pi/2)';
% t2 = linspace(pi/2,pi)';
% t3 = linspace(pi,3*pi/2)';
% t4 = linspace(3*pi/2,2*pi)';
% SP = [0.5 -1; 2.5 -1; 2.5+0.5*cos(t4), -0.5+0.5*sin(t4); 3 -delta_lfs/2; 4 -delta_lfs/2; 4 delta_lfs/2; 3 delta_lfs/2; 2.5+0.5*cos(t1), 0.5+0.5*sin(t1); 2.5 1; 0.5 1; 0.5+0.5*cos(t2), 0.5+0.5*sin(t2); 0 0.5; 0 -0.5; 0.5+0.5*cos(t3), -0.5+0.5*sin(t3);];
% % figure; plot(SP(:,1),SP(:,2),'m');

%% ========================================================================
% Compute medial axis
%==========================================================================
SPx = SP(:,1); SPy = SP(:,2);
% figure; plot(SPx,SPy,'m');
%
x = min(SPx) : dx : max(SPx);
y = min(SPy) : dx : max(SPy);
[BX,BY] = meshgrid(x,y);

%--------------------------------------------------------------------------
% Compute distance from shoreline and gradient of it
%--------------------------------------------------------------------------
PTS.Poly(1).x = SP(:,1);
PTS.Poly(1).y = SP(:,2);
[D,gradD] = SignedDistanceFunction_v2(PTS,BX,BY,delta,hmax);
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
[Vx,Vy] = Compute8SSED(BX,BY,SPx,SPy);
Vx(id_out) = 0; Vy(id_out) = 0;
Vmag = sqrt(Vx.^2 + Vy.^2);
fprintf('Computing 8SSED map is done\n');

%--------------------------------------------------------------------------
% Find medial axis
%--------------------------------------------------------------------------
% BW = double(~id_out);
% id_MA = bwmorph(BW,'skel',inf);
% id_MA = logical(id_MA);

% id_MA = GradD < 0.9 & D <= 0;
% id_MA = Skeletonization_Choi_v2(BX,BY,Vx,Vy,2*dx,rho_pruning,dtheta);
id_MA = Skeletonization_DivSSED_v2(BX,BY,Vx,Vy);
id_MA(id_out) = 0;
fprintf('Skeletonization is done\n');

figure; hold on;
plot(SPx,SPy,'k');
myPlotOnTop(BX(id_MA),BY(id_MA),'.r');
view(2); axis equal tight; box on;
%% ========================================================================
% Thinning MA
%%=========================================================================
bw = id_MA;
%--------------------------------------------------------------------------
% (???) no idea why.. but bwskel looks better...
%--------------------------------------------------------------------------
id_MA_thinned = bwskel(bw);
% id_MA_thinned = bwmorph(bw,'skel',inf);

%--------------------------------------------------------------------------
% Remove isolated pixels
%--------------------------------------------------------------------------
id_MA_thinned = bwmorph(id_MA_thinned,'clean',inf);

%--------------------------------------------------------------------------
% I didn't mean it but this part fills holes
% (it's better to merge joints latter...)
%--------------------------------------------------------------------------
[B_MA_thinned,L_MA_thinned] = bwboundaries(id_MA_thinned);
nL = zeros(max(L_MA_thinned(:)),1);
for i = 1 : max(L_MA_thinned(:))
    nL(i) = nnz(L_MA_thinned == i);
end
A_MA_thinned = nL;
id_MA_thinned = logical(L_MA_thinned > 0);
fprintf('Thinning is done\n');

% figure; hold on;
% plot(SPx,SPy,'k');
% myPlotOnTop(BX(id_MA_thinned),BY(id_MA_thinned),'.r');
% view(2); axis equal tight; box on;

%% ========================================================================
% Construct branches (with redundant joints)
%==========================================================================
MA = ConstructMedialAxisBranch(id_MA_thinned);
fprintf('Branch construction is done\n');
%--------------------------------------------------------------------------
% Plot..
%--------------------------------------------------------------------------
figure; hold on;
plot(SPx,SPy,'k');
for i = 1 : MA.nBranch
    plot(BX(MA.ID{i}),BY(MA.ID{i}),'linewidth',1.5);
end
% myPlotOnTop(BX(id_MA_joint),BY(id_MA_joint),'ob');
axis equal; box on;

% myScatter3(BX(id_MA_thinned),BY(id_MA_thinned),num_nghb_MA1D(id_MA_thinned))
% mycmap('rainbow',[-.5 4.5],5);

%% ========================================================================
% Merge redundant joints
%==========================================================================
MA_merged = MergeJoints2ClosestGridPoint(MA);
id_MA_joint = MA_merged.ID_Joint;
fprintf('Merging joints is done\n');

figure; hold on;
plot(SPx,SPy,'k');
for i = 1:MA_merged.nBranch
%     plot(MA_jointmerged.XY{i}(:,1),MA_jointmerged.XY{i}(:,2),'linewidth',1.5);
    plot(BX(MA_merged.ID{i}),BY(MA_merged.ID{i}),'r','linewidth',1.5);
%     plot(BX(MA_merged.ID{i}),BY(MA_merged.ID{i}),'linewidth',1.5);
end
% myPlotOnTop(BX(id_MA_joint),BY(id_MA_joint),'ob');
axis equal; box on;


%% ========================================================================
% Pruning corners
%==========================================================================
MA_pruned = SkeletonPruning_Choi_level1branch(BX,BY,Vx,Vy,D,MA_merged,rho_pruning,dtheta);
fprintf('Pruning corner branches is done\n');
ID = MA_pruned.ID;

figure; hold on;
plot(SPx,SPy,'k');
for i = 1 : length(ID)
    plot(BX(ID{i}),BY(ID{i}),'.r','linewidth',1.5);
%     plot(BX(ID{i}),BY(ID{i}),'linewidth',1.5);
end
% myPlotOnTop(BX(id_MA_joint),BY(id_MA_joint),'ob');
axis equal; box on;

%% ========================================================================
% Separate 1D & 2D area based on local feature size
%==========================================================================
id_MA_pruned = vertcat(MA_pruned.ID{:});

bw_MA_pruned = zeros(size(BX));
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
% Transform boundary points to actual domain
%--------------------------------------------------------------------------
% for i = 1 : length(L1D_boundary)
%     x = L1D_boundary{i}(:,2);
%     y = L1D_boundary{i}(:,1);
%     
%     b = size(BX,2); a = 1;
%     d = max(BX(:)); c = min(BX(:));
%     x = (d-c)/(b-a)*(x - a) + c;
%     b = size(BY,1); a = 1;
%     d = max(BY(:)); c = min(BY(:));
%     y = (d-c)/(b-a)*(y - a) + c;
%     
%     L1D_boundary{i} = [x,y];
% end
% 
% 
% for i = 1 : length(L2D_boundary)
%     x = L2D_boundary{i}(:,2);
%     y = L2D_boundary{i}(:,1);
%     
%     b = size(BX,2); a = 1;
%     d = max(BX(:)); c = min(BX(:));
%     x = (d-c)/(b-a)*(x - a) + c;
%     b = size(BY,1); a = 1;
%     d = max(BY(:)); c = min(BY(:));
%     y = (d-c)/(b-a)*(y - a) + c;
%     
%     L2D_boundary{i} = [x,y];
% end

%--------------------------------------------------------------------------
% put NaN for points out of domain (for visualization)
%--------------------------------------------------------------------------
% D(id_out) = nan;
% dMA_pruned(id_out) = nan;
% % GradD(id_out) = nan;
% lfs(id_out) = nan;
% lfs_2D(id_out) = nan;



%--------------------------------------------------------------------------
% Get indices for MA1D and MA2D
%--------------------------------------------------------------------------
id_MA1D = intersect(id_MA_pruned,find(L1D > 0));
id_MA2D = intersect(id_MA_pruned,find(L2D > 0));
XY_MA1D = [BX(id_MA1D),BY(id_MA1D)];
XY_MA2D = [BX(id_MA2D),BY(id_MA2D)];

%--------------------------------------------------------------------------
% Plot result
%--------------------------------------------------------------------------
% figure; hold on;
% plot(SPx,SPy,'k');
% plot(BX(id_MA),BY(id_MA),'.r');
% axis equal tight;

% figure; hold on;
% plot(SPx,SPy,'k');
% plot(BX(id_MA_pruned),BY(id_MA_pruned),'.r');
% axis equal tight;

% figure; hold on;
% surf(BX,BY,lfs,'edgecolor','interp','facecolor','interp');
% % myPlotOnTop(PtsMA_pruned(:,1),PtsMA_pruned(:,2),'.b');
% myPlotOnTop(SPx,SPy,'k');
% view(2); colorbar; axis equal; box on;
% mycmap('jet',[0 delta_lfs]);

% figure; hold on;
% myPlotOnTop(SPx,SPy,'k');
% myPlotOnTop(XY_MA1D(:,1),XY_MA1D(:,2),'.r');
% myPlotOnTop(XY_MA2D(:,1),XY_MA2D(:,2),'.b');
% view(2); colorbar; axis equal tight;
% mycmap('jet',[0 delta_lfs]);



% figure; hold on;
% % surf(BX,BY,L2D,'edgecolor','none','facecolor','interp');
% % surf(BX,BY,L1D,'edgecolor','none','facecolor',.5*[1 1 1]); view(2);
% % plot(BX(L2D > 0),BY(L2D > 0),'sw','MarkerFaceColor','w');
% plot(BX(L1D > 0),BY(L1D > 0),'s','color',.5*[1 1 1],'MarkerFaceColor',.5*[1 1 1]); view(2);
% myPlotOnTop(XY_MA1D(:,1),XY_MA1D(:,2),'.r');
% myPlotOnTop(XY_MA2D(:,1),XY_MA2D(:,2),'.b');
% myPlotOnTop(SPx,SPy,'k');
% view(2); axis equal tight; mycmap('binary',[0 1]);


fprintf('1D/2D area seperation is done\n');

%% ========================================================================
% Fill 2D area
%==========================================================================
L2D_Filled = L2D;
L1D_Filled = L1D;
for i = 1 : length(id_MA2D)
    k = id_MA2D(i);
    iL2D = L2D(k);
    
    x = BX(k); y = BY(k);
    iVmag = Vmag(k);
    
    num_nghb_buffer = ceil(iVmag/dx)+1;
    
    [I,J] = ind2sub(size(BX),k);
    
    I1 = max(I - num_nghb_buffer,1);
    I2 = min(I + num_nghb_buffer,size(BX,1));
    
    J1 = max(J - num_nghb_buffer,1);
    J2 = min(J + num_nghb_buffer,size(BX,2));
    
    I3 = I1:I2;
    J3 = J1:J2;
    [I3,J3] = meshgrid(I3,J3);
    I3 = I3(:);
    J3 = J3(:);
    
    k_nghb = sub2ind(size(BX),I3,J3);
    
%     k_nghb = k + [-N-1, -N, -N+1, -1, 1, N-1, N, N+1];
    
    tempx = BX(k_nghb);
    tempy = BY(k_nghb);
    
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

id_MA1D_Filled = intersect(id_MA_pruned, find(L1D_Filled > 0));
id_MA2D_Filled = intersect(id_MA_pruned, find(L2D_Filled > 0));
%--------------------------------------------------------------------------
% Plot result
%--------------------------------------------------------------------------
% temp = nan(size(BX));
% temp(L1D > 0) = 1; temp(L2D > 0) = 2;
% figure; hold on;
% surf(BX,BY,temp,'EdgeColor','none','FaceColor','interp');
% myPlotOnTop(SPx,SPy,'k');
% % myPlotOnTop(BX(id_MA1D),BY(id_MA1D),'.r');
% % myPlotOnTop(BX(id_MA2D),BY(id_MA2D),'.b');
% axis equal; box on;
% mycmap('gray',[0 3]);
% 
% temp = nan(size(BX));
% temp(L1D_Filled > 0) = 1; temp(L2D_Filled > 0) = 2;
% figure; hold on;
% surf(BX,BY,temp,'EdgeColor','none','FaceColor','interp');
% myPlotOnTop(SPx,SPy,'k');
% % myPlotOnTop(BX(id_MA1D_Filled),BY(id_MA1D_Filled),'.r');
% % myPlotOnTop(BX(id_MA2D_Filled),BY(id_MA2D_Filled),'.b');
% axis equal; box on; 
% mycmap('gray',[0 3]);

fprintf('Filling 2D area is done\n');
%% ========================================================================
% Get boundary of L2D
%==========================================================================
L2D_new_boundary = GetLayerBoundary(L2D_Filled,BX,BY,SPx,SPy,'Smoothing');
figure; hold on;
plot(SPx,SPy,'k');
for i = 1 : length(L2D_new_boundary)
    plot(L2D_new_boundary{i}(:,1),L2D_new_boundary{i}(:,2),'b'); 
end
axis equal; box on;

fprintf('Get boundary of 2D layer is done\n');
%% ========================================================================
% Make MA1D is connected to 2D area
%==========================================================================
xy_L2D_boundary = vertcat(L2D_new_boundary{:});

nBranch = 0;
ID = [];
Conn = [];
XY = [];
for i = 1 : length(MA_pruned.ID)
    id = MA_pruned.ID{i};
    is_MA1D = L1D_Filled(id) > 0;
    if any(is_MA1D)
        nBranch = nBranch + 1;
        id = id(is_MA1D);
        ID{nBranch} = id;
        Conn(nBranch,:) = [id(1), id(end)];
        XY{nBranch} = [BX(id),BY(id)];
    end
end
MA_connected.nBranch = nBranch;
MA_connected.ID = ID;
MA_connected.Conn = Conn;
MA_connected.XY = XY;

for i = 1 : length(MA_connected.ID)
    id = MA_connected.ID{i};
    k = L1D_Filled(id) > 0;
    id = id(k);
    if isempty(id)
        continue;
    end
    x = BX(id);
    y = BY(id);
    
    dist1 = sqrt((xy_L2D_boundary(:,1) - x(1)).^2 + (xy_L2D_boundary(:,2) - y(1)).^2);
    dist2 = sqrt((xy_L2D_boundary(:,1) - x(end)).^2 + (xy_L2D_boundary(:,2) - y(end)).^2);
    
    if any(dist1 < (2)*dx)
        [~,iend] = min(dist1);
%         MA_connected.ID{i} = vertcat(iend,MA_connected.ID{i});
        x = vertcat(xy_L2D_boundary(iend,1),x);
        y = vertcat(xy_L2D_boundary(iend,2),y);
        MA_connected.Conn(i,1) = 0;
    end
    
    if any(dist2 < (2)*dx)
        [~,iend] = min(dist2);
%         MA_connected.ID{i} = vertcat(MA_connected.ID{i},iend);
         x = vertcat(x,xy_L2D_boundary(iend,1));
         y = vertcat(y,xy_L2D_boundary(iend,2));
         MA_connected.Conn(i,2) = 0;
    end
    MA_connected.XY{i} = [x,y];    
end

%--------------------------------------------------------------------------
% Plot result..
%--------------------------------------------------------------------------
figure; hold on;
plot(SPx,SPy,'k');
for i = 1 : length(L2D_new_boundary)
    plot(L2D_new_boundary{i}(:,1),L2D_new_boundary{i}(:,2),'b'); 
end

for i = 1 : length(MA_pruned.ID)
    id = MA_pruned.ID{i};
    k = L1D_Filled(id) > 0;
    id = id(k);
    plot(BX(id),BY(id),'linewidth',1.5);
end
axis equal; box on;


figure; hold on;
plot(SPx,SPy,'k');
for i = 1 : length(L2D_new_boundary)
    myPlotOnTop(L2D_new_boundary{i}(:,1),L2D_new_boundary{i}(:,2),'b','linewidth',1.5); 
end
for i = 1 : length(MA_connected.ID)
    x = MA_connected.XY{i}(:,1);
    y = MA_connected.XY{i}(:,2);
    myPlotOnTop(x,y,'r','linewidth',1.5);
end
axis equal; box on;

fprintf('Connecting MA1D to 2D area is done\n');
%% ========================================================================
% Construct main stream
%==========================================================================
ID = MA_connected.ID;
ID1 = MA_connected.ID;
Conn = MA_connected.Conn;
Conn1 = MA_connected.Conn;
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
plot(SPx,SPy,'k');
axis equal; box on;
for ii = 1 : length(L2D_new_boundary)
    plot(L2D_new_boundary{ii}(:,1),L2D_new_boundary{ii}(:,2),'b'); 
end
for ii = 1 : length(ID_new)
    plot(BX(ID_new{ii}),BY(ID_new{ii}),'linewidth',1.5);
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
dx = mean(mean(diff(BX,1,2)));
dy = mean(mean(diff(BY,1,1)));
dx = max(dx,dy);
for i = 1 : max(L2D(:))
    id = (L2D(id_MA2D) == i);
    id_MA2D_i = id_MA2D(id);
    for j = 1 : length(id_MA2D_i)
        k = id_MA2D_i(j);
        x = BX(k); y = BY(k);
        iVmag = Vmag(k);
        
        num_nghb_buffer = ceil(iVmag/dx)+1;
        
        [I,J] = ind2sub(size(BX),k);
        
        I1 = max(I - num_nghb_buffer,1);
        I2 = min(I + num_nghb_buffer,size(BX,1));
        
        J1 = max(J - num_nghb_buffer,1);
        J2 = min(J + num_nghb_buffer,size(BX,2));
        
        I3 = I1:I2;
        J3 = J1:J2;
        [I3,J3] = meshgrid(I3,J3);
        I3 = I3(:);
        J3 = J3(:);
        
        k_nghb = sub2ind(size(BX),I3,J3);
        
        %     k_nghb = k + [-N-1, -N, -N+1, -1, 1, N-1, N, N+1];
        
        tempx = BX(k_nghb);
        tempy = BY(k_nghb);
        
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
    
    x = BX(k); y = BY(k);
    iVmag = Vmag(k);
    
    num_nghb_buffer = ceil(iVmag/dx)+1;
    
    [I,J] = ind2sub(size(BX),k);
    
    I1 = max(I - num_nghb_buffer,1);
    I2 = min(I + num_nghb_buffer,size(BX,1));
    
    J1 = max(J - num_nghb_buffer,1);
    J2 = min(J + num_nghb_buffer,size(BX,2));
    
    I3 = I1:I2;
    J3 = J1:J2;
    [I3,J3] = meshgrid(I3,J3);
    I3 = I3(:);
    J3 = J3(:);
    
    k_nghb = sub2ind(size(BX),I3,J3);
    
%     k_nghb = k + [-N-1, -N, -N+1, -1, 1, N-1, N, N+1];
    
    tempx = BX(k_nghb);
    tempy = BY(k_nghb);
    
    dist = sqrt( (tempx - x).^2 + (tempy - y).^2);
    
    id = (dist < iVmag);
    id = k_nghb(id);
    id(id_out(id)) = [];
%     I4 = I3(id);
%     J4 = J3(id);
    
    L2D_Filled(id) = iL2D;
    L1D_Filled(id) = 0;
end

temp = nan(size(BX));
temp(L1D > 0) = 1; temp(L2D > 0) = 2;
figure; hold on;
surf(BX,BY,temp,'EdgeColor','none','FaceColor','interp');
myPlotOnTop(SPx,SPy,'k');
axis equal; box on; 
mycmap('gray',[0 3]);



temp = nan(size(BX));
temp(L1D_Filled > 0) = 1; temp(L2D_Filled > 0) = 2;
figure; hold on;
surf(BX,BY,L2D_Filled,'EdgeColor','none','FaceColor','interp');
myPlotOnTop(SPx,SPy,'k');
axis equal; box on;  axis(ax);
mycmap('gray',[0 3]);

%% Find boundary
[I,J] = ind2sub(size(BX),find(id_2D_fill));
temp1 = sub2ind(size(BX),min(I+1,size(BX,1)),J);
temp2 = sub2ind(size(BX),max(I-1,1),J);
temp3 = sub2ind(size(BX),I,min(J+1,size(BX,2)));
temp4 = sub2ind(size(BX),I,max(J-1,1));

temp = vertcat(temp1,temp2,temp3,temp4);
id_2D_expand = false(size(BX));
id_2D_expand(temp) = true;
figure; plot(BX(id_2D_expand),BY(id_2D_expand),'.');

id = inpolygon(SPx,SPy,BX(temp),BY(temp));
figure; plot(SPx(id),SPy(id),'.r');

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
plot(SPx,SPy,'k');
for i = 1:nnz(iBranchLong)
    plot(MA.XY_new{i}(:,1),MA.XY_new{i}(:,2),'linewidth',1.5);
end
myPlotOnTop(BX(id_MA_joint),BY(id_MA_joint),'ob','linewidth',2);
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
plot(SPx,SPy,'k');
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

