function L2D_new_boundary = GetLayerBoundary_v1(L,X,Y,BX,BY,option)
% L : layer
% X,Y background grid
% BX,BY : boundary points
%

Debugging = 0;
nL = max(L(:));
L2D_new_boundary = cell(nL,1);
dx = mean(mean(diff(X,1,2)));
dy = mean(mean(diff(Y,1,1)));
dx = max(dx,dy);

for iL = 1 : nL
    iL2D = L == iL;
    if strcmpi(option,'Fill-Cut')
        t = 1 : length(BX);
        t2 = 1;
        dist = sqrt((BX(1:end-1)-BX(2:end)).^2 + (BY(1:end-1)-BY(2:end)).^2);
        
        for i = 1 : length(BX)-1
            n = max(2,round(dist(i)/dx));
            t2 = [t2, linspace(t2(end),t2(end)+1,n)];
        end
        t2 = unique(t2);
        t2 = t2(:);
        
        SPx2 = interp1(t,BX,t2);
        SPy2 = interp1(t,BY,t2);
        
        x = X(iL2D); y = Y(iL2D);
        [~,dist] = knnsearch([x,y],[SPx2,SPy2]);
        I = dist < sqrt(2)*dx;
        x1 = vertcat(X(iL2D),SPx2((I)));
        y1 = vertcat(Y(iL2D),SPy2((I)));
        iB = boundary(x1,y1,1);
        bx = x1(iB); by = y1(iB);
    elseif strcmpi(option,'Smoothing')
        B = bwboundaries(double(iL2D));
        if length(B) > 1
            warning('it might be a big issue... just ignore now');
%             B = B{1};
        end
        
        x = B{1}(:,2);
        y = B{1}(:,1);
        
        b = size(X,2); a = 1;
        d = max(X(:)); c = min(X(:));
        x = (d-c)/(b-a)*(x - a) + c;
        b = size(Y,1); a = 1;
        d = max(Y(:)); c = min(Y(:));
        y = (d-c)/(b-a)*(y - a) + c;


        SmoothingParam = 0.1;
        id_fixedPoints = [1,length(x)];
        t = 1 : length(x);
        w = ones(length(x),1);
        w(id_fixedPoints) = 1e8;
        Smooth_x = csaps(t,x,SmoothingParam,t,w);
        Smooth_y = csaps(t,y,SmoothingParam,t,w);
        Smooth_x(id_fixedPoints) = x(id_fixedPoints);
        Smooth_y(id_fixedPoints) = y(id_fixedPoints);
        bx = Smooth_x; by = Smooth_y;
    end
    
    L2D_new_boundary{iL} = [bx(:),by(:)];
end


if Debugging == 1
    figure; hold on;
    plot(SPx2,SPy2,'k');
    % myPlotOnTop(x,y,'.k');
    plot(L2D_new_boundary(:,1),L2D_new_boundary(:,2),'b')
    axis equal; box on;
end














