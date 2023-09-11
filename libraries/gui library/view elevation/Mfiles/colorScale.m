function colorScale(x,y,z,cmap)

% Limits
xlim = [min(x(:)),max(x(:))];
ylim = [min(y(:)),max(y(:))];
zlim = [min(z(:)),max(z(:))];
dz = abs(diff(zlim));

% Border width
bw = 1;
xyr = 1;

bwy = bw*diff(ylim)/100; % Y border width = 1%
bwx = bwy/xyr; % border width (in degree of longitude)

% Scale Legend
wsc = bwx;
xsc = xlim(2) + wsc*8;

zscale = linspace(zlim(1),zlim(2),length(cmap));

yscale = linspace(0,diff(ylim)/2,length(cmap));

%ysc = ylim(2)/2;

%wsc = bwy;
ysc = ylim(1) + wsc*8;

ddz = dtick(dz*max(0.5*xyr*diff(xlim)/yscale(end),1));

ztick = (ddz*ceil(zscale(1)./ddz)):ddz:zscale(end);

h = patch(xsc + repmat(wsc*[-1;1;1;-1],[1,length(cmap)]), ...
    ysc + [repmat(yscale,[2,1]);repmat(yscale + diff(yscale(1:2)),[2,1])], ...
    repmat(zscale,[4,1]), ...
    'EdgeColor','flat','LineWidth',.1,'FaceColor','flat','clipping','off');

set(h,'tag','Colorbar')
uistack(h,'top')

colormap(cmap);

caxis([zlim(1),zlim(2)]);

h = patch(xsc + wsc*[-1,1,1,-1],ysc + yscale(end)*[0,0,1,1],'k','FaceColor','none',...
    'Clipping','off');

set(h,'tag','Colorbar');
uistack(h,'top');

h = text(xsc + 2*wsc + zeros(size(ztick)),ysc + (ztick - zscale(1))*0.5*diff(ylim)/...
    diff(zscale([1,end])),num2str(ztick'), ...
    'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',10);

set(h,'tag','Colorbar');
uistack(h,'top');

% indicates min and max Z values
h = text(xsc,ysc - bwy/2,sprintf('%g m',roundsd(zlim(1),3)),'FontWeight','bold', ...
    'HorizontalAlignment','left','VerticalAlignment','top','FontSize',10);

set(h,'tag','Colorbar');
uistack(h,'top');

h = text(xsc,ysc + .5*diff(ylim) + bwy/2,sprintf('%g m',roundsd(zlim(2),3)),'FontWeight','bold', ...
    'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',10);

set(h,'tag','Colorbar');
uistack(h,'top');

    function dd = dtick(dlim,deg)
        %DTICK Tick intervals
        
        if nargin < 2
            deg = 0;
        end
        
        if deg && dlim <= 2/60
            % less than 2 minutes: base 36
            m = 10.^floor(log10(dlim*36))/36;
        elseif deg && dlim <= 2
            % less than 2 degrees: base 6
            m = 10.^floor(log10(dlim*6))/6;
        else
            % more than few degrees or not degrees: decimal rules
            m = 10.^floor(log10(dlim));
        end
        p = ceil(dlim/m);
        if p <= 1
            dd = .1*m;
        elseif p == 2
            dd = .2*m;
        elseif p <= 5
            dd = .5*m;
        else
            dd = m;
        end
    end

    function y=roundsd(x,n)
        
        og = 10.^(floor(log10(abs(x)) - n + 1));
        y = round(x./og).*og;
        y(x==0) = 0;
    end

end