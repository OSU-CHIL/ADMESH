function ZoomToBadElements(app)

%------------------------------------------------------------------------------
% Get GUI data
%------------------------------------------------------------------------------
% fig         = varargin{1};
% gui         = % guidata(fig);
% pH          = findobj('Tag','ADMESH Plot Window');
% meshPatch   = findobj('Tag','Mesh');
pH = app.UIAxes;
meshPatch = app.MESH;

%------------------------------------------------------------------------------
% Check if there is anything currently in the plot window
%------------------------------------------------------------------------------
plotItems = get(pH, 'Children');

if isempty(plotItems) || isempty(meshPatch)

    return
    
end

% Turn off plot hover function temporarily
% set(app.Window,'WindowButtonMotionFcn' , ''); 


%------------------------------------------------------------------------------
% Check for low quality elements, zoom to the first one we find
%------------------------------------------------------------------------------

% Get connectivity
t = app.MESH.ConnectivityList;

% Get vertices
p = app.MESH.Points;

% Get the element quality
[~, ~, EQ] = MeshQuality(p,t,0,'Triangle');

% Find the first bad element
EQi = find(EQ < app.MinEQ,1,'first');

% Get axis settings
xLIM = app.xLimits; yLIM = app.yLimits;

% Get current axis limits
pHxLIM = app.UIAxes.XLim;
pHyLIM = app.UIAxes.YLim;

% If there are no more low qualiy elements, inform user
if isempty(EQi)
    msgbox(['No more elements fall under the minimum element criteria set at ' num2str(app.MinEQ)],'ADMESH');
    return;
end

% Do we need to zoom out?

if ~isequal([xLIM yLIM], [pHxLIM, pHyLIM]) % Zoom out
    
    bbox = [pHxLIM pHyLIM];
    
    if size(t,1) > 100000
        dt = 1;
    else
        dt = 10;
    end

    dxmin = linspace(bbox(1),xLIM(1),dt);
    dxmax = linspace(bbox(2),xLIM(2),dt);
    dymin = linspace(bbox(3),yLIM(1),dt);
    dymax = linspace(bbox(4),yLIM(2),dt);
    
    for i = 1:dt
        app.UIAxes.XLim = [dxmin(i) dxmax(i)];
        app.UIAxes.YLim = [dymin(i) dymax(i)];
        drawnow
    end
    
end

% Zoom in

% Get element coordinates
x = p(t(EQi,:),1);
y = p(t(EQi,:),2);

% Compute bounding box to zoom into
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);

offsetx = (xmax-xmin)*5;
offsety = (ymax-ymin)*5;

offset = max(offsetx, offsety );

bbox = [xmin-offset xmax+offset ymin-offset ymax+offset];

    if size(t,1) > 100000
        dt = 1;
    else
        dt = 10;
    end

dxmin = linspace(xLIM(1),bbox(1),dt);
dxmax = linspace(xLIM(2),bbox(2),dt);
dymin = linspace(yLIM(1),bbox(3),dt);
dymax = linspace(yLIM(2),bbox(4),dt);


for i = 1:dt
    app.UIAxes.XLim = [dxmin(i) dxmax(i)];
    app.UIAxes.YLim = [dymin(i) dymax(i)];
    drawnow
end

drawnow

% Turn on coordinate display. 
% set(app.Window,'WindowButtonMotionFcn' , @CoordDisplay);

end