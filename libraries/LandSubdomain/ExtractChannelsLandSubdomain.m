function ExtractChannelsLandSubdomain(app,event)

% First call
if strcmpi(event,'firstcall')
    filename = app.ElevationDataFilename;
    if isempty(filename) || isempty(app.xyzFun)
        warndlg('No elevation data is read.','Error');
        return;
    end
    [~,~,ext] = fileparts(filename);
    
    if ~any(strcmpi(ext,{'.tif','.tiff'}))
        warndlg('The elevation data should be .tif or .tiff format.','Error');
    end
    
    PTS = app.PTS;
    if ~isempty(PTS.Constraints) && any([PTS.Constraints(:).num] == 18)
        choice = questdlg(...
            'Open channel constraints are found. Do you want to replace them?',...
            'ADMESH',...
            'Yes','No','No');
        
        if strcmpi(choice,'no')
            return;
        end
        
        I = [PTS.Constraints(:).num] == 18;
        PTS.Constraints(I) = [];
        app.PTS = PTS;
        
        PlotEdgeStructure(app,.1);
    end
    
    app.ProgressBarButton.Text = 'Extracting open channels from DEM...'; drawnow;
    [FD,A] = ExtractOpenChannels(filename);
    app.TTC.FD = FD;
    app.TTC.A = A;
end

if isempty(app.TTC)
    % Do nothing and return if not called with button first
    return;
end
    
FD = app.TTC.FD;
A = app.TTC.A;

app.ProgressBarButton.Text = 'Extracting open channels from DEM...'; drawnow;
warnStruct = warning;
warning('off');
FL1 = klargestconncomps(STREAMobj(FD,A > app.kEditField.Value));
warning(warnStruct);

if isempty(FL1)
    warndlg(['There are only 0 connected components in the stream network. ',...
        'Use smaller threshold to get stream network.'],'Error');
    app.TTC.FL = {};
    
    h = findobj(app.ViewAxes,'tag','landchannel');
    delete(h);
    return;
end

BP = [app.PTS.Poly.x(:), app.PTS.Poly.y(:)];
[sx,sy] = STREAMobj2XY(FL1);
FL = NaNdlm2struct([sx,sy],'Boundary',BP);

PTS = app.PTS;
n = length(PTS.Constraints);
for i = 1 : length(FL)
    PTS.Constraints(n+i).num = -18;
    PTS.Constraints(n+i).xy = FL{i};
    PTS.Constraints(n+i).type = 'line';
    PTS.Constraints(n+i).data = [];
    PTS.Constraints(n+i).Kappa = [];
end
app.PTS = PTS;

% Visualize
app.ProgressBarButton.Text = 'Draw extracted open channels...'; drawnow;
FL2 = cellfun(@(x) [x; nan(1,2)],FL,'UniformOutput',0);
FL2 = vertcat(FL2{:});

h = findobj(app.ViewAxes,'tag','landchannel');
delete(h);
h = plot(app.ViewAxes,FL2(:,1),FL2(:,2),'b');

set(h,'tag','landchannel');

app.ProgressBarButton.Text = 'Ready'; drawnow;

