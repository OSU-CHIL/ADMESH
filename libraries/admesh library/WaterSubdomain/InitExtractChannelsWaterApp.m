function InitExtractChannelsWaterApp(app, mainapp)

% Assign main app variables
app.MainApp = mainapp;
app.UIAxes = mainapp.UIAxes;
app.ProgressBarButton = mainapp.ProgressBarButton;
app.ElevationDataFilename = mainapp.ElevationDataFilename;
app.PTS = mainapp.PTS;
app.xyzFun = mainapp.xyzFun;

% Move subapp to center of mainapp
mainP = mainapp.UIFigure.Position;
subP = app.UIFigure.Position;
subP(1) = mainP(1) + 0.5*(mainP(3) - subP(3));
subP(2) = mainP(2) + 0.5*(mainP(4) - subP(4));
app.UIFigure.Position = subP;

% Check if edge structure data is loaded
if isempty(app.PTS)
    msg = 'You must first load an edge structure file (.shp, .mat, or .14) file.';
    [~] = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');
    delete(app);
    return;
end

% Set buttons on the app
CheckButtonsExtractChannelsWaterApp(app);

