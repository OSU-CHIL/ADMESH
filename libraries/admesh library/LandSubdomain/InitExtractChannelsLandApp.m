function InitExtractChannelsLandApp(app, mainapp)

% Assign main app variables
app.MainApp = mainapp;
app.PTS = mainapp.PTS;
app.UIAxes = mainapp.UIAxes;
app.ProgressBarButton = mainapp.ProgressBarButton;
app.ElevationDataFilename = mainapp.ElevationDataFilename;
app.xyzFun = mainapp.xyzFun;
app.MinDrainageAreaEditField.Value = mainapp.MinDrainageArea;

% Move subapp to center of mainapp
mainP = mainapp.UIFigure.Position;
subP = app.UIFigure.Position;
subP(1) = mainP(1) + 0.5*(mainP(3) - subP(3));
subP(2) = mainP(2) + 0.5*(mainP(4) - subP(4));
app.UIFigure.Position = subP;

% Check if elevatio data file is assigned
CheckElevationDataFile(app);