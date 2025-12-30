function CheckButtonsExtractChannelsWaterApp(app)

% Check if edge structure data is loaded
if isempty(app.PTS)
    msg = 'You must first load an edge structure file (.shp, .mat, or .14) file.';
    [~] = uiconfirm(app.UIFigure,msg,'ADMESH',...
        'Options',{'OK'},'DefaultOption',1,'Icon','Error');
    delete(app);
    return;
end

if ~isempty(app.pgon) && ~isempty(app.pgon_ext)
    app.SwitchLandWaterMaskButton.Enable = 'on';
    app.ExtractInternalConstraintsButton.Enable = 'on';
    if ~isempty(app.Constraints)
        app.WriteADMESHinputfileButton.Enable = 'on';
    end
end

