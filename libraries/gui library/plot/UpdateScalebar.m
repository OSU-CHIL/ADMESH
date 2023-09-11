function UpdateScalebar(app)

if strcmpi(app.CoordinateSystemDropDown.Value,'Projected (m)')
    UnitScale = 1;
elseif strcmpi(app.CoordinateSystemDropDown.Value,'Unprojected (decimal degree)')
    UnitScale = deg2km(1e3); % deg2m
end
scalebar('hAxes',app.UIAxes,'Location','southeast','Bold',1,'Unit','m','UnitScale',UnitScale);