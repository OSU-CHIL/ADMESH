function Z = CreateElevationGrid(X,Y,xyzFun,Settings,ProgressBar)


%------------------------------------------------------------------------------
% Check for On/Off & Available data
%------------------------------------------------------------------------------
OnOff = strcmp(Settings.B.Status,'On') || strcmp(Settings.T.Status,'On');

if OnOff && ~isempty(xyzFun) % If true, interpolate bathymetry to grid
    
    ProgressBar.Text = 'Interpolating bathymetry onto background mesh...'; drawnow;
    
    %---------------------------------------------------------------------------
    % Interpolate bathymetry (xyzFun) onto background mesh (X,Y)
    %---------------------------------------------------------------------------   
    Z = xyzFun(X',Y')';
        
    % Interpolate NaN's if they exist
    if any(isnan(Z(:)))
        Z = inpaint_nans(Z);
    end

else
    
    Z = [];

end


end



