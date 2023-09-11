function [res,hmax,hmin,K,R,s,C,T,g] = DealVariables(app)
% AssignAdmeshVariables - Function that assigns variables from the GUI
% workspace
%
% Syntax:  [res,hmax,hmin,K,R,xyz,s,C,T,g] = AssignAdmeshVariables(guiFig)
%
% Inputs:
%    guiFig - handle that identifies the figure
%
% Outputs:
%    res    - a factor multiple, either 1,2,or 3, of hmin controlling the grid spacing.
%    hmin   - Minimum element size
%    hmax   - Maximum element size
%    K      - parameter specifying the number of elements per radian of curvature
%    R      - parameter specifying the number of elements that span half a channel/tributary
%    s      - parameter specifying the density of elements in areas with large changes in elevation
%    C      - parameter specifying the number of elements per tidal wavelength
%    T      - frequency of wave propogation
%    g      - grading limiting
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% August 2013; Last revision: 08-August-2013

%------------- BEGIN CODE --------------

%--------------------------------------------------------------------------
% Assign maximum & minimum values
%--------------------------------------------------------------------------
hmax = app.MaxElementSizeEditField.Value;
hmin = app.MinElementSizeEditField.Value;

%--------------------------------------------------------------------------
% Assign Resolution Parameter
%--------------------------------------------------------------------------
res = app.SpatialResolutionEditField.Value;

%--------------------------------------------------------------------------
% Assign Curvature Parameters
%--------------------------------------------------------------------------
if strcmpi(app.BoundaryCurvatureDropDown.Value,'on')
    K   = app.BoundaryCurvatureEditField.Value;
else
    K   = [];
end

%--------------------------------------------------------------------------
% Assign Local Feature Size Parameters
%--------------------------------------------------------------------------
if strcmpi(app.LocalFeatureSizeDropDown.Value,'on')
    R   = app.LocalFeatureSizeEditField.Value;
else
    R   = [];
end

%--------------------------------------------------------------------------
% Assign Bathymetry Parameters
%--------------------------------------------------------------------------
% Check for bathymetry on and ensure bathymetry data exists
GI = isa(app.xyzFun, 'griddedInterpolant');
SI = isa(app.xyzFun, 'scatteredInterpolant');

if strcmpi(app.ElevationGradientsDropDown.Value,'on') && (GI || SI)
    s       = app.ElevationGradientsEditField.Value;
else
    s       = nan;
end

%--------------------------------------------------------------------------
% Assign Tidal Wavelength Parameters
%--------------------------------------------------------------------------
% Check for tide on and ensure bathymetry data exists
if ~strcmpi(app.DominateTideDropDown.Value,'off') && GI
    
    C = app.DominateTideEditField.Value;
    
    switch app.DominateTideDropDown.Value
        
        case 'M2'
            
            T = 12.42*3600;
            
        case 'S2'
            
            T = 12.00*3600;
            
        case 'N2'
            
            T = 12.658*3600;
            
        case 'K2'
            
            T = 11.967*3600;
            
    end
    
else
    C = [];
    T = [];
end

%--------------------------------------------------------------------------
% Assign Grading Parameter
%--------------------------------------------------------------------------
g = app.MeshGradingEditField.Value;

end