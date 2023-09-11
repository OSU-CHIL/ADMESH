function [res,hmax,hmin,K,R,s,C,T,g] = DealVariables(fig)
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
% Get GUI Data
%--------------------------------------------------------------------------
gui = % guidata(fig);

%--------------------------------------------------------------------------
% Assign maximum & minimum values
%--------------------------------------------------------------------------
hmax = str2double(strtrim(get(gui.MaxElementSize,'string')));
hmin = str2double(strtrim(get(gui.MinElementSize,'string')));

%--------------------------------------------------------------------------
% Assign Resolution Parameter
%--------------------------------------------------------------------------
res = get(gui.ResolutionStatus,'value');

%--------------------------------------------------------------------------
% Assign Curvature Parameters
%--------------------------------------------------------------------------
if get(gui.CurvatureStatus,'value') == 1
    K   = str2double(strtrim(get(gui.CurvatureValue,'string')));
else
    K   = [];
end

%--------------------------------------------------------------------------
% Assign Local Feature Size Parameters
%--------------------------------------------------------------------------
if get(gui.LFSStatus,'value') == 2
    R   = str2double(strtrim(get(gui.LFSValue,'string')));
else
    R   = [];
end

%--------------------------------------------------------------------------
% Assign Bathymetry Parameters
%--------------------------------------------------------------------------
% Check for bathymetry on and ensure bathymetry data exists
GI = isa(gui.xyzFun, 'griddedInterpolant');
SI = isa(gui.xyzFun, 'scatteredInterpolant');


if (get(gui.ElevStatus,'value') == 2) && (GI || SI)
    s       = str2double(strtrim(get(gui.ElevValue,'string')));
else
    s       = nan;
end

%--------------------------------------------------------------------------
% Assign Tidal Wavelength Parameters
%--------------------------------------------------------------------------
% Check for tide on and ensure bathymetry data exists
if (get(gui.TidalStatus,'value') > 1) && GI
    
    C = str2double(strtrim(get(gui.TidalValue,'string')));
    
    switch get(gui.TidalStatus,'value')
        
        case 2 % M2
            
            T = 12.42*3600;
            
        case 3 % S2
            
            T = 12.00*3600;
            
        case 4 % N2
            
            T = 12.658*3600;
            
        case 5 % K2
            
            T = 11.967*3600;
            
    end
    
else
    C = [];
    T = [];
end

%--------------------------------------------------------------------------
% Assign Grading Parameter
%--------------------------------------------------------------------------
g = str2double(strtrim(get(gui.GradingValue,'string')));

end