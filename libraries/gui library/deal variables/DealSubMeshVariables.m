function [PTS,xyzFun,res,hmax,hmin,K,R,s,C,T,g] = DealSubMeshVariables(varargin)
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
guiFig = findobj('Tag','ADmesh Figure'); guiH = % guidata(guiFig);

%--------------------------------------------------------------------------
% Assign PTS data structure
%--------------------------------------------------------------------------
PTS = guiH.subPTS;

%--------------------------------------------------------------------------
% Assign maximum & minimum values
%--------------------------------------------------------------------------
hmax = str2double(strtrim(get(guiH.MaxElementSize,'string')));

hmin = str2double(strtrim(get(guiH.MinElementSize,'string')));

%--------------------------------------------------------------------------
% Assign Resolution Parameter
%--------------------------------------------------------------------------

if get(guiH.HighRes,'value')
    res = 5;
elseif get(guiH.MedRes,'value')
    res = 3;
else
    res = 1;
end

%--------------------------------------------------------------------------
% Assign Curvature Parameters
%--------------------------------------------------------------------------

if get(guiH.CurvatureOn,'value')
    K   = str2double(strtrim(get(guiH.Curvature,'string')));
else
    K   = [];
end

%--------------------------------------------------------------------------
% Assign Local Feature Size Parameters
%--------------------------------------------------------------------------

if get(guiH.LFSOn,'value')
    R   = str2double(strtrim(get(guiH.LFS,'string')));
else
    R   = [];
end

%--------------------------------------------------------------------------
% Assign Bathymetry Parameters
%--------------------------------------------------------------------------
% Check for bathymetry on and ensure bathymetry data exists
SI = isa(guiH.xyzFun, 'scatteredInterpolant');
GI = isa(guiH.xyzFun, 'griddedInterpolant');
if get(guiH.BathymetryOn,'value') && (SI || GI)
    
    s       = str2double(strtrim(get(guiH.BathymetryValue,'string')));
    xyzFun  = guiH.xyzFun;
    
else
    s       = nan;
    xyzFun  = guiH.xyzFun;
end

%--------------------------------------------------------------------------
% Assign Tidal Wavelength Parameters
%--------------------------------------------------------------------------
% Check for tide on and ensure bathymetry data exists
if ~get(guiH.TideOff,'value') && (SI || GI)
    
    C = str2double(strtrim(get(guiH.TideValue,'string')));
    
    if get(guiH.TideM2,'value')
        
        T = 12.42*3600;
        
    elseif get(guiH.TideS2,'value')
        
        T = 12.00*3600;
        
    elseif get(guiH.TideN2,'value')
        
        T = 12.658*3600;
        
    elseif get(guiH.TideK2,'value')
        
        T = 11.967*3600;
        
    end
    
else
    C = [];
    T = [];
end

%--------------------------------------------------------------------------
% Assign Grading Parameter
%--------------------------------------------------------------------------
g    = str2double(strtrim(get(guiH.GradingValue,'string')));

end