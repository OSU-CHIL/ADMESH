function h0 = Dominate_tide(h0,Tide_Period, tide_value,Z,sz,hmax,hmin,Settings)
% Dominate_tide - Compute the mesh size based on dominate tides
%
% Syntax:  [h_tide] = Dominate_tide(Tide_Period, tide_value,B)
%
% Inputs:
%    Tide_Period - 
%    tide_value - number of elements per tidal wavelength
%    B - bathymetry interpolated onto background grid
%    delta - Grid spacing
%
% Outputs:
%    h_tide - mesh size based on tidal wavelength
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Colton Conroy Dustin West
% The Ohio State University
% email address: conroy.51@osu.edu, dww.425@gmail.com
% August 2013; Last revision: 08-August-2013

%------------- BEGIN CODE --------------

%------------------------------------------------------------------------------
% Elevation Status
%------------------------------------------------------------------------------
Status = Settings.T.Status;

if strcmp(Status,'On') && ~isempty(Z)
    
    g = 9.81; % gravity (m/s^2)
    
    h_tide = (Tide_Period/tide_value).*sqrt(g.*abs(Z));
    
    h_tide(h_tide == 0) = hmax;
    
    % Enforce boundary conditions
    h_tide(h_tide < hmin) = hmin;
    h_tide(h_tide > hmax) = hmax;

    % Compare initial conditions and save
    h0 = min(h_tide, h0);
    
end

end