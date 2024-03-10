function varargout = Cart2Geo(varargin)
% Meters2Geo - Converts given XY in Spherical Mercator 
% EPSG:900913 to lat/lon in WGS84 Datum
%
% Syntax:  [lon,lat] = Meters2Geo(x,y)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Reference: http://www.maptiler.org/google-maps-coordinates-tile-bounds-projection/
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
%---------------------------- BEGIN CODE --------------------------------------

R = 6378137; % Radius of the Earth

if nargin == 1 % Convert PTS or MESH data structure
    
    input = varargin{1}; % Assign input

    if isfield(input,'Points') % Mesh structure input
        
        input.Points(:,1) = (input.cpplon + input.Points(:,1)./(R*cos(input.cpplat))).*180/pi;
        input.Points(:,2) = (input.Points(:,2)./R).*(180/pi);
        
    else % PTS structure input
        
        % Convert edge structure
        for i = 1:length(input.Poly)
            input.Poly(i).x = (input.cpplon + input.Poly(i).x./(R*cos(input.cpplat))).*180/pi;
            input.Poly(i).y = (input.Poly(i).y./R).*(180/pi);
            
        end
        
        % Convert constraints
        if isfield(input,'Constraints')
            
            for i = 1:length(input.Constraints)
                input.Constraints(i).xy(:,1) = (input.cpplon + input.Constraints(i).xy(:,1)./(R*cos(input.cpplat))).*180/pi;
                input.Constraints(i).xy(:,2) = (input.Constraints(i).xy(:,2)./R).*(180/pi);
            end
            
        end
        
    end
    
    varargout{1} = input;
    
elseif nargin == 3 % Convert (x,y)
    
    [PTS,x,y] = deal(varargin{:});
    
    % Convert points
    lon = (PTS.cpplon + x./(R*cos(PTS.cpplat))).*180/pi;
    lat = (y./R).*(180/pi);
    
    varargout{1} = lon;
    varargout{2} = lat;
    
end