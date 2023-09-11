function [varargout] = GetElevationRGB(varargin)
%
% Modified from DEM Shaded relief image plot created by
%
%	Author: François Beauducel <beauducel@ipgp.fr>
%	Created: 2007-05-17
%	Updated: 2013-01-08

%	Copyright (c) 2013, François Beauducel, covered by BSD License.
%	All rights reserved.
%
%	Redistribution and use in source and binary forms, with or without
%	modification, are permitted provided that the following conditions are
%	met:
%
%	   * Redistributions of source code must retain the above copyright
%	     notice, this list of conditions and the following disclaimer.
%	   * Redistributions in binary form must reproduce the above copyright
%	     notice, this list of conditions and the following disclaimer in
%	     the documentation and/or other materials provided with the distribution
%
%	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%	ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%	POSSIBILITY OF SUCH DAMAGE.


%-----------------------------------------------------------------------
% Display elevation for mesh
%-----------------------------------------------------------------------
if nargin == 4
    
    % Distribute inputs
    [x,y,xyzFun,ColorMapping] = deal(varargin{:});

    % If no colormapping is selected or no bathymetry exists
    if (isnumeric(xyzFun) || strcmp(ColorMapping,'None'))
        
        % Assign coloring for element faces
        FaceVertexCData = ones(numel(x),1)*[0 0 1];
        
        % Create empty variable for colormap
        cmap = [];
        
        varargout{1} = FaceVertexCData;
        varargout{2} = cmap;
        
        return
        
    end
    
    % Obtain elevation at each of the new nodal points
    z = xyzFun(x,y);
    
    % Somtimes NaN's occur. Interpolate those.
    z = inpaint_nans(z);
    
    if strcmp(ColorMapping,'Jet')
        
        cmap = jet;
        
        % Determine zmin & zmax
        zcut = 0.5;
        zmin = nmedian(z(~isnan(z)),zcut/100);
        zmax = nmedian(z(~isnan(z)),1 - zcut/100);
        dz = zmax - zmin;
        
    elseif strcmp(ColorMapping,'Land & Sea')
        
        % load grayscale values
        csea = seacolor(256);
        cmap = landcolor(256).^1.3;
        
        % Determine zmin & zmax
        zcut = 0.5;
        zmin = nmedian(z(~isnan(z)),zcut/100);
        zmax = nmedian(z(~isnan(z)),1 - zcut/100);
        dz = zmax - zmin;
        
        % builds the colormap: concatenates seacolor and landcolor around 0
        if zmin < 0 && zmax > 0
            l = size(csea,1);
            r = size(cmap,1)*abs(zmin)/zmax/l;
            cmap = cat(1,interp1(1:l,csea,linspace(1,l,round(l*r)),'*linear'),cmap);
            
        elseif zmin < 0 && zmax < 0
            
            cmap = csea;
        end
        
    end
    
    % normalisation of Z using CMAP and convertion to RGB
    FaceVertexCData = ind2rgb(uint16((z - zmin)*(length(cmap)/dz)),cmap);
    
    % Reshape into an mx3 matrix
    FaceVertexCData = reshape(FaceVertexCData,size(FaceVertexCData,1),3);
    
    varargout{1} = FaceVertexCData;
    varargout{2} = cmap;
    
end

%-----------------------------------------------------------------------
% Display elevation with edge structure
%-----------------------------------------------------------------------
if nargin == 3
    
    [xyzFun,PTS,ColorMapping] = deal(varargin{:});
    
    % If no colormapping is selected or no bathymetry exists
    if (isnumeric(xyzFun) || strcmp(ColorMapping,'None'))
        
        % Assign coloring for element faces
        FaceVertexCData = ones(numel(x),1)*[0 0 1];
        
        % Create empty variable for colormap
        cmap = [];
        
        varargout{1} = FaceVertexCData;
        varargout{2} = cmap;
        
        return
        
    end
    
    %-----------------------------------------------------------------------
    % Return colormapping for gridded or scattered data set
    %-----------------------------------------------------------------------
    
    % Get Polygon limits for boundaing box
    xmin = min(vertcat(PTS.Poly(:).x));
    xmax = max(vertcat(PTS.Poly(:).x));
    ymin = min(vertcat(PTS.Poly(:).y));
    ymax = max(vertcat(PTS.Poly(:).y));
    
    xv = linspace(xmin,xmax,round(1000*(xmax-xmin)/(ymax-ymin)));
    yv = linspace(ymin,ymax,round(1000*(ymax-ymin)/(xmax-xmin)));
    
    [xq,yq] = meshgrid(xv,yv);

    z = xyzFun(xq',yq')';

    AlphaData = PointsInDomain(xq,yq,PTS);
  
    clear xq yq
    
    if strcmp(ColorMapping,'Jet')
        
        cmap = jet;
        
        % Determine zmin & zmax
        zcut = 0.5;
        zmin = nmedian(z(~isnan(z)),zcut/100);
        zmax = nmedian(z(~isnan(z)),1 - zcut/100);
        dz = zmax - zmin;
        
    elseif strcmp(ColorMapping,'Land & Sea')
        
        % load grayscale values
        csea = seacolor(256);
        cmap = landcolor(256).^1.3;
        
        % Determine zmin & zmax
        zcut = 0.5;
        zmin = nmedian(z(~isnan(z)),zcut/100);
        zmax = nmedian(z(~isnan(z)),1 - zcut/100);
        dz = zmax - zmin;
        
        % builds the colormap: concatenates seacolor and landcolor around 0
        if zmin < 0 && zmax > 0
            l = size(csea,1);
            r = size(cmap,1)*abs(zmin)/zmax/l;
            cmap = cat(1,interp1(1:l,csea,linspace(1,l,round(l*r)),'*linear'),cmap);
            
        elseif zmin < 0 && zmax < 0
            
            cmap = csea;
            
        end
        
    end
    
    % normalisation of Z using CMAP and convertion to RGB
    FaceVertexCData = ind2rgb(uint16((z - zmin)*(length(cmap)/dz)),cmap);
        
    varargout{1} = xv;
    varargout{2} = yv;
    varargout{3} = FaceVertexCData;
    varargout{4} = AlphaData;
    varargout{5} = cmap;
    
end

%-----------------------------------------------------------------------
% End Code
%-----------------------------------------------------------------------
end

%-----------------------------------------------------------------------
% nmedian
%-----------------------------------------------------------------------
function y = nmedian(x,n)
%NMEDIAN Generalized median filter
%	NMEDIAN(X,N) sorts elemets of X and returns N-th value (N normalized).
%	So:
%	   N = 0 is minimum value
%	   N = 0.5 is median value
%	   N = 1 is maximum value

if nargin < 2
    n = 0.5;
end
y = sort(x(:));
y = interp1(sort(y),n*(length(y)-1) + 1);
end

%-----------------------------------------------------------------------
% seacolor
%-----------------------------------------------------------------------
function y = seacolor(n)
%SEACOLOR Sea colormap adapted from NGDC ETOPO1
%
%	Author: Francois Beauducel <beauducel@ipgp.fr>

J = [ ...
    0.0392         0    0.4745
    0.1020         0    0.5373
    0.1020         0    0.5373
    0.1490         0    0.5961
    0.1490         0    0.5961
    0.1059    0.0118    0.6510
    0.1059    0.0118    0.6510
    0.0627    0.0235    0.7059
    0.0627    0.0235    0.7059
    0.0196    0.0353    0.7569
    0.0196    0.0353    0.7569
    0    0.0549    0.7961
    0    0.0549    0.7961
    0    0.0863    0.8235
    0    0.0863    0.8235
    0    0.1176    0.8471
    0    0.1176    0.8471
    0    0.1529    0.8745
    0    0.1529    0.8745
    0.0471    0.2667    0.9059
    0.0471    0.2667    0.9059
    0.1020    0.4000    0.9412
    0.1020    0.4000    0.9412
    0.0745    0.4588    0.9569
    0.0745    0.4588    0.9569
    0.0549    0.5216    0.9765
    0.0549    0.5216    0.9765
    0.0824    0.6196    0.9882
    0.0824    0.6196    0.9882
    0.1176    0.6980    1.0000
    0.1176    0.6980    1.0000
    0.1686    0.7294    1.0000
    0.1686    0.7294    1.0000
    0.2157    0.7569    1.0000
    0.2157    0.7569    1.0000
    0.2549    0.7843    1.0000
    0.2549    0.7843    1.0000
    0.3098    0.8235    1.0000
    0.3098    0.8235    1.0000
    0.3686    0.8745    1.0000
    0.3686    0.8745    1.0000
    0.5412    0.8902    1.0000
    0.5412    0.8902    1.0000
    0.7373    0.9020    1.0000
    ];

l = length(J);
if nargin < 1
    n = 256;
end
y = interp1(1:l,J,linspace(1,l,n),'*linear');


end

%-----------------------------------------------------------------------
% landcolor
%-----------------------------------------------------------------------
function y = landcolor(n)
%LANDCOLOR Land colormap
%
%	Author: Francois Beauducel <beauducel@ipgp.fr>
%	$Revision: 1.0.0 $   $Date: 2012/05/17 11:22:44 $

J = [ ...
    0.095678 0.53427 0.21682
    0.15785 0.5979 0.23274
    0.21286 0.64673 0.2514
    0.26411 0.68789 0.27268
    0.32959 0.72416 0.31308
    0.39794 0.75695 0.36038
    0.46153 0.7871 0.40624
    0.52108 0.81516 0.45135
    0.57702 0.84152 0.49547
    0.62973 0.86645 0.53891
    0.67946 0.89016 0.58187
    0.72647 0.91282 0.62427
    0.77095 0.93455 0.66619
    0.81306 0.95546 0.70772
    0.85292 0.97563 0.7489
    0.89066 0.99514 0.78976
    0.88379 0.98595 0.77038
    0.86389 0.96758 0.73236
    0.84615 0.94972 0.69623
    0.8303 0.93233 0.66186
    0.81612 0.91536 0.6291
    0.80341 0.8988 0.59784
    0.79201 0.8826 0.56795
    0.78191 0.86676 0.53946
    0.7729 0.85123 0.51224
    0.76479 0.83602 0.48615
    0.75747 0.8211 0.46111
    0.75084 0.80645 0.43704
    0.74506 0.79206 0.41414
    0.73981 0.77792 0.39211
    0.73501 0.76401 0.37089
    0.73068 0.75033 0.35052
    0.72683 0.73685 0.33106
    0.72042 0.72074 0.31228
    0.71032 0.70085 0.29417
    0.69761 0.67821 0.27694
    0.68489 0.65558 0.26026
    0.67235 0.63313 0.24418
    0.65997 0.61082 0.22889
    0.64775 0.58874 0.21406
    0.63568 0.56689 0.19983
    0.62376 0.54527 0.18622
    0.61197 0.52391 0.17299
    0.60033 0.50283 0.16046
    0.58881 0.48203 0.14832
    0.57742 0.46151 0.13667
    0.56616 0.44133 0.12555
    0.55502 0.4214 0.11472
    0.54398 0.4019 0.10456
    0.53306 0.38266 0.094633
    0.52226 0.36382 0.085242
    0.51155 0.3453 0.076179
    0.50095 0.32714 0.067515
    0.49045 0.30938 0.059259
    0.48005 0.29193 0.051294
    0.46973 0.27495 0.043796
    0.45951 0.25823 0.0365
    0.44938 0.24206 0.029715
    0.43934 0.22609 0.023063
    0.42938 0.21074 0.016949
    0.41951 0.19556 0.010917
    0.40971 0.18105 0.0054326
    0.4 0.16667 0
    ];

l = length(J);
if nargin < 1
    n = 256;
end
y = interp1(1:l,J,linspace(1,l,n),'*linear');

end

