function [PTS,status,xyzIncluded,xyzFun] = ImportTemporaryFile(fileType,mesh)
% ImportTemporaryFile - Import file for meshing sub domain
%
% Syntax:  [PTS,status] = ImportTemporaryFile(fileType)
%
% Inputs:
%    mesh
%
% Outputs:
%
% See also: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% October 2013; Last revision: 21-October-2013

%------------------------------- BEGIN CODE -----------------------------------

xyzIncluded = 0;
status = 1;
PTS = [];

switch fileType
    
    case '.mat' 
        
        [PTS,status,xyzIncluded,xyzFun] = importMAT;

    case '.14'
        
        [PTS,status,xyzIncluded,xyzFun] = import14;
        
    case '.kml'
        
        [PTS,status,xyzIncluded,xyzFun] = importKMLv2;
        
end

% Reconvert coordinate conversion to have same reference
if isfield(mesh,'cpplon') && isfield(PTS,'cpplon')
    
    PTS = Meters2Geo(PTS);
    
    for k = 1:length(PTS.Poly)
        
        [PTS.Poly(k).x,PTS.Poly(k).y] = Geo2Meters(...
            PTS.Poly(k).x,PTS.Poly(k).y,mesh.cpplon,mesh.cpplat);
        
    end

end










end