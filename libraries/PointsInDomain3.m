function IN = PointsInDomain3(X,Y,PTS)

% Loop over each polygon

% Test first polygon
IN = insidepoly(X,Y,PTS.Poly(1).x,PTS.Poly(1).y);

% Test Remaining polygons
for k = 2:length(PTS.Poly)
    
    in = insidepoly(X,Y,PTS.Poly(k).x,PTS.Poly(k).y);
    
    % If (x,y) are in an interior polygon we must remove it
    % If (x,y) are in an exterior polygon we must add it
    IN = xor(IN,in);
    
end

% Younghun: bankline constraints
if isfield(PTS,'Constraints')
for k = 1:length(PTS.Constraints)
    if PTS.Constraints(k).num == 18
        if isequal(PTS.Constraints(k).xy(1,:),PTS.Constraints(k).xy(end,:))
            in = insidepoly(X,Y,PTS.Constraints(k).xy(:,1),PTS.Constraints(k).xy(:,2));
            
            % If (x,y) are in an interior polygon we must remove it
%             IN(in) = 0;
%             IN = xor(IN,in);
        end
    end
    
end
end

end