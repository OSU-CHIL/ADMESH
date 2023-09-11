function p = rejectionMethod(p,DistFun,MeshFun,geps)

% Keep only d<0 points
p = p(DistFun(p) < geps, : );

% Probability to keep point
r0 = 1./(MeshFun(p)).^2;

% Rejection method
p = p(rand(size(p,1),1)<r0./max(r0),:);

end
