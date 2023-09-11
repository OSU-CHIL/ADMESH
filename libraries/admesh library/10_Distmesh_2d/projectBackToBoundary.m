function p = projectBackToBoundary(phi,p)

pdist = phi.f(p); 

%ox = find(pdist > 0);
ox = pdist > 0;

% pgradx = phi.fx(p(ox,:));
% pgrady = phi.fy(p(ox,:));
% 
% p(ox,:)= p(ox,:)-[pdist(ox).*pgradx,pdist(ox).*pgrady];

p(ox,:)= p(ox,:)-[pdist(ox).*phi.fx(p(ox,:)),pdist(ox).*phi.fy(p(ox,:))];
