function [p,t,pix]=fixmesh(p,t,ptol)
%FIXMESH  Remove duplicated/unused nodes and fix element orientation.
%   [P,T]=FIXMESH(P,T)

%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

if nargin<3, ptol=1024*eps; end
if nargin>=2 && (isempty(p) | isempty(t)), pix=1:size(p,1); return; end %#ok<*OR2>

snap=max(max(p,[],1)-min(p,[],1),[],2)*ptol;

[~,ix,jx]=unique(round(p/snap)*snap,'rows','stable');

p=p(ix,:);

if nargin>=2
    
    t=reshape(jx(t),size(t));
    
    [pix,~,jx1]=unique(t);
    
    t=reshape(jx1,size(t));
    
    p=p(pix,:);
    
    pix=ix(pix);
    
    if size(t,2)==size(p,2)+1
        flip=simpvol(p,t)<0;
        t(flip,[1,2])=t(flip,[2,1]);
    end
    
end

    function v=simpvol(p,t)
        %SIMPVOL Simplex volume.
        %   V=SIMPVOL(P,T)
        
        %   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.
        
        switch size(p,2)
            case 1
                d12=p(t(:,2),:)-p(t(:,1),:);
                v=d12;
            case 2
                d12=p(t(:,2),:)-p(t(:,1),:);
                d13=p(t(:,3),:)-p(t(:,1),:);
                v=(d12(:,1).*d13(:,2)-d12(:,2).*d13(:,1))/2;
            case 3
                d12=p(t(:,2),:)-p(t(:,1),:);
                d13=p(t(:,3),:)-p(t(:,1),:);
                d14=p(t(:,4),:)-p(t(:,1),:);
                v=dot(cross(d12,d13,2),d14,2)/6;
            otherwise
                v=zeros(size(t,1),1);
                for ii=1:size(t,1)
                    A=zeros(size(p,2)+1);
                    A(:,1)=1;
                    for jj=1:size(p,2)+1
                        A(jj,2:end)=p(t(ii,jj),:);
                    end
                    v(ii)=det(A);
                end
                v=v/factorial(size(p,2));
        end
    end


end