function [new_coast]=join_cst(coast,tol)
% JOIN_CST  Makes continuous coastline from fragmented
%           segments by joining endpoints of segments less than
%           TOL distance apart.   The resulting coastline can be
%           used for blanking or filling by adding in a few extra
%           points, and the islands will also be polygons.
%           The output coastline is sorted largest
%           segment first, followed by sucessively smaller segments.
%
%   USAGE: [new_coast,slen]=join_cst(coast,tol);
%
%           Where COAST contains the x and y coastline positions
%           in the first two columns, and coastine segments are
%           separated by a row containing [NaN NaN] or [-99999. -99999.],
%           and TOL is the distance within which coastline segments will
%           be joined.  SLEN is a vector containing the length of the
%           resulting segments.

%
% Rich Signell               |  E-mail: rsignell@crusty.er.usgs.gov
% U.S. Geological Survey     |  Voice : 508-457-2229
% 384 Woods Hole Road        |  FAX   : 508-457-2310
% Woods Hole, MA  02543-1598 |  WWW   : http://crusty.er.usgs.gov/rsignell.html
%
% requires FIXCOAST.M

coast=fixcoast(coast);   %replace -99999. with NaN, eliminate dup NaNs, etc.

iseg=find(isnan(coast(:,1)));
nseg=length(iseg)-1;
irem=ones(nseg,1);
z=coast(:,1)+sqrt(-1)*coast(:,2);
zstart=z(iseg(1:nseg)+1);
zstop=z(iseg(2:(nseg+1))-1);
zends=[zstart(:) zstop(:)];

znew=nan+sqrt(-1)*nan;
k=0;

%while there are remaining segments, process
while(~isempty(find(irem, 1))),
    %start at 1st remaining segment
    ii=find(irem);
    id0=ii(1);
    id=id0;
    z0=z(iseg(id)+1);   %first point on segment
    zind=(iseg(id)+1):(iseg(id+1)-1);
    %line(real(z(zind)),imag(z(zind)),'erasemode','none','color','white')
    zc=z(iseg(id+1)-1); %last point on segment
    
    irem(id)=0;    %add current segment to list of segments used
    
    found_next=1;
    tried_reverse=0;
    while (found_next),
        ii=find(irem);     %indices of remaining segments
        nrem=length(ii);
        dgap=abs((ones(size(zends(ii,:))))*zc-zends(ii,:));
        id=find(dgap<tol, 1);
        if(~isempty(id)),
            idt=find(dgap==min(dgap(:)));                       % find closest points
            idt=idt(1);                                %take first if  equal dist
            if(idt>nrem),                              %found stop point of segment
                id=ii(idt-nrem);                   %find which segment we matched
                zi=(iseg(id+1)-1):-1:(iseg(id)+1);  %segment coordinates
            else                                       %found start point of segment
                id=ii(idt);                        %find which segment we matched
                zi=(iseg(id)+1):1:(iseg(id+1)-1);   %segment coordinates
            end
            ni=length(zi);
            %line(real(z(zi)),imag(z(zi)),'erasemode','none','color','white') %draw seg
            irem(id)=0;             %mask out segment
            zc=z(zi(ni));              %next connecting point
            zind=[zind zi];         %add segment to index list
            
        else
            if(~tried_reverse),
                zc=z0;
                nzind=length(zind);
                zind=zind(nzind:-1:1);
                tried_reverse=1;
            else
                found_next=0;
            end
        end
    end
    
    k=k+1;
    znew=[znew; nan+sqrt(-1)*nan; z(zind)];  %#ok<*AGROW> %add on new contatenated segment
end

znew=[znew; nan+sqrt(-1)*nan];  %add nan on end
iseg=find(isnan(znew));
nseg=length(iseg)-1;
slen=(iseg(2:(nseg+1))-iseg(1:nseg)-2);    %data pts in segments
[~,isort]=sort(slen);
isort=flipud(isort);
zind=[];
%sort largest to smallest
for i=1:nseg;
    zind=[zind (iseg(isort(i)):iseg(isort(i)+1))];
end
znew=znew(zind);
new_coast=[real(znew) imag(znew)];
new_coast=fixcoast(new_coast);

end