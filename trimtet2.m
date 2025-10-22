% 2024: 6/23  Trim away tets that are outside the domain
% Uses TRI of boundary

function [TET,x,y,z,t1] = trimtet2(xt,yt,zt,TRI,x,y,z,TET,xr,yr,zr);

to = cputime;
ip = inpolyhedron(TRI,[xt,yt,zt],[x,y,z]);
A = ip(TET); a = find(sum(A') > 0);  % all vertices in
keep = find(sum(A') > 0);
it = tsearchn([x,y,z],TET,[xr,yr,zr]);  a = unique(it)';
b = unique([keep,a]);
TET = TET(b,:);

np=length(x);
rpts=[];
for i=1:np
    if sum(sum(ismember(TET,i)))==0
        rpts=[rpts;i];
    end
end

pts=[x y z];
pts(rpts,:)=[];

for i=1:length(TET)
    for j=1:4
        TET(i,j)=find(ismember(pts, [x(TET(i,j)),y(TET(i,j)),z(TET(i,j))],'rows'));
    end
end

x=pts(:,1); y=pts(:,2); z=pts(:,3);

t1 = cputime - to;
